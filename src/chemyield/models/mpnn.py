"""Directed message-passing network (Yang et al. 2019 style) for multi-component reactions.

Each component molecule is embedded, embeddings are sum-pooled within each role, and the role
vectors predict yield. Two switches make this a *diagnostic* model, not just a predictor:

  encoder     'graph' : D-MPNN over the molecular graph (default)
              'id'    : a learned vector per molecule, NO chemistry. A molecule never seen in
                        training gets a zero vector. If 'graph' cannot beat 'id' on held-out
                        molecules, the graph is not extracting transferable chemistry.
  interaction 'concat'      : role vectors are concatenated, an MLP combines them (default)
              'transformer' : role vectors are tokens (+ role embedding) that attend to each other
                              (self-attention) before the MLP, so halide x ligand x base
                              interactions are modelled explicitly.

A deep ensemble (n_models > 1) provides sigma. Only unique molecules in a minibatch are encoded.
"""

from __future__ import annotations

import copy

import numpy as np
import torch
from torch import nn

from ..features.graph import ATOM_DIM, BOND_DIM
from .base import MODELS, Model


class DMPNNEncoder(nn.Module):
    def __init__(self, atom_dim=ATOM_DIM, bond_dim=BOND_DIM, hidden=128, depth=3, dropout=0.0):
        super().__init__()
        self.depth = depth
        self.W_i = nn.Linear(atom_dim + bond_dim, hidden, bias=False)
        self.W_h = nn.Linear(hidden, hidden, bias=False)
        self.W_o = nn.Linear(atom_dim + hidden, hidden)
        self.act, self.drop = nn.ReLU(), nn.Dropout(dropout)

    def forward(self, x, ei, ea, rev, batch, n_mols, mol_ids=None):
        src, dst = ei
        h0 = self.act(self.W_i(torch.cat([x[src], ea], 1)))
        h = h0
        for _ in range(self.depth - 1):
            incoming = torch.zeros(x.size(0), h.size(1)).index_add(0, dst, h)
            # message for directed edge u->v: everything arriving at u EXCEPT from v (no backtracking)
            m = incoming[src] - h[rev]
            h = self.drop(self.act(h0 + self.W_h(m)))
        incoming = torch.zeros(x.size(0), h.size(1)).index_add(0, dst, h)
        atoms = self.drop(self.act(self.W_o(torch.cat([x, incoming], 1))))
        pooled = torch.zeros(n_mols, atoms.size(1)).index_add(0, batch, atoms)
        counts = torch.bincount(batch, minlength=n_mols).clamp(min=1).unsqueeze(1)
        return pooled / counts


class ReactionNet(nn.Module):
    def __init__(self, n_roles, extra_dim=0, hidden=128, depth=3, ffn_hidden=256, dropout=0.1,
                 encoder="graph", interaction="concat", n_layers=1, heads=4, id_lut=None):
        super().__init__()
        self.n_roles, self.hidden = n_roles, hidden
        self.encoder_kind, self.interaction = encoder, interaction
        if encoder == "graph":
            self.enc = DMPNNEncoder(hidden=hidden, depth=depth, dropout=dropout)
        elif encoder == "id":
            self.register_buffer("lut", id_lut)       # global molecule id -> compact index (0 = unknown)
            self.id_emb = nn.Embedding(int(id_lut.max()) + 1, hidden, padding_idx=0)
        else:
            raise ValueError(f"encoder must be 'graph' or 'id', got '{encoder}'")
        if interaction == "transformer":
            self.role_emb = nn.Parameter(torch.randn(n_roles, hidden) * 0.02)
            layer = nn.TransformerEncoderLayer(d_model=hidden, nhead=heads, dim_feedforward=2 * hidden,
                                               dropout=dropout, batch_first=True, norm_first=True)
            self.tr = nn.TransformerEncoder(layer, n_layers, enable_nested_tensor=False)
        elif interaction != "concat":
            raise ValueError(f"interaction must be 'concat' or 'transformer', got '{interaction}'")
        self.ffn = nn.Sequential(
            nn.Linear(n_roles * hidden + extra_dim, ffn_hidden), nn.ReLU(), nn.Dropout(dropout),
            nn.Linear(ffn_hidden, ffn_hidden // 2), nn.ReLU(), nn.Linear(ffn_hidden // 2, 1))

    def embed(self, graph):
        if self.encoder_kind == "graph":
            return self.enc(*graph)
        ids = graph[6]
        idx = torch.zeros_like(ids)
        inb = ids < self.lut.numel()
        idx[inb] = self.lut[ids[inb]]
        return self.id_emb(idx)

    def forward(self, graph, p_rxn, p_role, p_mol, n_rxn, extra=None):
        emb = self.embed(graph)
        flat = torch.zeros(n_rxn * self.n_roles, self.hidden).index_add(
            0, p_rxn * self.n_roles + p_role, emb[p_mol])
        tokens = flat.view(n_rxn, self.n_roles, self.hidden)
        if self.interaction == "transformer":
            tokens = self.tr(tokens + self.role_emb)
        z = tokens.reshape(n_rxn, self.n_roles * self.hidden)
        if extra is not None:
            z = torch.cat([z, extra], 1)
        return self.ffn(z).squeeze(1)


@MODELS.register("mpnn")
class MPNNRegressor(Model):
    input_kind = "graph"

    def __init__(self, hidden=128, depth=3, ffn_hidden=256, dropout=0.1, lr=1e-3, weight_decay=1e-5,
                 batch_size=128, epochs=60, patience=10, n_models=1, seed=0, roles="auto",
                 product=True, encoder="graph", interaction="concat", n_layers=1, heads=4,
                 threads=None, verbose=False):
        self.cfg = dict(hidden=hidden, depth=depth, ffn_hidden=ffn_hidden, dropout=dropout,
                        encoder=encoder, interaction=interaction, n_layers=n_layers, heads=heads)
        self.lr, self.wd, self.bs, self.epochs, self.patience = lr, weight_decay, batch_size, epochs, patience
        self.n_models, self.seed, self.verbose = n_models, seed, verbose
        self.roles, self.product = roles, product   # read by the runner to build GraphFeaturizer
        if threads:
            torch.set_num_threads(int(threads))

    # ---- extra (tabular side) features: impute + standardise with TRAIN stats ----------------
    def _prep_extra(self, gd, fit=False):
        if gd.extra is None:
            return gd
        e = gd.extra.astype(np.float32)
        if fit:
            self.e_mean_ = np.nan_to_num(np.nanmean(e, 0), nan=0.0)
            filled = np.where(np.isnan(e), self.e_mean_, e)
            self.e_std_ = np.where(filled.std(0) > 0, filled.std(0), 1.0).astype(np.float32)
        e = np.where(np.isnan(e), self.e_mean_, e)
        return gd.subset_extra((e - self.e_mean_) / self.e_std_)

    def _predict_one(self, net, gd, bs=512):
        net.eval()
        out = []
        with torch.no_grad():
            for s in range(0, len(gd), bs):
                out.append(net(*gd.batch(np.arange(s, min(s + bs, len(gd))))).numpy())
        return np.concatenate(out)

    @staticmethod
    def _id_lut(gd):
        """global molecule id -> 1..K for molecules present in TRAINING reactions, else 0 (unknown)."""
        lut = torch.zeros(len(gd.table.graphs), dtype=torch.long)
        for j, m in enumerate(sorted({m for it in gd.items for _, m in it})):
            lut[m] = j + 1
        return lut

    def fit(self, gd, y, gd_val=None, y_val=None):
        y = np.asarray(y, np.float32)
        self.y_mean_, self.y_std_ = float(y.mean()), float(y.std() or 1.0)
        gd = self._prep_extra(gd, fit=True)
        gd_val = self._prep_extra(gd_val) if gd_val is not None else None
        extra_dim = 0 if gd.extra is None else gd.extra.shape[1]
        id_lut = self._id_lut(gd) if self.cfg["encoder"] == "id" else None
        yt = torch.from_numpy((y - self.y_mean_) / self.y_std_)
        self.nets_ = []
        for k in range(self.n_models):
            torch.manual_seed(self.seed + k)
            rng = np.random.default_rng(self.seed + k)
            net = ReactionNet(gd.n_roles, extra_dim, id_lut=id_lut, **self.cfg)
            opt = torch.optim.AdamW(net.parameters(), lr=self.lr, weight_decay=self.wd)
            steps = self.epochs * int(np.ceil(len(gd) / self.bs))
            sched = torch.optim.lr_scheduler.OneCycleLR(opt, max_lr=self.lr, total_steps=steps,
                                                        pct_start=0.1)
            best, best_state, bad = np.inf, None, 0
            for ep in range(self.epochs):
                net.train()
                perm = rng.permutation(len(gd))
                tot = 0.0
                for s in range(0, len(perm), self.bs):
                    idx = perm[s:s + self.bs]
                    if len(idx) < 2:      # skip a trailing singleton batch
                        sched.step()
                        continue
                    pred = net(*gd.batch(idx))
                    loss = nn.functional.mse_loss(pred, yt[idx])
                    opt.zero_grad()
                    loss.backward()
                    nn.utils.clip_grad_norm_(net.parameters(), 5.0)
                    opt.step()
                    sched.step()
                    tot += loss.item() * len(idx)
                if gd_val is not None:
                    pv = self._predict_one(net, gd_val) * self.y_std_ + self.y_mean_
                    score = float(np.sqrt(np.mean((pv - np.asarray(y_val)) ** 2)))
                else:
                    score = tot / len(gd)
                if self.verbose:
                    print(f"  member {k} epoch {ep:3d} train_mse {tot/len(gd):.4f} val_rmse {score:.3f}")
                if score < best - 1e-4:
                    best, best_state, bad = score, copy.deepcopy(net.state_dict()), 0
                else:
                    bad += 1
                    if bad >= self.patience:
                        break
            net.load_state_dict(best_state)
            self.nets_.append(net)
        return self

    def predict(self, gd):
        gd = self._prep_extra(gd)
        P = np.stack([self._predict_one(n, gd) for n in self.nets_]) * self.y_std_ + self.y_mean_
        return P.mean(0), (P.std(0) if len(self.nets_) > 1 else None)
