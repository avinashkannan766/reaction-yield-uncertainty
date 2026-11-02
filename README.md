# chemyield

Config-driven reaction **yield prediction with uncertainty**. Built to be dataset-agnostic: swapping
the dataset means editing one block of a YAML file, not code.

```
raw source ──adapter──▶ standard table ──clean──▶ clean.csv ──split──▶ featurize ──▶ model ──▶ metrics
 (any format)           (components+product       (+ report)   (OOD)    (fit on       (GBM / D-MPNN   + conformal
                         +yield)                                         train only)    ensembles)     intervals
```

## Results at a glance

Two very different datasets, same pipeline (details and caveats in the sections below; 3 seeds, untuned
models, so differences below ~0.1 R² on group splits are within noise).

| | TDC Buchwald-Hartwig (HTE grid; 3,955 reactions, 46 distinct molecules) | TDC USPTO yields (patent-mined; 50k raw -> 46,929 reactions) |
|---|---|---|
| best test R², random split | 0.94 (LightGBM) | 0.17 (5-model LightGBM) |
| best test R², unseen groups | 0.47-0.52 unseen scaffold, 0.59-0.65 unseen aryl halide, 0.72-0.76 unseen additive, 0.27-0.41 unseen ligand (leave-one-out) | 0.13 (unseen k-means cluster of products) |
| 90 % interval coverage with ordinary calibration, unseen groups | 0.42-0.77: **under-covers** | 0.89-0.90: fine |
| with group-aware calibration | 0.84-1.00, at 2-4x the width | no change (0.906) |
| are the intervals informative? | yes except on unseen scaffolds (width 45-68 vs 75-83 for a constant predictor; 85-94 on scaffolds, no better than the constant) | barely (width 70-73 vs 75) |

What this shows:
* **A random split can measure memorisation, not chemistry.** On Buchwald a model with *no chemistry at all*
  (a learned vector per molecule) matches the best model on a random split (R² 0.957) and collapses on
  unseen aryl halides (0.26 vs ~0.6 for the graph model).
* **Tree models on fingerprints are the reference to beat.** The D-MPNN is behind them on the random,
  scaffold, additive and ligand splits and roughly on par on unseen aryl halides and bases; it only
  clearly beats the no-chemistry baseline on unseen aryl halides. All graph models untuned.
* **Conformal coverage fails under group shift unless calibration is group-aware** (Buchwald), and that
  fix costs interval width. Whether it is needed depends on how much real signal the model has (USPTO:
  not needed, but the intervals say little).
* **Always compare interval width with a constant predictor's**: coverage alone is not usefulness.

## Quick start

```bash
python -m venv .venv && source .venv/bin/activate     # Windows: .venv\Scripts\Activate.ps1
pip install torch --index-url https://download.pytorch.org/whl/cpu
pip install -e .[dev]        # rdkit, lightgbm, scikit-learn, pandas, pyyaml
pytest                       # 51 tests, ~1.5 min

# 1) pipeline smoke test on SYNTHETIC data (yields are invented; proves the plumbing only)
python -m chemyield.synthetic data/raw/synthetic/synthetic_bh.csv
python -m chemyield.benchmark --config configs/synthetic.yaml

# 2) the real thing: TDC Buchwald-Hartwig (Yields task, Harvard Dataverse file 6175640)
mkdir -p data/raw/tdc_buchwald
#    download buchwald-hartwig.pkl into that folder. Dataverse returns 403 to Python's default
#    User-Agent, so send a browser-style one (or use curl -L -A "Mozilla/5.0", or your browser):
python -c "import urllib.request as u; r=u.Request('https://dataverse.harvard.edu/api/access/datafile/6175640', headers={'User-Agent':'Mozilla/5.0'}); open('data/raw/tdc_buchwald/buchwald-hartwig.pkl','wb').write(u.urlopen(r).read())"

# clean only (no training), then READ the report before trusting any number
python -c "from chemyield.config import load_config; from chemyield.run import prepare_data; prepare_data(load_config('configs/tdc_buchwald.yaml'), force=True)"
cat data/processed/tdc_buchwald/cleaning_report.md

python -m chemyield.run --config configs/tdc_buchwald.yaml          # one model, one split
python -m chemyield.benchmark --config configs/tdc_buchwald.yaml --only mean_baseline ridge_morgan rf_morgan lgbm_morgan lgbm_morgan_desc
python -m chemyield.benchmark --config configs/tdc_buchwald.yaml --only mpnn_ensemble3   # slow on CPU
```

Outputs: `data/processed/<name>/clean.csv` + `cleaning_report.md`, `runs/<name>/<exp>/…`
(config, per-seed metrics and predictions), `results/<name>_benchmark.md` (the ablation table).
Note: `benchmark` writes its table only at the very end, and each call rewrites
`results/<name>_benchmark.md` with just the experiments it ran. Use `--suffix _x` to keep calls
separate, `--only <experiments>` / `--only-splits <substrings>` to run a subset, and, if a long run is
interrupted, rebuild the table from the per-experiment results already saved under `runs/`:
`python -m chemyield.collect --dataset <name> [--contains gkf] --out results/<name>_collected.md`.

### Second dataset: USPTO yields (optional, large)

```bash
mkdir -p data/raw/uspto      # 432 MB; Windows: curl.exe -L -A "Mozilla/5.0" -o data\raw\uspto\uspto_yields.pkl <url>
curl -L -A "Mozilla/5.0" -o data/raw/uspto/uspto_yields.pkl "https://dataverse.harvard.edu/api/access/datafile/4186956"

# vet it BEFORE trusting it: cleaning report + yield distribution + can group splits be built?
python -m chemyield.profile --config configs/uspto_yields.yaml --sample 50000 --group-by cluster:100 cluster:300 --out results/uspto_clusters.md
# tree models only (graph models are too slow on a CPU at this size); mean/ridge/LightGBM on the 3 cluster splits took ~30 min, random forest / 5-model / group-K-fold take far longer
python -m chemyield.benchmark --config configs/uspto_yields.yaml --only mean_baseline ridge_morgan lgbm_morgan --suffix _uspto
```

`profile` is the "investigate before committing" step for any new dataset: it runs the adapter and
cleaning and reports what was dropped, yield pile-ups, component counts and whether each `group_by` can
actually produce a test set.

## Standard schema (what every adapter must produce)

| column | meaning |
|---|---|
| `reaction_id` | string |
| `components` | list of `(role, smiles)` – reactants, reagents, ligand, base, solvent … |
| `product` | `.`-joined product SMILES |
| `yield` | float (cleaning normalises the scale to 0–100) |
| extras | e.g. `temperature`, listed in `keep_cols` |

Roles are what make the representations meaningful: fingerprints, descriptors and the MPNN all pool
**per role**, so "which molecule played which part" is explicit.

## Swapping datasets

1. Pick an adapter in the config: `csv` (generic; reaction-SMILES column, `role_cols`, or
   reactants/reagents columns), `tdc_pickle`, or write a 20-line class and `@ADAPTERS.register("name")`.
2. Point `dataset.params` at the file and name the columns (`dataset.sample: {n: ..., seed: ...}` takes a random subset of the raw rows for big files).
3. Run. Cleaning, splits, features, models and the benchmark table are unchanged.
   `tests/test_end_to_end.py::test_swapping_dataset_needs_only_the_dataset_block` proves this.

### The real TDC layout

`buchwald-hartwig.pkl` is a DataFrame with columns `ID`, `X`, `Y`: `X` is a dict
`{product, catalyst, reactant}` and `Y` is the yield **as a fraction (0–1)**. The `reactant` string is
one dot-joined bag containing the aryl halide, amine, Pd precatalyst, phosphine ligand, base and
additive. `configs/tdc_buchwald.yaml` therefore carries `cleaning.role_rules` (SMARTS → role, first
match wins) that split the bag into `catalyst / ligand / base / aryl_halide / amine / additive`.
`tests/test_adapters.py` pins this against the first real row.

Check the **role signatures** section of `cleaning_report.md`: each reaction should have one component
per role. Known ambiguity: a molecule with both a C–halogen bond and an NH₂ (e.g. 4-chloroaniline) is
classed as the aryl halide by these rules; it would show up as a missing amine in the signatures.

## Cleaning (auditable)

Missing values → yield scale/range (clip small overshoots, drop impossible) → canonical SMILES
(atom maps stripped; metal complexes kept via unsanitized fallback and counted) → role rules →
element-balance check → de-duplication (mean yield; replicates that disagree are flagged via
`yield_std`). Every step reports rows removed.

## Splits (`split:` in config)

* `random` – optimistic reference only.
* `group` – the **test set contains only groups never seen in training**: `group_by: scaffold`
  (product Murcko scaffold), `generic_scaffold` (atoms->C, bonds->single: topology only),
  `cluster:<k>` (k-means on product fingerprints; for DIVERSE data where scaffolds / products are
  near-unique and median group size is 1), `role:<name>` (e.g. `role:aryl_halide`, `role:additive`),
  `product`, `column:<name>`. Val (early stopping) and cal (conformal) are random draws from the remaining rows.
* `loo` – **leave-one-group-out**: `group_by: role:ligand` with seeds `[0,1,2,3]` holds out each of the
  4 ligands in turn (fold = `seed % n_groups`); the whole group is the test set. Use it for roles with
  too few members for a 20 % test set (ligand, base). The benchmark lists per-fold R² in
  `results/<name>_benchmark.md` and `results/<name>_loo_folds.csv`; a split may carry its own `seeds:`.
* `val_by_group: true` (on `group` / `loo` splits): validation (early stopping) and calibration
  (conformal) sets are **also whole groups unseen in training**, instead of random rows from seen
  groups. This is what a new-chemistry deployment looks like. Effect on SYNTHETIC data (6 seeds, single
  LightGBM, unseen aryl halides): 90 % coverage rose 0.60 → 0.79, but intervals roughly doubled in width
  and R² fell 0.61 → 0.46 because fewer groups remain for training; coverage stayed below 0.9 and noisy
  (2 calibration groups). On the real Buchwald file the effect is larger and cleaner - see "Further
  findings" below. Needs enough groups: skipped with
  "too few usable groups" for ligand/base leave-one-out and for scaffold (5 groups). Appears as
  `<split>+gval` in the results.
* A `group` split that cannot be built (too few groups) is skipped by the benchmark with a message.
* Featurizers fit on train only (vocabularies, standardisation), so held-out components look
  exactly as they would in deployment.

## Representations (`features:` in config) – each is one switch

| name | what it adds |
|---|---|
| `morgan` | per-role Morgan count fingerprints (+ product block) |
| `rdkit_desc` | per-role mean physchem descriptors |
| `conditions` | numeric columns (temperature, time…) standardised; multi-hot component identity per role |
| `lookup` | **your own descriptor table** joined by SMILES per role: Kamlet–Taft solvent params, base pKaH, ligand %Vbur/Tolman, Boltzmann-averaged 3D descriptors … |
| (model `mpnn`) | learned D-MPNN graph embeddings per component; the features above can be fed in as side inputs |

The `lookup` featurizer is the intended hook for the physical-condition / 3D-conformer ideas: compute
descriptors once per unique molecule (cheap – HTE data has few unique reagents), save a CSV, add a
`lookup` entry, and the ablation table tells you whether it helped.

## Models (`model:` in config)

`mean` (floor), `ridge`, `rf` (σ = tree std), `lightgbm` (`n_models>1` → bagged ensemble σ),
`mpnn` (D-MPNN, `n_models>1` → deep ensemble σ). Add one with `@MODELS.register`.

### Diagnostic switches on `mpnn`

| option | values | purpose |
|---|---|---|
| `encoder` | `graph` (default) / `id` | `id` = same network but a learned vector per molecule and **no chemistry**; molecules unseen in training get a zero vector. If `graph` cannot beat `id` on held-out molecules, the graph is not giving transferable chemistry. |
| `interaction` | `concat` (default) / `transformer` | `transformer`: role vectors become tokens (+ role embedding) and attend to each other (`n_layers`, `heads`) before the MLP, modelling halide × ligand × base interactions. |
| `product` | `true` / `false` | include the product graph. Redundant here (fixed by halide + amine). |

Buchwald experiments ready in `configs/tdc_buchwald.yaml`: `mpnn_noproduct`, `id_embed_noproduct`,
`mpnn_transformer_noproduct` (plus the LOO ligand/base splits):
```bash
python -m chemyield.benchmark --config configs/tdc_buchwald.yaml --only lgbm_morgan mpnn_ensemble3 mpnn_noproduct id_embed_noproduct mpnn_transformer_noproduct
```

### The graph model (`mpnn`)

Own PyTorch implementation of a Chemprop-style directed MPNN (not the `chemprop` package, and not
validated against its numbers). Every component, and the product, is a separate molecular graph
(56 atom features, 11 bond features). Messages travel along directed bonds without backtracking;
atoms are mean-pooled per molecule, molecules are sum-pooled per role, role vectors are concatenated
(fixed slot per role) and an MLP predicts the standardised yield. σ comes from the ensemble.

| hyperparameter | default | value in `tdc_buchwald.yaml` (`mpnn_ensemble3`) |
|---|---|---|
| `hidden` / `depth` | 128 / 3 | default |
| `ffn_hidden` / `dropout` | 256 / 0.1 | default |
| `lr` (OneCycle peak, 10 % warm-up) / `weight_decay` | 1e-3 / 1e-5 (AdamW) | default |
| `batch_size` | 128 | default |
| `epochs` / `patience` | 60 / 10 | 40 / 8 |
| `n_models` | 1 | 3 |
| `product` (include product graph) | true | default |
| loss / grad-clip | MSE / 5.0 | – |

None of these were tuned. Not implemented: GIN/GINE, AttentiveFP, attention readout, FiLM
conditioning, 3D/equivariant models, CGR (no atom mapping).

## Uncertainty

Split-conformal 90 % intervals (σ-normalised when the model provides σ). Coverage is reported on the
test set. **Under group (OOD) splits, coverage drops below the target** – conformal's guarantee needs
exchangeability. That gap is a result worth reporting, not a bug.

## Results so far: TDC Buchwald-Hartwig (3 seeds, mean ± std, test R²)

Measured with the configs in this repo; hyperparameters untuned for every model.

| experiment | random | unseen scaffold | unseen aryl halide | unseen additive |
|---|---|---|---|---|
| mean baseline | -0.003 ± 0.003 | -0.017 ± 0.001 | -0.029 ± 0.038 | -0.066 ± 0.058 |
| ridge + Morgan | 0.675 ± 0.005 | 0.501 ± 0.003 | 0.504 ± 0.141 | 0.550 ± 0.046 |
| random forest + Morgan | 0.920 ± 0.007 | 0.522 ± 0.014 | 0.590 ± 0.301 | 0.763 ± 0.069 |
| LightGBM + Morgan | 0.941 ± 0.001 | 0.474 ± 0.025 | 0.588 ± 0.283 | 0.721 ± 0.052 |
| LightGBM + Morgan + descriptors | 0.939 ± 0.002 | 0.452 ± 0.063 | 0.650 ± 0.229 | 0.723 ± 0.053 |
| D-MPNN ensemble (3) | 0.818 ± 0.009 | 0.163 ± 0.046 | 0.597 ± 0.203 | 0.627 ± 0.080 |

Reading it: random splits overstate performance (0.94 vs ~0.5 for unseen scaffolds); the graph model
is *behind* the tree models, most clearly on unseen scaffolds. Both follow-up hypotheses were then
tested (next section): the product graph is not the cause, and an ID-embedding baseline with no
chemistry matches the best model on random splits. Full tables incl. coverage are in
`results/tdc_buchwald_benchmark.md`.

## Further findings: TDC Buchwald-Hartwig (3 seeds; leave-one-out = every group once)

Measured on the real 3,955-reaction file. Everything is untuned. Differences below ~0.1 R² are inside
the seed/group noise on the group splits. "Hypothesis" marks explanations that were not tested.

**R² for the graph-model diagnostics** (all graph variants use 3-member ensembles; ID = learned vector
per molecule, no chemistry, unseen molecule -> zero vector; "-np" = product graph off):

| model | random | unseen scaffold | unseen aryl halide | unseen additive | leave-one-ligand-out | leave-one-base-out |
|---|---|---|---|---|---|---|
| LightGBM + Morgan | 0.941 | 0.474 | 0.588 | 0.721 | 0.406 | 0.520 |
| D-MPNN | 0.818 | 0.164 | 0.598 | 0.621 | 0.128 | 0.557 |
| D-MPNN -np | 0.817 | 0.218 | 0.612 | 0.603 | 0.093 | 0.553 |
| ID embedding -np | **0.957** | 0.320 | 0.258 | 0.592 | 0.309 | 0.503 |
| transformer-interaction -np | 0.879 | 0.139 | **0.690** | 0.708 | 0.210 | 0.319 |

* A model with **no chemistry** (ID) matches the best model on a random split: on this ~46-molecule
  grid a random split measures memorisation of combinations, not chemistry.
* Where chemistry must transfer (unseen aryl halide) the graph beats the ID baseline by ~0.35 R², so
  the graph does carry transferable information there; it does not beat the tree models.
* The product graph made no consistent difference (the original "product graph causes the scaffold
  collapse" hypothesis is mostly not supported).
* Transformer interaction helped with many groups (+0.06 to +0.10 on random / halide / additive) and
  hurt with few (-0.08 scaffold, -0.23 base leave-one-out, driven by one fold).

**Leave-one-ligand-out, R² per held-out ligand** (named from structure):

| held-out ligand | LightGBM | D-MPNN | ID baseline |
|---|---|---|---|
| t-BuXPhos | 0.77 | 0.72 | 0.71 |
| XPhos | -0.98 | -1.03 | -1.21 |
| t-BuBrettPhos | 0.93 | 0.79 | 0.86 |
| AdBrettPhos | 0.90 | **0.03** | 0.88 |

XPhos fails for every model (nothing in the other three ligands predicts its yield level). On
AdBrettPhos "ignore the ligand" works (0.88) while the graph model extrapolates badly (0.03) -
hypothesis: an unseen molecule gives the graph a confident but wrong embedding, whereas fingerprint
bits never seen in training simply contribute zeros.

**Uncertainty (90 % target)** - 5-model LightGBM vs single LightGBM, test coverage (interval width):
unseen scaffold 0.71 (36.6) vs 0.42 (19.3); unseen aryl halide 0.76 (27.9) vs 0.60 (18.5); unseen
additive 0.74 (24.2) vs 0.64 (18.6); on the random split both are ~0.91 at ~19-20, so the ensemble's
sigma widens intervals only on unfamiliar chemistry, at no accuracy cost (R² 0.942 vs 0.941).

**Coverage of the 90 % intervals on every split** (mean over seeds / folds; width in yield points;
random-row calibration):

| split | random forest | LightGBM | 5-model LightGBM |
|---|---|---|---|
| random | 0.90 (19.9) | 0.92 (20.2) | 0.91 (18.8) |
| unseen scaffold | 0.67 (37.1) | 0.42 (19.3) | 0.71 (36.6) |
| unseen aryl halide | 0.77 (28.9) | 0.60 (18.5) | 0.76 (27.9) |
| unseen additive | 0.76 (23.6) | 0.64 (18.6) | 0.74 (24.2) |
| leave-one-ligand-out | 0.65 (23.1) | 0.67 (18.4) | 0.75 (23.9) |
| leave-one-base-out | 0.57 (20.7) | 0.58 (18.8) | 0.67 (24.1) |

Coverage reaches the target only on the random split; every group split under-covers.

**Group-aware calibration (`val_by_group`)**, coverage (width) with R² for unseen aryl halide / additive:

| split | calibration | random forest | LightGBM | 5-model LightGBM |
|---|---|---|---|---|
| aryl halide | random rows | 0.77 (28.9), R² 0.590 | 0.60 (18.5), 0.588 | 0.76 (27.9), 0.642 |
| aryl halide | whole groups | 0.84 (39.2), 0.614 | 0.90 (56.4), 0.541 | 0.89 (55.1), 0.536 |
| additive | random rows | 0.76 (23.6), 0.763 | 0.64 (18.6), 0.721 | 0.74 (24.2), 0.731 |
| additive | whole groups | 0.87 (37.8), 0.741 | 0.86 (43.5), 0.729 | 0.86 (42.2), 0.720 |

Calibrating on whole unseen groups brings coverage to 0.84-0.90 (near the target) at the cost of much
wider intervals (38-56 yield points), and LightGBM loses some R² because fewer groups remain for
training; random forest keeps its R² and gives the narrowest intervals, so it is the reference model
for new-chemistry prediction here. Caveat: 3 seeds, coverage std 0.03-0.09, and the guarantee still
assumes a new test group resembles the calibration groups.

## Group-K-fold calibration (`conformal.method: group_kfold`)

Calibration residuals come from K refits of the TRAIN data: the train groups are split into K folds,
the model is refit on K-1 folds and predicts the held-out fold, so every train row gets a prediction
from a model that never saw its group. The conformal quantile is taken from those residuals; the final
model is trained on all training rows (the cal rows go back into training, so no groups are given up,
unlike `val_by_group`). It also works with few groups (ligand/base leave-one-out, scaffold), where
`val_by_group` is skipped. Cost: K+1 fits per seed. Config (per benchmark split):
`conformal: {method: group_kfold, k: 5}`; it shows up as `<split>+gkf`. On synthetic data it gives
higher coverage and wider intervals than random-row calibration (unit test).

**Measured on TDC Buchwald-Hartwig** (3 seeds / folds; 90 % target; coverage (width in yield points);
rebuilt from `runs/` with `python -m chemyield.collect --dataset tdc_buchwald --contains gkf`).
The **mean baseline** row is the yardstick: a model that predicts the training mean for everything
also reaches ~0.90 coverage, at a width of ~75-83. An interval is only informative if it is clearly
narrower than that.

| split | mean baseline | random forest | LightGBM | 5-model LightGBM | D-MPNN (3) |
|---|---|---|---|---|---|
| unseen scaffold | 0.92 (83.4) | 1.00 (94.2) | 0.97 (84.7) | 0.99 (91.4) | 0.87 (63.8) |
| unseen aryl halide | 0.92 (76.4) | 0.99 (61.7) | 0.92 (50.2) | 0.98 (55.1) | 0.91 (48.3) |
| unseen additive | 0.89 (75.0) | 0.95 (44.7) | 0.93 (47.0) | 0.93 (44.5) | 0.92 (52.4) |
| leave-one-ligand-out | 0.90 (76.7) | 0.88 (53.7) | 0.87 (46.1) | 0.89 (50.8) | 0.80 (45.0) |
| leave-one-base-out | 0.90 (77.6) | 0.96 (67.5) | 0.91 (54.7) | 0.94 (58.9) | 0.88 (47.4) |

Findings:
* Coverage is restored on every group split (0.4-0.7 with random-row calibration -> ~0.87-1.00), and
  it works where `val_by_group` cannot run (scaffold, ligand, base).
* On **unseen scaffolds the intervals are useless**: width 85-94 is no narrower than the constant
  predictor (83). Hypothesis (not tested): each refit trains on only 3 of the 4 training groups, so its
  residuals are larger than the final model's real error.
* Unseen additive / aryl halide / ligand / base give intervals clearly narrower than the baseline
  (about 45-68 vs 75-78) at near-target coverage. Random forest on additive: 0.95 (44.7), better
  coverage than 5-model LightGBM at the same width and a higher R² (0.770 vs 0.731).
* Compared with `val_by_group` (aryl halide / additive, same seeds): coverage is higher but intervals
  are wider (RF aryl halide 0.99 (61.7) vs 0.84 (39.2); LightGBM 0.92 (50.2) vs 0.90 (56.4), the one
  case where it is narrower). So it is not a general replacement; it is the safeguard for few groups.
* Point accuracy is unchanged (R² within seed noise of the random-row-calibration run).
* Coverage std across seeds / folds is 0.04-0.19; with 3 seeds and 3-15 groups these differences
  are indicative, not significant. No coverage guarantee holds for an arbitrary new chemical group.

## Second dataset: TDC USPTO yields (patent-mined, diverse)

Chosen to be the opposite of the Buchwald HTE grid: real but text-mined yields from US patents
(TDC `uspto_yields`, 853,638 reactions, Harvard Dataverse file id 4186956, 432 MB pickle; TDC's page
lists it as CC0 - not checked at the source). Config: `configs/uspto_yields.yaml`. **Vetted before use**
with `python -m chemyield.profile` (adapter + cleaning + group-structure report):

* **Layout** `[ID, X, Y]`, `X = {catalyst, reactant, product}`, atom-mapped SMILES. The `catalyst` key
  is not just catalysts (first row: dichloromethane), so it is renamed to role `agent`; inside
  `reactant`, fragments WITH atom maps are real reactants and unmapped ones (AlCl3 ions, ...) become
  `agent` too (`atom_map_roles: true`). Ionic reagents are split into ions, which loses a little.
* **Yields** are fractions with junk (min -1.0, max 693.65); cleaning drops values outside ~[-5, 105] %.
* **Sample**: the benchmark uses a random 50,000 of the RAW rows (`dataset.sample`) -> **46,929** clean
  reactions (96 out-of-range, 2,975 duplicates merged). 48,875 distinct reactants, 2,598 distinct
  agents, 46,309 distinct products. Profiles at 20k / 50k / 60k raw rows gave the same yield
  distribution (mean 64.0, std 26.0, ~20 % >= 90). The duplicate share grew with sample size (2.8 % at
  20k, 7.1 % at 60k), so it would be larger on the full file; the full file was never cleaned.
* **What it does NOT have**: temperature, time, amounts or a patent id - only structures and an
  unstructured agent bag (6.8 % of reactions have no agent). Condition-embedding ideas cannot be
  tested on it.
* **Why cluster splits**: structural groups are near-singletons (median group size 1; at 60k rows
  `generic_scaffold` had 13,825 groups with one holding 12.3 %), so a "held-out scaffold" is nearly a
  random split. `cluster:100` (k-means on product fingerprints) gives 100 groups, largest 3.1 %,
  93 with >= 20 reactions, median 374: a 20 % test set is ~20 whole clusters.

**Results** (test, mean ± std over 3 seeds unless marked; 90 % target; coverage (width in yield
points); random-row calibration unless noted):

| model | random: R² | random: coverage (width) | unseen cluster: R² | unseen cluster: coverage (width) |
|---|---|---|---|---|
| mean baseline | 0.000 | 0.91 (75.3) | -0.005 ± 0.006 | 0.91 (75.0) |
| ridge + Morgan | 0.092 ± 0.010 | 0.89 (71.3) | 0.071 ± 0.007 | 0.90 (71.9) |
| random forest | 0.153 ± 0.008 | 0.90 (70.1) | 0.110 (1 seed) | 0.90 (69.8) |
| LightGBM | 0.160 ± 0.008 | 0.89 (69.5) | 0.131 ± 0.015 | 0.90 (70.0) |
| 5-model LightGBM | 0.172 ± 0.008 | 0.90 (73.1) | 0.134 (1 seed) | 0.89 (72.3) |

LightGBM on unseen clusters, other calibrations: group-aware val/cal 0.906 (71.7), group-K-fold
0.906 (71.3); R² 0.131 / 0.136. Spearman: 0.36 (LightGBM), 0.32 (ridge) on clusters. Only tree models
and ridge were run: the graph models are too slow on a CPU for this size, and group-K-fold was not run
for the random forest / 5-model LightGBM (also too slow).

Findings:
* **Yield is only weakly predictable from structure here**: best R² 0.13 on unseen clusters (RMSE
  23.8 vs 25.7 for a constant), consistent with the literature on patent yields. This is a property of
  the data (no conditions, text-mined, selection toward successful reactions), not evidence the
  pipeline is wrong. Nonlinear models beat ridge about 2x.
* **Random vs cluster gap is small** (LightGBM 0.160 -> 0.131) - unlike Buchwald (0.94 -> ~0.5) -
  so memorisation of near-duplicates is not what is driving the random score.
* **Coverage is right (~0.90) with plain random-row calibration**; group-aware calibration adds
  nothing (0.895 -> 0.906). Unlike Buchwald (0.4-0.7 under-coverage). Hypothesis (not tested): the
  error here is mostly irreducible noise that looks the same inside and outside a cluster.
* **But the intervals are nearly uninformative**: width 70-73 vs 75 for the constant predictor.
  Coverage alone is not usefulness - compare the width to the `mean baseline` row.
* The 5-model ensemble gains ~0.01 R² and has wider intervals - not clearly worth its cost.
* Caveats: one 50k sample (seed 0), 3 seeds (1 for two cells), ~20 held-out clusters per seed,
  untuned models, no graph models.

## Roadmap (GNN deep-dive)

1. Diagnostics - **built and run on the real Buchwald file** (results above): `product: false`,
   ID-embedding baseline, leave-one-ligand/base-out splits. Not built: longer
   training / lr sweep / Huber loss.
2. Interaction-aware architecture - transformer over component tokens **built and run** (mixed: helps
   with many groups, hurts with few). Not
   built: FiLM conditioning, attention readout, Gasteiger/xtb atom features, heteroscedastic head.
   Component dropout was dropped from the plan: an unseen molecule gives the graph a confident but
   wrong embedding, not a missing one (leave-one-ligand-out: AdBrettPhos R² 0.03 for the graph vs 0.88
   for the ID baseline), so hiding components in training is unlikely to fix it.
3. Uncertainty - group-aware val/cal **built and run** (`val_by_group`: coverage 0.84-0.90 at ~2x width).
   Group-K-fold conformal **built and run** (see "Group-K-fold calibration" above): restores coverage
   but is too wide on scaffolds.
4. Stretch: cached 3D / Boltzmann-averaged embeddings, CGR from substructure-derived atom mapping,
   self-supervised pretraining (Colab).

## Known limitations / honest notes

* No atom mapping is used (the USPTO maps are stripped in cleaning; they would allow a CGR model there), so there is no CGR model yet. The D-MPNN encodes components
  independently and pools by role.
* `.pkl` files are Python pickles and can execute code on load – use only the official TDC file.
* Role rules: on the real file all 3,955 reactions have exactly one catalyst, ligand, base,
  aryl_halide, amine and additive (cleaning report, role signatures). Re-check this after any rule change.
* The TDC file has only 3,955 reactions, 46 distinct input molecules (15 halides, 4 ligands, 3 bases, 22
  additives, 1 amine, 1 catalyst) and 5 distinct products. It is a small, heavily benchmarked HTE set:
  `group_by: scaffold` there is effectively leave-one-aryl-group-out (5 groups), and a graph model has
  little structure to learn beyond ~46 molecule identities.
* The synthetic dataset has invented yields; it validates code, not chemistry.
