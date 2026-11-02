"""Tabular baselines. They exist so every fancier model has an honest floor to beat."""

from __future__ import annotations

import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.impute import SimpleImputer
from sklearn.linear_model import Ridge
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from .base import MODELS, Model


@MODELS.register("mean")
class MeanBaseline(Model):
    """Predicts the training mean. Any model with R2 <= 0 is not learning anything."""

    def __init__(self, seed: int = 0):
        pass

    def fit(self, X, y, X_val=None, y_val=None):
        self.mu_ = float(np.mean(y))
        return self

    def predict(self, X):
        return np.full(len(X), self.mu_), None


@MODELS.register("ridge")
class RidgeModel(Model):
    def __init__(self, alpha: float = 10.0, seed: int = 0):
        self.alpha = alpha

    def fit(self, X, y, X_val=None, y_val=None):
        self.pipe_ = make_pipeline(SimpleImputer(strategy="mean"), StandardScaler(),
                                   Ridge(alpha=self.alpha))
        self.pipe_.fit(X, y)
        return self

    def predict(self, X):
        return self.pipe_.predict(X), None


@MODELS.register("rf")
class RandomForestModel(Model):
    """sigma = std of the individual trees' predictions."""

    def __init__(self, n_estimators: int = 300, max_features: float = 0.3, min_samples_leaf: int = 1,
                 seed: int = 0, n_jobs: int = -1):
        self.kw = dict(n_estimators=n_estimators, max_features=max_features,
                       min_samples_leaf=min_samples_leaf, random_state=seed, n_jobs=n_jobs)

    def fit(self, X, y, X_val=None, y_val=None):
        self.rf_ = RandomForestRegressor(**self.kw).fit(np.nan_to_num(X, nan=-1.0), y)
        return self

    def predict(self, X):
        X = np.nan_to_num(X, nan=-1.0)
        per_tree = np.stack([t.predict(X) for t in self.rf_.estimators_])
        return per_tree.mean(0), per_tree.std(0)


@MODELS.register("lightgbm")
class LightGBMModel(Model):
    """Gradient boosting with early stopping on the validation set.

    n_models > 1 trains a bagged ensemble (different seeds + row subsampling); sigma is the
    ensemble std. With n_models == 1 there is no sigma (conformal then uses plain residuals).
    """

    def __init__(self, n_models: int = 1, n_estimators: int = 2000, learning_rate: float = 0.05,
                 num_leaves: int = 31, colsample_bytree: float = 0.3, subsample: float = 0.8,
                 min_child_samples: int = 10, reg_lambda: float = 1.0, early_stopping: int = 100,
                 seed: int = 0, n_jobs: int = -1):
        self.n_models, self.early_stopping, self.seed = n_models, early_stopping, seed
        self.kw = dict(n_estimators=n_estimators, learning_rate=learning_rate,
                       num_leaves=num_leaves, colsample_bytree=colsample_bytree,
                       subsample=subsample, subsample_freq=1, min_child_samples=min_child_samples,
                       reg_lambda=reg_lambda, n_jobs=n_jobs, verbosity=-1)

    def fit(self, X, y, X_val=None, y_val=None):
        import warnings

        import lightgbm as lgb

        # lightgbm >=4.6 renamed eval_set; the old name still works on every version we support
        warnings.filterwarnings("ignore", message=".*eval_set.*")

        self.models_ = []
        for k in range(self.n_models):
            m = lgb.LGBMRegressor(random_state=self.seed + k, **self.kw)
            if X_val is not None and self.early_stopping:
                m.fit(X, y, eval_set=[(X_val, y_val)],
                      callbacks=[lgb.early_stopping(self.early_stopping, verbose=False)])
            else:
                m.fit(X, y)
            self.models_.append(m)
        return self

    def predict(self, X):
        P = np.stack([m.predict(X) for m in self.models_])
        return P.mean(0), (P.std(0) if len(self.models_) > 1 else None)
