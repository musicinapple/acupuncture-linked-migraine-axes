#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import logging
import os
import warnings
from dataclasses import dataclass
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np
import pandas as pd
from sklearn.base import clone
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import RFECV
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import RocCurveDisplay, roc_auc_score
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import LinearSVC


LOGGER = logging.getLogger("human_subject_aware_biomarkers")


@dataclass(frozen=True)
class RunConfig:
    top_n_genes: int
    cv_splits: int
    random_seed: int
    test_splits: int
    test_fold: int


def configure_logging(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    # Keep third-party libraries quiet even when verbose is enabled.
    logging.getLogger("matplotlib").setLevel(logging.WARNING)
    logging.getLogger("matplotlib.font_manager").setLevel(logging.WARNING)
    logging.getLogger("matplotlib.backends").setLevel(logging.WARNING)
    logging.getLogger("fontTools").setLevel(logging.WARNING)

def _apply_mpl_style() -> None:
    import matplotlib as mpl
    import seaborn as sns
    sns.set_style("whitegrid")
    mpl.rcParams.update(
        {
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 10,
            "axes.titlesize": 12,
            "axes.labelsize": 10,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "legend.fontsize": 9,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 1.2, # SOTA: Thicker borders
        }
    )


def read_expression_matrix(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", index_col=0)
    if df.empty:
        raise ValueError(f"Empty matrix: {path}")
    return df


def read_meta(path: str, group_col: str, subject_col: str) -> pd.DataFrame:
    meta = pd.read_csv(path, sep="\t")
    required = {"sample_id", group_col, subject_col}
    missing = required - set(meta.columns)
    if missing:
        raise ValueError(f"Metadata missing columns: {sorted(missing)}")
    meta = meta.copy()
    meta["sample_id"] = meta["sample_id"].astype(str)
    meta[group_col] = meta[group_col].astype(str)
    meta[subject_col] = meta[subject_col].astype(str)
    return meta


def load_ranked_genes(path: str) -> List[str]:
    df = pd.read_csv(path)
    if "gene_symbol" not in df.columns:
        raise ValueError("Ranked gene table must contain gene_symbol column.")
    if "adj.P.Val" in df.columns:
        df = df.sort_values(["adj.P.Val", "P.Value"] if "P.Value" in df.columns else ["adj.P.Val"])
    elif "P.Value" in df.columns:
        df = df.sort_values(["P.Value"])
    return [str(g) for g in df["gene_symbol"].tolist() if str(g).strip() and str(g).lower() != "nan"]


def rank_genes_by_abs_mean_diff(expr: pd.DataFrame, y: np.ndarray) -> List[str]:
    if expr.empty:
        raise ValueError("Empty matrix.")
    if expr.shape[1] != len(y):
        raise ValueError("Expression samples do not match y length.")
    unique = set(np.unique(y).tolist())
    if unique - {0, 1}:
        raise ValueError(f"y must be binary {{0,1}} (got {sorted(unique)})")

    case_mask = y == 1
    ctrl_mask = y == 0
    if int(case_mask.sum()) == 0 or int(ctrl_mask.sum()) == 0:
        raise ValueError("Both case and control samples are required for supervised ranking.")

    case_mean = expr.loc[:, case_mask].mean(axis=1)
    ctrl_mean = expr.loc[:, ctrl_mask].mean(axis=1)
    score = (case_mean - ctrl_mean).abs()
    return score.sort_values(ascending=False).index.astype(str).tolist()


def rank_genes_by_variance(expr: pd.DataFrame) -> List[str]:
    if expr.empty:
        raise ValueError("Empty matrix.")
    score = expr.var(axis=1)
    return score.sort_values(ascending=False).index.astype(str).tolist()


def align_samples(expr: pd.DataFrame, meta: pd.DataFrame, group_col: str, subject_col: str) -> Tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    meta = meta.set_index("sample_id")
    missing = [s for s in expr.columns if s not in meta.index]
    if missing:
        raise ValueError(f"Metadata missing {len(missing)} sample(s) from matrix (first 5: {missing[:5]})")

    meta = meta.loc[expr.columns].copy()
    y = meta[group_col].map({"control": 0, "case": 1})
    if y.isna().any():
        bad = meta[group_col][y.isna()].unique().tolist()
        raise ValueError(f"Unexpected {group_col} labels: {bad}. Expected control/case.")

    groups = meta[subject_col].astype(str).to_numpy()
    return expr, y.astype(int).to_numpy(), groups


def verify_subject_label_consistency(y: np.ndarray, groups: np.ndarray) -> None:
    df = pd.DataFrame({"y": y.astype(int), "subject": groups.astype(str)})
    per_subject = df.groupby("subject")["y"].nunique()
    bad = per_subject[per_subject > 1]
    if not bad.empty:
        raise ValueError(f"Group labels vary within subject(s): {bad.index.tolist()[:10]}")


def get_positive_proba(model: Pipeline, X: np.ndarray, positive_label: int = 1) -> np.ndarray:
    proba = model.predict_proba(X)
    classes = list(model.classes_)
    if positive_label not in classes:
        raise ValueError(f"Positive label {positive_label} not in model.classes_: {classes}")
    return proba[:, classes.index(positive_label)]


def select_top_genes(expr: pd.DataFrame, ranked: Sequence[str], top_n: int) -> pd.DataFrame:
    chosen = [g for g in ranked if g in expr.index]
    if not chosen:
        raise ValueError("No ranked genes found in expression matrix.")
    chosen = chosen[: min(int(top_n), len(chosen))]
    out = expr.loc[chosen].copy()
    variances = out.var(axis=1)
    out = out.loc[variances > 0].copy()
    if out.empty:
        raise ValueError("All selected genes had zero variance.")
    return out


def make_holdout_split(
    y: np.ndarray,
    groups: np.ndarray,
    n_splits: int,
    test_fold: int,
    seed: int,
) -> Tuple[np.ndarray, np.ndarray]:
    splitter = StratifiedGroupKFold(n_splits=int(n_splits), shuffle=True, random_state=int(seed))
    splits = list(splitter.split(np.zeros_like(y), y, groups))
    if not (0 <= int(test_fold) < len(splits)):
        raise ValueError(f"--test-fold must be in [0, {len(splits)-1}] (got {test_fold})")
    train_idx, test_idx = splits[int(test_fold)]
    return train_idx, test_idx


def make_cv_splits(y: np.ndarray, groups: np.ndarray, n_splits: int, seed: int) -> List[Tuple[np.ndarray, np.ndarray]]:
    splitter = StratifiedGroupKFold(n_splits=int(n_splits), shuffle=True, random_state=int(seed))
    return list(splitter.split(np.zeros_like(y), y, groups))


def tune_lasso_c(
    X: np.ndarray,
    y: np.ndarray,
    groups: np.ndarray,
    splits: List[Tuple[np.ndarray, np.ndarray]],
    random_seed: int,
    auc_tolerance: float,
) -> Tuple[float, List[Dict[str, float]]]:
    pipeline = Pipeline(
        steps=[
            ("scaler", StandardScaler(with_mean=True, with_std=True)),
            ("clf", LogisticRegression(penalty="l1", solver="liblinear", max_iter=5000, random_state=random_seed)),
        ]
    )

    candidate_cs = np.logspace(-3, 2, 20)
    scores: List[Dict[str, float]] = []
    best_score = -1.0

    for c in candidate_cs:
        model = clone(pipeline)
        model.set_params(clf__C=float(c))
        fold_scores: List[float] = []
        for train_idx, test_idx in splits:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                model.fit(X[train_idx], y[train_idx])
                proba = get_positive_proba(model, X[test_idx])
                fold_scores.append(roc_auc_score(y[test_idx], proba))
        mean_auc = float(np.mean(fold_scores))
        scores.append({"C": float(c), "mean_auc": mean_auc})
        if mean_auc > best_score:
            best_score = mean_auc

    # Prefer sparsity: choose the smallest C whose AUC is within tolerance of the best.
    tol = max(0.0, float(auc_tolerance))
    candidates = [d for d in scores if d["mean_auc"] >= best_score - tol]
    if not candidates:
        best_c = float(max(scores, key=lambda d: d["mean_auc"])["C"])
    else:
        best_c = float(min(candidates, key=lambda d: d["C"])["C"])

    return best_c, scores


def fit_lasso(X: np.ndarray, y: np.ndarray, c_value: float, random_seed: int) -> Pipeline:
    model = Pipeline(
        steps=[
            ("scaler", StandardScaler(with_mean=True, with_std=True)),
            ("clf", LogisticRegression(penalty="l1", solver="liblinear", C=c_value, max_iter=5000, random_state=random_seed)),
        ]
    )
    model.fit(X, y)
    return model


def extract_lasso_selected_genes(model: Pipeline, feature_names: Sequence[str]) -> List[str]:
    clf: LogisticRegression = model.named_steps["clf"]
    coefs = clf.coef_.ravel()
    return [feature_names[i] for i, coef in enumerate(coefs) if abs(float(coef)) > 1e-9]


def run_svm_rfe(
    X: np.ndarray,
    y: np.ndarray,
    splits: List[Tuple[np.ndarray, np.ndarray]],
    n_jobs: int,
) -> Tuple[List[int], List[Dict[str, float]]]:
    scaler = StandardScaler(with_mean=True, with_std=True)
    Xs = scaler.fit_transform(X)
    estimator = LinearSVC(dual=False, max_iter=30000)
    rfecv = RFECV(
        estimator=estimator,
        step=1,
        cv=splits,
        scoring="roc_auc",
        min_features_to_select=5,
        n_jobs=int(n_jobs),
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        rfecv.fit(Xs, y)
    selected_idx = [i for i, ok in enumerate(rfecv.support_.tolist()) if ok]

    curve: List[Dict[str, float]] = []
    if hasattr(rfecv, "cv_results_") and "mean_test_score" in rfecv.cv_results_:
        scores = rfecv.cv_results_["mean_test_score"]
        n_features = rfecv.cv_results_["n_features"]
        for n, s in zip(n_features, scores):
            curve.append({"n_features": int(n), "mean_auc": float(s)})
    return selected_idx, curve


def fit_random_forest(X: np.ndarray, y: np.ndarray, random_seed: int) -> RandomForestClassifier:
    model = RandomForestClassifier(
        n_estimators=800,
        random_state=random_seed,
        n_jobs=-1,
        class_weight="balanced_subsample",
    )
    model.fit(X, y)
    return model


def save_list(path: str, items: Iterable[str]) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["gene_symbol"])
        for item in items:
            w.writerow([item])


def plot_lasso_cv(scores: List[Dict[str, float]], out_pdf: str) -> None:
    import matplotlib.pyplot as plt

    os.makedirs(os.path.dirname(out_pdf), exist_ok=True)
    _apply_mpl_style()
    xs = [d["C"] for d in scores]
    ys = [d["mean_auc"] for d in scores]
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.plot(xs, ys, marker="o", linewidth=2.0, markersize=6, color="#377EB8", alpha=0.8)
    ax.set_xscale("log")
    ax.set_xlabel("C (inverse regularization strength)", fontweight='bold')
    ax.set_ylabel("Mean AUC (subject CV)", fontweight='bold')
    ax.grid(axis="both", linestyle="--", alpha=0.3)
    ax.set_axisbelow(True)
    fig.tight_layout()
    fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_svm_rfe_curve(curve: List[Dict[str, float]], out_pdf: str) -> None:
    import matplotlib.pyplot as plt

    os.makedirs(os.path.dirname(out_pdf), exist_ok=True)
    _apply_mpl_style()
    if not curve:
        return
    xs = [d["n_features"] for d in curve]
    ys = [d["mean_auc"] for d in curve]
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    ax.plot(xs, ys, marker="o", linewidth=1.6, markersize=4.0, color="#4C72B0")
    ax.set_xlabel("Number of selected features")
    ax.set_ylabel("Mean AUC (subject CV)")
    ax.set_title("Human SVM-RFE (subject-aware CV)")
    ax.grid(axis="both", linestyle="--", alpha=0.18)
    ax.set_axisbelow(True)
    fig.tight_layout()
    fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_rf_importance(feature_names: Sequence[str], importances: np.ndarray, out_pdf: str, top_k: int = 30) -> None:
    import matplotlib.pyplot as plt
    import seaborn as sns

    os.makedirs(os.path.dirname(out_pdf), exist_ok=True)
    _apply_mpl_style()
    order = np.argsort(importances)[::-1][: int(top_k)]
    names = [feature_names[i] for i in order][::-1]
    vals = importances[order][::-1]
    
    fig, ax = plt.subplots(figsize=(7, 8))
    sns.barplot(x=vals, y=names, palette="viridis", ax=ax, edgecolor='white', linewidth=0.5)
    ax.set_xlabel("Gini importance", fontweight='bold')
    ax.grid(axis="x", linestyle="--", alpha=0.3)
    ax.set_axisbelow(True)
    fig.tight_layout()
    fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_venn(a: set[str], b: set[str], c: set[str], out_pdf: str) -> None:
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle

    os.makedirs(os.path.dirname(out_pdf), exist_ok=True)
    _apply_mpl_style()
    fig, ax = plt.subplots(figsize=(7.0, 5.0))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    circles = [
        (0.43, 0.55, "#4C72B0", "LASSO"),
        (0.60, 0.55, "#55A868", "SVM-RFE"),
        (0.515, 0.40, "#C44E52", "RF"),
    ]
    for x, y, color, label in circles:
        ax.add_patch(Circle((x, y), 0.23, alpha=0.28, color=color))
        ax.text(x, y + 0.26, label, ha="center", va="center", fontsize=11, fontweight="bold")

    only_a = len(a - b - c)
    only_b = len(b - a - c)
    only_c = len(c - a - b)
    ab = len((a & b) - c)
    ac = len((a & c) - b)
    bc = len((b & c) - a)
    abc = len(a & b & c)

    ax.text(0.33, 0.58, str(only_a), ha="center", va="center", fontsize=12)
    ax.text(0.70, 0.58, str(only_b), ha="center", va="center", fontsize=12)
    ax.text(0.515, 0.20, str(only_c), ha="center", va="center", fontsize=12)
    ax.text(0.515, 0.62, str(ab), ha="center", va="center", fontsize=12)
    ax.text(0.42, 0.37, str(ac), ha="center", va="center", fontsize=12)
    ax.text(0.61, 0.37, str(bc), ha="center", va="center", fontsize=12)
    ax.text(0.515, 0.46, str(abc), ha="center", va="center", fontsize=12, fontweight="bold")

    ax.set_title("Human Feature Selection Overlap", fontsize=12)
    fig.tight_layout()
    fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_roc(y_true: np.ndarray, y_score: np.ndarray, out_pdf: str, title: str, n_boot: int = 200, seed: int = 12345) -> float:
    import matplotlib.pyplot as plt
    from sklearn.metrics import roc_curve

    os.makedirs(os.path.dirname(out_pdf), exist_ok=True)
    _apply_mpl_style()
    auc = float(roc_auc_score(y_true, y_score))
    ci_text = ""
    fpr, tpr, _ = roc_curve(y_true, y_score)
    grid = np.linspace(0, 1, 101)
    tpr_grid = np.interp(grid, fpr, tpr)
    tpr_ci = None
    if int(n_boot) > 0:
        rng = np.random.default_rng(int(seed))
        idx_case = np.where(np.asarray(y_true) == 1)[0]
        idx_ctrl = np.where(np.asarray(y_true) == 0)[0]
        boot = []
        tpr_boot = []
        for _ in range(int(n_boot)):
            s_case = rng.choice(idx_case, size=len(idx_case), replace=True)
            s_ctrl = rng.choice(idx_ctrl, size=len(idx_ctrl), replace=True)
            s = np.concatenate([s_case, s_ctrl])
            if len(np.unique(np.asarray(y_true)[s])) < 2:
                continue
            y_b = np.asarray(y_true)[s]
            s_b = np.asarray(y_score)[s]
            boot.append(float(roc_auc_score(y_b, s_b)))
            fpr_b, tpr_b, _ = roc_curve(y_b, s_b)
            tpr_boot.append(np.interp(grid, fpr_b, tpr_b))
        if boot:
            lo, hi = np.percentile(np.asarray(boot), [2.5, 97.5]).tolist()
            ci_text = f" (95% CI {lo:.3f}-{hi:.3f})"
        if tpr_boot:
            arr = np.vstack(tpr_boot)
            tpr_ci = (
                np.percentile(arr, 2.5, axis=0),
                np.percentile(arr, 97.5, axis=0),
            )
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.plot(grid, tpr_grid, color="#377EB8", linewidth=2.5, label=f"ROC (AUC = {auc:.3f})")
    if tpr_ci is not None:
        ax.fill_between(grid, tpr_ci[0], tpr_ci[1], color="#377EB8", alpha=0.15, label="95% CI")
    ax.plot([0, 1], [0, 1], linestyle="--", color="#777777", linewidth=1.5, alpha=0.8)
    
    ax.set_xlabel("False Positive Rate", fontweight='bold')
    ax.set_ylabel("True Positive Rate", fontweight='bold')
    ax.grid(axis="both", linestyle="--", alpha=0.3)
    ax.set_axisbelow(True)
    ax.legend(frameon=True, framealpha=0.9, loc="lower right")
    
    # Text annotation for AUC summary
    ax.text(0.95, 0.05, f"AUC: {auc:.3f}{ci_text}", transform=ax.transAxes, 
            ha='right', va='bottom', fontsize=10, fontweight='bold',
            bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))

    fig.tight_layout()
    fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return auc


def main() -> int:
    parser = argparse.ArgumentParser(description="Task 3.3: Human biomarker selection with subject-aware (no leakage) splitting.")
    parser.add_argument("--matrix", default="data/processed/human_expression_matrix.tsv", help="Human expression matrix (genes x samples).")
    parser.add_argument("--meta", default="data/processed/human_sample_metadata.tsv", help="Human sample metadata TSV.")
    parser.add_argument(
        "--prefilter",
        choices=("train_effect", "variance", "ranked_file"),
        default="train_effect",
        help=(
            "How to prefilter genes before feature selection. "
            "train_effect = rank genes by |mean(case)-mean(control)| computed on TRAIN subjects only (recommended). "
            "variance = rank genes by variance on TRAIN subjects only (unsupervised). "
            "ranked_file = use --ranked-genes (can leak if computed on full dataset)."
        ),
    )
    parser.add_argument(
        "--ranked-genes",
        default="results/tables/T2.6_human_deg_emtab13397.csv",
        help="Ranked gene table used only when --prefilter ranked_file is selected.",
    )
    parser.add_argument("--group-col", default="group", help="Metadata column for class labels (control/case).")
    parser.add_argument("--subject-col", default="subject_id", help="Metadata column used as subject grouping.")
    parser.add_argument("--top-n", type=int, default=500, help="Top N ranked genes to consider as features.")
    parser.add_argument("--cv-splits", type=int, default=5, help="Subject-aware CV splits for hyperparameter tuning.")
    parser.add_argument("--test-splits", type=int, default=5, help="Subject-aware splits used to create a held-out test fold.")
    parser.add_argument("--test-fold", type=int, default=0, help="Which fold index to hold out as test (0-based).")
    parser.add_argument(
        "--lasso-auc-tolerance",
        type=float,
        default=0.01,
        help="Choose the smallest LASSO C within (best_auc - tolerance) to favor sparse signatures.",
    )
    parser.add_argument(
        "--svm-rfe-jobs",
        type=int,
        default=1,
        help="Parallel jobs for SVM-RFE (use 1 to avoid noisy warnings from worker processes).",
    )
    parser.add_argument("--seed", type=int, default=20251222, help="Random seed for splitting/models.")
    parser.add_argument("--out-dir", default="results/tables", help="Output directory for tables.")
    parser.add_argument("--fig-dir", default="figures/raw_plots", help="Output directory for figures.")
    parser.add_argument("--tag", default="T3.3", help="Prefix tag for outputs (default: T3.3).")
    parser.add_argument("--verbose", action="store_true", help="Verbose logging.")
    args = parser.parse_args()

    configure_logging(args.verbose)
    cfg = RunConfig(
        top_n_genes=int(args.top_n),
        cv_splits=int(args.cv_splits),
        random_seed=int(args.seed),
        test_splits=int(args.test_splits),
        test_fold=int(args.test_fold),
    )
    run_at = dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

    warnings.filterwarnings(
        "ignore",
        message=r".*encountered in matmul.*",
        category=RuntimeWarning,
        module=r"sklearn\\.utils\\.extmath",
    )

    expr = read_expression_matrix(args.matrix)
    meta = read_meta(args.meta, group_col=args.group_col, subject_col=args.subject_col)

    expr, y, groups = align_samples(expr, meta, group_col=args.group_col, subject_col=args.subject_col)
    verify_subject_label_consistency(y, groups)

    train_idx, test_idx = make_holdout_split(y, groups, cfg.test_splits, cfg.test_fold, cfg.random_seed)
    y_train, y_test = y[train_idx], y[test_idx]
    groups_train, groups_test = groups[train_idx], groups[test_idx]

    expr_train = expr.iloc[:, train_idx]
    if args.prefilter == "ranked_file":
        ranked = load_ranked_genes(args.ranked_genes)
    elif args.prefilter == "variance":
        ranked = rank_genes_by_variance(expr_train)
    else:
        ranked = rank_genes_by_abs_mean_diff(expr_train, y_train)

    expr_small = select_top_genes(expr, ranked, cfg.top_n_genes)
    feature_names = expr_small.index.astype(str).tolist()
    X_all = expr_small.T.to_numpy(dtype=float)

    X_train, X_test = X_all[train_idx], X_all[test_idx]

    def _subject_counts(y_arr: np.ndarray, grp_arr: np.ndarray) -> Dict[str, int]:
        df = pd.DataFrame({"y": y_arr.astype(int), "subject": grp_arr.astype(str)}).drop_duplicates("subject")
        return {"n_subjects": int(len(df)), "n_case_subjects": int((df["y"] == 1).sum()), "n_control_subjects": int((df["y"] == 0).sum())}

    train_subject_counts = _subject_counts(y_train, groups_train)
    test_subject_counts = _subject_counts(y_test, groups_test)

    LOGGER.info(
        "Holdout split: train=%s samples (%s subjects), test=%s samples (%s subjects)",
        len(train_idx),
        train_subject_counts["n_subjects"],
        len(test_idx),
        test_subject_counts["n_subjects"],
    )
    LOGGER.info("Train subjects: %s", train_subject_counts)
    LOGGER.info("Test subjects:  %s", test_subject_counts)

    cv_splits = make_cv_splits(y_train, groups_train, cfg.cv_splits, cfg.random_seed)
    LOGGER.info("Inner CV splits: %s", len(cv_splits))

    # LASSO (train only)
    best_c, lasso_scores = tune_lasso_c(
        X_train,
        y_train,
        groups_train,
        cv_splits,
        cfg.random_seed,
        auc_tolerance=float(args.lasso_auc_tolerance),
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        lasso_model = fit_lasso(X_train, y_train, best_c, cfg.random_seed)
    lasso_selected = extract_lasso_selected_genes(lasso_model, feature_names)
    if not lasso_selected:
        candidates = sorted([d["C"] for d in lasso_scores], reverse=True)
        for c in candidates:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                trial = fit_lasso(X_train, y_train, float(c), cfg.random_seed)
            selected = extract_lasso_selected_genes(trial, feature_names)
            if selected:
                best_c = float(c)
                lasso_model = trial
                lasso_selected = selected
                break
    LOGGER.info("LASSO selected %s / %s genes (C=%s)", len(lasso_selected), len(feature_names), best_c)

    # SVM-RFE (train only)
    svm_selected_indices, svm_curve = run_svm_rfe(X_train, y_train, cv_splits, n_jobs=int(args.svm_rfe_jobs))
    svm_selected = [feature_names[i] for i in svm_selected_indices]
    LOGGER.info("SVM-RFE selected %s / %s genes", len(svm_selected), len(feature_names))

    # RF (train only)
    rf_model = fit_random_forest(X_train, y_train, cfg.random_seed)
    rf_importances = rf_model.feature_importances_
    rf_top = np.argsort(rf_importances)[::-1][:50]
    rf_selected = [feature_names[i] for i in rf_top]
    LOGGER.info("RF selected top %s genes for intersection", len(rf_selected))

    inter = sorted(set(lasso_selected) & set(svm_selected) & set(rf_selected))
    LOGGER.info("Intersection genes: %s", len(inter))

    tag = str(args.tag).strip()
    out_lasso = os.path.join(args.out_dir, f"{tag}_human_lasso_selected.csv")
    out_svm = os.path.join(args.out_dir, f"{tag}_human_svm_rfe_selected.csv")
    out_rf = os.path.join(args.out_dir, f"{tag}_human_rf_top50.csv")
    out_intersection = os.path.join(args.out_dir, f"{tag}_human_intersection_genes.csv")
    out_config = os.path.join(args.out_dir, f"{tag}_human_run_config.json")
    out_perf = os.path.join(args.out_dir, f"{tag}_human_performance.json")

    save_list(out_lasso, lasso_selected)
    save_list(out_svm, svm_selected)
    save_list(out_rf, rf_selected)
    save_list(out_intersection, inter)

    os.makedirs(os.path.dirname(out_config), exist_ok=True)
    with open(out_config, "w", encoding="utf-8") as f:
        json.dump(
            {
                "run_at_utc": run_at,
                "matrix": args.matrix,
                "meta": args.meta,
                "prefilter": args.prefilter,
                "ranked_genes": args.ranked_genes,
                "group_col": args.group_col,
                "subject_col": args.subject_col,
                "top_n": cfg.top_n_genes,
                "seed": cfg.random_seed,
                "test_splits": cfg.test_splits,
                "test_fold": cfg.test_fold,
                "cv_splits": cfg.cv_splits,
                "best_lasso_C": best_c,
                "n_features_input": len(feature_names),
                "n_lasso_selected": len(lasso_selected),
                "n_svm_selected": len(svm_selected),
                "n_rf_selected": len(rf_selected),
                "n_intersection": len(inter),
                "train_samples": int(len(train_idx)),
                "test_samples": int(len(test_idx)),
                "train_subjects": int(len(set(groups_train.tolist()))),
                "test_subjects": int(len(set(groups_test.tolist()))),
            },
            f,
            ensure_ascii=False,
            indent=2,
        )
        f.write("\n")

    fig_lasso = os.path.join(args.fig_dir, "Fig4A.png")
    fig_svm = os.path.join(args.fig_dir, "Fig4_supp_svm.png")
    fig_rf = os.path.join(args.fig_dir, "Fig4B.png")
    fig_venn = os.path.join(args.fig_dir, "Fig4_supp_venn.png")
    fig_roc = os.path.join(args.fig_dir, "Fig4C.png")

    plot_lasso_cv(lasso_scores, fig_lasso)
    plot_svm_rfe_curve(svm_curve, fig_svm)
    plot_rf_importance(feature_names, rf_importances, fig_rf, top_k=30)
    plot_venn(set(lasso_selected), set(svm_selected), set(rf_selected), fig_venn)

    # Evaluate on held-out subjects using LASSO model trained on train only.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        test_proba = get_positive_proba(lasso_model, X_test)
    auc_test = plot_roc(y_test, test_proba, fig_roc, title="E-MTAB-13397 subject-held-out ROC (LASSO)")

    with open(out_perf, "w", encoding="utf-8") as f:
        json.dump(
            {
                "run_at_utc": run_at,
                "model": "lasso_logistic_regression",
                "best_lasso_C": best_c,
                "holdout_auc": auc_test,
                "train_subject_counts": train_subject_counts,
                "test_subject_counts": test_subject_counts,
                "notes": [
                    "Subject-aware split prevents leakage across repeated measures.",
                    "Gene prefiltering is computed on TRAIN subjects only unless --prefilter ranked_file is used.",
                ],
            },
            f,
            ensure_ascii=False,
            indent=2,
        )
        f.write("\n")

    LOGGER.info("Held-out subject ROC AUC: %.4f", auc_test)
    LOGGER.info("Saved tables: %s %s %s %s %s %s", out_lasso, out_svm, out_rf, out_intersection, out_config, out_perf)
    LOGGER.info("Saved figures: %s %s %s %s %s", fig_lasso, fig_svm, fig_rf, fig_venn, fig_roc)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
