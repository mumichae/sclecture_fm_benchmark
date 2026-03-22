"""
Compute kNN-based Jaccard dissimilarity for a single pair of embeddings
across a wide range of k values (k sweep), with curve and UMAP visualization.

Author: Tengxiao Gao

Description:
This script compares two embedding spaces by measuring how their local
neighbourhood structures differ. It computes per-cell Jaccard dissimilarity
between k-nearest neighbour (kNN) graphs over a range of k values.

Main features:
- Compare one pair of embeddings (e.g. scVI vs Harmony)
- Sweep k over a wide range (e.g. 50 → 2000)
- Efficient computation by reusing kNN up to max_k
- Outputs:
    1. Summary CSV (mean/median dissimilarity vs k)
    2. Per-cell CSV for a selected k
    3. Line plot (dissimilarity vs k)
    4. UMAP plot (cell-wise dissimilarity)

Inputs:
- Zarr file with embeddings stored in `adata.obsm`

Example usage:
python script.py \
    --zarr_path /vol/data/data/output/HRCA/RGC.zarr \
    --emb_a scvi_lineage \
    --emb_b harmony_lineage \
    --k_start 50 \
    --k_end 2000 \
    --k_step 50 \
    --umap_k 500 \
    --out_prefix results/scvi_vs_harmony

Notes:
- Jaccard dissimilarity = 1 - (intersection / union) of kNN sets
- Larger values indicate greater structural differences
- Results can be affected by cell density and choice of k
- Optional subsampling is available for large datasets

"""
import os
import argparse
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
from sklearn.neighbors import NearestNeighbors


def find_obsm_key(adata, kind, name):
    """
    kind: 'emb' or 'umap'
    name: e.g. 'scvi_lineage'
    """
    candidates = [
        f"X_{kind}--metrics:{name}",
        f"X_{kind}--metrics/{name}",
        f"X_{kind}:{name}",
        f"X_{kind}/{name}",
    ]
    for k in candidates:
        if k in adata.obsm_keys():
            return k

    hits = [k for k in adata.obsm_keys() if (f"X_{kind}" in str(k) and name in str(k))]
    if len(hits) == 1:
        return hits[0]
    elif len(hits) > 1:
        raise KeyError(
            f"Multiple possible {kind} keys found for '{name}': {hits}\n"
            f"Please inspect adata.obsm_keys() and choose one."
        )
    else:
        raise KeyError(
            f"Could not find {kind} key for '{name}'.\n"
            f"Available obsm keys:\n{list(adata.obsm_keys())}"
        )


def validate_embedding(X, key):
    X = np.asarray(X)
    if X.ndim != 2:
        raise ValueError(f"Embedding {key} must be 2D, got shape={X.shape}")
    if not np.isfinite(X).all():
        raise ValueError(f"Embedding {key} contains NaN or inf")
    return X


def compute_knn_indices_maxk(X, max_k, metric="euclidean"):
    """
    Compute neighbors once up to max_k, excluding self.
    Then smaller k can be obtained by slicing [:, :k].
    """
    nn = NearestNeighbors(n_neighbors=max_k + 1, metric=metric, algorithm="auto")
    nn.fit(X)
    indices = nn.kneighbors(X, return_distance=False)
    return indices[:, 1:]


def jaccard_dissimilarity_from_knn(knn_a, knn_b):
    """
    Per-cell Jaccard dissimilarity:
        1 - |A ∩ B| / |A ∪ B|
    """
    n = knn_a.shape[0]
    out = np.empty(n, dtype=np.float64)

    for i in range(n):
        sa = set(knn_a[i])
        sb = set(knn_b[i])
        inter = len(sa & sb)
        union = len(sa | sb)
        out[i] = 1.0 - (inter / union if union > 0 else 0.0)

    return out


def main():
    parser = argparse.ArgumentParser(
        description="Compute single-pair kNN Jaccard dissimilarity over a wide k range and plot curve + one UMAP example."
    )
    parser.add_argument(
        "--zarr_path",
        type=str,
        required=True,
        help="Path to zarr, e.g. /vol/data/data/output/HRCA/RGC.zarr",
    )
    parser.add_argument(
        "--emb_a",
        type=str,
        default="scvi_lineage",
        help="Embedding A. Default: scvi_lineage",
    )
    parser.add_argument(
        "--emb_b",
        type=str,
        default="harmony_lineage",
        help="Embedding B. Default: harmony_lineage",
    )
    parser.add_argument(
        "--umap_name",
        type=str,
        default="scvi_lineage",
        help="Which UMAP to use as 2D carrier. Default: scvi_lineage",
    )
    parser.add_argument(
        "--k_start",
        type=int,
        default=50,
        help="Start of k range. Default: 50",
    )
    parser.add_argument(
        "--k_end",
        type=int,
        default=2000,
        help="End of k range (inclusive). Default: 2000",
    )
    parser.add_argument(
        "--k_step",
        type=int,
        default=50,
        help="Step size of k. Default: 50",
    )
    parser.add_argument(
        "--umap_k",
        type=int,
        default=500,
        help="Which k to display on UMAP. Default: 500",
    )
    parser.add_argument(
        "--metric",
        type=str,
        default="euclidean",
        help="Distance metric for kNN. Default: euclidean",
    )
    parser.add_argument(
        "--max_cells",
        type=int,
        default=0,
        help="Optional subsampling. 0 means use all cells.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=0,
        help="Random seed for subsampling.",
    )
    parser.add_argument(
        "--out_prefix",
        type=str,
        default="RGC_scvi_lineage_vs_harmony_lineage_k50_to_2000",
        help="Output file prefix.",
    )

    args = parser.parse_args()

    if not os.path.exists(args.zarr_path):
        raise FileNotFoundError(f"Zarr not found: {args.zarr_path}")

    k_list = list(range(args.k_start, args.k_end + 1, args.k_step))
    if len(k_list) == 0:
        raise ValueError("k_list is empty. Check k_start, k_end, k_step.")

    if args.umap_k not in k_list:
        raise ValueError(f"--umap_k={args.umap_k} must be in generated k list: {k_list}")

    print(f"[INFO] k_list: {k_list[:5]} ... {k_list[-5:]}  (n={len(k_list)})")

    print(f"[INFO] Reading zarr: {args.zarr_path}")
    adata = ad.read_zarr(args.zarr_path)

    emb_key_a = find_obsm_key(adata, "emb", args.emb_a)
    emb_key_b = find_obsm_key(adata, "emb", args.emb_b)
    umap_key = find_obsm_key(adata, "umap", args.umap_name)

    print(f"[INFO] emb_a key : {emb_key_a}")
    print(f"[INFO] emb_b key : {emb_key_b}")
    print(f"[INFO] umap key  : {umap_key}")

    X_a = validate_embedding(adata.obsm[emb_key_a], emb_key_a)
    X_b = validate_embedding(adata.obsm[emb_key_b], emb_key_b)
    U = validate_embedding(adata.obsm[umap_key], umap_key)

    if U.shape[1] != 2:
        raise ValueError(f"UMAP must be 2D, got shape={U.shape}")

    n_cells = X_a.shape[0]
    if X_b.shape[0] != n_cells or U.shape[0] != n_cells:
        raise ValueError("Cell number mismatch across embeddings/UMAP")

    if args.max_cells > 0 and n_cells > args.max_cells:
        rng = np.random.default_rng(args.seed)
        keep_idx = np.sort(rng.choice(n_cells, size=args.max_cells, replace=False))
        print(f"[INFO] Subsampling {args.max_cells}/{n_cells} cells")
        X_a = X_a[keep_idx]
        X_b = X_b[keep_idx]
        U = U[keep_idx]
        cell_indices = keep_idx
    else:
        cell_indices = np.arange(n_cells)

    n_used = len(cell_indices)
    max_k = max(k_list)

    if max_k >= n_used:
        raise ValueError(
            f"Max k={max_k} must be smaller than number of used cells={n_used}"
        )

    print(f"[INFO] Number of cells used: {n_used}")
    print(f"[INFO] Computing neighbors once up to max_k={max_k} for both embeddings...")

    knn_a_all = compute_knn_indices_maxk(X_a, max_k=max_k, metric=args.metric)
    knn_b_all = compute_knn_indices_maxk(X_b, max_k=max_k, metric=args.metric)

    summary_rows = []
    umap_percell_df = None

    for k in k_list:
        print(f"[INFO] Computing dissimilarity for k={k}")
        knn_a = knn_a_all[:, :k]
        knn_b = knn_b_all[:, :k]

        dissim = jaccard_dissimilarity_from_knn(knn_a, knn_b)

        summary_rows.append(
            {
                "embedding_a": args.emb_a,
                "embedding_b": args.emb_b,
                "pair": f"{args.emb_a} vs {args.emb_b}",
                "k": k,
                "mean_dissimilarity": float(np.mean(dissim)),
                "median_dissimilarity": float(np.median(dissim)),
                "n_cells": len(dissim),
            }
        )

        if k == args.umap_k:
            umap_percell_df = pd.DataFrame(
                {
                    "embedding_a": args.emb_a,
                    "embedding_b": args.emb_b,
                    "pair": f"{args.emb_a} vs {args.emb_b}",
                    "k": k,
                    "cell_index": cell_indices,
                    "dissimilarity": dissim,
                }
            )

    if umap_percell_df is None:
        raise RuntimeError("Failed to create per-cell dataframe for UMAP.")

    summary_df = pd.DataFrame(summary_rows)

    # save CSVs
    summary_csv = f"{args.out_prefix}_summary.csv"
    percell_csv = f"{args.out_prefix}_percell_k{args.umap_k}.csv"

    summary_df.to_csv(summary_csv, index=False)
    umap_percell_df.to_csv(percell_csv, index=False)

    print(f"[INFO] Saved summary CSV : {summary_csv}")
    print(f"[INFO] Saved per-cell CSV: {percell_csv}")

    # line plot
    line_png = f"{args.out_prefix}_lineplot.png"
    plt.figure(figsize=(9, 6))
    plt.plot(
        summary_df["k"],
        summary_df["mean_dissimilarity"],
        marker="o",
        linewidth=2.2,
    )
    plt.xlabel("k (number of neighbors)")
    plt.ylabel("Mean Jaccard dissimilarity")
    plt.title(f"Graph dissimilarity vs k\n{args.emb_a} vs {args.emb_b}")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(line_png, dpi=250)
    plt.close()

    print(f"[INFO] Saved line plot   : {line_png}")

    # UMAP plot
    umap_png = f"{args.out_prefix}_UMAP_k{args.umap_k}.png"
    plt.figure(figsize=(7, 6))
    sc = plt.scatter(
        U[:, 0],
        U[:, 1],
        c=umap_percell_df["dissimilarity"].to_numpy(),
        s=4,
        alpha=0.9,
        linewidths=0,
    )
    plt.colorbar(sc, label="cell-wise Jaccard dissimilarity")
    plt.title(
        f"Dissimilarity on UMAP\n{args.emb_a} vs {args.emb_b}, k={args.umap_k}\nUMAP: {args.umap_name}"
    )
    plt.xlabel("UMAP1")
    plt.ylabel("UMAP2")
    plt.tight_layout()
    plt.savefig(umap_png, dpi=250)
    plt.close()

    print(f"[INFO] Saved UMAP plot   : {umap_png}")
    print("[DONE]")


if __name__ == "__main__":
    main()