"""
Compute kNN-based Jaccard dissimilarity across multiple embeddings.

Author: Tengxiao Gao

Description:
This script computes pairwise cell-wise kNN Jaccard dissimilarity
between multiple embedding spaces (e.g. scVI, Harmony, scPoli, DRVI).
It evaluates how similar the local neighbourhood structure is across embeddings.

Main features:
- Supports multiple embeddings (default: 8 → 28 unordered pairs)
- Computes kNN graphs for different k values
- Calculates per-cell Jaccard dissimilarity
- Outputs:
    1. Summary CSV (mean/median dissimilarity per pair & k)
    2. Per-cell dissimilarity CSV (for selected pair)
    3. Line plot (dissimilarity vs k for all pairs)
    4. UMAP visualization (cell-wise dissimilarity)

Inputs:
- Zarr file containing embeddings stored in `adata.obsm`

Example usage:
python script.py \
    --zarr_path /vol/data/data/output/HRCA/RGC.zarr \
    --k_list 5 10 20 50 100 200 \
    --out_prefix results/RGC

Notes:
- Jaccard dissimilarity = 1 - (intersection / union) of kNN sets
- Larger values indicate more different neighbourhood structures
- Optional subsampling is available for large datasets
- Designed for embedding comparison in single-cell data integration tasks
"""

import os
import argparse
import itertools
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


def compute_knn_indices(X, k, metric="euclidean"):
    """
    Return kNN indices excluding self.
    """
    nn = NearestNeighbors(n_neighbors=k + 1, metric=metric, algorithm="auto")
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
        description="Compute pairwise kNN Jaccard dissimilarity for 8 embeddings (28 unordered pairs), plot multi-line figure, and one UMAP example."
    )
    parser.add_argument(
        "--zarr_path",
        type=str,
        required=True,
        help="Path to zarr, e.g. /vol/data/data/output/HRCA/RGC.zarr",
    )
    parser.add_argument(
        "--embeddings",
        type=str,
        nargs="+",
        default=[
            "drvi_global",
            "drvi_lineage",
            "harmony_global",
            "harmony_lineage",
            "scpoli_global",
            "scpoli_lineage",
            "scvi_global",
            "scvi_lineage",
        ],
        help="Embedding names to compare. Default: 8 embeddings (4 global + 4 lineage).",
    )
    parser.add_argument(
        "--k_list",
        type=int,
        nargs="+",
        default=[5, 10, 20, 50, 100, 200],
        help="List of k values. Example: --k_list 5 10 20 50 100 200",
    )
    parser.add_argument(
        "--metric",
        type=str,
        default="euclidean",
        help="Distance metric for kNN. Default: euclidean",
    )
    parser.add_argument(
        "--umap_pair_a",
        type=str,
        default="scvi_lineage",
        help="Embedding A for the UMAP example pair.",
    )
    parser.add_argument(
        "--umap_pair_b",
        type=str,
        default="harmony_lineage",
        help="Embedding B for the UMAP example pair.",
    )
    parser.add_argument(
        "--umap_name",
        type=str,
        default="scvi_lineage",
        help="Which UMAP to use as 2D carrier for coloring.",
    )
    parser.add_argument(
        "--umap_k",
        type=int,
        default=20,
        help="Which k to show on UMAP.",
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
        default="RGC_all8_28pairs",
        help="Output file prefix.",
    )

    args = parser.parse_args()

    if not os.path.exists(args.zarr_path):
        raise FileNotFoundError(f"Zarr not found: {args.zarr_path}")

    print(f"[INFO] Reading zarr: {args.zarr_path}")
    adata = ad.read_zarr(args.zarr_path)

    emb_names = args.embeddings
    if len(emb_names) < 2:
        raise ValueError("Need at least 2 embeddings.")

    # load embeddings
    emb_key_map = {}
    emb_data_map = {}
    for name in emb_names:
        key = find_obsm_key(adata, "emb", name)
        X = validate_embedding(adata.obsm[key], key)
        emb_key_map[name] = key
        emb_data_map[name] = X
        print(f"[INFO] Embedding loaded: {name} -> {key}, shape={X.shape}")

    # load one UMAP carrier
    umap_key = find_obsm_key(adata, "umap", args.umap_name)
    U = validate_embedding(adata.obsm[umap_key], umap_key)
    if U.shape[1] != 2:
        raise ValueError(f"UMAP must be 2D, got shape={U.shape}")
    print(f"[INFO] UMAP loaded: {args.umap_name} -> {umap_key}, shape={U.shape}")

    n_cells = U.shape[0]
    for name in emb_names:
        if emb_data_map[name].shape[0] != n_cells:
            raise ValueError(f"Cell number mismatch for {name}")

    # optional subsampling
    if args.max_cells > 0 and n_cells > args.max_cells:
        rng = np.random.default_rng(args.seed)
        keep_idx = np.sort(rng.choice(n_cells, size=args.max_cells, replace=False))
        print(f"[INFO] Subsampling {args.max_cells}/{n_cells} cells")
        U = U[keep_idx]
        for name in emb_names:
            emb_data_map[name] = emb_data_map[name][keep_idx]
        cell_indices = keep_idx
    else:
        cell_indices = np.arange(n_cells)

    if max(args.k_list) >= len(cell_indices):
        raise ValueError(
            f"Max k={max(args.k_list)} must be smaller than number of cells={len(cell_indices)}"
        )

    # precompute kNN for each embedding and each k
    knn_cache = {}
    for name in emb_names:
        X = emb_data_map[name]
        for k in args.k_list:
            print(f"[INFO] Computing kNN: embedding={name}, k={k}")
            knn_cache[(name, k)] = compute_knn_indices(X, k=k, metric=args.metric)

    # unordered pairs only: C(n,2)
    pair_list = list(itertools.combinations(emb_names, 2))
    print(f"[INFO] Number of unordered pairs = {len(pair_list)}")

    summary_rows = []
    per_cell_rows = []

    for emb_a, emb_b in pair_list:
        print(f"[INFO] Pair: {emb_a} vs {emb_b}")
        for k in args.k_list:
            knn_a = knn_cache[(emb_a, k)]
            knn_b = knn_cache[(emb_b, k)]
            dissim = jaccard_dissimilarity_from_knn(knn_a, knn_b)

            summary_rows.append(
                {
                    "embedding_a": emb_a,
                    "embedding_b": emb_b,
                    "pair": f"{emb_a} vs {emb_b}",
                    "k": k,
                    "mean_dissimilarity": float(np.mean(dissim)),
                    "median_dissimilarity": float(np.median(dissim)),
                    "n_cells": len(dissim),
                }
            )

            if emb_a == args.umap_pair_a and emb_b == args.umap_pair_b:
                tmp = pd.DataFrame(
                    {
                        "embedding_a": emb_a,
                        "embedding_b": emb_b,
                        "pair": f"{emb_a} vs {emb_b}",
                        "k": k,
                        "cell_index": cell_indices,
                        "dissimilarity": dissim,
                    }
                )
                per_cell_rows.append(tmp)

            # also allow reversed order for UMAP example safety
            elif emb_a == args.umap_pair_b and emb_b == args.umap_pair_a:
                tmp = pd.DataFrame(
                    {
                        "embedding_a": args.umap_pair_a,
                        "embedding_b": args.umap_pair_b,
                        "pair": f"{args.umap_pair_a} vs {args.umap_pair_b}",
                        "k": k,
                        "cell_index": cell_indices,
                        "dissimilarity": dissim,
                    }
                )
                per_cell_rows.append(tmp)

    summary_df = pd.DataFrame(summary_rows)

    if len(per_cell_rows) == 0:
        raise ValueError(
            f"UMAP example pair not found: {args.umap_pair_a} vs {args.umap_pair_b}\n"
            f"Available embeddings: {emb_names}"
        )

    per_cell_df = pd.concat(per_cell_rows, axis=0, ignore_index=True)

    # save CSV
    summary_csv = f"{args.out_prefix}_summary_28pairs.csv"
    per_cell_csv = f"{args.out_prefix}_{args.umap_pair_a}_vs_{args.umap_pair_b}_percell.csv"

    summary_df.to_csv(summary_csv, index=False)
    per_cell_df.to_csv(per_cell_csv, index=False)

    print(f"[INFO] Saved summary CSV : {summary_csv}")
    print(f"[INFO] Saved per-cell CSV: {per_cell_csv}")

    # line plot: 28 lines
    line_png = f"{args.out_prefix}_lineplot_28pairs.png"

    plt.figure(figsize=(16, 9))
    for pair_name, df_pair in summary_df.groupby("pair"):
        df_pair = df_pair.sort_values("k")
        plt.plot(
            df_pair["k"],
            df_pair["mean_dissimilarity"],
            marker="o",
            linewidth=2,
            label=pair_name,
        )

    plt.xlabel("k (number of neighbors)", fontsize=14)
    plt.ylabel("Mean Jaccard dissimilarity", fontsize=14)
    plt.title("Graph dissimilarity vs k", fontsize=18)
    plt.xticks(args.k_list)
    plt.grid(True, alpha=0.3)
    plt.legend(bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0.)
    plt.tight_layout()
    plt.savefig(line_png, dpi=250, bbox_inches="tight")
    plt.close()

    print(f"[INFO] Saved line plot   : {line_png}")

    # one UMAP example
    if args.umap_k not in set(args.k_list):
        raise ValueError(f"--umap_k={args.umap_k} must be one of {args.k_list}")

    df_umap = per_cell_df.loc[per_cell_df["k"] == args.umap_k].copy()
    df_umap = df_umap.sort_values("cell_index")

    umap_png = f"{args.out_prefix}_UMAP_{args.umap_pair_a}_vs_{args.umap_pair_b}_k{args.umap_k}.png"

    plt.figure(figsize=(7, 6))
    sc = plt.scatter(
        U[:, 0],
        U[:, 1],
        c=df_umap["dissimilarity"].to_numpy(),
        s=4,
        alpha=0.9,
        linewidths=0,
    )
    plt.colorbar(sc, label="cell-wise Jaccard dissimilarity")
    plt.title(
        f"Dissimilarity on UMAP\n{args.umap_pair_a} vs {args.umap_pair_b}, k={args.umap_k}\nUMAP: {args.umap_name}"
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