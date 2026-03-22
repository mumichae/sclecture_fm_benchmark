"""
Compute kNN-based Jaccard dissimilarity between two embeddings
(scVI vs Harmony) across a large k range, with curve and UMAP visualization.

Author: Tengxiao Gao

Description:
This script compares two embedding spaces by measuring how their
k-nearest neighbour (kNN) structures differ. It computes per-cell
Jaccard dissimilarity for a range of k values (k = 1000–8000).

Outputs:
- Summary CSV (mean/median dissimilarity vs k)
- Per-cell CSV (for selected k)
- Line plot (dissimilarity vs k)
- UMAP plot (cell-wise dissimilarity)

Input:
- Zarr file containing embeddings in `adata.obsm`

Usage example:
python script.py \
    --zarr_path /vol/data/data/output/HRCA/RGC.zarr \
    --umap_k 5000

Notes:
- Jaccard dissimilarity = 1 - (intersection / union) of kNN sets
- Larger values indicate greater differences between embeddings
- This script requires substantial memory for large k (e.g. k > 1000)
- Due to VM memory limitations, full runs (k up to 8000) may not complete
"""
import os
import argparse
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
from sklearn.neighbors import NearestNeighbors


def find_obsm_key(adata, kind, name):

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

    raise KeyError(f"Could not find embedding for {name}")


def compute_knn_indices_maxk(X, max_k):

    nn = NearestNeighbors(
        n_neighbors=max_k + 1,
        metric="euclidean",
        algorithm="auto"
    )

    nn.fit(X)

    indices = nn.kneighbors(X, return_distance=False)

    return indices[:, 1:]


def jaccard_dissimilarity_from_knn(knn_a, knn_b):

    n = knn_a.shape[0]

    out = np.empty(n)

    for i in range(n):

        sa = set(knn_a[i])
        sb = set(knn_b[i])

        inter = len(sa & sb)
        union = len(sa | sb)

        out[i] = 1 - (inter / union if union > 0 else 0)

    return out


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--zarr_path",
        type=str,
        required=True
    )

    parser.add_argument(
        "--umap_k",
        type=int,
        default=5000
    )

    args = parser.parse_args()

    # k list
    k_list = list(range(1000, 8001, 500))

    print("[INFO] k_list:", k_list)

    print("[INFO] Reading zarr")

    adata = ad.read_zarr(args.zarr_path)

    emb_a = find_obsm_key(adata, "emb", "scvi_lineage")
    emb_b = find_obsm_key(adata, "emb", "harmony_lineage")

    umap_key = find_obsm_key(adata, "umap", "scvi_lineage")

    X_a = np.asarray(adata.obsm[emb_a])
    X_b = np.asarray(adata.obsm[emb_b])

    U = np.asarray(adata.obsm[umap_key])

    n_cells = X_a.shape[0]

    max_k = max(k_list)

    print("[INFO] cells:", n_cells)
    print("[INFO] max_k:", max_k)

    print("[INFO] computing knn")

    knn_a = compute_knn_indices_maxk(X_a, max_k)
    knn_b = compute_knn_indices_maxk(X_b, max_k)

    summary_rows = []

    umap_df = None

    for k in k_list:

        print("[INFO] k =", k)

        neigh_a = knn_a[:, :k]
        neigh_b = knn_b[:, :k]

        dissim = jaccard_dissimilarity_from_knn(neigh_a, neigh_b)

        summary_rows.append({
            "k": k,
            "mean_dissimilarity": float(np.mean(dissim)),
            "median_dissimilarity": float(np.median(dissim))
        })

        if k == args.umap_k:

            umap_df = pd.DataFrame({
                "cell_index": np.arange(n_cells),
                "dissimilarity": dissim
            })

    summary_df = pd.DataFrame(summary_rows)

    prefix = "RGC_scvi_vs_harmony_k1000_to_8000"

    summary_csv = f"{prefix}_summary.csv"
    umap_csv = f"{prefix}_percell_k{args.umap_k}.csv"

    summary_df.to_csv(summary_csv, index=False)
    umap_df.to_csv(umap_csv, index=False)

    print("[INFO] saved csv")

    # line plot

    plt.figure(figsize=(9,6))

    plt.plot(
        summary_df["k"],
        summary_df["mean_dissimilarity"],
        marker="o"
    )

    plt.xlabel("k")
    plt.ylabel("Mean Jaccard dissimilarity")

    plt.title("Graph dissimilarity vs k\nscvi_lineage vs harmony_lineage")

    plt.grid(True)

    line_png = f"{prefix}_lineplot.png"

    plt.savefig(line_png, dpi=300)

    plt.close()

    print("[INFO] saved line plot")

    # umap

    if umap_df is not None:

        plt.figure(figsize=(7,6))

        sc = plt.scatter(
            U[:,0],
            U[:,1],
            c=umap_df["dissimilarity"],
            s=4
        )

        plt.colorbar(sc)

        plt.title(f"dissimilarity on UMAP k={args.umap_k}")

        umap_png = f"{prefix}_UMAP_k{args.umap_k}.png"

        plt.savefig(umap_png, dpi=300)

        plt.close()

        print("[INFO] saved umap")

    print("[DONE]")


if __name__ == "__main__":
    main()