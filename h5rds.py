import hdf5plugin
import scanpy as sc
from scipy import io
import pandas as pd
import subprocess
import os
import argparse
import gzip


def gzip_file(gzipfilename, filename):
    """Gzip compress file"""
    buffer_size = 1024 * 1024
    with open(filename, "rb") as f_in, gzip.open(gzipfilename, "wb") as f_out:
        buffer = f_in.read(buffer_size)
        while buffer:
            f_out.write(buffer)
            buffer = f_in.read(buffer_size)


def create_10x_data(h5ad, path, umap: bool):
    """Create 10X data"""
    try:
        adata = sc.read_h5ad(h5ad)

        if os.path.exists(path):
            print(f"path {path} already exists")
            return False

        os.mkdir(path)

        print("Creating barcodes.tsv ...")
        with open(f"{path}/barcodes.tsv", "w") as f:
            for item in adata.obs_names:
                f.write(item + "\n")

        print("Creating features.tsv ...")
        with open(f"{path}/features.tsv", "w") as f:
            for item in ["\t".join([x, x, "Gene Expression"]) for x in adata.var_names]:
                f.write(item + "\n")

        print("Creating matrix.mtx ...")
        io.mmwrite(f"{path}/matrix", adata.X.T)

        print("Creating metadata.csv ...")
        adata.obs.to_csv(f"{path}/metadata.csv")

        if umap is True:
            umap_df = adata.obsm["X_umap"]
            umap_df = pd.DataFrame(
                umap_df, columns=["UMAP_1", "UMAP_2"], index=adata.obs.index
            )
            umap_df.to_csv(f"{path}/umap_coordinates.csv", index=True)
        return True
    except Exception as exp:
        print(exp)
        return False


def creating_gzips(path):
    try:
        print("Gzipping matrix.mtx... This can take a long time :-(")
        gzip_file(f"{path}/matrix.mtx.gz", f"{path}/matrix.mtx")

        print("Gzipping barcodes.tsv...")
        gzip_file(f"{path}/barcodes.tsv.gz", f"{path}/barcodes.tsv")

        print("Gzipping features.tsv...")
        gzip_file(f"{path}/features.tsv.gz", f"{path}/features.tsv")

        print("Finished gzipping")
    except Exception as exp:
        print(exp)
        return False
    return True


def delete_raw_files(path):
    try:
        print(f"deleting {path}/matrix.mtx")
        os.remove(f"{path}/matrix.mtx")

        print(f"deleting {path}/barcodes.tsv")
        os.remove(f"{path}/barcodes.tsv")

        print(f"deleting {path}/features.tsv")
        os.remove(f"{path}/features.tsv")

        print("finished deleting")

        return True

    except Exception as exp:
        print(exp)
        return False


def parse_commad_line():
    parser = argparse.ArgumentParser(description="Process h5ad for rds extraction")
    parser.add_argument("h5ad", help="h5ad file")
    parser.add_argument("path", help="directory path")
    parser.add_argument(
        "--umap",
        type=bool,
        const=True,
        default=False,
        dest="umap",
        help="if you want umap",
        nargs="?",
    )
    args = parser.parse_args()
    return args


def r_run(path, has_umap=""):
    """Script to run R script"""
    # command = "C:/Path/To/bin/Rscript.exe"
    command = "Rscript"
    arg = "--vanilla"

    print("Passing command to R")
    try:
        prc = subprocess.Popen(
            [command, arg, "seurat_object_creation.R", path, has_umap],
            cwd=os.getcwd(),
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )

        output, error = prc.communicate()

        if prc.returncode == 0:
            print(f"R OUTPUT:\n {output.decode('utf-8')}")
        else:
            print(f"R ERROR:\n {error.decode('utf-8')}")

        return True

    except Exception as exc:
        print(f"r_run error: path={path}, has_umap={has_umap}")
        print(exc)

        return False


def main(h5ad, path, umap: bool):
    print(f"running with parameters:\n  h5ad: {h5ad}\n  path: {path}\n  umap: {umap}")
    if not create_10x_data(h5ad, path, umap):
        return
    if not creating_gzips(path):
        return
    delete_raw_files(path)
    print("All good ;-)")


def main_from_cli():
    args = parse_commad_line()
    main(args.h5ad, args.path, args.umap)
    if args.umap:
        ump = "dummy"
    ump = "_" if args.umap else ""
    r_run(args.path, ump)


if __name__ == "__main__":
    main_from_cli()
