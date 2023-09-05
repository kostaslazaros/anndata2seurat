# Anndata .h5ad to Seurat .rds conversion

## Powered by: BiHELab, Ionian University.

This is a combined python/R script that can be used to convert .h5ad anndata files to .rds seurat files.

To use, follow the steps below:

## Python requirements

- Python version >= 3.9

- Create virtual environment
   - for windows:
       ``````
       python -m venv venv

       venv\Scripts\activate

       pip install -r requirements.txt
       ``````

    - for linux:

        ``````
        python -m venv venv

        venv/bin/activate

        pip install -r requirements.txt
        ``````

    - for conda:

        ``````
        conda create -n <your_env_name>

        conda activate <your_env_name>

        conda install --file requirements.txt
        ``````

## R requirements

- R version >= 4.2

- Install Seurat:

    ``````
    install.packages('Seurat')
    ``````

## How to run script

- Open console.

- Activate virtual environment.

- Run the following command:

    ``````
    python h5rds.py <h5ad_file> <dest_folder> [--umap]
    ``````
- umap argument is optional.

### Acknowledgement

- This script has been inspired by  [Sandbomics](https://github.com/mousepixels/sanbomics_scripts)

