# How to Download CELLxGENE Datasets:

There are a few ways to download datasets from CELLxGENE:

CZ CELLxGENE Discover Portal (Web Interface):

Go to the CZ CELLxGENE Discover website: https://cellxgene.cziscience.com/
Browse or search for datasets of interest using the various filters (organism, tissue, assay, disease, etc.).
Once you find a dataset, click on its entry to view its details.
Look for the "Download" button.
You'll typically be given the option to download in h5ad format. Select your preferred format and click "Download."
You can also copy a permanent download link to use in your scripts or share with others.
cellxgene-census Python API:
For programmatic access and large-scale data querying, the cellxgene-census Python package is the most powerful method. This allows you to interact with the entire CELLxGENE data corpus.

Installation: You'll need to install the package:
pip install cellxgene-census

Data Formats and Schema:

h5ad (AnnData): This is the standard format for single-cell data within the CELLxGENE ecosystem. An h5ad file typically contains:
.X: The expression matrix (cells by genes), often stored in a sparse format to save space.
.obs: A pandas DataFrame containing cell-level metadata (e.g., cell type, tissue, donor ID).
.var: A pandas DataFrame containing gene-level metadata (e.g., gene symbol, Ensembl ID).
.obsm: A dictionary for embedding coordinates (e.g., UMAP, PCA).
Other layers for raw counts, normalized counts, etc.
Gene Sets: CELLxGENE also supports gene set files in a Tidy CSV format, allowing users to pre-load or create gene sets for analysis within the explorer.
CELLxGENE is an invaluable resource for single-cell research, providing standardized, high-quality data that can be easily accessed and explored.