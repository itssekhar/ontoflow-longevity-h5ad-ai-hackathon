import cellxgene_census
import anndata
import pandas as pd

# You can specify a Census version. The latest stable version is usually a good choice.
CENSUS_VERSION = "2023-12-15" # Or use cellxgene_census.get_census_version_directory() to see available versions

with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as census:
    # Get all human cell metadata
    human_obs = cellxgene_census.get_obs(
        census,
        organism="homo_sapiens",
        column_names=["assay", "cell_type", "tissue", "sex", "disease"]
    )
    print(f"Total human cells: {len(human_obs)}")
    print(human_obs.head())

    # Get a specific slice of data as an AnnData object
    # For example, all cells from 'lung' tissue in 'homo_sapiens'
    adata = cellxgene_census.get_anndata(
        census,
        organism="homo_sapiens",
        measurement_name="RNA",
        obs_value_filter="tissue == 'lung'",
        # You can also filter by genes, e.g., var_value_filter="feature_name in ['CD19', 'CD3G']"
        # And specify observation and variable column names to retrieve
        obs_column_names=["cell_type", "tissue", "disease", "assay"],
        var_column_names=["feature_name", "feature_id"]
    )

    print(f"AnnData object loaded with {adata.n_obs} cells and {adata.n_vars} genes.")
    print(adata)
    print("\nObservation metadata:")
    print(adata.obs.head())
    print("\nVariable metadata (genes):")
    print(adata.var.head())
    print("\nExpression matrix (first 5x5):")
    # The expression matrix is usually sparse
    print(adata.X[:5, :5].toarray())

    # You can also download the *source* H5AD file for a specific dataset ID
    # (if you know the dataset ID from the CELLxGENE Discover portal)
    dataset_id = "f01d9db1-a76d-41a2-8890-dd3b039a5d83.h5ad"
    cellxgene_census.download_source_h5ad(dataset_id, "downloaded_data.h5ad", census_version=CENSUS_VERSION)
    print(f"Downloaded source H5AD for dataset {dataset_id} to downloaded_data.h5ad")