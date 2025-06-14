import anndata
import numpy as np
import pandas as pd
import scanpy as sc

# Create a dummy AnnData object
n_obs = 100
n_vars = 50
X = np.random.rand(n_obs, n_vars) * 10
obs_data = pd.DataFrame({
    'cell_type': np.random.choice(['T cell', 'B cell', 'Macrophage'], n_obs),
    'tissue': np.random.choice(['Lung', 'Spleen'], n_obs),
    'patient_age': np.random.randint(20, 80, n_obs)
})
var_data = pd.DataFrame(index=[f'gene_{i}' for i in range(n_vars)])

adata = anndata.AnnData(X, obs=obs_data, var=var_data)

# Add a simple embedding for visualization (e.g., UMAP or PCA)
# cellxgene expects embeddings to be in .obsm with keys prefixed by 'X_'
sc.pp.pca(adata)
sc.tl.umap(adata, random_state=42)
adata.obsm['X_umap'] = adata.obsm['X_umap'] # Ensure the key is X_umap

# Save the AnnData object to an h5ad file
output_file = "my_single_cell_data.h5ad"
adata.write(output_file)
print(f"Dummy AnnData saved to {output_file}")

# Now, launch cellxgene from your terminal:
# Open your terminal or command prompt in the directory where you saved `my_single_cell_data.h5ad`
# and run the following command:
# cellxgene launch my_single_cell_data.h5ad --open

# To see all options:
# cellxgene launch --help