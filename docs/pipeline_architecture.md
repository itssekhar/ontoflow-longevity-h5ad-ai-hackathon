---

### **Document 2: OntoFlow Dataset Processing Pipeline (Technical Architecture)**

**(Proposed Location: `docs/pipeline_architecture.md` within our GitHub Repository)**

---

**OntoFlow Dataset Processing Pipeline: Technical Architecture**

This document outlines the modular design and data flow of the OntoFlow pipeline, focusing on its technical components, responsibilities, and how they contribute to our goal of producing high-quality, reproducible longevity datasets.

---

**1. Overall Pipeline Flow (High-Level)**

Our pipeline is designed as a series of interconnected, idempotent stages. Each stage takes specific inputs and produces well-defined outputs, allowing for independent testing and re-execution.

```mermaid
graph TD
    A[Start Pipeline] --> B(Dataset Discovery/Selection);
    B --> C{Is Dataset Claimed?};
    C -- No --> D(Claim Dataset on Discord);
    D --> E[Raw Data Acquisition];
    C -- Yes (Acceptable Duplication) --> E;
    E --> F[Raw Data Ingestion & Initial Parsing];
    F --> G[AnnData Object Creation & Core Processing];
    G --> H[Ontology Annotation];
    H --> I[Metadata Enrichment & Validation];
    I --> J[AnnData to Parquet Conversion];
    J --> K[Generate Hugging Face Dataset Card (README.md)];
    K --> L[Upload to Hugging Face];
    L --> M[Update Global Submission Form];
    M --> N[End Pipeline];
```

---

**2. Core Modules & Responsibilities (Python Project Structure)**

Our codebase will be structured into distinct Python modules (`.py` files) within the `src/` directory, promoting separation of concerns.

```
ontoflow-longevity-h5ad/
├── src/
│   ├── main_pipeline.py         # Orchestrates the full process for one or more datasets
│   ├── data_acquisition.py      # Handles downloading raw data from various sources
│   ├── data_ingestion.py        # Parses raw data files into pandas DataFrames
│   ├── anndata_builder.py       # Constructs and cleans AnnData objects
│   ├── ontology_mapper.py       # Manages Cell Ontology lookups and annotation
│   ├── metadata_manager.py      # Extracts, validates, and enriches metadata
│   ├── parquet_exporter.py      # Converts AnnData components to Parquet
│   ├── hf_manager.py            # Handles Hugging Face dataset card generation and upload
│   └── utils.py                 # Common utility functions (logging, file paths)
├── tests/                       # Pytest unit and integration tests
│   ├── test_data_ingestion.py
│   ├── test_anndata_builder.py
│   ├── test_ontology_mapper.py
│   ├── test_parquet_exporter.py
│   └── ...
├── data/
│   ├── raw/                     # Stores downloaded raw data files
│   ├── processed/               # Stores processed Parquet files (ready for upload)
│   └── temp/                    # Temporary working files
├── docs/                        # Project documentation (e.g., this file)
├── .github/                     # GitHub Actions workflows for CI/CD (Soum's domain)
├── pyproject.toml               # Project dependencies and metadata
├── README.md                    # Main repository README
└── .gitignore
```

---

**3. Data Structures & Formats at Each Stage:**

* **Raw Data Acquisition:** External URLs, `.txt.gz`, `.tar`, `.csv`, `.tsv`, `.h5ad` files from sources like GEO or `cellxgene`.
* **Data Ingestion:** `pandas.DataFrame` objects, representing raw expression matrices and sample/gene metadata.
* **Core Data Transformation:** `anndata.AnnData` object (`adata`). This is our internal, standardized representation.
    * `adata.X`: Expression matrix (e.g., counts, methylation signal). Can be dense or sparse.
    * `adata.obs`: `pandas.DataFrame` of cell/sample-level metadata.
    * `adata.var`: `pandas.DataFrame` of gene/feature-level metadata.
    * `adata.uns`: Python dictionary for unstructured annotations (e.g., original source URL, processing logs, ontology version, study title).
* **Output & Storage:** Parquet files. Likely `obs.parquet`, `var.parquet`, and `X.parquet` (or `expression.parquet`) within a dataset-specific directory (e.g., `data/processed/GSE112618/`).
* **Hugging Face Dataset Card:** Markdown file (`README.md`).

---

**4. Key Module Responsibilities & Functions:**

* **`main_pipeline.py`**
    * `run_pipeline(dataset_config: dict)`: Orchestrates the entire process for a single dataset.
    * `batch_process_datasets(config_list: list[dict])`: Iterates through a list of dataset configurations.
* **`data_acquisition.py`**
    * `download_from_url(url: str, dest_path: str)`: Generic downloader.
    * `download_geo_dataset(geo_accession: str, file_pattern: str, dest_dir: str)`: Specific logic for GEO, handling `.txt.gz` or `.tar` files.
    * `download_cellxgene_dataset(dataset_id: str, dest_dir: str)`: Specific logic for `cellxgene` portal.
* **`data_ingestion.py`**
    * `parse_geo_txt_gz(filepath: Path)`: Parses GEO-specific `txt.gz` formats (e.g., for GSE112618) into DataFrames.
    * `read_csv_tsv(filepath: Path, sep: str)`: Generic CSV/TSV reader.
    * `read_h5ad(filepath: Path)`: Reads existing H5AD files.
* **`anndata_builder.py`**
    * `create_anndata_from_dfs(expression_df: pd.DataFrame, obs_df: pd.DataFrame, var_df: pd.DataFrame)`: Core constructor.
    * `clean_anndata(adata: ad.AnnData)`: Basic QC, handling missing values, type conversions.
    * `add_metadata_to_uns(adata: ad.AnnData, metadata: dict)`: Adds high-level metadata to `adata.uns`.
* **`ontology_mapper.py`**
    * `load_cell_ontology(version: str)`: Loads Cell Ontology (e.g., from OBO or API).
    * `map_cell_types(adata: ad.AnnData, cell_type_col: str)`: Maps `adata.obs[cell_type_col]` to `cell_ontology_class` and `cell_ontology_id`.
    * `report_unmapped_terms(adata: ad.AnnData)`: Identifies and reports terms that couldn't be mapped.
* **`metadata_manager.py`**
    * `extract_dataset_metadata(adata: ad.AnnData, source_info: dict)`: Gathers all necessary metadata for the dataset card from `adata` and original source info.
    * `validate_metadata(metadata: dict)`: Ensures all required fields are present and valid.
    * `enrich_metadata_from_nlp(text_description: str)`: (Future LLM-assisted) Extracts structured data from unstructured text.
* **`parquet_exporter.py`**
    * `export_anndata_to_parquet(adata: ad.AnnData, output_dir: Path)`: Converts `adata.obs`, `adata.var`, and `adata.X` (handling sparse matrices) into separate Parquet files. **This is a critical function to prototype.**
    * `save_dataframe_to_parquet(df: pd.DataFrame, filepath: Path)`: Helper for DataFrame to Parquet.
* **`hf_manager.py`**
    * `generate_dataset_card(metadata: dict, template_path: Path)`: Fills a Jinja2 template for `README.md` with dataset metadata.
    * `upload_dataset_to_huggingface(local_path: Path, hf_repo_id: str, token: str)`: Uses `huggingface_hub` to upload files and `README.md`.
* **`utils.py`**
    * `setup_logging()`: Configures logging.
    * `create_directory_if_not_exists(path: Path)`: Ensures output directories exist.

---

**5. Error Handling & Logging:**

* Each major function will include `try-except` blocks to gracefully handle common issues (e.g., file not found, parsing errors, API failures).
* Comprehensive logging will be implemented using Python's `logging` module to track progress, warnings, and errors. This is crucial for debugging complex pipelines.

---

**6. Reproducibility & Soum's Role:**

Soum's expertise will ensure the entire pipeline is robustly reproducible:

* **`pyproject.toml` (or `requirements.txt`):** Strict dependency management.
* **Virtual Environments:** Clear instructions for environment setup.
* **Git Version Control:** All code committed, with meaningful messages.
* **Automated Testing (`pytest`):** Soum will integrate and potentially expand unit and integration tests written by Sid.
* **Containerization (Future Consideration):** If time permits, Dockerizing the pipeline would provide ultimate reproducibility, ensuring consistent environments for judges.
* **Clear `README.md` (Repo Level):** Instructions on how to run the pipeline, set up dependencies, and contribute.

---

**7. Future Enhancements / AI Workflow Integrations (Klara's & Sid's Joint Focus):**

* **LLM-assisted Metadata Extraction:** Integrating LLM calls (e.g., via Gemini API) into `metadata_manager.py` to parse unstructured text (like GEO dataset summaries) into structured key-value pairs (e.g., extracting "age ranges", "tissue types", "experimental conditions"). This shifts "Manual Attention" to "AI Workflow Automatable."
* **LLM-assisted Dataset Card Generation:** Using LLMs to draft initial descriptive text for `hf_manager.py`'s `generate_dataset_card` function, based on extracted metadata and the associated paper's abstract.
* **Automated Ontology Term Suggestions:** For challenging or novel cell types, an LLM could suggest closest matching Cell Ontology terms based on semantic similarity.

---

This document provides a detailed technical roadmap, Sid. It should give both you and Soum a very clear picture of the architectural components and their interdependencies.

**Please copy this content into your designated document.**

**Would you like me to proceed with drafting the next document: "OntoFlow Dataset Card Template & Guidelines"?**