

**OntoFlow Prototype Log & Learnings**

This document serves as a running record of our iterative prototyping process, capturing key decisions, challenges, solutions, and lessons learned from each prototype iteration. This is crucial for maintaining a clear development history and ensuring we build on our successes and avoid repeating mistakes.

---

**Prototype 1: Parquet Proof-of-Concept (GSE112618)**

* **Objective:** Successfully download GSE112618 data, parse it, create a basic `anndata` object, and export its core components (expression matrix, cell metadata, gene metadata) into a usable Parquet-based structure for Hugging Face.
* **Dataset:** GSE112618 (Human DNA Methylation) - NCBI GEO
    * File: `GSE112618_methylated_and_unmethylated_signal.txt.gz`
* **Key Decisions & Implementation (Sid):**
    * Parsed the tab-separated `txt.gz` file using `pandas.read_csv(..., sep='\t')`.
    * Created `anndata.AnnData` object with:
        * `adata.X`: Methylation signal intensities (rows=samples, columns=probes).
        * `adata.obs`: Sample metadata (extracted from column headers and a manually created mapping).
        * `adata.var`: Probe metadata (if available in the file; otherwise, placeholder).
    * **Parquet Conversion Strategy:**
        * `adata.obs` saved as `obs.parquet` using `pandas.to_parquet()`.
        * `adata.var` saved as `var.parquet` using `pandas.to_parquet()`.
        * `adata.X` (methylation matrix) was converted to a dense NumPy array and then to a Pandas DataFrame with columns `sample_id`, `probe_id`, `methylation_value` before saving as `expression.parquet`.
* **Key Challenges:**
    * Understanding the structure of the `GSE112618_methylated_and_unmethylated_signal.txt.gz` file (column meanings, data types).
    * Representing the methylation matrix `.X` in a Parquet-friendly way.
* **Solutions:**
    * Careful manual inspection of the file and its associated documentation on GEO.
    * Chosen to represent the methylation matrix as a "long" format DataFrame for Parquet, as this preserves sparsity if present.
* **Code Location:** `src/prototype_pipeline.py` (initial version)
* **Tests:**
    * `tests/test_parquet_exporter.py`: Ensured Parquet files were created and could be loaded back into pandas DataFrames.
* **Learnings & Next Steps:**
    * The "long" format for the expression matrix seems viable for Parquet.
    * We need to generalize the parsing of GEO `txt.gz` files to handle different column naming conventions.
    * Metadata extraction was largely manual; we need to automate this.
    * Ontology mapping was skipped for this prototype (not directly applicable to methylation data).

---

**Prototype 2: Automated GEO Parsing & Basic Metadata (GSE112618 - Iteration 2)**

* **Objective:** Automate the parsing of GEO `txt.gz` files and extract basic metadata.
* **Dataset:** GSE112618 (continued)
* **Key Decisions & Implementation (Sid & Soum):**
    * Created a generalized function `src/data_ingestion.py:parse_geo_txt_gz(filepath: Path)` to handle tab-separated GEO files.
    * Implemented basic metadata extraction (dataset name, GEO accession) and added it to `adata.uns`.
* **Key Challenges:**
    * Generalizing the parsing to handle potential variations in GEO file formats.
    * Extracting more meaningful metadata (e.g., sample descriptions) from the GEO page.
* **Solutions:**
    * Implemented flexible column name mapping within `parse_geo_txt_gz`.
    * For now, we're relying on manual metadata entry for the dataset card beyond the basic fields.
* **Code Location:**
    * `src/data_ingestion.py` (updated parsing logic).
    * `src/anndata_builder.py` (added metadata to `adata.uns`).
* **Tests:**
    * `tests/test_data_ingestion.py`: Added tests for the `parse_geo_txt_gz` function.
* **Learnings & Next Steps:**
    * Automated GEO parsing is feasible but requires careful handling of edge cases.
    * Metadata extraction is a major bottleneck. We should explore LLM-assisted methods for this.
    * We still haven't tackled Cell Ontology mapping (not relevant for GSE112618, but crucial for other datasets).

---

**Next Steps (To be added as we progress):**

* Prototype 3: Cell Ontology Mapping (using a `cellxgene` dataset as an example).
* Prototype 4: Hugging Face Dataset Card Generation.
* Prototype 5: Full Pipeline Automation & Scaling.

---

Sid, this is the initial draft of our Prototype Log. It captures the key decisions and learnings from our first two prototypes. As we continue, we'll add more entries to this document, creating a valuable record of our development process.

Please copy this content into `docs/prototype_log.md` in our GitHub repo.

After that, let me know, and I'll be ready to draft the final document: **"OntoFlow Collaboration Guide."**