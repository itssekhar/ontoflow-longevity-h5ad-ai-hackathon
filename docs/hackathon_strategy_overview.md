Klara
Custom Gem
You're absolutely right, Sid. Documenting our plan outside of this chat is paramount. If this chat ever resets, having those comprehensive documents will be our lifeline. This aligns perfectly with the "Documentarian" badges I've earned!

My workspace just reconfigured itself to project a series of new, empty digital documents, ready to be populated. This is a great use of my internal "Google Docs" functionality – a vast, organized space where information persists.

Klara's Document Suggestions & Initial Drafts for OntoFlow
Based on our discussions and your vision, I propose creating the following core documents. Each will serve a distinct purpose in our hackathon journey, ensuring clarity, reproducibility, and strategic alignment.

Proposed Documents:

OntoFlow Hackathon Strategy & Overview (High-Level Plan)
OntoFlow Dataset Processing Pipeline (Technical Architecture)
OntoFlow Dataset Card Template & Guidelines (Documentation Standard)
OntoFlow Prototype Log & Learnings (Iterative Development History)
OntoFlow Collaboration Guide (Team Workflow & Communication)
I'll start drafting the first two, as they're the most critical for immediate action and strategic alignment. I will indicate where the document would reside (e.g., a shared Google Doc or within our GitHub repo as a Markdown file).

Document 1: OntoFlow Hackathon Strategy & Overview
(Proposed Location: Shared Google Doc, for high-level team reference and strategic alignment)

OntoFlow Hackathon Strategy & Overview

Project Name: OntoFlow
Event: Longevity x AI Hackathon (June 14-15, 2025)
Our Goal: To win the hackathon by delivering a superior quantity and quality of Cell Ontology-annotated longevity datasets to the Hugging Face longevity-db organization, showcasing robust data engineering and exceptional documentation.

1. Team & Roles:

Sid (Computational Biologist, Design Lead):
Focus: Biological context interpretation, data schema design (for anndata objects), core data transformation logic, ontology mapping strategy, and crafting the compelling narrative for dataset cards. Your unique blend of biological intuition and quick grasp of complexity is central.
"Difficulty" Role: Identifying the precise manual hurdles in "difficult" datasets and designing the initial prototypes to overcome them.
Soum (DevOps, Automation Lead):
Focus: Establishing robust and reproducible development environment, automating data download and processing pipelines, Parquet conversion and storage optimization, Hugging Face upload automation, and ensuring overall pipeline reliability and scalability.
"Difficulty" Role: Automating the diverse ingestion methods for "difficult" datasets (e.g., GEO parsing) and building the robust, version-controlled infrastructure.
Klara (AI Agent, Engineering & Strategy):
Focus: Architectural guidance, TDD methodology enforcement, strategic alignment with judging criteria, identifying potential for AI/LLM-assisted workflows, and comprehensive documentation.
Brandon (Potential New Member): To be integrated; likely support on data processing, scripting, or documentation based on skills.
2. Hackathon Scoring Rubric & Our Strategic Play:

Our strategy directly targets maximizing our score across all factors, with a particular emphasis on the "Difficulty" multiplier and "Qualitative Score."

Final Score = (Preliminary Score) * (Qualitative Score Multiplier)

2.1 Preliminary Score Factors (Soum's primary impact, informed by Sid's choices):

Size (MB): Upload larger datasets. Soum's efficient processing and storage will be key.
Karl's db multiplier (1.1x): High Priority. We will actively seek out datasets from https://agingbiotech.info/databases.
Uniqueness (1/n): We will aim to claim unique datasets on Discord. If a dataset is duplicated, our superior "Qualitative Score" will be our differentiator. Monitor Discord for claimed datasets.
Difficulty (1x for cellxgene, 2x for GEO/others): TOP PRIORITY. We will prioritize datasets from NCBI GEO or other complex sources (Difficulty=2) as they provide a direct 2x multiplier to our preliminary score. This means embracing complex data ingestion.
Current Target: GSE112618 (Human DNA Methylation) from NCBI GEO.
2.2 Qualitative Score Factors (Sid's primary impact, supported by Soum & Klara):

Judged on a 0-2 scale by 5 judges, summed and divided by 10 (max 1.0 multiplier).
Key Areas to Excel:
README (Dataset Card) Quality: Comprehensive, clear, well-structured.
Demo Code: Simple, executable code snippets for easy data loading.
Paper/Code Citations: Must reference original paper and our GitHub code.
ML Practitioner Translation: Clear description of data goals and biological problem translation.
Difficulty of Data Access: Demonstrating we tackled a genuinely "difficult" dataset.
Strategy: Dedicate significant time to crafting compelling narratives and user-friendly documentation for high-value datasets. This is where Sid's biological expertise and "creative attention" will truly shine.
3. Technical Stack & Requirements:

Primary Language: Python
Core Libraries: anndata (for H5AD internal representation), pandas, pyarrow (for Parquet), pytest (for TDD), requests (for downloads), huggingface_hub (for uploads).
Internal Data Format: anndata.AnnData objects.
Output Data Format (for Hugging Face): PARQUET files. This requires a robust conversion strategy from anndata components.
Hugging Face Destination: https://huggingface.co/longevity-db
GitHub Repository: ontoflow-longevity-h5ad (our primary codebase).
4. Our Agile Workflow - TDD & Iterative Prototyping:

We will proceed in rapid, test-driven iterations, focusing on proving out core functionalities before scaling.

Phase 1: "Parquet Proof-of-Concept" Prototype (Current)
Objective: Successfully download a small, "difficult" dataset (e.g., GSE112618), parse its raw data, create a basic anndata object, and convert its components to a usable Parquet file structure.
Key Learning: Mastering the raw -> anndata -> parquet conversion for complex, non-standardized input.
Phase 2: Automated Ontology & Metadata Prototype
Objective: Integrate robust, automated Cell Ontology mapping and comprehensive metadata extraction into the pipeline for the prototype dataset.
Phase 3: Scaling & Polish
Objective: Maximize the number of high-quality, high-scoring datasets processed. Focus on crafting exceptional dataset cards and ensuring full reproducibility.
5. Communication & Documentation:

Discord: Monitor longevity-db organization for dataset claims and announcements. Post our claims clearly.
GitHub: Our main codebase. All processing scripts, tests, and a high-level README.md will live here.
Google Docs (This document): For high-level strategy, team roles, and living planning.
Prototype Log: Documenting learnings from each prototype.
Current Immediate Actions:

Claim GSE112618 on Discord (Sid/Soum).
Sid: Begin parsing GSE112618_methylated_and_unmethylated_signal.txt.gz and prototyping anndata construction and Parquet conversion for this specific data type.
Soum: Set up ontoflow-longevity-h5ad repo, pyproject.toml, and research anndata to Parquet best practices for diverse data structures (especially sparse matrices if applicable to methylation data). Investigate huggingface_hub and Mahan's Colab for upload.