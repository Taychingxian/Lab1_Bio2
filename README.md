# Lab1_Bio2
# 🧬 Lab 1: UniProt Protein Analyzer

**Student Name:** TAY CHING XIAN  
**Matric ID:** A23CS0307  
**Course:** Bioinformatics Lab 1

## 📖 Description
This is a web-based bioinformatics tool developed using **Python** and **Streamlit**. It automates the retrieval and physicochemical analysis of protein sequences from the **UniProt (Swiss-Prot)** database.

The application allows users to input a UniProt Accession ID, fetches the corresponding biological data via the ExPASy/Entrez API, and performs an immediate analysis of the protein's properties using Biopython.

## ✨ Key Features
* **Automated Data Retrieval:** Fetches live data including Protein Name, Description, Organism, and Sequence directly from UniProt.
* **Physicochemical Analysis:**
    * **Sequence Length:** Total number of amino acids.
    * **Molecular Weight:** Calculated in Daltons (Da).
    * **Isoelectric Point (pI):** The pH at which the protein carries no net electrical charge.
* **Visualization:** Generates an interactive bar chart displaying the percentage composition of amino acids.
