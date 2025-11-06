"""
Lab 1: UniProt Protein Analyzer
Description: This app retrieves and analyzes protein sequences from UniProt.
"""

# Import necessary libraries
import streamlit as st
from Bio import ExPASy  # To connect to UniProt/ExPASy
from Bio import SeqIO    # To parse sequence files
from Bio.SeqUtils.ProtParam import ProteinAnalysis # For protein analysis
import pandas as pd      # To display analysis results neatly

# ==============================================================================
# 1. CUSTOM FUNCTION: retrieve_data(id)
# ==============================================================================

def retrieve_data(id):
    """
    Retrieves a protein record from the UniProt database using the ExPASy server.

    Parameter:
        id (str): The UniProt ID (e.g., "P0DTC2").
    Return:
        record (SeqRecord): A Biopython SeqRecord object, or None if not found.
    """
    
    # !!! IMPORTANT: Replace with your actual email address !!!
    # UniProt requires an email address for database access.
    ExPASy.email = "YOUR_EMAIL@example.com"
    
    try:
        # st.cache_data can be used here to cache results for repeated queries
        @st.cache_data
        def get_record(uniprot_id):
            handle = ExPASy.get_sprot_raw(uniprot_id)
            record = SeqIO.read(handle, "swiss")
            handle.close()
            return record
        
        record = get_record(id)
        return record

    except Exception as e:
        # Handle cases where the ID is not found or a network error occurs
        st.error(f"Error retrieving data for ID '{id}': {e}")
        st.info("Please check the UniProt ID (e.g., P0DTC2) and your internet connection.")
        return None

# ==============================================================================
# 2. CUSTOM FUNCTION: get_basic_analysis(sequence)
# ==============================================================================

def get_basic_analysis(sequence):
    """
    Performs basic protein sequence analysis.

    Parameter:
        sequence (Seq): A Biopython Seq object.
    Return:
        analysis (dict): A dictionary containing analysis results.
    """
    
    # Convert Seq object to a string for ProteinAnalysis
    # We clean the sequence to remove ambiguous 'X' or 'U' residues
    # which can cause errors in molecular_weight or pI calculations.
    clean_sequence = str(sequence).replace("X", "").replace("U", "")
    
    if not clean_sequence:
        return {
            "length": len(sequence),
            "composition": {},
            "molecular_weight": 0,
            "isoelectric_point": 0,
            "error": "Sequence contains only ambiguous residues."
        }

    # Create the analysis object
    analysis_obj = ProteinAnalysis(clean_sequence)
    
    # 1. Get sequence length (from original sequence)
    seq_length = len(sequence)
    
    # 2. Get amino acid composition (returns a dict)
    aa_composition = analysis_obj.get_amino_acids_percent()
    
    # 3. Get molecular weight
    mw = analysis_obj.molecular_weight()
    
    # 4. Get isoelectric point
    pi = analysis_obj.isoelectric_point()
    
    # Store results in a dictionary
    analysis_outcome = {
        "length": seq_length,
        "composition": aa_composition,
        "molecular_weight": mw,
        "isoelectric_point": pi
    }
    
    return analysis_outcome

# ==============================================================================
# 3. STREAMLIT APP INTERFACE
# (This section follows the template provided in 1.4)
# ==============================================================================

# --- Page Configuration ---
st.set_page_config(
    page_title="UniProt Protein Analyzer",
    page_icon="🧬",
    layout="wide"
)

# --- Header Section ---
st.title("🧬 UniProt Protein Sequence Analyzer")
st.markdown("""
Welcome to the simple protein analyzer. This tool uses Biopython to retrieve
protein data from the [UniProt](https://www.uniprot.org/) database.

**Enter a UniProt ID below to start.** (Example: `P0DTC2` for Spike Glycoprotein - COVID-19)
""")

# --- Input Section (Using a sidebar for a clean layout) ---
st.sidebar.header("Input Controls")
uniprot_id = st.sidebar.text_input("Enter UniProt ID:", "P0DTC2")
submit_button = st.sidebar.button("Analyze Protein")

# --- Main Content Area (Processing and Display) ---
if submit_button and uniprot_id:
    # 1. Retrieve the data
    with st.spinner(f"Retrieving record for {uniprot_id}..."):
        record = retrieve_data(uniprot_id)
    
    if record:
        st.success(f"Successfully retrieved data for **{record.id}**")
        
        # 2. Analyze the sequence
        with st.spinner("Analyzing protein sequence..."):
            sequence = record.seq
            analysis_results = get_basic_analysis(sequence)

        # 3. Display the outcomes creatively
        
        # --- Display Basic Info ---
        st.subheader("1. Basic Protein Information")
        st.markdown(f"**Description:** {record.description}")
        st.markdown(f"**Organism:** *{record.annotations.get('organism', 'N/A')}*")
        
        # Use columns for other metadata
        col1, col2 = st.columns(2)
        col1.markdown(f"**Gene Name:** `{record.annotations.get('gene_name', 'N/A')}`")
        col2.markdown(f"**Accessions:** {', '.join(record.annotations.get('accessions', []))}")
        
        # Display sequence in an expandable box
        with st.expander("View Full Protein Sequence (FASTA format)"):
            st.code(record.format("fasta"))

        # --- Display Analysis Outcome ---
        st.subheader("2. Basic Sequence Analysis")
        
        if "error" in analysis_results:
            st.warning(f"Analysis Error: {analysis_results['error']}")
        else:
            # Use columns for key metrics
            col_len, col_mw, col_pi = st.columns(3)
            col_len.metric("Sequence Length (AA)", analysis_results['length'])
            col_mw.metric("Molecular Weight (Da)", f"{analysis_results['molecular_weight']:.2f}")
            col_pi.metric("Isoelectric Point (pI)", f"{analysis_results['isoelectric_point']:.2f}")

            # --- Display Amino Acid Composition ---
            st.subheader("3. Amino Acid Composition (%)")
            
            # Convert the composition dictionary to a Pandas DataFrame for a nice table
            comp_dict = analysis_results['composition']
            df_comp = pd.DataFrame(list(comp_dict.items()), columns=['Amino Acid', 'Percentage'])
            
            # Convert percentage from 0.xx to xx.x%
            df_comp['Percentage'] = df_comp['Percentage'] * 100
            
            # Sort by percentage
            df_comp = df_comp.sort_values(by="Percentage", ascending=False).set_index("Amino Acid")
            
            # Display a bar chart and the data table
            st.bar_chart(df_comp)
            with st.expander("View Composition Data Table"):
                st.dataframe(df_comp.style.format("{:.2f}%"))

elif submit_button and not uniprot_id:
    st.warning("Please enter a UniProt ID in the sidebar.")

# --- Footer Section ---
st.markdown("---")
st.markdown(
    "**Lab 1 Submission by: [TAY CHING XIAN]** | "
)