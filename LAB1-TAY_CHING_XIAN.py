# LAB1-TAY_CHING_XIAN.py
# A simple Streamlit app to retrieve and analyze protein sequences from UniProt.

import streamlit as st
from Bio import ExPASy       # To connect to UniProt/ExPASy
from Bio import SeqIO        # To parse sequence files
from Bio import Entrez       # To set email for NCBI/UniProt access
from Bio.SeqUtils.ProtParam import ProteinAnalysis # For protein analysis
import pandas as pd          # To display analysis results neatly
import urllib.error          # For handling connection errors

# --- Set Email (As required by Biopython for API access) ---
# ** PLEASE CHANGE THIS TO YOUR OWN EMAIL ADDRESS **
Entrez.email = "tayxian04@gamil.com"
# **********************************************************


# --- 1. Required Custom Function: retrieve_data ---

def retrieve_data(uniprot_id):
    """
    Retrieves protein sequence data from UniProt based on the given ID.
    
    Parameter:
    uniprot_id (str): The UniProt ID (e.g., "P01308").
    
    Return:
    retrieved record (Bio.SeqRecord.SeqRecord): The full record object from UniProt.
    """
    st.write(f"Attempting to retrieve data for UniProt ID: {uniprot_id}...")
    try:
        # Use get_sprot_raw to get the Swiss-Prot text record
        handle = ExPASy.get_sprot_raw(uniprot_id)
        # Parse the record using SeqIO
        record = SeqIO.read(handle, "swiss")
        handle.close()
        st.success("Data retrieval successful!")
        return record
        
    except urllib.error.HTTPError:
        st.error(f"Could not find UniProt ID: {uniprot_id}. Please check the ID.")
        return None
    except Exception as e:
        st.error(f"An error occurred during data retrieval: {e}")
        return None

# --- 2. Required Custom Function: get_basic_analysis ---

def get_basic_analysis(sequence):
    """
    Analyzes a given protein sequence and returns basic analysis outcomes.
    
    Parameter:
    sequence (Bio.Seq.Seq): The protein sequence object (e.g., record.seq).
    
    Return:
    analysis_outcome (dict): A dictionary containing:
        - Sequence Length
        - Molecular Weight
        - Isoelectric Point
        - Amino Acid Composition
    """
    
    # ProtParam analysis can fail on non-standard AAs (X, B, Z, *)
    # We must clean the sequence string first for accurate analysis
    
    seq_str = str(sequence).upper()
    valid_aas = "ACDEFGHIKLMNPQRSTVWY"
    cleaned_seq_str = "".join([aa for aa in seq_str if aa in valid_aas])

    # Check if the sequence is empty after cleaning
    if not cleaned_seq_str:
        return {"error": "Sequence contains no standard amino acids for analysis."}
        
    try:
        # Create the analysis object from the cleaned sequence
        analyst = ProteinAnalysis(cleaned_seq_str)
        
        # 1. Get sequence length (of the cleaned sequence)
        length = len(cleaned_seq_str)
        
        # 2. Get amino acid composition
        # This returns a dict like {'A': 8.7, 'C': 1.3, ...}
        aa_comp_percent = analyst.get_amino_acids_percent()
        
        # 3. Get molecular weight
        mol_weight = analyst.molecular_weight()
        
        # 4. Get isoelectric point
        iso_point = analyst.isoelectric_point()
        
        # Organize all outcomes into a single dictionary to return
        analysis_outcome = {
            "Sequence Length": length,
            "Molecular Weight": mol_weight,
            "Isoelectric Point": iso_point,
            "Amino Acid Composition": aa_comp_percent
        }
        
        return analysis_outcome

    except ValueError as e:
        st.error(f"Error during sequence analysis: {e}")
        return None

# --- 3. Streamlit Application Interface ---

# Set the page title and layout
st.set_page_config(page_title="Lab 1 Protein Analyzer", layout="wide")

# Main title
st.title("🧬 Lab 1: UniProt Protein Analyzer")
st.write("By: TAY CHING XIAN")

# Sidebar for user input
st.sidebar.header("Input")
uniprot_id_input = st.sidebar.text_input("Enter a UniProt ID:", "P01308")
analyze_button = st.sidebar.button("Analyze Protein")

# Main content area
if analyze_button:
    # Check if input is provided
    if uniprot_id_input:
        
        # --- Step 1: Retrieve Data ---
        with st.spinner(f"Retrieving data for {uniprot_id_input} from UniProt..."):
            protein_record = retrieve_data(uniprot_id_input)
        
        # Check if retrieval was successful
        if protein_record:
            
            # --- Step 2: Display Basic Info ---
            st.header("1. Protein Information")
            
            # Display key info from the record
            st.markdown(f"**ID:** `{protein_record.id}`")
            st.markdown(f"**Name:** `{protein_record.name}`")
            st.markdown(f"**Description:** {protein_record.description}")
            
            # Get organism from annotations
            if 'organism' in protein_record.annotations:
                st.markdown(f"**Organism:** *{protein_record.annotations['organism']}*")
            
            # Show a snippet of the sequence
            st.subheader("Sequence Snippet")
            st.code(f"{str(protein_record.seq[:70])}...")
            
            st.divider() # Visual separator

            # --- Step 3: Analyze and Display Outcome ---
            st.header("2. Basic Analysis")
            
            with st.spinner("Analyzing protein sequence..."):
                analysis_results = get_basic_analysis(protein_record.seq)
            
            if analysis_results:
                # Use columns for a neat metric display
                col1, col2, col3 = st.columns(3)
                
                # Display metrics from the analysis outcome
                col1.metric("Sequence Length (aa)", 
                            f"{analysis_results['Sequence Length']}")
                
                col2.metric("Molecular Weight (Da)", 
                            f"{analysis_results['Molecular Weight']:.2f}")
                
                col3.metric("Isoelectric Point (pI)", 
                            f"{analysis_results['Isoelectric Point']:.2f}")

                st.divider()

                # --- Step 4: Display Amino Acid Composition Creatively ---
                st.subheader("Amino Acid Composition (%)")
                
                # Convert the composition dictionary to a DataFrame for plotting
                aa_comp_dict = analysis_results['Amino Acid Composition']
                
                # Convert dict to DataFrame, orient="index" makes AAs the rows
                aa_df = pd.DataFrame.from_dict(
                    aa_comp_dict, 
                    orient="index", 
                    columns=["Percentage"]
                )
                
                # Convert from ratio (0.08) to percentage (8.0)
                aa_df["Percentage"] = aa_df["Percentage"] * 100
                
                # Display a bar chart
                st.bar_chart(aa_df)

    else:
        st.sidebar.warning("Please enter a UniProt ID to start.")

else:
    # Default message on load
    st.info("Enter a UniProt ID in the sidebar (e.g., P01308) and click 'Analyze Protein'.")
