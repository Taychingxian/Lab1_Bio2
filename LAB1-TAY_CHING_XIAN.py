# LAB1-TAY_CHING_XIAN.py 
# TAY CHING XIAN (A23CS0307)

import streamlit as st
from Bio import ExPASy       
from Bio import SeqIO        
from Bio import Entrez      
from Bio.SeqUtils.ProtParam import ProteinAnalysis 
import pandas as pd          
import urllib.error         

# --- Set Email (As required by Biopython for API access) ---
Entrez.email = "tayxian04@gamil.com"



def retrieve_data(uniprot_id):
    st.write(f"Attempting to retrieve data for UniProt ID: {uniprot_id}...")
    try:
        handle = ExPASy.get_sprot_raw(uniprot_id)
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


def get_basic_analysis(sequence):

    seq_str = str(sequence).upper()
    valid_aas = "ACDEFGHIKLMNPQRSTVWY"
    cleaned_seq_str = "".join([aa for aa in seq_str if aa in valid_aas])

    if not cleaned_seq_str:
        return {"error": "Sequence contains no standard amino acids for analysis."}
        
    try:
        analyst = ProteinAnalysis(cleaned_seq_str)
   
        length = len(cleaned_seq_str)

        aa_comp_percent = analyst.get_amino_acids_percent()

        mol_weight = analyst.molecular_weight()

        iso_point = analyst.isoelectric_point()

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

st.set_page_config(page_title="Lab 1 Protein Analyzer", layout="wide")

st.title("🧬 Lab 1: UniProt Protein Analyzer")
st.write("By: TAY CHING XIAN")


st.sidebar.header("Input")
uniprot_id_input = st.sidebar.text_input("Enter a UniProt ID:", "P01308")
analyze_button = st.sidebar.button("Analyze Protein")

# Main content area
if analyze_button:
    if uniprot_id_input:

        with st.spinner(f"Retrieving data for {uniprot_id_input} from UniProt..."):
            protein_record = retrieve_data(uniprot_id_input)

        if protein_record:
            
            st.header("1. Protein Information")

            st.markdown(f"**ID:** `{protein_record.id}`")
            st.markdown(f"**Name:** `{protein_record.name}`")
            st.markdown(f"**Description:** {protein_record.description}")
        
            if 'organism' in protein_record.annotations:
                st.markdown(f"**Organism:** *{protein_record.annotations['organism']}*")
      
            st.subheader("Sequence Snippet")
            st.code(f"{str(protein_record.seq[:70])}...")
            
            st.divider() 

            st.header("2. Basic Analysis")
            
            with st.spinner("Analyzing protein sequence..."):
                analysis_results = get_basic_analysis(protein_record.seq)
            
            if analysis_results:
         
                col1, col2, col3 = st.columns(3)
                
             
                col1.metric("Sequence Length (aa)", 
                            f"{analysis_results['Sequence Length']}")
                
                col2.metric("Molecular Weight (Da)", 
                            f"{analysis_results['Molecular Weight']:.2f}")
                
                col3.metric("Isoelectric Point (pI)", 
                            f"{analysis_results['Isoelectric Point']:.2f}")

                st.divider()

                
                st.subheader("Amino Acid Composition (%)")
                
                
                aa_comp_dict = analysis_results['Amino Acid Composition']

                aa_df = pd.DataFrame.from_dict(
                    aa_comp_dict, 
                    orient="index", 
                    columns=["Percentage"]
                )

                aa_df["Percentage"] = aa_df["Percentage"] * 100

                st.bar_chart(aa_df)

    else:
        st.sidebar.warning("Please enter a UniProt ID to start.")

else:
    st.info("Enter a UniProt ID in the sidebar (e.g., P01308) and click 'Analyze Protein'.")
