import streamlit as st
import time
import requests
import re
from bio_engine import (
    parse_sequence, detect_sequence_type,
    transcribe, translate, build_polypeptide, characterize_protein
)

# Streamlit Page Config
st.set_page_config(
    page_title="DNA & RNA Sequence Analyzer",
    page_icon="🧬",
    layout="wide"
)

st.title("🧬 DNA & RNA Sequence Analyzer")
st.write("From nucleotide sequence to protein — the complete molecular biology pipeline.")

# ─── API Search Functions ───
def _search_uniprot(seq):
    url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": seq[:50],
        "format": "json",
        "size": 5,
        "fields": "accession,protein_name,organism_name,cc_function"
    }
    resp = requests.get(url, params=params, timeout=15)
    if resp.status_code != 200: return []
    data = resp.json()
    results = []
    for entry in data.get("results", [])[:5]:
        protein_name = "Unknown"
        if entry.get("proteinDescription", {}).get("recommendedName"):
            protein_name = entry["proteinDescription"]["recommendedName"].get("fullName", {}).get("value", "Unknown")
        elif entry.get("proteinDescription", {}).get("submissionNames"):
            protein_name = entry["proteinDescription"]["submissionNames"][0].get("fullName", {}).get("value", "Unknown")
        
        organism = entry.get("organism", {}).get("scientificName", "Unknown")
        function_text = "Not available"
        for comment in entry.get("comments", []):
            if comment.get("commentType") == "FUNCTION":
                texts = comment.get("texts", [])
                if texts: function_text = texts[0].get("value", "Not available")
        results.append({"accession": entry.get("primaryAccession", ""), "protein_name": protein_name, "organism": organism, "function": function_text})
    return results

def _search_blast(seq):
    put_url = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"
    put_params = {"CMD": "Put", "PROGRAM": "blastp", "DATABASE": "nr", "QUERY": seq, "FORMAT_TYPE": "JSON2"}
    resp = requests.post(put_url, data=put_params, timeout=30)
    rid_match = re.search(r"RID = (\S+)", resp.text)
    if not rid_match: return []
    rid = rid_match.group(1)
    
    get_params = {"CMD": "Get", "RID": rid, "FORMAT_TYPE": "JSON2"}
    for _ in range(12):
        time.sleep(5)
        r = requests.get(put_url, params=get_params, timeout=30)
        if "WAITING" not in r.text:
            try:
                data = r.json()
                hits = data.get("BlastOutput2", [{}])[0].get("report", {}).get("results", {}).get("search", {}).get("hits", [])
                results = []
                for hit in hits[:5]:
                    desc = hit.get("description", [{}])[0]
                    results.append({"accession": desc.get("accession", ""), "protein_name": desc.get("title", "Unknown"), "organism": desc.get("sciname", "Unknown"), "function": f"Score: {hit.get('hsps', [{}])[0].get('bit_score', 'N/A')}"})
                return results
            except Exception:
                return []
    return []

# ─── UI ───
st.header("1. 🔤 Sequence Input")
input_method = st.radio("Choose Input Method", ["Type/Paste Sequence", "Upload File"])

raw_text = ""
if input_method == "Type/Paste Sequence":
    raw_text = st.text_area("Enter DNA/RNA sequence", height=150, placeholder="ATGCGATCG...")
else:
    uploaded_file = st.file_uploader("Upload Sequence File", type=["txt", "fasta", "fa"])
    if uploaded_file is not None:
        raw_text = uploaded_file.getvalue().decode("utf-8")
        st.text_area("File Contents Preview", raw_text, height=100, disabled=True)

if raw_text:
    sequence = parse_sequence(raw_text)
    if sequence:
        st.header("2. 🔍 Sequence Detection")
        seq_type, detect_explanation = detect_sequence_type(sequence)
        
        if seq_type == "INVALID":
            st.error(detect_explanation)
        else:
            st.success(f"**{seq_type} Detected**")
            st.info(detect_explanation)
            
            strand_type = "non-template"
            if seq_type == "DNA":
                strand_type = st.radio("Select DNA Strand Type", ["non-template", "template"], 
                                       format_func=lambda x: "Non-template (Coding/Sense)" if x == "non-template" else "Template (Antisense)")
            
            if st.button("🔬 Analyze Pipeline", type="primary"):
                # Transcription
                st.header("3. 📝 Transcription")
                mrna, transcription_explanation = transcribe(sequence, seq_type, strand_type)
                st.code(f"mRNA: {mrna}", language="text")
                st.info(transcription_explanation)
                
                # Translation
                st.header("4. 🔄 Translation")
                codon_list, translation_explanation = translate(mrna)
                st.write("**Codons:**")
                st.write(", ".join([f"`{c}` ({aa[1]})" for c, aa in codon_list]))
                st.info(translation_explanation)
                
                # Polypeptide
                st.header("5. 🧱 Amino Acids (Polypeptide Chain)")
                chain, polypeptide_explanation = build_polypeptide(codon_list)
                if not chain:
                    st.warning("No start codon (AUG) found, or sequence ended before start.")
                else:
                    chain_str = "-".join([aa[2] for aa in chain])
                    st.code(chain_str, language="text")
                    st.info(polypeptide_explanation)
                    
                    # Protein Characterization
                    st.header("6. 🧫 Protein Characterisation")
                    props, protein_explanation = characterize_protein(chain)
                    
                    col1, col2, col3, col4 = st.columns(4)
                    col1.metric("Amino Acids", props["length"])
                    col2.metric("Molecular Weight", f"{props['molecular_weight_da']} Da")
                    col3.metric("Net Charge", props["net_charge"])
                    col4.metric("Hydrophobic", f"{props['hydrophobic_pct']}%")
                    
                    st.info(protein_explanation)
                    
                    st.subheader("🌐 Database Lookup")
                    if st.button("Search BLAST / UniProt"):
                        with st.spinner("Searching databases..."):
                            protein_seq = props["sequence"]
                            results = _search_uniprot(protein_seq)
                            source = "UniProt"
                            if not results:
                                results = _search_blast(protein_seq)
                                source = "NCBI BLAST"
                            
                            if results:
                                st.success(f"Matches found via {source}:")
                                for r in results:
                                    with st.container(border=True):
                                        st.write(f"**{r['protein_name']}**")
                                        st.write(f"*Organism:* {r['organism']}")
                                        st.write(f"*Accession:* {r['accession']}")
                                        st.write(f"*Function:* {r['function']}")
                            else:
                                st.warning("No matches found. Sequence may be too short or synthetic.")
