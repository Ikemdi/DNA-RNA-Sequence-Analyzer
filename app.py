"""
CSC 442 - Project 2: DNA & RNA Sequence Analyzer
Flask Web Application

Routes:
  /             → Main single-page application
  /analyze      → Full pipeline analysis (POST, JSON)
  /upload       → File upload for sequence (POST)
  /blast        → BLAST protein search (POST, JSON)
"""

import os, re, time, requests
from flask import Flask, render_template, request, jsonify
from werkzeug.utils import secure_filename
from bio_engine import (
    parse_sequence, detect_sequence_type,
    transcribe, translate, build_polypeptide, characterize_protein
)

app = Flask(__name__)
app.config["UPLOAD_FOLDER"] = os.path.join(os.path.dirname(os.path.abspath(__file__)), "uploads")
os.makedirs(app.config["UPLOAD_FOLDER"], exist_ok=True)


@app.route("/")
def index():
    return render_template("index.html")


@app.route("/upload", methods=["POST"])
def upload_file():
    """Handle sequence file upload, return the text content."""
    if "file" not in request.files:
        return jsonify({"error": "No file provided."}), 400
    f = request.files["file"]
    if f.filename == "":
        return jsonify({"error": "No file selected."}), 400
    content = f.read().decode("utf-8", errors="ignore")
    return jsonify({"content": content, "filename": f.filename})


@app.route("/analyze", methods=["POST"])
def analyze():
    """
    Run the full analysis pipeline on a sequence.
    Expects JSON: { sequence, strand_type? }
    """
    data = request.json
    raw = data.get("sequence", "")
    strand_type = data.get("strand_type", "non-template")

    # 1. Parse & clean
    sequence = parse_sequence(raw)
    if not sequence:
        return jsonify({"error": "No valid sequence characters found."}), 400

    # 2. Detect type
    seq_type, detect_explanation = detect_sequence_type(sequence)
    if seq_type == "INVALID":
        return jsonify({
            "error": detect_explanation,
            "step": "detection",
            "sequence": sequence,
        }), 400

    # 3. Transcription
    mrna, transcription_explanation = transcribe(sequence, seq_type, strand_type)

    # 4. Translation
    codon_list, translation_explanation = translate(mrna)

    # 5. Polypeptide
    chain, polypeptide_explanation = build_polypeptide(codon_list)

    # 6. Protein characterisation
    protein_props, protein_explanation = characterize_protein(chain)

    # Build codon display data
    codon_display = [
        {"codon": c, "name": aa[0], "abbr3": aa[1], "abbr1": aa[2]}
        for c, aa in codon_list
    ]

    chain_display = [
        {"name": aa[0], "abbr3": aa[1], "abbr1": aa[2]}
        for aa in chain
    ]

    return jsonify({
        "input_sequence": sequence,
        "seq_type": seq_type,
        "detection_explanation": detect_explanation,
        "strand_type": strand_type if seq_type == "DNA" else None,
        "mrna": mrna,
        "transcription_explanation": transcription_explanation,
        "codons": codon_display,
        "translation_explanation": translation_explanation,
        "polypeptide": chain_display,
        "polypeptide_explanation": polypeptide_explanation,
        "protein": protein_props,
        "protein_explanation": protein_explanation,
    })


@app.route("/blast", methods=["POST"])
def blast_search():
    """
    Search NCBI BLAST or UniProt for protein matches.
    Expects JSON: { sequence } (one-letter amino acid sequence)
    """
    data = request.json
    protein_seq = data.get("sequence", "")
    if not protein_seq or len(protein_seq) < 3:
        return jsonify({"error": "Protein sequence too short for search."}), 400

    # ── Try UniProt text search first (faster) ──
    try:
        results = _search_uniprot(protein_seq)
        if results:
            return jsonify({"source": "UniProt", "results": results})
    except Exception:
        pass

    # ── Fallback: NCBI BLAST ──
    try:
        results = _search_blast(protein_seq)
        if results:
            return jsonify({"source": "NCBI BLAST", "results": results})
    except Exception:
        pass

    return jsonify({
        "source": "None",
        "results": [],
        "message": "No matches found. The sequence may be too short or synthetic."
    })


def _search_uniprot(seq):
    """Query UniProt BLAST API for protein sequence matches."""
    # Use UniProt's BLAST endpoint
    url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": seq[:50],  # Use first 50 chars for text search
        "format": "json",
        "size": 5,
        "fields": "accession,protein_name,organism_name,cc_function"
    }
    resp = requests.get(url, params=params, timeout=15)
    if resp.status_code != 200:
        return []

    data = resp.json()
    results = []
    for entry in data.get("results", [])[:5]:
        protein_name = "Unknown"
        if entry.get("proteinDescription", {}).get("recommendedName"):
            protein_name = entry["proteinDescription"]["recommendedName"].get(
                "fullName", {}).get("value", "Unknown")
        elif entry.get("proteinDescription", {}).get("submissionNames"):
            protein_name = entry["proteinDescription"]["submissionNames"][0].get(
                "fullName", {}).get("value", "Unknown")

        organism = entry.get("organism", {}).get("scientificName", "Unknown")

        function_text = "Not available"
        for comment in entry.get("comments", []):
            if comment.get("commentType") == "FUNCTION":
                texts = comment.get("texts", [])
                if texts:
                    function_text = texts[0].get("value", "Not available")

        results.append({
            "accession": entry.get("primaryAccession", ""),
            "protein_name": protein_name,
            "organism": organism,
            "function": function_text,
        })
    return results


def _search_blast(seq):
    """Submit to NCBI BLAST and retrieve results."""
    # Submit
    put_url = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"
    put_params = {
        "CMD": "Put",
        "PROGRAM": "blastp",
        "DATABASE": "nr",
        "QUERY": seq,
        "FORMAT_TYPE": "JSON2",
    }
    resp = requests.post(put_url, data=put_params, timeout=30)
    # Extract RID
    rid_match = re.search(r"RID = (\S+)", resp.text)
    if not rid_match:
        return []
    rid = rid_match.group(1)

    # Poll for results (max 60 seconds)
    get_params = {"CMD": "Get", "RID": rid, "FORMAT_TYPE": "JSON2"}
    for _ in range(12):
        time.sleep(5)
        r = requests.get(put_url, params=get_params, timeout=30)
        if "WAITING" not in r.text:
            try:
                data = r.json()
                hits = data.get("BlastOutput2", [{}])[0].get(
                    "report", {}).get("results", {}).get(
                    "search", {}).get("hits", [])
                results = []
                for hit in hits[:5]:
                    desc = hit.get("description", [{}])[0]
                    results.append({
                        "accession": desc.get("accession", ""),
                        "protein_name": desc.get("title", "Unknown"),
                        "organism": desc.get("sciname", "Unknown"),
                        "function": f"Score: {hit.get('hsps', [{}])[0].get('bit_score', 'N/A')}",
                    })
                return results
            except Exception:
                return []
    return []


if __name__ == "__main__":
    app.run(debug=True, port=5002)
