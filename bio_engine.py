"""
CSC 442 - Project 2: DNA & RNA Sequence Analyzer
Core Biology Engine

Handles sequence detection, transcription, translation,
amino acid chain building, and protein characterization.
"""

# ─── Standard Codon Table (mRNA codon → amino acid) ──────────────────────────
CODON_TABLE = {
    "UUU": ("Phenylalanine", "Phe", "F"), "UUC": ("Phenylalanine", "Phe", "F"),
    "UUA": ("Leucine", "Leu", "L"),       "UUG": ("Leucine", "Leu", "L"),
    "CUU": ("Leucine", "Leu", "L"),       "CUC": ("Leucine", "Leu", "L"),
    "CUA": ("Leucine", "Leu", "L"),       "CUG": ("Leucine", "Leu", "L"),
    "AUU": ("Isoleucine", "Ile", "I"),    "AUC": ("Isoleucine", "Ile", "I"),
    "AUA": ("Isoleucine", "Ile", "I"),    "AUG": ("Methionine", "Met", "M"),  # START
    "GUU": ("Valine", "Val", "V"),        "GUC": ("Valine", "Val", "V"),
    "GUA": ("Valine", "Val", "V"),        "GUG": ("Valine", "Val", "V"),
    "UCU": ("Serine", "Ser", "S"),        "UCC": ("Serine", "Ser", "S"),
    "UCA": ("Serine", "Ser", "S"),        "UCG": ("Serine", "Ser", "S"),
    "CCU": ("Proline", "Pro", "P"),       "CCC": ("Proline", "Pro", "P"),
    "CCA": ("Proline", "Pro", "P"),       "CCG": ("Proline", "Pro", "P"),
    "ACU": ("Threonine", "Thr", "T"),     "ACC": ("Threonine", "Thr", "T"),
    "ACA": ("Threonine", "Thr", "T"),     "ACG": ("Threonine", "Thr", "T"),
    "GCU": ("Alanine", "Ala", "A"),       "GCC": ("Alanine", "Ala", "A"),
    "GCA": ("Alanine", "Ala", "A"),       "GCG": ("Alanine", "Ala", "A"),
    "UAU": ("Tyrosine", "Tyr", "Y"),      "UAC": ("Tyrosine", "Tyr", "Y"),
    "UAA": ("Stop", "***", "*"),          "UAG": ("Stop", "***", "*"),
    "CAU": ("Histidine", "His", "H"),     "CAC": ("Histidine", "His", "H"),
    "CAA": ("Glutamine", "Gln", "Q"),     "CAG": ("Glutamine", "Gln", "Q"),
    "AAU": ("Asparagine", "Asn", "N"),    "AAC": ("Asparagine", "Asn", "N"),
    "AAA": ("Lysine", "Lys", "K"),        "AAG": ("Lysine", "Lys", "K"),
    "GAU": ("Aspartic acid", "Asp", "D"), "GAC": ("Aspartic acid", "Asp", "D"),
    "GAA": ("Glutamic acid", "Glu", "E"), "GAG": ("Glutamic acid", "Glu", "E"),
    "UGU": ("Cysteine", "Cys", "C"),     "UGC": ("Cysteine", "Cys", "C"),
    "UGA": ("Stop", "***", "*"),
    "UGG": ("Tryptophan", "Trp", "W"),
    "CGU": ("Arginine", "Arg", "R"),      "CGC": ("Arginine", "Arg", "R"),
    "CGA": ("Arginine", "Arg", "R"),      "CGG": ("Arginine", "Arg", "R"),
    "AGU": ("Serine", "Ser", "S"),        "AGC": ("Serine", "Ser", "S"),
    "AGA": ("Arginine", "Arg", "R"),      "AGG": ("Arginine", "Arg", "R"),
    "GGU": ("Glycine", "Gly", "G"),      "GGC": ("Glycine", "Gly", "G"),
    "GGA": ("Glycine", "Gly", "G"),      "GGG": ("Glycine", "Gly", "G"),
}

# Average molecular weights of amino acids (Daltons)
AA_WEIGHTS = {
    "G": 57.05, "A": 71.08, "V": 99.13, "L": 113.16, "I": 113.16,
    "P": 97.12, "F": 147.18, "W": 186.21, "M": 131.20, "S": 87.08,
    "T": 101.10, "C": 103.14, "Y": 163.18, "H": 137.14, "D": 115.09,
    "E": 129.12, "N": 114.10, "Q": 128.13, "K": 128.17, "R": 156.19,
}

DNA_BASES = set("ATCG")
RNA_BASES = set("AUCG")


def parse_sequence(raw_text):
    """Clean raw input: strip FASTA headers, whitespace, numbers."""
    lines = raw_text.strip().split("\n")
    cleaned = []
    for line in lines:
        line = line.strip()
        if line.startswith(">"):
            continue  # FASTA header
        line = "".join(ch for ch in line if ch.isalpha())
        cleaned.append(line.upper())
    return "".join(cleaned)


def detect_sequence_type(sequence):
    """
    Detect if a sequence is DNA, RNA, or invalid.
    Returns: ("DNA"|"RNA"|"INVALID", explanation_text)
    """
    bases = set(sequence.upper())
    has_t = "T" in bases
    has_u = "U" in bases

    if has_t and has_u:
        return "INVALID", (
            "The sequence contains both Thymine (T) and Uracil (U). "
            "DNA uses T while RNA uses U — a valid biological sequence "
            "cannot contain both. Please check your input."
        )

    invalid_chars = bases - (DNA_BASES | RNA_BASES)
    if invalid_chars:
        return "INVALID", (
            f"The sequence contains invalid character(s): {', '.join(sorted(invalid_chars))}. "
            "DNA sequences should only contain A, T, C, G and RNA sequences "
            "should only contain A, U, C, G."
        )

    if len(bases) == 0:
        return "INVALID", "The sequence is empty."

    if has_t and not has_u:
        return "DNA", (
            "This sequence was identified as DNA because it contains "
            "Thymine (T), which is only found in DNA. DNA (Deoxyribonucleic Acid) "
            "uses the bases Adenine (A), Thymine (T), Cytosine (C), and Guanine (G)."
        )

    if has_u and not has_t:
        return "RNA", (
            "This sequence was identified as RNA because it contains "
            "Uracil (U), which is only found in RNA. RNA (Ribonucleic Acid) "
            "uses the bases Adenine (A), Uracil (U), Cytosine (C), and Guanine (G)."
        )

    # Only A, C, G — ambiguous, default to DNA
    return "DNA", (
        "The sequence contains only A, C, and G bases (no T or U). "
        "Since it does not contain Uracil, it has been assumed to be DNA. "
        "If this is RNA, the result will still be valid since no T→U conversion is needed."
    )


def transcribe(sequence, seq_type, strand_type="non-template"):
    """
    Perform transcription: DNA → mRNA.
    If already RNA, returns as-is.

    strand_type: "non-template" (coding/sense) or "template" (antisense)
    """
    if seq_type == "RNA":
        return sequence, (
            "Since the input is already an RNA sequence, no transcription is needed. "
            "The sequence can be used directly as mRNA for translation."
        )

    complement_map = {"A": "U", "T": "A", "C": "G", "G": "C"}
    direct_map = {"A": "A", "T": "U", "C": "C", "G": "G"}

    if strand_type == "template":
        # Template (antisense) → complement each base, T→U
        mrna = "".join(complement_map[b] for b in sequence)
        explanation = (
            "Transcription was performed on the template (antisense) strand. "
            "During transcription, RNA polymerase reads the template strand in the "
            "3' to 5' direction and builds the mRNA by pairing complementary bases: "
            "A pairs with U, T pairs with A, C pairs with G, and G pairs with C. "
            f"The template DNA strand '{sequence[:20]}...' was converted to mRNA "
            f"'{mrna[:20]}...' by reading each base and writing its RNA complement."
        )
    else:
        # Non-template (coding/sense) → simply replace T with U
        mrna = "".join(direct_map[b] for b in sequence)
        explanation = (
            "Transcription was performed on the non-template (coding/sense) strand. "
            "The coding strand has the same sequence as the mRNA, except DNA uses "
            "Thymine (T) where RNA uses Uracil (U). So to get the mRNA, we simply "
            f"replace every T with U. The DNA sequence '{sequence[:20]}...' becomes "
            f"mRNA '{mrna[:20]}...'."
        )

    return mrna, explanation


def translate(mrna):
    """
    Translate mRNA into a list of (codon, amino_acid_tuple) pairs.
    amino_acid_tuple = (full_name, three_letter, one_letter)
    """
    codons = []
    for i in range(0, len(mrna) - 2, 3):
        codon = mrna[i:i+3]
        if len(codon) == 3:
            aa = CODON_TABLE.get(codon, ("Unknown", "???", "?"))
            codons.append((codon, aa))

    explanation = (
        "Translation is the process of reading the mRNA sequence in groups of three "
        "bases called codons. Each codon specifies a particular amino acid. The ribosome "
        "reads from the start codon (AUG, which codes for Methionine) and continues "
        "until it reaches a stop codon (UAA, UAG, or UGA), which signals the end of "
        "the protein. The amino acids are joined together in order to form a chain "
        "called a polypeptide."
    )

    return codons, explanation


def build_polypeptide(codon_list):
    """
    Extract the polypeptide chain from translation results.
    Starts at first AUG, stops at first stop codon.
    Returns list of amino acid tuples and explanation.
    """
    chain = []
    started = False

    for codon, aa in codon_list:
        if not started and codon == "AUG":
            started = True
        if started:
            if aa[0] == "Stop":
                break
            chain.append(aa)

    explanation = (
        "Amino acids are the building blocks of proteins. There are 20 standard amino acids, "
        "each with a unique chemical structure. When amino acids are linked together by peptide "
        "bonds during translation, they form a chain called a polypeptide. This polypeptide "
        "chain then folds into a specific 3D shape to become a functional protein. "
        "The sequence of amino acids determines the protein's structure and function."
    )

    return chain, explanation


def characterize_protein(chain):
    """
    Provide basic characterisation of the protein from its amino acid chain.
    Returns a dict of properties and explanation.
    """
    if not chain:
        return {}, "No protein could be characterised because no amino acids were found."

    one_letter_seq = "".join(aa[2] for aa in chain)
    total_weight = sum(AA_WEIGHTS.get(aa[2], 110) for aa in chain) - (len(chain) - 1) * 18.02

    # Count amino acid composition
    composition = {}
    for aa in chain:
        name = aa[0]
        composition[name] = composition.get(name, 0) + 1

    # Charge at pH 7 (simplified)
    pos_charged = sum(1 for aa in chain if aa[2] in ("K", "R", "H"))
    neg_charged = sum(1 for aa in chain if aa[2] in ("D", "E"))
    hydrophobic = sum(1 for aa in chain if aa[2] in ("A", "V", "I", "L", "M", "F", "W", "P"))

    props = {
        "sequence": one_letter_seq,
        "length": len(chain),
        "molecular_weight_da": round(total_weight, 2),
        "positive_charged": pos_charged,
        "negative_charged": neg_charged,
        "net_charge": pos_charged - neg_charged,
        "hydrophobic_residues": hydrophobic,
        "hydrophobic_pct": round(100 * hydrophobic / len(chain), 1) if chain else 0,
        "composition": dict(sorted(composition.items(), key=lambda x: -x[1])),
    }

    explanation = (
        "A protein is a large molecule made up of one or more polypeptide chains that fold into "
        "a specific 3D shape. Proteins perform nearly every function in living cells — from "
        "catalysing chemical reactions (enzymes) to providing structural support. The properties "
        "shown above describe this protein's basic characteristics: its molecular weight, charge, "
        "and how many of its amino acids are hydrophobic (water-repelling). These properties "
        "help scientists predict how the protein might behave and what role it plays in the organism."
    )

    return props, explanation
