/* ──────────────────────────────────────────────────────────────
   CSC 442 – Project 2  ·  DNA & RNA Sequence Analyzer  ·  JS
   ────────────────────────────────────────────────────────────── */

document.addEventListener("DOMContentLoaded", () => {

    // ── Element references ───────────────────────────────────
    const seqTextarea   = document.getElementById("seq-textarea");
    const dropZone      = document.getElementById("drop-zone");
    const fileInput     = document.getElementById("file-input");
    const browseBtn     = document.getElementById("browse-btn");
    const dropContent   = document.getElementById("drop-content");
    const dropSuccess   = document.getElementById("drop-success");
    const fileNameDisp  = document.getElementById("file-name-display");
    const strandSel     = document.getElementById("strand-selector");
    const btnAnalyze    = document.getElementById("btn-analyze");
    const loader        = document.getElementById("loader");
    const results       = document.getElementById("results-container");

    // Amino acid color palette
    const AA_COLORS = {
        G:"#4ade80",A:"#22c55e",V:"#16a34a",L:"#15803d",I:"#166534",
        P:"#facc15",F:"#f97316",W:"#ea580c",M:"#eab308",S:"#38bdf8",
        T:"#0ea5e9",C:"#f472b6",Y:"#ec4899",H:"#a78bfa",D:"#ef4444",
        E:"#dc2626",N:"#8b5cf6",Q:"#7c3aed",K:"#3b82f6",R:"#2563eb",
    };

    // ── File upload (browse + drag/drop) ─────────────────────
    browseBtn.addEventListener("click", e => { e.preventDefault(); fileInput.click(); });
    dropZone.addEventListener("click", e => {
        if (e.target === dropZone || e.target.closest(".drop-zone-content")) fileInput.click();
    });

    fileInput.addEventListener("change", () => { if (fileInput.files[0]) handleFile(fileInput.files[0]); });

    ["dragenter","dragover"].forEach(e =>
        dropZone.addEventListener(e, ev => { ev.preventDefault(); dropZone.classList.add("drag-over"); }));
    ["dragleave","drop"].forEach(e =>
        dropZone.addEventListener(e, ev => { ev.preventDefault(); dropZone.classList.remove("drag-over"); }));

    dropZone.addEventListener("drop", e => {
        const file = e.dataTransfer.files[0];
        if (file) handleFile(file);
    });

    function handleFile(file) {
        const fd = new FormData();
        fd.append("file", file);
        fetch("/upload", { method: "POST", body: fd })
            .then(r => r.json())
            .then(data => {
                if (data.error) { alert(data.error); return; }
                seqTextarea.value = data.content;
                dropContent.classList.add("hidden");
                dropSuccess.classList.remove("hidden");
                fileNameDisp.textContent = data.filename;
                checkSequenceType();
            });
    }

    // ── Live sequence type detection for strand selector ─────
    seqTextarea.addEventListener("input", debounce(checkSequenceType, 400));

    function checkSequenceType() {
        const raw = seqTextarea.value.toUpperCase().replace(/[^A-Z]/g, "");
        const hasT = raw.includes("T");
        const hasU = raw.includes("U");
        if (hasT && !hasU && raw.length > 0) {
            strandSel.classList.remove("hidden");
        } else {
            strandSel.classList.add("hidden");
        }
    }

    // ── Radio card toggle ────────────────────────────────────
    document.querySelectorAll('.radio-card input').forEach(radio => {
        radio.addEventListener('change', () => {
            document.querySelectorAll('.radio-card').forEach(c => c.classList.remove('active'));
            radio.closest('.radio-card').classList.add('active');
        });
    });

    // ── Analyze ──────────────────────────────────────────────
    btnAnalyze.addEventListener("click", () => {
        const seq = seqTextarea.value.trim();
        if (!seq) { alert("Please enter a sequence."); return; }

        const strandRadio = document.querySelector('input[name="strand"]:checked');
        const strandType = strandRadio ? strandRadio.value : "non-template";

        loader.classList.remove("hidden");
        results.classList.add("hidden");

        fetch("/analyze", {
            method: "POST",
            headers: { "Content-Type": "application/json" },
            body: JSON.stringify({ sequence: seq, strand_type: strandType }),
        })
        .then(r => r.json())
        .then(data => {
            loader.classList.add("hidden");
            if (data.error) { alert(data.error); return; }
            renderResults(data);
            results.classList.remove("hidden");
            document.getElementById("step-detection").scrollIntoView({ behavior: "smooth", block: "start" });
        })
        .catch(err => { loader.classList.add("hidden"); alert("Analysis failed: " + err); });
    });

    // ── Render all results ───────────────────────────────────
    function renderResults(data) {
        // Detection
        const badge = document.getElementById("badge-type");
        badge.textContent = data.seq_type === "DNA" ? "🧬 DNA Detected" : "🔗 RNA Detected";
        badge.className = "result-badge " + data.seq_type.toLowerCase();
        document.getElementById("display-input-seq").textContent = formatSeq(data.input_sequence);
        document.getElementById("explain-detection").textContent = data.detection_explanation;

        // Transcription
        document.getElementById("display-input-seq2").textContent = formatSeq(data.input_sequence);
        document.getElementById("display-mrna").textContent = formatSeq(data.mrna);
        document.getElementById("explain-transcription").textContent = data.transcription_explanation;

        // Translation — codon grid
        const grid = document.getElementById("codon-grid");
        grid.innerHTML = data.codons.map(c => {
            let cls = "codon-cell";
            if (c.codon === "AUG") cls += " start";
            if (c.name === "Stop") cls += " stop";
            return `<div class="${cls}">
                <div class="codon-code">${c.codon}</div>
                <div class="codon-aa">${c.abbr3}</div>
            </div>`;
        }).join("");
        document.getElementById("explain-translation").textContent = data.translation_explanation;

        // Amino acids
        const chainDisp = document.getElementById("chain-display");
        chainDisp.innerHTML = data.polypeptide.map(aa => {
            const color = AA_COLORS[aa.abbr1] || "#666";
            return `<div class="chain-bead" style="background:${color}" title="${aa.name} (${aa.abbr3})">${aa.abbr1}</div>`;
        }).join("");

        const tbody = document.getElementById("chain-tbody");
        tbody.innerHTML = data.polypeptide.map((aa, i) =>
            `<tr><td>${i+1}</td><td>${aa.name}</td><td>${aa.abbr3}</td><td>${aa.abbr1}</td></tr>`
        ).join("");
        document.getElementById("explain-amino").textContent = data.polypeptide_explanation;

        // Protein stats
        const stats = document.getElementById("protein-stats");
        const p = data.protein;
        if (p && p.length) {
            stats.innerHTML = `
                <div class="stat-card"><div class="stat-value">${p.length}</div><div class="stat-label">Amino Acids</div></div>
                <div class="stat-card"><div class="stat-value">${p.molecular_weight_da.toLocaleString()}</div><div class="stat-label">Mol. Weight (Da)</div></div>
                <div class="stat-card"><div class="stat-value">${p.net_charge > 0 ? '+' : ''}${p.net_charge}</div><div class="stat-label">Net Charge</div></div>
                <div class="stat-card"><div class="stat-value">${p.hydrophobic_pct}%</div><div class="stat-label">Hydrophobic</div></div>
                <div class="stat-card"><div class="stat-value">${p.positive_charged}</div><div class="stat-label">+ve Residues</div></div>
                <div class="stat-card"><div class="stat-value">${p.negative_charged}</div><div class="stat-label">−ve Residues</div></div>
            `;
        } else {
            stats.innerHTML = '<p style="color:var(--text-dim)">No protein produced (no start codon found or sequence too short).</p>';
        }
        document.getElementById("explain-protein").textContent = data.protein_explanation;

        // Store protein sequence for BLAST
        window._proteinSeq = p ? p.sequence : "";

        // Reset BLAST section
        document.getElementById("blast-results").classList.add("hidden");
        document.getElementById("blast-loader").classList.add("hidden");
    }

    // ── BLAST search ─────────────────────────────────────────
    document.getElementById("btn-blast").addEventListener("click", () => {
        const seq = window._proteinSeq;
        if (!seq || seq.length < 3) { alert("Protein sequence too short for database search."); return; }

        document.getElementById("blast-loader").classList.remove("hidden");
        document.getElementById("blast-results").classList.add("hidden");

        fetch("/blast", {
            method: "POST",
            headers: { "Content-Type": "application/json" },
            body: JSON.stringify({ sequence: seq }),
        })
        .then(r => r.json())
        .then(data => {
            document.getElementById("blast-loader").classList.add("hidden");
            const container = document.getElementById("blast-results");
            container.classList.remove("hidden");

            document.getElementById("blast-source").textContent = "Source: " + (data.source || "N/A");

            const list = document.getElementById("blast-list");
            if (!data.results || data.results.length === 0) {
                list.innerHTML = `<div class="blast-card"><p>${data.message || "No matches found."}</p></div>`;
                return;
            }
            list.innerHTML = data.results.map(r => `
                <div class="blast-card">
                    <h4>${r.protein_name}</h4>
                    <p><strong>Organism:</strong> ${r.organism}</p>
                    <p><strong>Accession:</strong> ${r.accession}</p>
                    <p><strong>Function:</strong> ${r.function}</p>
                </div>
            `).join("");

            // Update explanation
            document.getElementById("explain-blast").textContent =
                "The system searched the " + data.source + " database using the protein " +
                "sequence generated from your input. The results above show real proteins " +
                "that match or closely resemble your sequence. Each result includes the " +
                "protein's name, the organism it was found in, and what function it performs. " +
                "This helps identify whether your sequence corresponds to a known protein in nature.";
        })
        .catch(() => {
            document.getElementById("blast-loader").classList.add("hidden");
            alert("Database search failed. Please try again.");
        });
    });

    // ── Helpers ───────────────────────────────────────────────
    function formatSeq(seq) {
        // Insert space every 10 chars for readability
        return seq.replace(/(.{10})/g, "$1 ").trim();
    }

    function debounce(fn, ms) {
        let timer;
        return (...args) => { clearTimeout(timer); timer = setTimeout(() => fn(...args), ms); };
    }
});
