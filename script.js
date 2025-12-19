const codonToAminoAcid = {
    "UUU": "Phe", "UUC": "Phe", "UUA": "Leu", "UUG": "Leu",
    "CUU": "Leu", "CUC": "Leu", "CUA": "Leu", "CUG": "Leu",
    "AUU": "Ile", "AUC": "Ile", "AUA": "Ile", "AUG": "Met (Start)",
    "GUU": "Val", "GUC": "Val", "GUA": "Val", "GUG": "Val",
    "UCU": "Ser", "UCC": "Ser", "UCA": "Ser", "UCG": "Ser",
    "CCU": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro",
    "ACU": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr",
    "GCU": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala",
    "UAU": "Tyr", "UAC": "Tyr", "UAA": "STOP", "UAG": "STOP", "UGA": "STOP",
    "CAU": "His", "CAC": "His", "CAA": "Gln", "CAG": "Gln",
    "AAU": "Asn", "AAC": "Asn", "AAA": "Lys", "AAG": "Lys",
    "GAU": "Asp", "GAC": "Asp", "GAA": "Glu", "GAG": "Glu",
    "UGU": "Cys", "UGC": "Cys", "UGG": "Trp",
    "CGU": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg",
    "AGU": "Ser", "AGC": "Ser", "AGA": "Arg", "AGG": "Arg",
    "GGU": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly"
};

const restrictionEnzymes = {
    "EcoRI": "GAATTC",
    "BamHI": "GGATCC",
    "HindIII": "AAGCTT",
    "PstI": "CTGCAG"
};

const dnaInput = document.getElementById("dnaInput");
const rnaOutput = document.getElementById("rnaOutput");
const codonList = document.getElementById("codonList");
const validationFeedback = document.getElementById("validationFeedback");
const sequenceStats = document.getElementById("sequenceStats");
const reverseMode = document.getElementById("reverseMode");
const outputLabel = document.getElementById("outputLabel");
const complementOutput = document.getElementById("complementOutput");
const showComplement = document.getElementById("showComplement");
const reverseComplement = document.getElementById("reverseComplement");
const proteinOutput = document.getElementById("proteinOutput");
const orfList = document.getElementById("orfList");
const restrictionList = document.getElementById("restrictionList");
const palindromeList = document.getElementById("palindromeList");
const baseFreqChart = document.getElementById("baseFreqChart");
const compareInput = document.getElementById("compareInput");
const alignmentOutput = document.getElementById("alignmentOutput");

function convertSequence() {
    let inputSequence = dnaInput.value.toUpperCase().trim();
    let isReverse = reverseMode.checked;
    let outputText = "";
    let dnaSequence = "";

    if (!inputSequence) {
        resetOutputs();
        return;
    }

    if (isReverse) {
        if (!/^[AUCG]+$/.test(inputSequence)) {
            validationFeedback.innerText = "⚠️ Invalid RNA Sequence! Use A, U, C, G.";
            resetOutputs();
            return;
        }
        outputText = inputSequence.replace(/U/g, "T");
        dnaSequence = outputText;
        outputLabel.innerText = "DNA Sequence:";
    } else {
        if (!/^[ATCG]+$/.test(inputSequence)) {
            validationFeedback.innerText = "⚠️ Invalid DNA Sequence! Use A, T, C, G.";
            resetOutputs();
            return;
        }
        outputText = inputSequence.replace(/T/g, "U");
        dnaSequence = inputSequence;
        outputLabel.innerText = "RNA Sequence:";
    }

    validationFeedback.innerText = "";
    rnaOutput.innerText = outputText;

    let complementText = "";
    if (showComplement.checked) {
        complementText = getComplement(outputText, !isReverse);
        if (reverseComplement.checked) {
            complementText = complementText.split("").reverse().join("");
        }
        complementOutput.innerText = `Complement: ${complementText}`;
    } else {
        complementOutput.innerText = "";
    }

    let codons = [];
    let aminoAcids = [];
    let proteinSequence = "";
    let orfs = [];
    let sequenceForORF = isReverse ? inputSequence : outputText;
    if (sequenceForORF) {
        for (let i = 0; i < sequenceForORF.length; i += 3) {
            let codon = sequenceForORF.substring(i, i + 3);
            if (codon.length === 3) {
                codons.push(codon);
                let amino = codonToAminoAcid[codon] || "Unknown";
                aminoAcids.push(amino);
                if (codon === "AUG") {
                    proteinSequence += "M";
                } else if (["UAA", "UAG", "UGA"].includes(codon)) {
                    proteinSequence += "*";
                } else if (proteinSequence) {
                    proteinSequence += amino.substring(0, 1);
                }
            }
        }
        codonList.innerHTML = codons.map((codon, index) => {
            let className = codon === "AUG" ? "codon-start" : (["UAA", "UAG", "UGA"].includes(codon) ? "codon-stop" : "");
            return `<li><span class="codon ${className}">${codon}</span> → <span class="amino">${aminoAcids[index]}</span></li>`;
        }).join("");

        proteinOutput.innerText = proteinSequence || "No protein sequence found";
        orfs = findORFs(sequenceForORF);
        orfList.innerHTML = orfs.length ? orfs.map(orf => `<li>${orf}</li>`).join("") : "<li>No ORFs found</li>";
    } else {
        codonList.innerHTML = "";
        proteinOutput.innerText = "";
        orfList.innerHTML = "";
    }

    let restrictionSites = findRestrictionSites(dnaSequence);
    restrictionList.innerHTML = restrictionSites.length ? restrictionSites.map(site => `<li>${site}</li>`).join("") : "<li>No restriction sites found</li>";

    let palindromes = findPalindromes(inputSequence);
    palindromeList.innerHTML = palindromes.length ? palindromes.map(p => `<li>${p}</li>`).join("") : "<li>No palindromes found</li>";

    let freq = calculateBaseFrequency(inputSequence);
    baseFreqChart.innerHTML = `
        <div class="bar" style="height: ${freq.A * 5}px;">A: ${freq.A}</div>
        <div class="bar" style="height: ${freq[isReverse ? "U" : "T"] * 5}px;">${isReverse ? "U" : "T"}: ${freq[isReverse ? "U" : "T"]}</div>
        <div class="bar" style="height: ${freq.C * 5}px;">C: ${freq.C}</div>
        <div class="bar" style="height: ${freq.G * 5}px;">G: ${freq.G}</div>
    `;

    let compareSeq = compareInput.value.toUpperCase().trim();
    if (compareSeq) {
        let alignmentScore = simpleAlignment(inputSequence, compareSeq);
        alignmentOutput.innerText = `Alignment Score: ${alignmentScore} / ${Math.min(inputSequence.length, compareSeq.length)}`;
    } else {
        alignmentOutput.innerText = "";
    }

    let gcContent = calculateGCContent(inputSequence);
    let tm = calculateTm(inputSequence);
    sequenceStats.innerText = `Length: ${inputSequence.length} | Codons: ${codons.length} | GC Content: ${gcContent}% | Tm: ${tm}°C`;
}

function resetOutputs() {
    rnaOutput.innerText = "";
    complementOutput.innerText = "";
    codonList.innerHTML = "";
    proteinOutput.innerText = "";
    orfList.innerHTML = "";
    restrictionList.innerHTML = "";
    palindromeList.innerHTML = "";
    baseFreqChart.innerHTML = "";
    alignmentOutput.innerText = "";
    sequenceStats.innerText = "Length: 0 | Codons: 0 | GC Content: 0% | Tm: 0°C";
}

function getComplement(sequence, isRNA) {
    const pairs = isRNA ? { "A": "U", "U": "A", "C": "G", "G": "C" } : { "A": "T", "T": "A", "C": "G", "G": "C" };
    return sequence.split("").map(base => pairs[base]).join("");
}

function calculateGCContent(sequence) {
    let gcCount = (sequence.match(/[GC]/g) || []).length;
    return sequence.length ? ((gcCount / sequence.length) * 100).toFixed(2) : 0;
}

function calculateTm(sequence) {
    let atCount = (sequence.match(/[ATU]/g) || []).length;
    let gcCount = (sequence.match(/[GC]/g) || []).length;
    return (atCount * 2 + gcCount * 4).toFixed(1);
}

function findORFs(sequence) {
    let orfs = [];
    for (let i = 0; i < sequence.length - 2; i++) {
        if (sequence.substring(i, i + 3) === "AUG") {
            let orf = "AUG";
            for (let j = i + 3; j < sequence.length - 2; j += 3) {
                let codon = sequence.substring(j, j + 3);
                if (codon.length === 3) {
                    orf += codon;
                    if (["UAA", "UAG", "UGA"].includes(codon)) {
                        orfs.push(orf);
                        break;
                    }
                }
            }
        }
    }
    return orfs;
}

function findRestrictionSites(sequence) {
    let sites = [];
    for (let enzyme in restrictionEnzymes) {
        let site = restrictionEnzymes[enzyme];
        if (sequence.includes(site)) {
            sites.push(`${enzyme}: ${site}`);
        }
    }
    return sites;
}

function findPalindromes(sequence) {
    let palindromes = [];
    for (let i = 0; i < sequence.length - 3; i++) {
        for (let len = 4; len <= 10 && i + len <= sequence.length; len += 2) {
            let sub = sequence.substring(i, i + len);
            if (isPalindrome(sub)) {
                palindromes.push(sub);
            }
        }
    }
    return palindromes;
}

function isPalindrome(str) {
    return str === str.split("").reverse().join("");
}

function calculateBaseFrequency(sequence) {
    return {
        A: (sequence.match(/A/g) || []).length,
        T: (sequence.match(/T/g) || []).length,
        U: (sequence.match(/U/g) || []).length,
        C: (sequence.match(/C/g) || []).length,
        G: (sequence.match(/G/g) || []).length
    };
}

function simpleAlignment(seq1, seq2) {
    let score = 0;
    for (let i = 0; i < Math.min(seq1.length, seq2.length); i++) {
        if (seq1[i] === seq2[i]) score++;
    }
    return score;
}

function generateRandomSequence(length = 30) {
    const isRNA = reverseMode.checked; // Check mode at runtime
    const bases = isRNA ? ["A", "U", "C", "G"] : ["A", "T", "C", "G"];
    return Array.from({ length }, () => bases[Math.floor(Math.random() * 4)]).join("");
}

function copyRNA() {
    let outputText = rnaOutput.innerText;
    if (!outputText || outputText.startsWith("⚠️")) return;

    navigator.clipboard.writeText(outputText);
    let copyBtn = document.getElementById("copyBtn");
    copyBtn.innerText = "✅ Copied!";
    copyBtn.classList.add("copied");

    setTimeout(() => {
        copyBtn.innerText = reverseMode.checked ? "📋 Copy DNA" : "📋 Copy RNA";
        copyBtn.classList.remove("copied");
    }, 2000);
}

function exportAll() {
    let data = [
        `${outputLabel.innerText}: ${rnaOutput.innerText}`,
        complementOutput.innerText,
        `Protein Sequence: ${proteinOutput.innerText}`,
        `ORFs: ${Array.from(orfList.children).map(li => li.innerText).join(", ")}`,
        `Restriction Sites: ${Array.from(restrictionList.children).map(li => li.innerText).join(", ")}`,
        `Palindromes: ${Array.from(palindromeList.children).map(li => li.innerText).join(", ")}`,
        sequenceStats.innerText
    ].filter(line => line && !line.includes("No ")).join("\n");

    let blob = new Blob([data], { type: "text/plain" });
    let link = document.createElement("a");
    link.href = URL.createObjectURL(blob);
    link.download = "sequence_analysis.txt";
    link.click();
}

document.getElementById("clearBtn").addEventListener("click", () => {
    dnaInput.value = "";
    compareInput.value = "";
    resetOutputs();
    validationFeedback.innerText = "";
});

document.getElementById("downloadBtn").addEventListener("click", () => {
    let outputText = rnaOutput.innerText;
    if (!outputText || outputText.startsWith("⚠️")) return;

    let blob = new Blob([outputText], { type: "text/plain" });
    let link = document.createElement("a");
    link.href = URL.createObjectURL(blob);
    link.download = reverseMode.checked ? "dna_sequence.txt" : "rna_sequence.txt";
    link.click();
});

document.getElementById("exportAllBtn").addEventListener("click", exportAll);

document.getElementById("themeToggle").addEventListener("click", () => {
    document.body.classList.toggle("light-mode");
    let themeBtn = document.getElementById("themeToggle");
    themeBtn.innerText = document.body.classList.contains("light-mode") ? "🌞" : "🌙";
});

reverseMode.addEventListener("change", () => {
    dnaInput.placeholder = reverseMode.checked ? "Enter RNA sequence (A, U, C, G)" : "Enter DNA sequence (A, T, C, G)";
    convertSequence();
});

showComplement.addEventListener("change", convertSequence);
reverseComplement.addEventListener("change", convertSequence);

document.getElementById("randomSeqBtn").addEventListener("click", () => {
    dnaInput.value = generateRandomSequence();
    convertSequence();
});

dnaInput.addEventListener("input", convertSequence);
compareInput.addEventListener("input", convertSequence);
