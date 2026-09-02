from pathlib import Path
import re
from collections import Counter

members_file = Path("work/codon/final_soluble_matched_members.tsv")
protein_fasta = Path("work/codon/final_soluble_matched_all.faa")
out_file = Path("work/codon/final_soluble_cds_protein_qc.tsv")

# Standard bacterial genetic code for internal codons.
CODON_TABLE = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L",
    "TCT":"S","TCC":"S","TCA":"S","TCG":"S",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*",
    "TGT":"C","TGC":"C","TGA":"*","TGG":"W",

    "CTT":"L","CTC":"L","CTA":"L","CTG":"L",
    "CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
    "CGT":"R","CGC":"R","CGA":"R","CGG":"R",

    "ATT":"I","ATC":"I","ATA":"I","ATG":"M",
    "ACT":"T","ACC":"T","ACA":"T","ACG":"T",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K",
    "AGT":"S","AGC":"S","AGA":"R","AGG":"R",

    "GTT":"V","GTC":"V","GTA":"V","GTG":"V",
    "GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "GAT":"D","GAC":"D","GAA":"E","GAG":"E",
    "GGT":"G","GGC":"G","GGA":"G","GGG":"G",
}


def read_fasta(path):
    records = {}
    current = None
    seq = []

    with path.open() as f:
        for line in f:
            line = line.rstrip()

            if line.startswith(">"):
                if current is not None:
                    records[current] = "".join(seq).upper()

                current = line[1:].split()[0]
                seq = []
            else:
                seq.append(line.strip())

        if current is not None:
            records[current] = "".join(seq).upper()

    return records


def translate_dna(seq):
    aa = []

    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        aa.append(CODON_TABLE.get(codon, "X"))

    return "".join(aa)


# Load target members.
targets = []

with members_file.open() as f:
    next(f)  # skip header

    for line in f:
        if not line.strip():
            continue

        group, full_id = line.rstrip("\n").split("\t")

        if "|" not in full_id:
            raise RuntimeError(f"Unexpected protein ID: {full_id}")

        genome, protein_id = full_id.split("|", 1)
        targets.append((group, genome, protein_id, full_id))


# Load protein sequences.
protein_sequences = read_fasta(protein_fasta)


# Load CDS sequences from all 10 assemblies.
cds_sequences = {}

cds_files = list(Path("data").rglob("cds_from_genomic.fna"))

print("CDS files found:", len(cds_files))

for path in cds_files:

    # Get assembly accession from the path.
    match = re.search(r"(GCF_\d+\.\d+)", str(path))

    if not match:
        continue

    genome = match.group(1)

    current_id = None
    seq = []

    def save_current():
        if current_id is not None:
            cds_sequences[(genome, current_id)] = "".join(seq).upper()

    with path.open() as f:
        for line in f:
            line = line.rstrip()

            if line.startswith(">"):

                save_current()

                header = line[1:]
                seq = []

                # Native NCBI CDS FASTA
                m = re.search(r"\[protein_id=([^\]]+)\]", header)

                if m:
                    current_id = m.group(1)

                else:
                    # Reconstructed Salmonella FASTA:
                    # >NP_459006.1
                    current_id = header.split()[0]

            else:
                seq.append(line.strip())

        save_current()


rows = []
status_counter = Counter()

for group, genome, protein_id, full_id in targets:

    protein = protein_sequences.get(full_id)
    cds = cds_sequences.get((genome, protein_id))

    if protein is None:
        status = "missing_protein"
        translated = ""
        cds_len = 0

    elif cds is None:
        status = "missing_cds"
        translated = ""
        cds_len = 0

    elif len(cds) % 3 != 0:
        status = "cds_not_multiple_of_3"
        translated = translate_dna(cds)
        cds_len = len(cds)

    else:
        translated = translate_dna(cds)
        cds_len = len(cds)

        # Remove terminal stop.
        if translated.endswith("*"):
            translated = translated[:-1]

        if translated == protein:
            status = "exact"

        # Bacterial alternative start codons may translate as M
        # in the annotated protein.
        elif (
            len(translated) == len(protein)
            and len(protein) > 1
            and protein[0] == "M"
            and translated[1:] == protein[1:]
        ):
            status = "start_codon_only"

        elif "*" in translated:
            status = "internal_stop"

        elif len(translated) != len(protein):
            status = "length_mismatch"

        else:
            status = "sequence_mismatch"

    status_counter[status] += 1

    rows.append(
        (
            group,
            genome,
            protein_id,
            full_id,
            len(protein) if protein else 0,
            cds_len,
            status,
        )
    )


with out_file.open("w") as out:
    out.write(
        "orthogroup\tgenome\tprotein_id\tfull_id\t"
        "protein_length\tcds_length\tstatus\n"
    )

    for row in rows:
        out.write("\t".join(map(str, row)) + "\n")


print("Targets:", len(targets))

for status, n in sorted(status_counter.items()):
    print(status, n)

print("Output:", out_file)