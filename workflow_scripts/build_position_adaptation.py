from pathlib import Path
from collections import defaultdict, Counter
import math
import re

members_file = Path("work/codon/topology_pass225_members.tsv")
weights_file = Path("work/codon/species_codon_weights.tsv")

positions_out = Path("work/codon/codon_position_adaptation.tsv")
thresholds_out = Path("work/codon/species_low_adaptation_thresholds.tsv")

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


# ------------------------------------------------------------
# Load target proteins
# ------------------------------------------------------------

targets = []

with members_file.open() as f:
    for line in f:
        if not line.strip():
            continue

        group, full_id = line.rstrip("\n").split("\t")
        genome, protein_id = full_id.split("|", 1)

        targets.append((group, genome, protein_id))


# ------------------------------------------------------------
# Load species-specific codon weights
# ------------------------------------------------------------

weights = {}

with weights_file.open() as f:
    next(f)

    for line in f:
        genome, aa, codon, count, weight = line.rstrip("\n").split("\t")
        weights[(genome, codon)] = float(weight)


# ------------------------------------------------------------
# Load all CDS sequences
# key = (genome, protein_id)
# ------------------------------------------------------------

cds_sequences = {}

for path in Path("data").rglob("cds_from_genomic.fna"):

    m = re.search(r"(GCF_\d+\.\d+)", str(path))

    if not m:
        continue

    genome = m.group(1)

    current_id = None
    seq_parts = []

    def save_record():
        if current_id is not None:
            cds_sequences[(genome, current_id)] = "".join(seq_parts).upper()

    with path.open() as f:

        for line in f:
            line = line.rstrip("\n")

            if line.startswith(">"):

                save_record()

                header = line[1:]
                seq_parts = []

                # Native NCBI CDS FASTA
                pm = re.search(r"\[protein_id=([^\]]+)\]", header)

                if pm:
                    current_id = pm.group(1)

                else:
                    # Reconstructed Salmonella CDS FASTA
                    current_id = header.split()[0]

            else:
                seq_parts.append(line.strip())

        save_record()


# ------------------------------------------------------------
# First pass:
# collect all eligible codon-position weights per species
# ------------------------------------------------------------

species_values = defaultdict(list)
records = []

missing = 0

for group, genome, protein_id in targets:

    cds = cds_sequences.get((genome, protein_id))

    if cds is None:
        missing += 1
        continue

    # Position 1 is initiation codon and is excluded.
    codon_position = 0

    for i in range(0, len(cds) - 2, 3):

        codon_position += 1
        codon = cds[i:i+3]

        aa = CODON_TABLE.get(codon)

        # Ignore unknown codons
        if aa is None:
            continue

        # Ignore stop codons
        if aa == "*":
            continue

        # Exclude initiating codon
        if codon_position == 1:
            continue

        weight = weights.get((genome, codon))

        if weight is None:
            continue

        species_values[genome].append(weight)

        records.append(
            (
                group,
                genome,
                protein_id,
                codon_position,
                codon,
                aa,
                weight,
            )
        )


# ------------------------------------------------------------
# Calculate empirical 10th percentile threshold
#
# nearest-rank definition:
# rank = ceil(0.10 * N)
#
# All values <= threshold are classified low.
# Ties are retained, so actual fraction may exceed 10%.
# ------------------------------------------------------------

thresholds = {}

for genome, values in sorted(species_values.items()):

    values = sorted(values)

    n = len(values)

    rank = math.ceil(0.10 * n)

    threshold = values[rank - 1]

    low_n = sum(v <= threshold for v in values)

    thresholds[genome] = threshold

    print(
        genome,
        "positions=", n,
        "threshold=", f"{threshold:.6f}",
        "low_positions=", low_n,
        "actual_fraction=", f"{low_n/n:.4f}",
    )


# ------------------------------------------------------------
# Write thresholds
# ------------------------------------------------------------

with thresholds_out.open("w") as out:

    out.write(
        "genome\tpositions\tthreshold\tlow_positions\tactual_low_fraction\n"
    )

    for genome in sorted(species_values):

        values = species_values[genome]
        threshold = thresholds[genome]

        low_n = sum(v <= threshold for v in values)

        out.write(
            f"{genome}\t"
            f"{len(values)}\t"
            f"{threshold:.6f}\t"
            f"{low_n}\t"
            f"{low_n/len(values):.6f}\n"
        )


# ------------------------------------------------------------
# Write position-level table
# ------------------------------------------------------------

with positions_out.open("w") as out:

    out.write(
        "orthogroup\tgenome\tprotein_id\tcodon_position\t"
        "codon\tamino_acid\trelative_weight\tlow_adaptation\n"
    )

    for (
        group,
        genome,
        protein_id,
        codon_position,
        codon,
        aa,
        weight,
    ) in records:

        low = weight <= thresholds[genome]

        out.write(
            f"{group}\t"
            f"{genome}\t"
            f"{protein_id}\t"
            f"{codon_position}\t"
            f"{codon}\t"
            f"{aa}\t"
            f"{weight:.6f}\t"
            f"{low}\n"
        )


print()
print("Targets:", len(targets))
print("Missing CDS:", missing)
print("Codon positions written:", len(records))
print("Position table:", positions_out)
print("Threshold table:", thresholds_out)