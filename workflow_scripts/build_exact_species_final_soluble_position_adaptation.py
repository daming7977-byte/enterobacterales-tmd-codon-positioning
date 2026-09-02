from pathlib import Path
import re

members_file = Path(
    "work/codon/exact_species_final_soluble_members.tsv"
)

weights_file = Path(
    "work/codon/species_codon_weights.tsv"
)

thresholds_file = Path(
    "work/codon/species_low_adaptation_thresholds.tsv"
)

positions_out = Path(
    "work/codon/exact_species_final_soluble_codon_position_adaptation.tsv"
)


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
# Load final 197 soluble target proteins
# ------------------------------------------------------------

targets = []

with members_file.open() as f:

    next(f)  # skip header

    for line in f:

        if not line.strip():
            continue

        parts = line.rstrip("\n").split("\t")

        if len(parts) != 4:
            raise RuntimeError(
                f"Expected 4 columns, got {len(parts)}: {line}"
            )

        group, genome, protein_id, full_id = parts

        if full_id != f"{genome}|{protein_id}":
            raise RuntimeError(
                f"Inconsistent member row: {line}"
            )

        targets.append(
            (group, genome, protein_id)
        )


# ------------------------------------------------------------
# Load frozen species-specific codon weights
# ------------------------------------------------------------

weights = {}

with weights_file.open() as f:

    next(f)

    for line in f:

        genome, aa, codon, count, weight = \
            line.rstrip("\n").split("\t")

        weights[(genome, codon)] = float(weight)


# ------------------------------------------------------------
# Load FROZEN low-adaptation thresholds
# from original membrane analysis
# ------------------------------------------------------------

thresholds = {}

with thresholds_file.open() as f:

    next(f)

    for line in f:

        fields = line.rstrip("\n").split("\t")

        genome = fields[0]
        threshold = float(fields[2])

        thresholds[genome] = threshold


print("Frozen thresholds loaded:", len(thresholds))

for genome in sorted(thresholds):

    print(
        genome,
        "threshold=",
        f"{thresholds[genome]:.6f}"
    )


# ------------------------------------------------------------
# Load CDS sequences
# ------------------------------------------------------------

cds_sequences = {}

for path in Path("data").rglob("cds_from_genomic.fna"):

    m = re.search(
        r"(GCF_\d+\.\d+)",
        str(path)
    )

    if not m:
        continue

    genome = m.group(1)

    current_id = None
    seq_parts = []

    def save_record():

        if current_id is not None:

            cds_sequences[
                (genome, current_id)
            ] = "".join(seq_parts).upper()

    with path.open() as f:

        for line in f:

            line = line.rstrip("\n")

            if line.startswith(">"):

                save_record()

                header = line[1:]
                seq_parts = []

                pm = re.search(
                    r"\[protein_id=([^\]]+)\]",
                    header
                )

                if pm:
                    current_id = pm.group(1)

                else:
                    current_id = header.split()[0]

            else:

                seq_parts.append(
                    line.strip()
                )

        save_record()


# ------------------------------------------------------------
# Build position-level table
# using FROZEN membrane-derived thresholds
# ------------------------------------------------------------

records = []

missing = 0

species_positions = {}
species_low = {}

for group, genome, protein_id in targets:

    cds = cds_sequences.get(
        (genome, protein_id)
    )

    if cds is None:
        missing += 1
        continue

    codon_position = 0

    for i in range(
        0,
        len(cds) - 2,
        3
    ):

        codon_position += 1

        codon = cds[i:i+3]

        aa = CODON_TABLE.get(codon)

        if aa is None:
            continue

        # exclude stop codon
        if aa == "*":
            continue

        # same initiation-codon exclusion
        # as frozen membrane analysis
        if codon_position == 1:
            continue

        weight = weights.get(
            (genome, codon)
        )

        if weight is None:
            continue

        threshold = thresholds.get(genome)

        if threshold is None:
            raise RuntimeError(
                f"No frozen threshold for genome: {genome}"
            )

        low = weight <= threshold

        species_positions[genome] = \
            species_positions.get(genome, 0) + 1

        if low:
            species_low[genome] = \
                species_low.get(genome, 0) + 1

        records.append(
            (
                group,
                genome,
                protein_id,
                codon_position,
                codon,
                aa,
                weight,
                low,
            )
        )


# ------------------------------------------------------------
# Output
# ------------------------------------------------------------

with positions_out.open("w") as out:

    out.write(
        "orthogroup\t"
        "genome\t"
        "protein_id\t"
        "codon_position\t"
        "codon\t"
        "amino_acid\t"
        "relative_weight\t"
        "low_adaptation\n"
    )

    for row in records:

        out.write(
            f"{row[0]}\t"
            f"{row[1]}\t"
            f"{row[2]}\t"
            f"{row[3]}\t"
            f"{row[4]}\t"
            f"{row[5]}\t"
            f"{row[6]:.6f}\t"
            f"{row[7]}\n"
        )


print()
print("Targets:", len(targets))
print("Missing CDS:", missing)
print("Codon positions written:", len(records))

print()
print("Observed final soluble low fractions:")

for genome in sorted(species_positions):

    n = species_positions[genome]
    low_n = species_low.get(genome, 0)

    print(
        genome,
        "positions=", n,
        "low=", low_n,
        "fraction=", f"{low_n/n:.4f}"
    )

print()
print("Output:", positions_out)