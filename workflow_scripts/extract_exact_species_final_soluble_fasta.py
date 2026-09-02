from pathlib import Path
import csv

MEMBERS = Path(
    "work/codon/exact_species_final_soluble_members.tsv"
)

ALL_FASTA = Path(
    "work/mmseqs_clean/all_proteins_chromosomal.faa"
)

OUT = Path(
    "work/codon/exact_species_final_soluble_all.faa"
)

wanted = set()

with MEMBERS.open() as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:
        wanted.add(row["full_id"])

print("Target proteins:", len(wanted))


found = set()
written = 0

with ALL_FASTA.open() as inp, OUT.open("w") as out:

    keep = False
    current_id = None

    for line in inp:

        if line.startswith(">"):

            current_id = line[1:].strip().split()[0]

            keep = current_id in wanted

            if keep:
                found.add(current_id)
                written += 1
                out.write(line)

        elif keep:
            out.write(line)


missing = wanted - found

print("Proteins written:", written)
print("Unique proteins found:", len(found))
print("Missing proteins:", len(missing))

if missing:
    print()
    print("Missing IDs:")
    for x in sorted(missing):
        print(x)

print()
print("Output:", OUT)