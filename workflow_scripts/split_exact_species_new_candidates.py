from pathlib import Path

INPUT = Path(
    "work/codon/exact_species_new_candidates_all.faa"
)

OUTDIR = Path(
    "work/codon/exact_species_new_candidate_batches"
)

BATCH_SIZE = 500

OUTDIR.mkdir(parents=True, exist_ok=True)

records = []

header = None
seq = []

with INPUT.open() as f:
    for line in f:
        line = line.rstrip("\n")

        if line.startswith(">"):
            if header is not None:
                records.append(
                    (header, "".join(seq))
                )

            header = line
            seq = []

        else:
            seq.append(line.strip())

    if header is not None:
        records.append(
            (header, "".join(seq))
        )

print("Total sequences:", len(records))

batch_no = 0

for start in range(
    0,
    len(records),
    BATCH_SIZE
):

    batch_no += 1

    batch = records[
        start:start + BATCH_SIZE
    ]

    out = OUTDIR / (
        f"exact_species_batch_{batch_no:02d}.faa"
    )

    with out.open("w") as f:

        for header, sequence in batch:

            f.write(header + "\n")

            for i in range(
                0,
                len(sequence),
                80
            ):
                f.write(
                    sequence[i:i+80]
                    + "\n"
                )

    print(
        out,
        "sequences=",
        len(batch)
    )

print("Batches:", batch_no)