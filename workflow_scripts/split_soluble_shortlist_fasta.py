from pathlib import Path

INPUT = Path("work/codon/soluble_control_shortlist_all.faa")
OUTDIR = Path("work/codon/soluble_control_batches")

BATCH_SIZE = 500

OUTDIR.mkdir(parents=True, exist_ok=True)

records = []

current_header = None
current_seq = []

with INPUT.open() as f:
    for line in f:
        line = line.rstrip("\n")

        if line.startswith(">"):
            if current_header is not None:
                records.append(
                    (current_header, "".join(current_seq))
                )

            current_header = line
            current_seq = []

        else:
            current_seq.append(line.strip())

    if current_header is not None:
        records.append(
            (current_header, "".join(current_seq))
        )

print("Total sequences:", len(records))

batch_number = 0

for start in range(0, len(records), BATCH_SIZE):

    batch_number += 1

    batch = records[start:start + BATCH_SIZE]

    out = OUTDIR / f"soluble_batch_{batch_number:02d}.faa"

    with out.open("w") as f:
        for header, seq in batch:
            f.write(header + "\n")

            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + "\n")

    print(
        out,
        "sequences=",
        len(batch)
    )

print("Batches:", batch_number)