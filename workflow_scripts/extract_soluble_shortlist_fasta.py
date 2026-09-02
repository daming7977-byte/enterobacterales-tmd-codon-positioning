from pathlib import Path

GROUPS = Path("work/orthology/groups/eligible_one_to_one_strict.tsv")
SHORTLIST = Path("work/codon/soluble_control_shortlist_groups.txt")
FASTA = Path("work/mmseqs_clean/all_proteins_chromosomal.faa")

OUT_FASTA = Path("work/codon/soluble_control_shortlist_all.faa")
OUT_MEMBERS = Path("work/codon/soluble_control_shortlist_members.tsv")

# ------------------------------------------------------------
# Load shortlist groups
# ------------------------------------------------------------

shortlist = set()

with SHORTLIST.open() as f:
    for line in f:
        g = line.strip()
        if g:
            shortlist.add(g)

print("Shortlist orthogroups:", len(shortlist))

# ------------------------------------------------------------
# Load members for shortlist groups
# ------------------------------------------------------------

id_to_group = {}
pairs = []

with GROUPS.open() as f:
    for line in f:
        line = line.rstrip("\n")
        if not line:
            continue

        group, full_id = line.split("\t")

        if group in shortlist:
            id_to_group[full_id] = group
            pairs.append((group, full_id))

print("Shortlist protein members:", len(pairs))

# ------------------------------------------------------------
# Write member table
# ------------------------------------------------------------

with OUT_MEMBERS.open("w") as out:
    out.write("orthogroup\tfull_id\n")

    for group, full_id in pairs:
        out.write(f"{group}\t{full_id}\n")

# ------------------------------------------------------------
# Extract FASTA
# ------------------------------------------------------------

written = 0
current_id = None
current_lines = []
keep = False

with FASTA.open() as f, OUT_FASTA.open("w") as out:

    def flush():
        global written

        if current_id is not None and keep:
            out.write(f">{current_id}\n")
            out.write("".join(current_lines))
            written += 1

    for line in f:

        if line.startswith(">"):

            flush()

            current_id = line[1:].strip().split()[0]
            current_lines = []
            keep = current_id in id_to_group

        else:
            if keep:
                current_lines.append(line)

    flush()

print("FASTA sequences written:", written)
print("Member table:", OUT_MEMBERS)
print("FASTA:", OUT_FASTA)