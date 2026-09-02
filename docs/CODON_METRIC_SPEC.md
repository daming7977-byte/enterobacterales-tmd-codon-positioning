# Codon Metric Specification

Date: 2026-08-25

## Analysis set

225 topology-qualified alpha-helical orthogroups.

Total protein/CDS pairs: 1982.

CDS-protein QC:
- exact: 1691
- start-codon-only: 291
- failed: 0

## Codon adaptation metric

Codon adaptation is calculated independently for each species.

For each amino acid, synonymous codon frequencies are calculated from the
complete chromosomal CDS set of that species.

For codon c encoding amino acid a:

    w(c) = frequency(c) / max frequency among synonymous codons for a

Thus the most frequently used synonymous codon has w = 1.

Stop codons are excluded.

The annotated initiating codon is excluded from low-adaptation calls because
bacterial alternative initiation codons are translated as methionine and are
subject to initiation-specific constraints.

## Normalization

Codon adaptation is species-specific and amino-acid-specific.

No codon-frequency information is pooled across species.

## Low-adaptation definition

For each species, codon positions in the topology-qualified analysis set are
ranked by relative adaptation.

The bottom decile is classified as low adaptation.

The threshold is defined before examining TMD-relative positions.

## Low-adaptation segment

A low-adaptation segment is defined as at least 3 consecutive low-adaptation
codons.

## Outcome blindness

No TMD-relative codon-spacing result was inspected before fixing this
implementation.