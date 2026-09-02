# Matched Non-TMD Control Specification

Date: 2026-08-25

## Purpose

Test whether cross-species conservation of low-adaptation codon segment
position is specifically associated with homologous transmembrane helices,
rather than representing a general tendency for conserved codon features to
occur at similar protein-relative positions.

## Control anchor definition

For each eligible protein-TMD unit:

1. Generate one matched non-TMD anchor within the same protein.
2. Match the anchor as closely as possible to the relative sequence position
   of the corresponding TMD center.
3. The control anchor must lie outside all predicted TM helices.
4. The control anchor must be at least 30 amino acids from the nearest TM
   helix boundary.
5. If no valid matched anchor exists, that protein-TMD unit is excluded from
   the matched-control comparison.
6. No low-adaptation feature information is used when choosing the control
   anchor.

## Analysis

The same low-adaptation segments, presence criteria, nearest-feature
assignment, and cross-species variance statistics used for the real TMD
analysis will be applied unchanged to matched non-TMD anchors.

The control definition is frozen before examining matched-control outcomes.## Distance-threshold feasibility rule

Before examining any matched-control codon outcome, candidate minimum
distances of 5, 10, 15, 20, 25, and 30 aa were evaluated using topology
geometry alone.

The strictest threshold retaining at least 80% of protein-TMD units was
selected.

Observed feasibility:

- 5 aa: 100.0%
- 10 aa: 97.17%
- 15 aa: 87.27%
- 20 aa: 69.66%
- 25 aa: 56.90%
- 30 aa: 46.44%

Therefore the frozen minimum distance is 15 aa.

No low-adaptation codon outcome was examined during threshold selection.