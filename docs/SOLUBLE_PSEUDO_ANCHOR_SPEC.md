Matched soluble pseudo-anchor specification

1. Matched control pairs
   - Final control set: exact species-set matched, median protein length within ±5%, strict all-GLOB DeepTMHMM classification, one-to-one orthogroup matching, minimum-cost by relative length difference.
   - Final matched pairs: 197.

2. Alignment
   - Each matched soluble orthogroup is aligned using MAFFT.
   - All 197 soluble orthogroups have valid alignments.

3. Pseudo-anchor principle
   - Soluble proteins have no true TMDs.
   - For each membrane orthogroup, each modal TMD provides a relative protein-position target.
   - The matched soluble orthogroup is assigned homologous pseudo-anchor columns at corresponding relative positions.

4. Occupancy requirement
   - Only alignment columns with 100% residue occupancy across all matched species are eligible.
   - 197/197 matched soluble orthogroups have sufficient 100%-occupancy columns for all required pseudo-anchors.

5. Anchor ordering
   - Pseudo-anchors are ordered along the soluble alignment to preserve the order of membrane TMDs.
   - Each membrane TMD index maps to one unique soluble pseudo-anchor.

6. Codon independence
   - No soluble low-adaptation segment information is used when constructing or selecting pseudo-anchors.

7. Downstream comparison
   - Soluble low-adaptation segments will be assigned to their nearest pseudo-anchor using the same one-segment-to-one-anchor assignment rule as the membrane analysis.
   - Cross-species positional variance will then be compared with the membrane TMD-associated result.