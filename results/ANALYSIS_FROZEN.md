# Enterobacterales TMD-spacing pilot — Analysis Freeze

Freeze date: 2026-08-31

## Scope

This file records the frozen computational analysis for the Enterobacterales
TMD-relative low-adaptation codon segment project.

After this freeze, biological thresholds, matching criteria, topology QC rules,
feature definitions, and statistical inclusion criteria will not be modified
based on observed outcomes.

## Dataset

- Enterobacterales species: 10
- Strict one-to-one orthogroups with >=8 species: 1836
- Final membrane orthogroups passing topology QC: 225
- Final exact-species matched soluble pairs: 197

## Frozen feature definition

Low-adaptation codons were defined using species-specific relative synonymous
codon adaptation weights.

- Codon weights: species-specific
- Low-adaptation threshold: frozen from the membrane analysis
- Threshold applied unchanged to matched soluble proteins
- Initiation codon excluded
- Low-adaptation segment: >=3 consecutive low-adaptation codons
- Qualifying cross-species unit:
  - species_total >= 8
  - species_with_feature >= 3

## Frozen TMD analysis

Primary anchors:

- TMD start
- TMD end

Supplementary anchor:

- TMD center

Primary statistic:

- Cross-species variance of feature position relative to the homologous anchor

Clustered unit:

- Ortholog group

## Primary observed results

### TMD start

- Qualifying units: 93
- Median variance: 22.7 aa^2

### TMD end

- Qualifying units: 93
- Median variance: 34.55 aa^2

### TMD center

- Qualifying units: 93
- Median variance: 28.6667 aa^2

## Synonymous-composition-preserving null

1000 permutations.

### TMD start

- Observed median variance: 22.7
- Null median-of-medians: 154.3333
- Empirical one-sided p: 0.000999

### TMD end

- Observed median variance: 34.55
- Null median-of-medians: 147.0
- Empirical one-sided p: 0.000999

## Sensitivity: excluding TMD1

Center-based analysis.

- Observed qualifying units: 54
- Observed median variance: 17.208333
- Null permutations: 1000
- Null median-of-medians: 116.812500
- Null minimum: 20.2125
- Null maximum: 949.0
- Null <= observed: 0
- Empirical one-sided p: 0.000999

Conclusion:
The positional-conservation signal is not dependent on the first TMD.

## Fixed-presence positional null

- Qualifying observed units: 80
- Observed median variance: 31.875
- Conditional-null median: 93.1958
- Empirical one-sided p: 0.000999

Conclusion:
The positional-conservation signal is not explained solely by cross-species
feature presence.

## Homologous non-TMD pseudo-anchor control

Frozen design:

- Pseudo-anchors are homologous non-TMD alignment columns
- >=80% occupancy
- >=15 aa from a TMD boundary
- Candidate pool >=10
- Symmetric 1000-permutation analysis

Final comparison:

- Usable membrane orthogroups: 147
- Qualifying observed TMD units in restricted comparison: 66
- Observed median variance: 47.2917
- Non-TMD null median: 979.2917
- Empirical one-sided p: 0.000999

Conclusion:
Low-adaptation segment positioning is substantially more constrained relative
to real TMD anchors than to homologous non-TMD pseudo-anchors within membrane
protein families.

## Family dominance and robustness

- Primary qualifying units: 93
- Orthogroups represented: 76
- Largest family contribution: 3 units (3.23%)
- Top 5 families: 13.98%
- Top 10 families: 24.73%

Leave-one-orthogroup-out START analysis:

- Full median variance: 22.7
- LOO minimum: 20.1667
- LOO maximum: 26.6190
- Largest absolute change: 3.9190

Conclusion:
The result is not driven by one or a small number of orthogroups.

## Final matched-soluble control

Final matching criteria:

- Exact species set
- Protein length within +/-5%
- DeepTMHMM all-GLOB
- One-to-one maximum-cardinality min-cost matching
- 197 matched family pairs
- 100% occupancy pseudo-anchor candidate columns
- Pseudo-anchor absolute relative-position error <=0.05

Final soluble codon dataset:

- Proteins: 1755
- CDS/protein usable: 1755/1755
- Codon positions: 620717
- Low-adaptation segments: 675
- Proteins with >=1 segment: 520

Final soluble conservation:

### START

- Pseudo-anchor units: 1474
- Qualifying units: 47
- Median variance: 24.3333

### END

- Pseudo-anchor units: 1478
- Qualifying units: 48
- Median variance: 22.7083

### Family-paired comparison

START:

- Shared qualifying matched families: 17
- Membrane paired median variance: 28.6667
- Soluble paired median variance: 169.0313
- Median paired difference (soluble - membrane): 140.3646
- Soluble > membrane: 12/17
- Exact one-sided sign-test p: 0.0717
- Family-level sign-flip p: 0.0102

END:

- Shared qualifying matched families: 17
- Membrane paired median variance: 48.375
- Soluble paired median variance: 75.0
- Median paired difference (soluble - membrane): 75.0
- Soluble > membrane: 10/17
- Exact one-sided sign-test p: 0.3145
- Family-level sign-flip p: 0.0738

Interpretation:
The matched-soluble comparison provides anchor-dependent support. The START
comparison favors stronger positional constraint in membrane families, whereas
the END comparison is weaker and not statistically robust. This control will
not be used to claim universal membrane-specificity.

## Frozen interpretation

Primary conclusion:

Low-adaptation synonymous codon segments show evolutionarily constrained
positioning relative to homologous transmembrane-domain boundaries in
Enterobacterales.

Supported statements:

1. Positional conservation is stronger than expected under synonymous
   composition-preserving randomization.
2. The signal persists after excluding TMD1.
3. The signal persists under a fixed-feature-presence positional null.
4. The signal is robust to removal of individual orthogroups.
5. Real TMD anchors show substantially stronger positional constraint than
   homologous non-TMD pseudo-anchors.
6. Matched soluble controls provide stronger support at TMD starts than at
   TMD ends.

Not supported:

- The analysis does not establish that positional constraint is unique to
  membrane proteins.
- The analysis does not demonstrate a causal translation-pause mechanism.
- The analysis does not by itself demonstrate effects on membrane insertion
  or protein function.

## Freeze rule

From 2026-08-31 onward:

- Do not change the low-adaptation threshold based on outcome.
- Do not change the >=3-species feature criterion based on outcome.
- Do not change the >=8-species orthogroup criterion based on outcome.
- Do not change the soluble matching criteria based on outcome.
- Do not change the pseudo-anchor QC threshold based on outcome.
- Do not select alternative null models solely because they produce more
  favorable significance.

Further analyses are limited to:

- visualization,
- descriptive summaries,
- prespecified sensitivity presentation,
- annotation of representative families,
- manuscript preparation,
- or genuinely new biological validation clearly labeled as such.