# Abstract

Synonymous codon usage can influence translation and cotranslational protein biogenesis, but whether local synonymous codon organization is evolutionarily constrained relative to membrane-protein topology remains unclear. Here, I developed a comparative framework to test whether low-adaptation synonymous codon segments show conserved positioning relative to homologous transmembrane-domain (TMD) boundaries across Enterobacterales.

I analyzed strict one-to-one orthologous groups spanning 10 representative genomes and identified 225 topology-qualified multi-pass membrane-protein families. Species-specific synonymous codon adaptation weights were used to define low-adaptation codons, and segments of at least three consecutive low-adaptation codons were mapped relative to homologous TMD starts and ends. Cross-species positional conservation was quantified as the variance of segment positions relative to each homologous TMD anchor.

Qualifying TMD units showed substantially lower positional variance than expected under 1,000 within-protein synonymous-codon permutations for both TMD starts and ends. The signal persisted after exclusion of first TMDs, under a fixed-presence conditional positional null, and in leave-one-orthogroup-out analyses. Homologous non-TMD pseudo-anchors within the same membrane-protein families produced markedly higher positional variance than bona fide TMD anchors. A separately constructed set of 197 one-to-one matched soluble-protein families provided additional but anchor-dependent evidence, with stronger support for TMD-start-relative constraint than for TMD-end-relative constraint.

Together, these results reveal evolutionarily constrained positioning of low-adaptation synonymous codon segments relative to homologous TMD boundaries in Enterobacterales. This pattern is not readily explained by synonymous composition, feature presence, family dominance, or generic conservation around homologous sequence positions. Although the underlying mechanism remains unresolved, the findings identify TMD-relative synonymous codon positioning as a testable feature of membrane-protein coding-sequence evolution.

# Introduction

Membrane-protein biogenesis requires the coordinated translation, targeting, insertion, and folding of nascent polypeptide chains at the cytoplasmic membrane. In bacteria, ribosomes translating inner-membrane proteins are cotranslationally targeted to the membrane through pathways involving the signal recognition particle (SRP), its receptor FtsY, the SecYEG translocon, and the YidC insertase (Mercier et al., 2022). For multi-pass membrane proteins, successive transmembrane domains (TMDs) emerge during translation and must acquire the appropriate topology and interactions within the membrane (Mercier et al., 2022; Smalinskaitė and Hegde, 2023). The timing of translation relative to the emergence of individual TMDs may therefore influence membrane-protein biogenesis.

Synonymous codon usage provides one potential layer of translational regulation. Although synonymous substitutions do not alter the encoded amino-acid sequence, codon choice can influence translation elongation, translational accuracy, mRNA abundance, and cotranslational protein folding (Liu et al., 2021; Moss et al., 2024). Differences in synonymous codon usage can generate nonuniform decoding behavior along an mRNA, and experimental and comparative studies have linked synonymous codon organization to cotranslational folding and protein biogenesis (Liu et al., 2021; Moss et al., 2024). At the same time, codon usage is shaped by multiple genome-wide and gene-specific factors, making it difficult to distinguish locally constrained codon organization from background synonymous variation.

Previous studies have reported position-specific evolutionary conservation of synonymous rare-codon clusters and other signatures consistent with selection on local translation dynamics across homologous coding sequences (Chartier et al., 2012; Chaney et al., 2017; Jacobs and Shakhnovich, 2017). These observations suggest that selection can act not only on overall codon usage but also on the positions of local synonymous-codon features. However, it remains unclear whether such features are evolutionarily constrained specifically relative to homologous individual TMD boundaries. In particular, this relationship has not been systematically quantified in a strict orthology framework while separating positional conservation from synonymous composition, feature presence, generic conservation around homologous non-TMD positions, and comparable patterns in soluble proteins.

Here, I developed a comparative framework to test whether locally low-adaptation synonymous codon segments show evolutionarily constrained positioning relative to homologous TMDs in Enterobacterales. Using strict one-to-one orthologous groups across 10 representative genomes, curated membrane-topology information, species-specific synonymous codon adaptation metrics, and multiple positional null models, I quantified cross-species variation in low-adaptation segment positions relative to homologous TMD starts and ends. I found that these segments showed substantially lower positional variance relative to bona fide TMD boundaries than expected under within-protein synonymous-codon permutation and homologous non-TMD pseudo-anchor controls. Matched soluble-family analyses provided additional but anchor-dependent evidence, with stronger support for TMD-start-relative constraint than for TMD-end-relative constraint. These results reveal evolutionarily constrained positioning of low-adaptation synonymous codon segments relative to homologous TMD boundaries, while leaving open whether this pattern directly reflects modulation of translation kinetics, membrane insertion, cotranslational folding, or another selective constraint.

# Results
Construction of a comparative membrane-protein dataset across Enterobacterales

To examine whether synonymous codon organization is evolutionarily constrained relative to transmembrane domains (TMDs), I assembled a comparative dataset spanning 10 representative Enterobacterales genomes. Protein orthology was defined using a strict reciprocal-best-hit procedure followed by a one-to-one family filter, with orthologous groups required to be represented in at least eight species. This procedure yielded 1,836 orthologous groups suitable for comparative analysis.

Membrane-protein families were subsequently subjected to topology quality control. Reviewed topology annotations were used as the primary source, supplemented by DeepTMHMM-assisted topology completion and consistency assessment. Families with ambiguous TMD order, inconsistent transmembrane architecture, beta-barrel topology, or insufficient numbers of TMDs were excluded. Following these filters, 225 membrane-protein orthogroups were retained for downstream analysis.

For each species, synonymous codon adaptation was quantified using species-specific relative codon weights calculated within synonymous amino-acid groups. Low-adaptation codons were defined using frozen species-specific thresholds, and low-adaptation segments were identified as runs of at least three consecutive low-adaptation codons. TMD starts and ends were treated as the primary positional anchors, with TMD centers examined as a supplementary reference. Cross-species positional conservation was quantified as the variance in low-adaptation segment position relative to homologous TMD anchors.

Low-adaptation codon segments show conserved positioning relative to TMD boundaries

Across the 225 topology-qualified membrane-protein families, 93 homologous orthogroup–TMD units satisfied the primary inclusion criterion of representation in at least eight species and presence of a low-adaptation segment in at least three species. The median cross-species positional variance was 22.70 aa² relative to TMD starts and 34.55 aa² relative to TMD ends. The corresponding median variance relative to TMD centers was 28.67 aa².

Representative families illustrated that this positional conservation was not confined to a single topological context. In the FocA-family orthogroup_540, low-adaptation segments occurred approximately 27 residues upstream of the first TMD in six of nine species, with a TMD-start-relative variance of only 0.075 aa². In the YghB/DedA-family orthogroup_36, low-adaptation segments overlapped the second TMD in five of eight species, whereas in the ClcA-family orthogroup_105 they clustered within the eighth TMD in five of nine species. These examples indicate that conserved low-adaptation segments can occur both upstream of and within TMDs and are not restricted to a particular TMD position within the protein.

To determine whether the observed positional conservation could arise from synonymous codon composition alone, I performed 1,000 within-protein synonymous-codon permutations. The observed median variance relative to TMD starts was 22.70 aa², compared with a null median-of-medians of 154.33 aa². None of the 1,000 permutations produced a value as low as the observed statistic, corresponding to an empirical one-sided P value of 9.99 × 10⁻⁴. A similar pattern was observed for TMD ends, for which the observed median variance was 34.55 aa² compared with a null median of 147.00 aa² (empirical one-sided P = 9.99 × 10⁻⁴).

Together, these results indicate that low-adaptation synonymous codon segments exhibit substantially stronger cross-species positional conservation relative to homologous TMD boundaries than expected from synonymous codon composition alone.

TMD-relative positional conservation is robust to positional, family-level, and homologous non-TMD controls

I next tested whether the observed positional conservation could be explained by factors unrelated to TMD-associated organization. Because N-terminal sequence context may impose distinct constraints on translation and targeting, I first repeated the analysis after excluding all first TMDs. Fifty-four homologous TMD units remained eligible under the frozen inclusion criteria. Their observed median positional variance was 17.21 aa², substantially below the median of the corresponding 1,000-permutation null distribution, 116.81 aa². None of the permutations produced a median variance as low as the observed value (empirical one-sided P = 9.99 × 10⁻⁴), indicating that the signal was not driven solely by first-TMD-associated effects.

I next asked whether the result could arise simply because low-adaptation segments were repeatedly present in the same orthologous families, irrespective of their precise position. To distinguish feature presence from positional conservation, I applied a fixed-presence positional null in which the observed cross-species pattern of feature occurrence was retained while segment positions were randomized. Eighty observed TMD units were eligible for this comparison. Their restricted observed median variance was 31.88 aa², whereas the corresponding conditional-null median was 93.20 aa² (empirical one-sided P = 9.99 × 10⁻⁴). Thus, repeated occurrence of low-adaptation segments across species was insufficient by itself to explain their constrained positioning relative to TMDs.

The signal was also broadly distributed across membrane-protein families rather than being dominated by a small number of orthogroups. The 93 primary qualifying TMD units were distributed across 76 orthogroups, with the largest individual family contributing only three units, corresponding to 3.23% of the total. The five and ten most highly represented families accounted for 13.98% and 24.73% of qualifying units, respectively. In a leave-one-orthogroup-out analysis of TMD-start-relative variance, the full-data median of 22.70 aa² remained stable, ranging from 20.17 to 26.62 aa² after removal of individual orthogroups.

Finally, I tested whether similarly conserved positioning would be observed relative to homologous sequence positions not associated with TMDs. Within the same membrane-protein families, homologous non-TMD alignment columns were selected as pseudo-anchors under frozen occupancy, distance, and candidate-pool criteria. In the restricted comparison, 66 observed TMD units had a median positional variance of 47.29 aa², whereas the median across 1,000 homologous non-TMD pseudo-anchor permutations was 979.29 aa². None of the permutations produced a median variance as low as the observed TMD-associated value (empirical one-sided P = 9.99 × 10⁻⁴).

Together, these analyses show that the observed positional conservation cannot be readily explained by first-TMD effects, feature-presence patterns, dominance by a small number of protein families, or generic conservation around arbitrary homologous sequence positions. Instead, low-adaptation segments are preferentially position-constrained relative to bona fide TMD boundaries.

Matched soluble-family controls reveal anchor-dependent differences in positional constraint

To further assess whether the observed positional constraint was associated with membrane-protein architecture rather than reflecting a general feature of protein-coding sequences, I constructed a matched soluble-protein control set. Each membrane-protein orthogroup was matched to a soluble orthogroup with the same species composition and a protein-length difference within 5%. Candidate soluble families were required to be classified as globular by DeepTMHMM, and one-to-one matching was obtained using a maximum-cardinality minimum-cost procedure. This yielded 197 matched membrane–soluble family pairs.

Pseudo-anchors in the soluble proteins were defined using fully occupied alignment columns positioned to approximate the relative locations of the corresponding membrane TMD boundaries. Pseudo-anchors with an absolute relative-position error greater than 0.05 were excluded. The same frozen species-specific codon weights, low-adaptation thresholds, and segment definition used for membrane proteins were then applied without modification to the matched soluble dataset. In total, 1,755 soluble proteins were analyzed, yielding 675 low-adaptation segments in 520 proteins.

At the level of all qualifying soluble pseudo-anchor units, median positional variances were 24.33 aa² for start-matched pseudo-anchors and 22.71 aa² for end-matched pseudo-anchors. Because qualifying units were sparse and did not necessarily occur at the same individual TMD indices in both members of a matched pair, I additionally performed a family-level paired analysis. For each membrane family and its matched soluble family, qualifying unit variances were summarized by their within-family median, and only matched pairs in which both members contained at least one qualifying unit were retained.

For TMD-start-relative positions, 17 matched family pairs were informative. The paired membrane-family median variance was 28.67 aa², compared with 169.03 aa² in the matched soluble families, corresponding to a median paired difference of 140.36 aa². Twelve of the 17 pairs showed higher variance in the soluble family. A family-level sign-flip test supported greater positional variance in the soluble controls (P = 0.0102), whereas the exact one-sided sign test provided weaker evidence (P = 0.0717).

For TMD-end-relative positions, the same 17 matched family pairs were informative. The paired membrane-family median variance was 48.38 aa², compared with 75.00 aa² in the matched soluble families. Ten of the 17 pairs showed higher variance in the soluble control. Statistical support was weaker than for TMD starts, with a family-level sign-flip P value of 0.0738 and an exact one-sided sign-test P value of 0.3145.

Thus, the matched soluble controls revealed an anchor-dependent pattern. The TMD-start comparison was consistent with stronger positional constraint in membrane-protein families, whereas the TMD-end comparison provided weaker evidence. Accordingly, these analyses do not support a universal membrane-protein-specific effect. Rather, they provide supplementary evidence that is consistent with the stronger positional constraint observed relative to bona fide TMD boundaries in the primary analysis and the homologous non-TMD control.

# Discussion

The comparative analysis identifies a reproducible evolutionary relationship between local synonymous codon organization and transmembrane-domain architecture in Enterobacterales. Low-adaptation synonymous codon segments showed substantially lower cross-species positional variance relative to homologous TMD boundaries than expected under within-protein synonymous-codon permutation. This pattern persisted after exclusion of first TMDs, after conditioning on the observed cross-species presence of low-adaptation segments, and after controlling for family-level dominance. Moreover, homologous non-TMD sequence positions within the same membrane-protein families did not reproduce the low positional variance observed relative to bona fide TMDs. Together, these results indicate that the location of low-adaptation synonymous codon segments is evolutionarily constrained with respect to membrane-protein topology and is not readily explained by synonymous composition, recurrent feature presence, or generic positional conservation within homologous proteins.

The use of homologous TMD boundaries as internal positional landmarks extends previous work on evolutionarily conserved synonymous codon patterns. Conserved rare-codon clusters have been identified at corresponding positions across homologous coding sequences, suggesting that evolutionary constraint can act on the local positioning of synonymous codon features rather than solely on specific codon identities (Chartier et al., 2012; Chaney et al., 2017). This analysis addresses a more specific form of positional conservation by asking whether such features retain their positions relative to recurrent structural and biogenetic landmarks. Because membrane-protein synthesis is coupled to targeting, translocon engagement, topology establishment, and membrane insertion, homologous TMD boundaries provide biologically meaningful positional anchors for this comparison (Mercier et al., 2022).

Several mechanisms could, in principle, generate selection on synonymous codon positioning around TMDs. Synonymous codon usage can influence translation elongation and cotranslational protein folding, and both experimental and evolutionary analyses have implicated local codon-usage patterns in the timing of nascent-chain folding events (Jacobs and Shakhnovich, 2017; Liu et al., 2021; Moss et al., 2024). One possible interpretation of the present pattern is therefore that local translation dynamics are evolutionarily constrained relative to particular stages of TMD emergence or insertion. However, the present analysis does not directly measure elongation rate, ribosome dwell time, targeting, or membrane insertion. Low adaptation should therefore not be equated with a translational pause. The observed positional constraint could also reflect other sequence-level factors associated with TMD architecture, including local mRNA structure, translational accuracy, nucleotide composition, or regulatory sequence constraints.

The matched soluble-family analysis further indicates that the strength of the membrane-associated signal depends on the positional anchor considered. At TMD starts, membrane families showed lower family-level positional variance than their matched soluble controls, with support from the family-level sign-flip analysis, although the exact sign test provided only suggestive evidence. In contrast, the corresponding TMD-end comparison showed weaker statistical support. This asymmetry argues against a simple model in which membrane proteins are uniformly more positionally constrained than soluble proteins at all comparable sequence positions. Instead, it is consistent with the possibility that positional constraints differ among stages of TMD emergence or insertion. At the same time, the modest number of informative matched family pairs and the weaker distribution-free evidence require caution. The start-versus-end pattern should therefore be regarded as hypothesis-generating rather than as evidence for a defined mechanistic asymmetry.

The present study has several limitations. First, the species-specific relative codon adaptation metric captures codon-usage preference but is not a direct measurement of translation speed. Relationships between synonymous codon usage, elongation kinetics, and protein folding are context dependent (Liu et al., 2021; Moss et al., 2024). Second, the analysis was restricted to Enterobacterales and to strict one-to-one orthologous families with sufficiently conserved multi-pass membrane topology. The extent to which the same positional constraint applies across more deeply diverged bacterial lineages remains unknown. Third, although the synonymous, fixed-presence, family-level, homologous non-TMD, and matched-soluble controls address several alternative explanations, other sequence-level covariates were not modeled explicitly. Finally, the comparative and observational nature of the analysis precludes causal inference about the biological effects of individual low-adaptation segments.

Future experimental studies could test whether the conserved segments identified here influence membrane-protein biogenesis. A direct test would be to synonymously recode selected conserved low-adaptation segments while preserving the amino-acid sequence and then measure effects on protein abundance, membrane insertion, topology, folding, or function. Ribosome profiling or other direct measurements of local elongation behavior could test whether these conserved segments correspond to reproducible changes in translation dynamics. Experimental synonymous recoding has been shown to perturb cotranslational folding and cellular fitness in a protein- and context-dependent manner, supporting the feasibility of direct functional tests of conserved low-adaptation segments (Walsh et al., 2020; Liu et al., 2021; Moss et al., 2024). Extending the comparative framework to broader phylogenetic groups and integrating codon organization with experimentally measured translation and membrane-insertion phenotypes should help determine whether TMD-relative synonymous organization represents a broader feature of membrane-protein evolution.

In summary, the results support evolutionary constraint on the positioning of low-adaptation synonymous codon segments relative to homologous TMD boundaries in Enterobacterales. The signal is robust to multiple composition- and position-based controls, whereas comparison with matched soluble proteins reveals a more nuanced, anchor-dependent pattern. These findings identify TMD-relative synonymous codon positioning as a testable feature of membrane-protein coding-sequence evolution and provide a framework for future mechanistic studies of its relationship to cotranslational membrane-protein biogenesis.

# Methods

## Genome selection and sequence resources

Ten representative Enterobacterales reference genomes were selected for comparative analysis. Protein-coding sequences, translated protein sequences, and genome annotations were obtained from NCBI RefSeq assemblies (O'Leary et al., 2016). The species panel comprised Escherichia coli, Shigella flexneri, Salmonella enterica, Yersinia pestis, Citrobacter koseri, Enterobacter cloacae, Proteus mirabilis, Klebsiella pneumoniae, Pantoea ananatis, and Serratia marcescens.

Chromosomal protein sequences were compiled into a combined protein FASTA file for orthology analysis. Coding sequences were retrieved from the corresponding RefSeq `cds_from_genomic.fna` files. For the Salmonella assembly, coding sequences that were not available in the same format were reconstructed from the genomic annotation before downstream codon-level analysis.

All CDS–protein pairs used for codon analysis were validated by translation using the standard bacterial genetic code. A CDS was considered usable when its translated sequence either exactly matched the annotated protein or differed only at the initiator residue because of an alternative bacterial start codon. Sequences containing internal stop codons, length mismatches, or other sequence mismatches were excluded.

## Ortholog identification and filtering

Orthologous relationships were inferred using MMseqs2 reciprocal best hits (Steinegger and Söding, 2017). Pairwise searches were conducted across the 10 proteomes using a minimum amino-acid sequence identity of 30% and minimum query and target alignment coverage of 70%.

Reciprocal-best-hit relationships were represented as an undirected graph, and connected components were used to define candidate orthologous groups. To minimize paralog ambiguity, only strict one-to-one groups containing no more than one protein from each species were retained. Orthologous groups were required to be represented in at least eight of the 10 species.

This procedure identified 1,836 strict one-to-one orthologous groups. Of these, 417 contained proteins from eight species, 994 from nine species, and 425 from all 10 species.

## Membrane-protein topology annotation and quality control

Reviewed UniProtKB/Swiss-Prot topology annotations were used as the primary source for identifying membrane-protein families and transmembrane-domain boundaries (The UniProt Consortium, 2025). Orthologous groups containing at least one reviewed membrane-protein annotation with three or more transmembrane domains (TMDs) were selected for further topology assessment.

Because reviewed topology annotations were not available for every member of each orthologous group, DeepTMHMM was used as a secondary topology-completion and consistency-assessment tool (Hallgren et al., 2022). DeepTMHMM predictions were generated for all proteins belonging to candidate membrane-protein families.

Protein sequences within each candidate orthologous group were aligned using MAFFT (Katoh and Standley, 2013). TMD architecture was considered consistent when the number and order of homologous transmembrane domains could be unambiguously matched across species. Families with inconsistent TMD architecture, ambiguous TMD correspondence, beta-barrel topology, or fewer than three homologous TMDs were excluded.

After topology quality control, 225 membrane-protein orthologous groups were retained for the primary analysis.

For each retained family, homologous TMDs were indexed according to their conserved order along the protein sequence. TMD start and end positions were used as the primary positional anchors, whereas TMD centers were analyzed as a supplementary anchor.

## Species-specific synonymous codon adaptation metric

Synonymous codon adaptation was quantified separately for each species to account for differences in codon usage among genomes.

For each amino acid, the relative weight of a synonymous codon c was defined as

w(c) = f(c) / max[f(c')]

where f(c) is the observed frequency of codon c among coding sequences from that species and the denominator is the frequency of the most frequently used synonymous codon encoding the same amino acid.

Stop codons were excluded from the calculation. Codon weights were derived independently for each of the 10 genomes using the corresponding complete chromosomal coding-sequence set.

For every codon position in each protein included in the analysis, the corresponding species-specific relative codon weight was assigned. The annotated initiation codon was excluded from positional analyses.

Within each species, low-adaptation codons were defined using the frozen bottom-decile threshold derived from the membrane-protein analysis set. Ties at the threshold were retained. These species-specific thresholds were fixed before analysis of TMD-relative positional conservation and were subsequently applied unchanged to matched soluble-protein controls.

## Definition of low-adaptation segments

Low-adaptation codons were identified independently within each species using the frozen species-specific thresholds described above. A low-adaptation segment was defined as a run of at least three consecutive codons whose relative codon weights were at or below the corresponding species-specific threshold.

Segments were called directly from consecutive codon positions after excluding the initiation codon and stop codons from low-adaptation classification. Runs shorter than three codons were discarded. For each retained segment, the start and end codon positions and the segment center were recorded.

Segment assignment was performed independently for each anchor type. For the TMD-start, TMD-end, and TMD-center analyses, each low-adaptation segment was assigned to the TMD whose corresponding start, end, or center position, respectively, had the smallest absolute distance from the segment center. If multiple segments from the same protein were assigned to the same TMD, only the segment with the smallest absolute anchor-relative distance was retained.

Anchor-relative position was defined as the segment-center position minus the corresponding TMD-anchor position, such that negative values indicated positions upstream of the anchor and positive values indicated positions downstream.

## Cross-species conservation of TMD-relative segment positioning

Homologous TMDs were indexed according to their conserved order within each topology-qualified orthologous group. Only proteins whose TMD number matched the modal TMD count of the corresponding orthogroup were included in the primary positional analysis.

For each anchor type, low-adaptation segments were assigned independently to the nearest homologous TMD anchor as described above. For each orthogroup–TMD unit, anchor-relative segment positions were then collected across species.

A unit was eligible for the primary analysis when the homologous TMD was represented in at least eight species. Among these units, positional conservation was evaluated only when at least three species contained an assigned low-adaptation segment.

Cross-species positional conservation was quantified using the sample variance of anchor-relative segment positions across species. Lower variance therefore indicated stronger conservation of segment position relative to the homologous TMD anchor.

TMD starts and ends were treated as the primary anchors, whereas TMD centers were analyzed separately as a supplementary anchor using the same anchor-specific assignment procedure. The primary summary statistic for each analysis was the median variance across all qualifying orthogroup–TMD units.

## Null models and robustness analyses

### Within-protein synonymous-codon permutation

To determine whether the observed positional conservation could arise from the synonymous codon composition of individual protein-coding sequences, I performed within-protein synonymous-codon permutations.

For each protein, codon positions were grouped according to the amino acid they encoded. Within each amino-acid class, the observed synonymous codons were randomly permuted among positions encoding that same amino acid. The initiation codon and stop codons were excluded from permutation. This procedure preserved the amino-acid sequence exactly and retained the synonymous-codon composition of each individual coding sequence, while randomizing the positions of synonymous codons along the protein.

Following each permutation, species-specific relative codon weights and the frozen low-adaptation thresholds were reapplied to the permuted coding sequences. Low-adaptation segments were recalled de novo using the same requirement of at least three consecutive low-adaptation codons as in the observed dataset.

Segment assignment was performed independently for each anchor type. In the TMD-start analysis, segment centers were assigned to the nearest TMD start; in the TMD-end analysis, they were assigned to the nearest TMD end. If multiple segments from the same protein were assigned to the same TMD, only the segment with the smallest absolute anchor-relative distance was retained.

For each homologous orthogroup–TMD unit represented in at least eight species, units containing an assigned segment in at least three species were retained. Cross-species positional variance was calculated using the sample variance of anchor-relative segment positions, and the median variance across all qualifying units was used as the summary statistic for each permutation.

Separate analyses were performed for TMD starts and TMD ends, each using 1,000 permutations and random seed 20260825. The empirical one-sided P value was calculated as

P = (k + 1) / (N + 1),

where k was the number of permutations with a median variance less than or equal to the observed median variance and N = 1,000.

### Sensitivity analysis excluding first TMDs

To assess whether the primary signal was driven by N-terminal or first-TMD-associated effects, the conservation analysis was repeated after excluding all TMDs with index 1.

The same inclusion criteria and within-protein synonymous-codon permutation procedure were retained. For each permutation, qualifying non-first-TMD units were identified and their median cross-species variance was recalculated.

### Fixed-presence positional null

To distinguish conservation of feature position from conservation of feature occurrence, I constructed a conditional fixed-presence positional null based on the TMD-center analysis.

First, 1,000 within-protein synonymous-codon permutations were generated using the same procedure as in the primary synonymous null. For each shuffled coding sequence, low-adaptation segments were recalled de novo using the frozen species-specific codon weights, low-adaptation thresholds, and minimum segment-length criterion. Segment centers were assigned to the nearest TMD center within the same protein, and when multiple shuffled segments were assigned to the same TMD, only the segment with the smallest absolute TMD-center-relative distance was retained.

For each orthogroup–TMD unit that qualified in the observed dataset, the exact set of species containing an observed low-adaptation segment was recorded. Null positional distributions were then collected separately for each observed-positive orthogroup–TMD–species combination, but only from synonymous-shuffled replicates in which a low-adaptation segment was assigned to that same TMD in that same species. Thus, the null distributions were explicitly conditioned on feature presence.

To ensure adequate sampling of the conditional positional distributions, an orthogroup–TMD unit was retained only when every observed-positive species in that unit had at least 20 conditional null positions available. This filtering yielded 80 usable observed TMD units.

A second-stage conditional resampling procedure was then performed for 1,000 replicates. In each replicate, one null TMD-center-relative position was sampled independently from the species-specific conditional pool for every observed-positive species in each usable unit. The sample variance across species was calculated for each unit, and the median variance across all usable units was used as the replicate-level summary statistic.

The observed statistic was recalculated using exactly the same set of 80 usable units. The empirical one-sided P value was calculated as

P = (k + 1) / (N + 1),

where k was the number of conditional-resampling replicates with a median variance less than or equal to the observed restricted median and N = 1,000. Random seed 20260825 was used throughout.

### Family-dominance and leave-one-orthogroup-out analyses

To evaluate whether the global signal was driven disproportionately by a small number of protein families, the number of qualifying TMD units contributed by each orthogroup was tabulated.

In addition, a leave-one-orthogroup-out analysis was performed for the TMD-start analysis. Each orthogroup was removed in turn, and the global median variance was recalculated from the remaining qualifying units. Stability of the median across these leave-one-orthogroup-out datasets was used to assess robustness to individual orthogroups.

### Homologous non-TMD pseudo-anchor control

To test whether similarly low positional variance could arise around arbitrary homologous sequence positions rather than bona fide transmembrane domains, I constructed a symmetric homologous non-TMD pseudo-anchor null within the same membrane-protein orthogroups.

Protein sequences within each orthogroup were aligned with MAFFT (Katoh and Standley, 2013). Candidate pseudo-anchor columns were defined directly from these multiple-sequence alignments. A candidate alignment column was required to contain a residue in at least 80% of eligible family members. In addition, at least 80% of eligible members had to map that column to a valid non-TMD residue located at least 15 amino acids from the nearest TMD boundary. Only orthogroups with at least 10 eligible candidate columns were retained.

For each retained orthogroup, the number of pseudo-anchors sampled in each permutation was set equal to the modal number of TMDs in that family. Candidate alignment columns were sampled without replacement and then ordered according to their alignment position, generating an ordinal series of pseudo-anchors analogous to the ordered TMD series.

For each family member, sampled alignment columns were mapped back to residue coordinates. Each observed low-adaptation segment was assigned to the nearest available pseudo-anchor within the same protein according to the absolute distance between the segment center and pseudo-anchor residue position. If multiple segments from the same protein were assigned to the same pseudo-anchor, only the segment with the smallest absolute distance was retained.

For each pseudo-anchor unit, the number of species in which the anchor residue was present was used as the species denominator. Units were retained when the pseudo-anchor was represented in at least eight species and at least three species contained an assigned low-adaptation segment. Cross-species positional variance was calculated using the sample variance of segment-center distances relative to the pseudo-anchor.

The observed TMD-based comparison statistic was restricted to the same set of orthogroups that passed the non-TMD candidate-pool criteria. A total of 1,000 pseudo-anchor permutations were performed using random seed 20260826. For each permutation, the median positional variance across all qualifying pseudo-anchor units was calculated.

The empirical one-sided P value was calculated as

P = (k + 1) / (N + 1),

where k was the number of pseudo-anchor permutations with a median variance less than or equal to the observed restricted TMD-associated median and N = 1,000.

## Matched soluble-protein controls

### Matched soluble-family selection

To provide an external control independent of the membrane-protein families themselves, each membrane orthogroup was matched to a soluble-protein orthogroup.

Candidate soluble orthogroups were drawn from the same strict one-to-one orthology dataset as the membrane-protein families. For each membrane orthogroup, soluble candidates were required to contain exactly the same set of species and to have a median protein length within 5% of the membrane-family median length.

Candidate soluble families were screened with DeepTMHMM, and only orthogroups for which all analyzed members were classified as globular were retained. Among the eligible membrane–soluble candidate pairs, one-to-one matching was optimized using a maximum-cardinality minimum-cost bipartite matching procedure. Relative median protein-length difference was used as the primary matching cost, whereas candidate rank was used only as a deterministic tie-breaker. Each membrane family and each soluble family could therefore contribute to at most one final pair.

This procedure produced 197 unique one-to-one matched membrane–soluble orthogroup pairs for downstream analysis.

### Soluble pseudo-anchor construction

For each matched membrane–soluble orthogroup pair, homologous membrane TMD boundaries were converted to relative sequence coordinates using

r = (p - 1) / (L - 1),

where p denotes the residue position of the TMD boundary and L the protein length. Only membrane proteins matching the modal TMD count of the corresponding orthogroup were used. For each homologous TMD index, the median relative position across membrane-family members was taken as the target position, separately for TMD starts and TMD ends.

Protein sequences from the matched soluble family were aligned using MAFFT (Katoh and Standley, 2013). Candidate pseudo-anchor positions were restricted to alignment columns occupied by a residue in 100% of soluble-family members. For each fully occupied column, the corresponding residue position was converted to a relative sequence coordinate in each soluble protein, and the median relative coordinate across family members was calculated.

Membrane TMD target positions were then mapped to soluble alignment columns using an ordered minimum-cost dynamic-programming procedure. Each TMD target was assigned to a unique soluble alignment column, target order along the protein was preserved, and the total absolute difference between membrane target relative positions and soluble candidate relative positions was minimized. TMD-start and TMD-end pseudo-anchors were constructed independently.

For each selected soluble pseudo-anchor column, the corresponding residue position was recorded separately for every soluble-family member. Matching quality was quantified as the absolute difference between the membrane-family target relative position and the median relative position of the selected soluble pseudo-anchor. Pseudo-anchors with an absolute relative-position difference greater than 0.05 were excluded from downstream conservation analysis.

### Soluble low-adaptation segment analysis

The same frozen species-specific codon weights and low-adaptation thresholds derived from the membrane-protein analysis were applied unchanged to the matched soluble dataset. Stop codons were excluded, and the annotated initiation codon was omitted from positional classification using the same rule as in the membrane-protein analysis. No codon weights or low-adaptation thresholds were re-estimated from the soluble control dataset.

Low-adaptation segments were assigned independently to the nearest eligible soluble pseudo-anchor for the start-matched and end-matched analyses. For each membrane orthogroup–TMD index, only the closest assigned low-adaptation segment per species was retained.

Soluble pseudo-anchor units were considered qualifying when the anchor was represented in at least eight species and at least three species contained an assigned low-adaptation segment. Cross-species positional conservation was quantified using the sample variance of pseudo-anchor-relative segment positions.

### Family-level membrane–soluble comparison

For the family-level paired comparison, qualifying unit variances within each membrane family were summarized by their median, and qualifying pseudo-anchor variances within the corresponding matched soluble family were summarized in the same manner.

Only matched membrane–soluble family pairs with qualifying statistics on both sides were retained. Paired differences were defined as soluble-family median variance minus membrane-family median variance.

Directional support for higher positional variance in soluble families was evaluated using a one-sided exact sign test and a family-level sign-flip permutation test with 100,000 permutations. The sign-flip analysis used random seed 20260831. The orthologous family was treated as the clustered statistical unit.

# Figure legends

## Figure 1. Comparative framework for testing TMD-relative conservation of low-adaptation synonymous codon segments

Overview of the comparative analysis workflow. Ten representative Enterobacterales genomes were used to construct 1,836 strict one-to-one orthologous groups represented in at least eight species. Reviewed membrane-topology annotations, supplemented by DeepTMHMM-assisted topology completion and consistency assessment, were used to define 225 topology-qualified multi-pass membrane-protein orthogroups. Species-specific synonymous codon adaptation weights were then calculated, and low-adaptation segments were defined as runs of at least three consecutive codons at or below the frozen species-specific threshold. Segment positions were quantified relative to homologous TMD starts and ends, with TMD centers used as a supplementary anchor. Cross-species positional conservation was evaluated using the variance of anchor-relative segment positions and tested using synonymous-codon permutation and additional positional controls.

## Figure 2. Representative membrane-protein families showing conserved low-adaptation segment positioning relative to homologous TMDs

Representative orthologous membrane-protein families illustrating distinct patterns of conserved low-adaptation segment positioning relative to homologous TMD starts. Each row represents one species within an orthologous family, with TMD positions and low-adaptation segments shown relative to the corresponding homologous TMD start. (A) orthogroup_540 (FocA family), in which low-adaptation segments occur approximately 27 residues upstream of TMD1 in six of nine species. (B) orthogroup_36 (YghB/DedA family), in which low-adaptation segments overlap TMD2 in five of eight species. (C) orthogroup_105 (ClcA family), in which low-adaptation segments cluster around TMD8 in five of nine species. These examples illustrate that conserved positioning can occur upstream of or within TMDs and is not limited to the first TMD.

## Figure 3. Low-adaptation segments show stronger TMD-relative positional conservation than expected under multiple null models

Distributions of permutation-derived median cross-species positional variances are shown together with the corresponding observed statistic. Lower variance indicates stronger conservation of segment position relative to the specified anchor. (A) TMD-start analysis using 1,000 within-protein synonymous-codon permutations. The observed median variance was 22.70 aa², compared with a null median of 154.33 aa² (empirical one-sided P = 9.99 × 10⁻⁴). (B) TMD-end analysis using the same permutation framework. The observed median variance was 34.55 aa², compared with a null median of 147.00 aa² (P = 9.99 × 10⁻⁴). (C) Sensitivity analysis excluding all first TMDs. The observed median variance was 17.21 aa², compared with a null median of 116.81 aa² (P = 9.99 × 10⁻⁴). (D) Symmetric homologous non-TMD pseudo-anchor control within the same membrane-protein families. The observed restricted TMD-associated median variance was 47.29 aa², whereas the median across 1,000 non-TMD pseudo-anchor permutations was 979.29 aa² (P = 9.99 × 10⁻⁴).

## Figure 4. Matched soluble-family controls reveal anchor-dependent differences in positional constraint

Family-level paired comparison between membrane-protein orthogroups and matched soluble-protein controls. Soluble families were matched one-to-one to membrane families using exact species composition, median protein length within 5%, and an all-globular DeepTMHMM classification. For each family, qualifying orthogroup–anchor unit variances were summarized by their median. Lines connect matched membrane and soluble family statistics. (A) TMD-start-relative comparison. Seventeen matched family pairs were informative; the paired membrane median was 28.67 aa² and the paired soluble median was 169.03 aa². Twelve of 17 pairs showed greater variance in the soluble family. The family-level sign-flip test gave P = 0.0102, whereas the exact one-sided sign test gave P = 0.0717. (B) TMD-end-relative comparison. Seventeen matched family pairs were informative; the paired membrane median was 48.38 aa² and the paired soluble median was 75.00 aa². Ten of 17 pairs showed greater variance in the soluble family. Statistical support was weaker than for TMD starts (sign-flip P = 0.0738; exact one-sided sign-test P = 0.3145).

# References

1. Mercier E, Wang X, Bögeholz LAK, Wintermeyer W, Rodnina MV. Cotranslational Biogenesis of Membrane Proteins in Bacteria. Frontiers in Molecular Biosciences. 2022;9:871121. doi:10.3389/fmolb.2022.871121.

2. Liu Y, Yang Q, Zhao F. Synonymous but Not Silent: The Codon Usage Code for Gene Expression and Protein Folding. Annual Review of Biochemistry. 2021;90:375–401. doi:10.1146/annurev-biochem-071320-112701.

3. Moss MJ, Chamness LM, Clark PL. The Effects of Codon Usage on Protein Structure and Folding. Annual Review of Biophysics. 2024;53:87–108. doi:10.1146/annurev-biophys-030722-020555.

4. Chartier M, Gaudreault F, Najmanovich R. Large-scale analysis of conserved rare codon clusters suggests an involvement in co-translational molecular recognition events. Bioinformatics. 2012;28:1438–1445. doi:10.1093/bioinformatics/bts149.

5. Chaney JL, Steele A, Carmichael R, et al. Widespread position-specific conservation of synonymous rare codons within coding sequences. PLoS Computational Biology. 2017;13:e1005531. doi:10.1371/journal.pcbi.1005531.

6. Smalinskaitė L, Hegde RS. The Biogenesis of Multipass Membrane Proteins. Cold Spring Harbor Perspectives in Biology. 2023;15(4):a041251. doi:10.1101/cshperspect.a041251.

7. Walsh IM, Bowman MA, Soto Santarriaga IF, Rodriguez A, Clark PL. Synonymous codon substitutions perturb cotranslational protein folding in vivo and impair cell fitness. Proceedings of the National Academy of Sciences of the United States of America. 2020;117(7):3528–3534. doi:10.1073/pnas.1907126117.

8. Jacobs WM, Shakhnovich EI. Evidence of evolutionary selection for cotranslational folding. Proceedings of the National Academy of Sciences of the United States of America. 2017;114(43):11434–11439. doi:10.1073/pnas.1705772114.

9. Steinegger M, Söding J. MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. Nat Biotechnol. 2017;35:1026–1028. doi:10.1038/nbt.3988.

10. Katoh K, Standley DM. MAFFT multiple sequence alignment software version 7: improvements in performance and usability. Mol Biol Evol. 2013;30(4):772–780. doi:10.1093/molbev/mst010.

11. Hallgren J, Tsirigos KD, Pedersen MD, et al. DeepTMHMM predicts alpha and beta transmembrane proteins using deep neural networks. bioRxiv. 2022. doi:10.1101/2022.04.08.487609.

12. O'Leary NA, Wright MW, Brister JR, et al. Reference sequence (RefSeq) database at NCBI: current status, taxonomic expansion, and functional annotation. Nucleic Acids Res. 2016;44(D1):D733–D745. doi:10.1093/nar/gkv1189.

13. The UniProt Consortium. UniProt: the Universal Protein Knowledgebase in 2025. Nucleic Acids Res. 2025;53(D1):D609–D617. doi:10.1093/nar/gkae1010.

# Data availability

All genome assemblies, coding sequences, protein sequences, and genome annotations used in this study were obtained from publicly available NCBI RefSeq resources. Reviewed membrane-protein topology annotations were obtained from UniProtKB/Swiss-Prot. Protein-family topology completion and consistency assessment were performed using DeepTMHMM.

The study analyzed 10 representative Enterobacterales genomes: Escherichia coli, Shigella flexneri, Salmonella enterica, Yersinia pestis, Citrobacter koseri, Enterobacter cloacae, Proteus mirabilis, Klebsiella pneumoniae, Pantoea ananatis, and Serratia marcescens. Assembly accession identifiers for all genomes are provided in Supplementary Table S1.

Processed analysis tables underlying the main figures and statistical analyses, including ortholog-group assignments, topology-qualified membrane-protein families, low-adaptation codon segments, TMD-relative positional measurements, permutation-null summaries, matched soluble-family controls, and family-level comparison statistics, are publicly available in the associated reproducibility repository and its archived Zenodo release (https://github.com/daming7977-byte/enterobacterales-tmd-codon-positioning; https://doi.org/10.5281/zenodo.22254908).

# Code availability

Custom scripts used for ortholog filtering, topology integration, codon-adaptation calculation, low-adaptation segment detection, TMD-relative positional analysis, permutation tests, homologous non-TMD pseudo-anchor controls, matched soluble-family construction, and family-level statistical analyses are publicly available at https://github.com/daming7977-byte/enterobacterales-tmd-codon-positioning and are archived at Zenodo under DOI: 10.5281/zenodo.22254908.

The repository includes frozen parameter files, software-version information, processed datasets, final statistical outputs, original workflow scripts, and portable core scripts for reproducing the released principal summary figures from the included processed data.

# Funding

This research received no specific external funding.

# Author contributions

Li Ming conceived the study, designed the analyses, performed the computational analyses, interpreted the results, prepared the figures, and wrote the manuscript.

# Acknowledgments

None.

# Conflict of interest

The author declares no competing interests.
