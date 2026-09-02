Matched soluble control specification

1. Matching universe
   - Strict one-to-one Enterobacterales orthogroups.
   - Membrane orthogroups: topology-pass set.
   - Soluble candidates are non-membrane orthogroups from the same strict orthology universe.

2. Matching variables
   - Identical species coverage.
   - Median protein length within ±5% of the membrane orthogroup.
   - The ±5% threshold was selected before inspecting codon-position outcomes.

3. Candidate ranking
   - Candidates are ranked by relative median-length difference.
   - Top 5 candidates per membrane orthogroup were screened by DeepTMHMM.

4. Soluble topology definition
   - A candidate orthogroup is considered strict soluble only if every member is classified as GLOB by DeepTMHMM.
   - TM, BETA, SP, and SP+TM predictions are excluded.

5. Final match
   - Among strict all-GLOB candidates, the highest-ranked candidate is selected as the unique matched soluble control.
   - No codon-position information is used during matching.

6. Matching result
   - 222 membrane orthogroups had candidates under the frozen species-coverage and ±5% length criteria.
   - 220/222 (99.1%) retained at least one strict all-GLOB candidate within the top 5.
   - No further expansion or relaxation of matching criteria is performed.