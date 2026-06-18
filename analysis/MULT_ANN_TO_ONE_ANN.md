# Strategy to reduce the complexity of the data
Through simplification of SnpEff ANN's annotation field.

## Current problem
- SnpEff ANN works by annotating a single variant to every transcript that are affected by that variation.
- However, SIFT4G only annotates the impact of every variation within *one* gene, instead of *transcript*.
- To make things simple, it is advised to reduce (aka eliminates) all ANN fields and retain only *one*, that corresponds to a selected transcript.

## Methodological approach
### Case 1 — SIFT annotation exists (apple-to-apple):
Match the ANN entry whose (allele, transcript_id) equals the SIFTINFO entry's key. SIFT already anchors to one specific transcript, so this pairs the SnpEff effect and SIFT score on the same transcript. No tiebreaking needed.
### Case 2 — SIFT absent (gene-level umbrella):
Collapse multiple ANN entries to one per gene using a single deterministic multi-key sort, then take the first row. Sort keys in order:

- Impact severity — HIGH > MODERATE > LOW > MODIFIER (most important; this is what scoring consumes)
- SnpEff's own ANN order — severity proxy within the same impact class (capture before any filtering, while order is intact)
- Biotype — prefer protein_coding over pseudogene/nc
