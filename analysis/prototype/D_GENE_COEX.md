# D_GENE_COEX.md — Sorghum D-gene Co-expression Positive Control

> Context file for Claude Code. Read this fully before writing any code.
> It defines the scientific objective, the **hard blindness constraints**, the data
> inventory, and the staged pipeline. The blindness constraints are guardrails, not
> suggestions — violating them invalidates the entire result.

---

## 1. Objective

Independently re-derive the sorghum **D-gene** using **co-expression** as the
detection modality, and use that re-derivation as a **blind positive control** that
validates a general-purpose co-expression pipeline.

- **Target gene (the answer we are testing for):** `D` = **Sobic.006g147400**, a NAC
  transcription factor on chromosome 6. Arabidopsis ortholog: **ANAC074**.
- **Original discovery (Fujimoto et al. 2018, PNAS):** D was found by **forward
  genetics + positional cloning** (dry × juicy F2 → chr6 → functional validation via
  ANAC074 knockout / ectopic expression). **It was NOT found by co-expression.**
- **D gene in EGI**: 8084661 

### Why this framing matters
We are **not** replicating Fujimoto's method and inheriting its blindness. We are
asking a *different* question: **would an independent co-expression pipeline have
nominated D, with no knowledge that D is the answer?** A "yes" predicts how the
pipeline behaves on master regulators nobody has cloned yet — but only if no
D-specific knowledge leaked into the pipeline. Sealing those leaks is the whole job.

---

## 2. Hard blindness constraints (GUARDRAILS — do not violate)

These are the leakage channels, ordered by how badly they rig the result. The agent
must refuse to implement any step that breaks one of these, and flag it to the user.

1. **Genetic-contrast leak (worst).** Do NOT build the network from a contrast that
   isolates D: no dry-stem vs `d-NIL`, no biparental population segregating at D as
   the network basis. Those make D the trivial top axis of variation — you'd be
   reading out the construct, not discovering anything. Anchor instead on the
   **progression of pith PCD within functional-D material** (stem internode
   developmental series, or pith-vs-rind), where co-expression reflects co-regulation.
2. **Locus leak.** Do NOT restrict candidates to the chr6 QTL interval. That interval
   *is* the published answer.
3. **Family leak.** Do NOT pre-filter candidates to the NAC family. Restricting to the
   ~1800-gene **TF universe is allowed** ("the master switch is a TF" is a generic
   prior). Restricting to NACs is not — recovering a NAC must be an *output*.
4. **Bait leak.**
   - Seeding with **D** itself = circular. Forbidden.
   - Seeding with **SbXCP1** = soft leak. Casto et al. 2018 published XCP1 as
     D-regulated, so a skeptic can dismiss it. Keep XCP1 **out** of the bait set. If
     ever included, disclose and report results with and without it.
   - Allowed bait: **SEN102 + a PCD-execution signature defined by orthology**, blind
     to the D papers (see §4).
5. **Parameter-tuning leak (the cardinal sin).** Do NOT adjust thresholds,
   module-merge cutoffs, MR decay, or the ranking metric until D appears. All such
   choices are frozen in the pre-registration (§7) **before** unblinding.

**Definition of "unblinding":** the first moment any code queries, prints, ranks, or
filters on `Sobic.006g147400` / chr6 interval / NAC family. Everything before that
must be hypothesis-agnostic. Keep unblinding in a single, clearly-marked, late stage.

---

## 3. Data inventory (what the user has)

| Asset | Description | Status | Remark |
|---|---|---|---|
| TF universe | ~1800 sorghum transcription factors | obtained | Located in `/Users/daffa/workspace/infobio/thesis/resources/TFDB` |
| Bait seed | **SEN102** — sorghum ortholog of Arabidopsis AtCEP1/2/3 (KDEL-tailed papain-like C1A cysteine protease) | obtained | EGI is 8058001 |
| Co-expression resource | **ATTED-II** (sorghum = `sbi`), Mutual-Rank based, *condition-independent*. v13 (Feb 2026). | available | Located in `/Users/daffa/workspace/infobio/Sbi-r.c1-0/Sbi-r.v25-12.G21627-S807.combat_pca.subagging.z.d` |
| Reference genome | **NCBI v3** sorghum (per project standard) | in use | Located in `/Users/daffa/workspace/infobio/thesis/resources/ref` |

---

## 4. About the anchor genes Bait / PCD-execution signature

Currently there are two families (or groups?) of bait genes:
1.  **TF** group
    Containing about 2000 genes (1800 after converted to EGI).
2.  **TF-regulated genes** group
    Currently, only includes thiol protease SEN102. I plan to expand the list after this iteration of attempt.

---

## 5. Current iteration

I want you to recommend me an analysis method to:
1.  Find which TFs among the 2000 TFs that are having co-expression value with the SEN102 gene, that is the sole PCD-related gene in this attempt.
2.  Report the analysis method in this markdown file (continue writing ## 6.)

---

## 6. Analysis method

### 6.0 What the data already gives us
ATTED-II's `…subagging.z.d` resource is a **precomputed, condition-independent
co-expression neighborhood per gene**: one file per gene, named by its EGI id. Each
file is a two-column, descending-sorted table:

```
<partner_EGI>\t<coex_score>
```

where `coex_score` is the **subagging z-transformed Mutual-Rank index** (higher =
tighter co-expression; the value is already symmetric/rank-based, so it is comparable
across baits). Concretely:

- **SEN102 = EGI `8058001`** → file `…subagging.z.d/8058001` already contains SEN102's
  genome-wide ranked co-expression partners.
- **TF universe** → `resources/TFDB/Sbi_TF_gene_ids_EGI.txt` (~1800 EGI ids).

Because the bait's neighborhood is precomputed, this iteration is **not** a network
rebuild (WGCNA/MR-from-scratch). It is a **single-bait guilt-by-association ranking**:
read SEN102's neighbor file, keep the rows whose partner is a TF, and rank them. This
is the cheapest method that still answers §5.1, and — critically — it touches no
D-specific channel (§2).

### 6.1 Method: bait-restricted reciprocal-rank ranking of the TF universe

**Inputs (frozen):**
- Bait set `B = { SEN102 (8058001) }` — single PCD-execution gene. XCP1 excluded (§2.4).
- Candidate set `C = ` the ~1800 EGI TFs in `Sbi_TF_gene_ids_EGI.txt`. No NAC/chr6/locus
  sub-filter (§2.2–2.3).
- Resource: `Sbi-r.v25-12.G21627-S807.combat_pca.subagging.z.d`.

**Procedure:**
1. **Forward score.** Load `…z.d/8058001`. For every TF `t ∈ C`, record
   `s_fwd(t)` = its coex_score in SEN102's file, and `r_fwd(t)` = its **1-based rank**
   within SEN102's *full* neighborhood (genome-wide, before restricting to TFs). TFs
   absent from the file get the worst rank / `s = NA`.
2. **Reciprocal score.** For each TF `t`, open `…z.d/<t>` and record SEN102's rank in
   *that* file → `r_rev(t)`. This guards against asymmetric/one-sided hits even though
   MR is nominally symmetric.
3. **Combined statistic (frozen ranking metric).** Mutual-Rank-style geometric mean of
   the two ranks:  `MR(t) = sqrt( r_fwd(t) · r_rev(t) )`. Lower `MR` = stronger,
   reciprocated co-expression. Report `s_fwd` alongside for interpretability. `MR` —
   not `s_fwd` — is the primary sort key, fixed before unblinding.
4. **Empirical significance (frozen, hypothesis-agnostic null).** SEN102's score is
   ranked against the *whole genome*, so calibrate per-TF significance against the
   genome-wide background of SEN102's own neighborhood: for each TF, its percentile =
   `r_fwd(t) / N_total`. Additionally build a **degree-matched permutation null** —
   draw 1000 random gene sets of size |C| from all genes (not just TFs) and recompute
   the rank distribution — to get an empirical FDR on "how enriched are TFs near SEN102
   vs. random genes of equal count." This null is defined without reference to D.
5. **Output (the pre-registered deliverable).** A single ranked table
   `tf_coex_with_SEN102.tsv` with columns:
   `EGI, Sobic_id, TF_family, s_fwd, r_fwd, r_rev, MR, percentile, emp_FDR`,
   sorted ascending by `MR`. Plus a fixed **call rule**: TFs with `emp_FDR < 0.05`
   (or, if that is empty, the **top-20 by MR**) are the "nominated" set.

**Tie/edge handling (frozen):** missing partner → rank = `N_total + 1`; ties in
coex_score → average rank; a TF that is its own bait is N/A (SEN102 is a protease, not
in C, so this does not arise).

### 6.2 Where the blindness wall sits
Steps 1–5 are **fully hypothesis-agnostic**: nothing queries `8084661` (D), the chr6
interval, or the NAC family. The `TF_family` column is populated from the generic TFDB
annotation for *all* TFs at once (not a NAC filter), so reading it does not unblind.

**Unblinding (separate, late, single step — §7 pre-registration must be signed first):**
only after `tf_coex_with_SEN102.tsv` is written and frozen do we ask the three sealed
questions, in this order, and report all three regardless of outcome:
1. What is D's (`8084661`) rank and `MR` in the frozen table?
2. Is the nominated set (emp_FDR<0.05 / top-20) enriched for NAC-family TFs?
3. Does D fall inside the nominated set under the *pre-committed* call rule — with **no**
   post-hoc threshold change (§2.5)?

### 6.3 Honest limitations of this single-bait iteration
- **One bait = low specificity.** A lone protease neighbor cannot distinguish D from any
  other PCD-co-regulated TF; expect a broad nominated set. This is acceptable for a
  *positive control* (we only ask "would D have been nominated?"), but the user's planned
  expansion of the TF-regulated bait set (§4.2: add SEN102-like C1A proteases, vacuolar
  processing enzymes, other orthology-defined PCD-execution genes) is what will sharpen
  specificity in the next iteration. When that happens, aggregate baits by **mean rank
  across baits** (keep the same `MR`-of-mean-rank machinery) rather than re-tuning.
- **Condition-independent resource.** ATTED-II MR averages over all conditions, so a
  pith-PCD-specific signal may be diluted. A Tier-2 confirmation on a pith/internode
  developmental-series GCN (MOROKOSHI / PlantNexus, §A) is the natural follow-up, but
  must use the *same* frozen metric and call rule to stay blind.

### 6.4 Concrete implementation sketch (read-only, no unblinding)
```python
# pseudo — operates only on 8058001 + the TF EGI list; never opens 8084661
import os, numpy as np, pandas as pd

ZD   = "Sbi-r.v25-12.G21627-S807.combat_pca.subagging.z.d"
TFS  = set(open(".../TFDB/Sbi_TF_gene_ids_EGI.txt").read().split())

def neighbors(egi):                      # -> {partner: (score, rank)}
    out = {}
    for i, line in enumerate(open(os.path.join(ZD, egi)), start=1):
        g, s = line.split()
        out[g] = (float(s), i)
    return out, i                        # i = N_total in this file

fwd, N = neighbors("8058001")            # SEN102 forward neighborhood
rows = []
for t in TFS:
    s_fwd, r_fwd = fwd.get(t, (np.nan, N + 1))
    try:
        rev, Nrev = neighbors(t)
        r_rev = rev.get("8058001", (np.nan, Nrev + 1))[1]
    except FileNotFoundError:
        r_rev = N + 1
    rows.append((t, s_fwd, r_fwd, r_rev, np.sqrt(r_fwd * r_rev), r_fwd / N))

df = (pd.DataFrame(rows,
        columns=["EGI","s_fwd","r_fwd","r_rev","MR","percentile"])
        .sort_values("MR"))
df.to_csv("tf_coex_with_SEN102.tsv", sep="\t", index=False)   # FROZEN deliverable
# emp_FDR + Sobic_id/TF_family annotation joined here, then STOP. Unblind only in §6.2.
```

---

## A. Key references

- **Fujimoto et al. 2018**, *PNAS* 115(37):E8783–E8792 — D / ANAC074, PCD of stem pith
  parenchyma. DOI: 10.1073/pnas.1807501115. *(target gene; positive control source)*
- **Casto et al. 2018**, *Plant Direct* — SbNAC_D = Sobic.006g147400; SbXCP1 link
  *(reason XCP1 is excluded from bait)*.
- **ATTED-II** (Obayashi et al.; 2018 MR-index paper; v13 2026) — `sbi`, Mutual Rank,
  condition-independent co-expression. https://atted.jp
- **MOROKOSHI** (Makita et al. 2015) — RIKEN sorghum transcriptome + co-expression.
- **PlantNexus** — tissue-specific GCNs for sorghum & barley *(Tier 2 / reproducibility precedent)*.

---

## B. Glossary

- **D / Dry gene** — Sobic.006g147400, NAC TF, induces PCD of stem pith parenchyma;
  functional allele → dry/pithy stems, non-functional → juicy stems.
- **SEN102** — sorghum AtCEP1/2/3 ortholog; KDEL-tailed papain-like (C1A) cysteine
  protease; PCD-execution bait.
- **MR (Mutual Rank)** — ATTED-II's symmetric, rank-based co-expression index.
- **kME** — intramodular connectivity / module membership (WGCNA), used to rank hubs.
- **Pith parenchyma** — central stem storage tissue whose programmed death drives the
  juicy→dry transition.