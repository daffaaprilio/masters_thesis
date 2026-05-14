# List of Analyses to Populate the Results & Discussions Section
Samples Phenotype
| Sample | TAA Production | Callus Formation |
|--------|-----------|-----------------------|
| SBC4 | ++ | Mid |
| SBC10 | +++ | Good |
| SBC11 | - | Mid |
| SBC23 | ++ | Good |

## Data flow for each analysis
```mermaid
flowchart LR
    R1[Raw reads] --> REF1[BTx623.fna\nreference genome]
    R2[Raw reads] --> REF1
    R3[Raw reads] --> REF1
    R4[Raw reads] --> REF1

    REF1 --> B4[SBC4.bam]
    REF1 --> B10[SBC10.bam]
    REF1 --> B11[SBC11.bam]
    REF1 --> B23[SBC23.bam]

    REF2[BTx623.fna\nreference genome] --> V4[SBC4.vcf]
    REF2 --> V10[SBC10.vcf]
    REF2 --> V11[SBC11.vcf]
    REF2 --> V23[SBC23.vcf]

    B4 --> REF2
    B10 --> REF2
    B11 --> REF2
    B23 --> REF2

    V4 --> VL["(1) Variant landscape"]
    V10 --> VL
    V11 --> VL
    V23 --> VL

    V4 --> MRG[merged.vcf]
    V10 --> MRG
    V11 --> MRG
    V23 --> MRG

    MRG --> ANN[annotate.merged.vcf]
    SNPEFF["SnpEFF variant annotation database\n(based on BTx623 reference genome\nwith Ensembl annotation)"] --> ANN

    ANN --> SKG["(4) Sorgoleone key genes exploration"]

    style VL fill:#3b2f6b,stroke:#7c6fcf,color:#fff
    style SKG fill:#3b2f6b,stroke:#7c6fcf,color:#fff
```

## Sequencing & data Quality

## Variant landscape

## Epigenomics landscape

## Identification of TAA biosynthetic pathway genes
`analysis/TAA/`

## Exploration of sorgoleone key genes
`analysis/sorgoleone/`

