## =========================================================
## PARAMETERS: Virtanen et al. 2018 Human BAT RNA-seq
## =========================================================
## Source: Cell Metabolism 27, 1-9 (2018)
## DOI: 10.1016/j.cmet.2018.05.020
## GEO: GSE113764
## =========================================================

## ---- Study design ----------------------------------------
SPECIES            <- "human"
GEO_ACCESSION      <- "GSE113764"
DATA_TYPE           <- "bulkRNAseq_paired_BAT_WAT"
STUDY_TITLE         <- "Virtanen_Human_BAT"
STUDY_DESCRIPTION   <- paste(
  "Paired supraclavicular BAT and subcutaneous WAT biopsies",
  "from 14 healthy adult subjects (10F/4M, age 21-52).",
  "Illumina HiSeq 2500, TruSeq RNA Access."
)

## ---- Sample information ----------------------------------
N_SUBJECTS         <- 14
N_SAMPLES_EXPECTED <- 28          # 14 BAT + 14 WAT (paired)
CONDITION_COL      <- "tissue"    # column name in colData
CONDITION_REF      <- "WAT"       # reference level (denominator)
CONDITION_TEST     <- "BAT"       # test level (numerator)
SUBJECT_COL        <- "subject"   # for paired design

## ---- DE thresholds (from paper) --------------------------
PADJ_CUTOFF        <- 0.05
LFC_CUTOFF         <- 1.0         # log2 fold-change >= 1

## ---- Expected DEG counts (paper Table S1 / Figure 1) ----
EXPECTED_UP_BAT    <- 847         # genes higher in BAT
EXPECTED_UP_WAT    <- 928         # genes higher in WAT
DEG_COUNT_TOLERANCE <- 0.20       # allow 20% deviation

## ---- BAT marker genes (from paper) ----------------------
BAT_MARKER_GENES <- c(
  "UCP1",                         # thermogenin (canonical BAT marker)
  "CKMT1A", "CKMT1B", "CKMT2",   # creatine kinase mitochondrial

  "ACTC1",                        # alpha-cardiac actin
  "MYH7",                         # myosin heavy chain 7
  "ACTA1",                        # alpha-skeletal actin
  "SLC25A20",                     # carnitine transporter
  "KCNK3",                        # potassium channel TASK-1
  "PM20D1"                        # N-acyl amino acid hydrolase
)

## ---- WAT marker genes (from paper) ----------------------
WAT_MARKER_GENES <- c(
  "HOXC8", "HOXC9",              # HOX cluster (WAT identity)
  "HOTAIR",                       # HOX antisense intergenic RNA
  "LEP",                          # leptin
  "ADH1B",                        # alcohol dehydrogenase 1B
  "COL6A3",                       # collagen VI alpha 3
  "FRZB",                         # frizzled-related protein
  "PLIN1"                         # perilipin 1
)

## ---- Lipid metabolism genes of interest ------------------
LIPID_METABOLISM_GENES <- c(
  "CPT1B", "CPT2", "ACADL", "ACADM", "ACADVL",  # fatty acid oxidation
  "SLC25A20", "HADHA", "HADHB",                   # beta-oxidation
  "FABP3", "FABP4",                                # fatty acid binding
  "LIPE", "PNPLA2", "MGLL",                       # lipolysis
  "DGAT1", "DGAT2",                                # triglyceride synthesis
  "UCP1", "UCP2", "UCP3",                          # uncoupling proteins
  "PPARGC1A", "PPARGC1B", "PPARG"                  # transcription factors
)

## ---- Plot dimensions -------------------------------------
VOLCANO_W          <- 8
VOLCANO_H          <- 6
HEATMAP_W          <- 2000       # pixels
HEATMAP_H          <- 3000
HEATMAP_RES        <- 200
PCA_W              <- 8
PCA_H              <- 6
PLOT_DPI           <- 300

## ---- Heatmap parameters ----------------------------------
HEATMAP_TOP_N      <- 50         # top N DEGs for heatmap
HEATMAP_SCALE      <- "row"      # z-score scaling

## ---- QC thresholds ---------------------------------------
MIN_COUNTS_PER_GENE  <- 10       # minimum total counts to keep gene
MIN_SAMPLES_DETECTED <- 3        # gene must be detected in >= N samples
LOW_COUNT_FILTER     <- TRUE     # apply low-count filtering

## ---- Annotation database ---------------------------------
ORG_DB             <- "org.Hs.eg.db"
GENOME             <- "hg38"
GENE_ID_TYPE       <- "ENSEMBL"  # GEO counts use ENSEMBL IDs
KEGG_ORGANISM      <- "hsa"

## ---- Reproducibility -------------------------------------
SEED               <- 42
GENERATE_PLOTS     <- TRUE
VERBOSE            <- TRUE

cat("[OK] parameters.R loaded for", STUDY_TITLE, "\n")
