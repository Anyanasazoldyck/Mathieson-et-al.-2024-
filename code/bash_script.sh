#!/bin/bash



DATA_DIR="/mnt/d/CAF_subtypes/data/late_nsclc"
OUT_DIR="/mnt/d/CAF_subtypes/data"

Rscript /mnt/d/CAF_subtypes/code/read_samples.R \
  --data_dir "$DATA_DIR" \
  --out_dir "$OUT_DIR"