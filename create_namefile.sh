#!/usr/bin/env bash

set -euo pipefail

if [ "$#" -ne 5 ]; then
    echo "Usage:"
    echo "  $0 <ref_mat.fa> <ref_pat.fa> <raw_dir> <chromosome_number> <output_pathfile>"
    exit 1
fi

RAW_DIR="$3"
CHR="$4"
OUTFILE="$5"
REF_MAT="$1"
REF_PAT="$2"

TMP=$(mktemp)

# --- Add references (chromosome-specific if needed) ---
echo -e "hg002_mat\t${REF_MAT}" >> "$TMP"
echo -e "hg002_pat\t${REF_PAT}" >> "$TMP"

# --- Find chromosome-specific fasta files ---
find "$RAW_DIR" -type f -name "*.fna.filtered.chr${CHR}" | sort | while read -r file; do

    # Extract GCA ID safely
    relative="${file#${RAW_DIR}/}"
    gca_id="${relative%%/*}"

    echo -e "${gca_id}\t${file}" >> "$TMP"

done

# --- Build header line ---
names=$(cut -f1 "$TMP" | paste -sd "," -)

{
    echo "(${names});"
    echo ""
    cat "$TMP"
} > "$OUTFILE"

rm "$TMP"

echo "Chromosome ${CHR} pathfile written to $OUTFILE"
