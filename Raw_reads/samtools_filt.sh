#!/usr/bin/env bash

INDIR="/media/labgenoma5/slbonatto_20b/tursiops/00_paleomix"
OUTDIR="/media/labgenoma5/slbonatto_20b/tursiops/00_masked_bam"
BED="/media/labgenoma5/slbonatto_20b/references/tursiops/ref_with_mito/GCA_011762595.2_mTurTru1.mat.Y_genomic.nonRepeat.autossomalChr.bed"

#### mkdir -p "$OUTDIR"

for BAM in "$INDIR"/*.bam; do
    SAMPLE=$(basename "$BAM" .bam)
    OUTBAM="$OUTDIR/${SAMPLE}.nonRepeat.bam"

    echo "==== Filtering $SAMPLE ===="

    if [ -e "$OUTBAM" ]; then
        echo "Skipping $SAMPLE (already filtered)"
        continue
    fi

    samtools view -@ 4 -b \
        -L "$BED" \
        "$BAM" \
        > "$OUTBAM"

    samtools index "$OUTBAM"

    echo "Finished $SAMPLE"
    echo
done
