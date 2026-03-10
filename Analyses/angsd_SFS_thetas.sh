#!/usr/bin/env bash

INDIR="/media/labgenoma5/slbonatto_20b/tursiops/00_masked_bam"
OUTDIR="/media/labgenoma5/slbonatto_20b/tursiops/03_het/autossom_het"
REF="/media/labgenoma5/slbonatto_20b/references/tursiops/ref_with_mito/GCA_011762595.2_mTurTru1.mat.Y_genomic.fasta"

# Creates output dir
mkdir -p "$OUTDIR"

for BAM in "$INDIR"/*.bam; do
    SAMPLE=$(basename "$BAM" .bam) #define sample name

    echo "==== $SAMPLE ===="
  
    #Check if final pestPG file exists. If it does, skip, otherwise performs the required steps
    if [ -e "$OUTDIR/${SAMPLE}.thetaStat.pestPG" ]; then
        echo "Skipping $SAMPLE (already done)"
        continue
    fi

    echo "Running ANGSD for $SAMPLE"
    #Run ansdg with default quality metrics
    angsd \
        -i "$BAM" \
        -anc "$REF" \
        -doSaf 1 \
        -GL 1 \
        -minMapQ 30 \
        -minQ 30 \
        -remove_bads 1 \
        -uniqueOnly 1 \
        -only_proper_pairs 1 \
        -out "$OUTDIR/${SAMPLE}" \
        -nThreads 6

    echo "Running realSFS"
    # Create SFS for the sample. This file will be used for heterozigosity estimates and as input for theta estimates
    realSFS "$OUTDIR/${SAMPLE}.saf.idx" -fold 1 \
        > "$OUTDIR/${SAMPLE}.sfs"

    echo "Running saf2theta"
    realSFS saf2theta \
        "$OUTDIR/${SAMPLE}.saf.idx" \
        -sfs "$OUTDIR/${SAMPLE}.sfs" \
        -outname "$OUTDIR/${SAMPLE}" \
        -anc "$REF"

    echo "Running thetaStat"
    # Run thetaStat to create theta estimates for the sample
    thetaStat do_stat \
        "$OUTDIR/${SAMPLE}.thetas.idx" \
        -outnames "$OUTDIR/${SAMPLE}.thetaStat"

    echo "Cleaning intermediate files"
    # Remove intermediary files to reduce storage usage, keeping only final sfs and pestPG text-like files
    rm -f \
        "$OUTDIR/${SAMPLE}.saf.idx" \
        "$OUTDIR/${SAMPLE}.saf.gz" \
        "$OUTDIR/${SAMPLE}.saf.pos.gz" \
        "$OUTDIR/${SAMPLE}.thetas.gz" \
        "$OUTDIR/${SAMPLE}.thetas.idx"

    echo "Finished $SAMPLE"
    echo
done
