#!/usr/bin/env bash
# =============================================================================
# Script: 01_extract_SLC30A8_variants.sh
# Project: Population and Evolutionary Analysis of SLC30A8 in Type 2 Diabetes
#
# Description:
#   Downloads 1000 Genomes Phase 1 VCF data for chromosome 8 and GENCODE
#   gene annotations, then extracts SLC30A8 variants and classifies them
#   as exonic or non-exonic using bedtools intersect.
#
# Requirements:
#   - wget
#   - gzip / awk / sed (standard Unix tools)
#   - bedtools (https://bedtools.readthedocs.io)
#     Expected at: ./bedtools2-master/bin/intersectBed
#
# Usage:
#   chmod +x 01_extract_SLC30A8_variants.sh
#   ./01_extract_SLC30A8_variants.sh
#
# Output files:
#   SLC30A8_variants.bed   — all SLC30A8 variants from 1000 Genomes
#   SLC30A8_exons.bed      — exon coordinates from GENCODE v49
#   exon_overlapping.bed   — variants overlapping exonic regions
#   non_exon.bed           — variants in non-exonic regions
# =============================================================================

set -euo pipefail   # exit on error, undefined variable, or pipe failure

# ----- Configuration ---------------------------------------------------------
SLC30A8_CHR=8
SLC30A8_START=117134995
SLC30A8_END=117176714

cd ~/Desktop/5330_Workspace_KimLuong/project
WORKDIR="$(pwd)"    # run from your project directory
BEDTOOLS="./bedtools2-master/bin/intersectBed"

echo "Working directory: ${WORKDIR}"
echo "SLC30A8 region: chr${SLC30A8_CHR}:${SLC30A8_START}-${SLC30A8_END}"
echo ""


# =============================================================================
# STEP 1: Download 1000 Genomes Phase 1 VCF (chr8)
# =============================================================================
echo "[1/4] Downloading 1000 Genomes chr8 VCF..."

VCF="ALL.chr8.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf.gz"
VCF_URL="https://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase1/analysis_results/integrated_call_sets/${VCF}"

wget -c "${VCF_URL}"
echo "      Done: ${VCF}"


# =============================================================================
# STEP 2: Extract SLC30A8 Variants -> BED format
# =============================================================================
echo "[2/4] Extracting SLC30A8 variants from VCF..."

gzip -cd "${VCF}" \
  | awk -v chr="${SLC30A8_CHR}" \
        -v start="${SLC30A8_START}" \
        -v end="${SLC30A8_END}" \
    'BEGIN {OFS="\t"} $1==chr && $2>=start && $2<=end {
        print $1, $2-1, $2, $3, $8, ".", "+"
    }' \
  > SLC30A8_variants.bed

echo "      Done: SLC30A8_variants.bed ($(wc -l < SLC30A8_variants.bed) variants)"


# =============================================================================
# STEP 3: Download GENCODE v49 Annotations and Extract SLC30A8 Exons
# =============================================================================
echo "[3/4] Downloading GENCODE v49 annotation and extracting SLC30A8 exons..."

GTF_GZ="gencode.v49.annotation.gtf.gz"
GTF_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz"

wget -c "${GTF_URL}"

# Extract exon entries for SLC30A8 region, strip "chr" prefix to match VCF
gzip -cd "${GTF_GZ}" \
  | awk -F'\t' -v start="${SLC30A8_START}" -v end="${SLC30A8_END}" \
    'BEGIN {OFS="\t"} $1=="chr8" && $3=="exon" && $4>=start && $5<=end {
        print $1, $4, $5, $9, ".", $7
    }' \
  | sed 's/^chr//' \
  > SLC30A8_exons.bed

echo "      Done: SLC30A8_exons.bed ($(wc -l < SLC30A8_exons.bed) exon entries)"


# =============================================================================
# STEP 4: Classify Variants — Exonic vs Non-Exonic (bedtools intersect)
# =============================================================================
echo "[4/4] Running bedtools intersect to classify variants..."

# Variants overlapping exons
"${BEDTOOLS}" \
  -u \
  -a "${WORKDIR}/SLC30A8_variants.bed" \
  -b "${WORKDIR}/SLC30A8_exons.bed" \
  > "${WORKDIR}/exon_overlapping.bed"

echo "      Exonic variants:     $(wc -l < exon_overlapping.bed)"

# Variants NOT overlapping exons (intronic/intergenic)
"${BEDTOOLS}" \
  -v \
  -a "${WORKDIR}/SLC30A8_variants.bed" \
  -b "${WORKDIR}/SLC30A8_exons.bed" \
  > "${WORKDIR}/non_exon.bed"

echo "      Non-exonic variants: $(wc -l < non_exon.bed)"


# =============================================================================
# DONE
# =============================================================================
echo ""
echo "Pipeline complete. Output files:"
echo "  SLC30A8_variants.bed   — all SLC30A8 variants"
echo "  SLC30A8_exons.bed      — GENCODE v49 exon coordinates"
echo "  exon_overlapping.bed   — exonic variants"
echo "  non_exon.bed           — non-exonic variants"
