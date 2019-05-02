#!/bin/bash
# Base quality score recalibration
usage() {
    echo "Usage: $0 [-b <string>] [-k <string>] [-o <string>] [-r <string>]"
    echo "        -b bam file"
    echo "        -k known sites (i.e dbsnp)"
    echo "        -o output name"
    echo "        -r reference"
    exit 1;
}
while getopts ":b:k:o:r:" x; do
    case "${x}" in
        b)
            b=${OPTARG}
            ;;
        k)
            k=${OPTARG}
            ;;
        o)
            o=${OPTARG}
	        ;;
        r)
            r=${OPTARG}
            ;;
        *)
            usage
            ;;
    esac
done
shift $((OPTIND-1))
if [ -z "${b}" ]; then
    usage
fi
if [ -z "${k}" ]; then
    k="common_all_w_chr.vcf.gz"
fi
if [ -z "${o}" ]; then
    usage
fi
if [ -z "${r}"]; then
    r="ref/hg19_k/hg19.fasta"
fi

## BaseRecalibrator
gatk BaseRecalibrator \
    -R $r \
    -I $b \
    --known-sites $k \
    -O ${b%/*}/${o}.table

## applyBQSR
gatk ApplyBQSR \
    -R $r \
    -I $b \
    --bqsr-recal-file ${b%/*}/${o}.table \
    -O ${b%.bam}.bqsr.bam

b2=${b%.bam}.bqsr.bam

gatk QualityScoreDistribution \
    -I ${b2} \
    -O ${b%/*}/qual_dist_score.txt \
    -CHART ${b%/*}/qual_dist_score.pdf

## second pass
gatk BaseRecalibrator \
    -R $r \
    -I $b2 \
    --known-sites $k \
    -o ${o}2.table
 
# gatk ApplyBQSR \
#     -R $r \
#     -I $b2 \
#     --bqsr-recal-file ${b%/*}/${o}2.table \
#     -O ${b%.bam}.pass2.bam
 
gatk AnalyzeCovariates \
     -before ${b%/*}/${o}.table \
     -after ${b%/*}/${o}2.table \
     -plots ${b%/*}/AnalyzieCovariates.pdf
