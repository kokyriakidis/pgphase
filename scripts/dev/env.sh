# Configuration shared by scripts/dev. Sourced, not executed.
#
# Every PGDEV_* value may come from the environment. Machine-specific values
# belong in scripts/dev/local.env (untracked), written as defaults so that the
# environment still wins, e.g.  : "${PGDEV_TRUTH_VCF:=/data/truth.vcf.gz}"

PGDEV_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
if [[ -f "${PGDEV_ROOT}/scripts/dev/local.env" ]]; then
    # shellcheck disable=SC1091
    source "${PGDEV_ROOT}/scripts/dev/local.env"
fi

: "${PGDEV_DATA:=${PGDEV_ROOT}/test_data}"
: "${PGDEV_REF:=${PGDEV_DATA}/chm13v2.0.chr20.renamed.fa}"
: "${PGDEV_BAM:=${PGDEV_DATA}/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam}"
: "${PGDEV_SITES:=${PGDEV_DATA}/chr20.sites.striped.vcf.gz}"
: "${PGDEV_GAF:=${PGDEV_DATA}/HG002.chr20.annotated.coord.gaf.gz}"
: "${PGDEV_CONTIG:=CHM13#0#chr20}"
# Parental read labels (read<TAB>MATERNAL|PATERNAL), from scripts/make_truth_hap_map.sh.
: "${PGDEV_TRUTH_READS:=${PGDEV_DATA}/derived/chr20_truth_hap.tsv}"
# Phased truth VCF for variant-level scores and the phaseable-read split.
: "${PGDEV_TRUTH_VCF:=}"
# Space-separated NAME=DIR pairs; each DIR holds phased.bam and phased.vcf[.gz].
: "${PGDEV_COMPETITORS:=}"
# Centromere/satellite interval on the contig, 0-based half-open.
: "${PGDEV_CEN:=26000000-32000000}"
# A Python with pysam.
: "${PGDEV_PYTHON:=python3}"
: "${PGDEV_CACHE:=${PGDEV_ROOT}/.dev-cache}"
: "${PGDEV_THREADS:=$(nproc)}"

export PGDEV_ROOT PGDEV_DATA PGDEV_REF PGDEV_BAM PGDEV_SITES PGDEV_GAF PGDEV_CONTIG \
       PGDEV_TRUTH_READS PGDEV_TRUTH_VCF PGDEV_COMPETITORS PGDEV_CEN PGDEV_PYTHON \
       PGDEV_CACHE PGDEV_THREADS
