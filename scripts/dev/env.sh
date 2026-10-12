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
# PGDEV_PROFILE=chr12 (etc.) selects scripts/dev/profiles/chr12.env; the
# default is HG002 chr20.
if [[ -n "${PGDEV_PROFILE:-}" ]]; then
    [[ -f "${PGDEV_ROOT}/scripts/dev/profiles/${PGDEV_PROFILE}.env" ]] ||
        { echo "unknown PGDEV_PROFILE ${PGDEV_PROFILE}" >&2; return 1 2>/dev/null || exit 1; }
    # shellcheck disable=SC1090
    source "${PGDEV_ROOT}/scripts/dev/profiles/${PGDEV_PROFILE}.env"
fi
: "${PGDEV_CHROM:=chr20}"

: "${PGDEV_REF:=${PGDEV_DATA}/chm13v2.0.chr20.renamed.fa}"
: "${PGDEV_BAM:=${PGDEV_DATA}/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam}"
: "${PGDEV_SITES:=${PGDEV_DATA}/chr20.sites.striped.vcf.gz}"
: "${PGDEV_GAF:=${PGDEV_DATA}/HG002.chr20.annotated.coord.gaf.gz}"
: "${PGDEV_CONTIG:=CHM13#0#chr20}"
# Parental read labels (read<TAB>MATERNAL|PATERNAL), from scripts/make_truth_hap_map.sh.
: "${PGDEV_TRUTH_READS:=${PGDEV_DATA}/derived/chr20_truth_hap.tsv}"
# Evaluation root holding <chrom>/truth.vcf.gz and <chrom>/<tool>/ competitor
# outputs (phased.bam, phased.vcf[.gz]); set it in local.env.
: "${PGDEV_EVAL_ROOT:=}"
# Phased truth VCF for variant-level scores and the phaseable-read split.
: "${PGDEV_TRUTH_VCF:=${PGDEV_EVAL_ROOT:+${PGDEV_EVAL_ROOT}/${PGDEV_CHROM}/truth.vcf.gz}}"
# Space-separated NAME=DIR pairs; each DIR holds phased.bam and phased.vcf[.gz].
: "${PGDEV_COMPETITORS:=${PGDEV_EVAL_ROOT:+hiphase=${PGDEV_EVAL_ROOT}/${PGDEV_CHROM}/hiphase longphase=${PGDEV_EVAL_ROOT}/${PGDEV_CHROM}/longphase}}"
# Regions excluded from scoring because the truth contradicts itself (BED);
# the default file is for chr20, a profile sets its own (empty = none).
: "${PGDEV_EXCLUDE=${PGDEV_ROOT}/scripts/dev/truth_exclude.bed}"
# Centromere/satellite interval on the contig, 0-based half-open (chr20 default).
: "${PGDEV_CEN:=26000000-32000000}"
# A Python with pysam.
: "${PGDEV_PYTHON:=python3}"
: "${PGDEV_CACHE:=${PGDEV_ROOT}/.dev-cache}"
: "${PGDEV_THREADS:=$(nproc)}"

# Runs and baselines live per profile.
PGDEV_STATE="${PGDEV_CACHE}${PGDEV_PROFILE:+/${PGDEV_PROFILE}}"

export PGDEV_ROOT PGDEV_DATA PGDEV_PROFILE PGDEV_CHROM PGDEV_EVAL_ROOT PGDEV_STATE PGDEV_REF PGDEV_BAM PGDEV_SITES PGDEV_GAF PGDEV_CONTIG \
       PGDEV_TRUTH_READS PGDEV_TRUTH_VCF PGDEV_COMPETITORS PGDEV_EXCLUDE PGDEV_CEN PGDEV_PYTHON \
       PGDEV_CACHE PGDEV_THREADS
