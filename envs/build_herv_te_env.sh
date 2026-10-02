#!/bin/bash
# Build the TE-analysis environment (herv_te). Run on a compute node (srun) or the login node.
#   conda/bioconda: python 3.11, samtools, bedtools, pysam
#   PyPI:           IRescue (subfamily TE counts for 10x, UMI-aware EM), pychromvar, scanpy
#   GitHub:         scTE (JiekaiLab/scTE), installed from source
# Stellarscope (locus level) is added in a second step, once the subfamily analysis is running.
set -euo pipefail
source "$(dirname "$0")/../config.sh"
MM=/projects/hzhu10@xsede.org/software/micromamba_bootstrap/bin/micromamba
export MAMBA_ROOT_PREFIX=/projects/hzhu10@xsede.org/software/micromamba_root

"$MM" create -y -p "$TE_ENV" -c conda-forge -c bioconda \
  python=3.11 samtools bedtools pysam pip git
"$TE_ENV/bin/pip" install --no-cache-dir irescue pychromvar scanpy anndata pandas numpy scipy
"$TE_ENV/bin/pip" install --no-cache-dir "git+https://github.com/JiekaiLab/scTE.git"

echo "--- versions"
"$TE_ENV/bin/irescue" --version || true
"$TE_ENV/bin/python" -c "import pychromvar, scanpy; print('pychromvar', pychromvar.__version__, 'scanpy', scanpy.__version__)"
"$TE_ENV/bin/scTE" -h | head -3 || true
"$TE_ENV/bin/samtools" --version | head -1
