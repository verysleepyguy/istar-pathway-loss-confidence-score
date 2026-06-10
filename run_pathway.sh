#!/bin/bash
set -e

prefix=$1      # e.g. data/demo/
pathways=$2    # e.g. data/pathways.tsv

device="${DEVICE:-cuda}"      # set DEVICE=cpu to run without GPU
pixel_size="${PIXEL_SIZE:-0.5}"
n_genes="${N_GENES:-1000}"    # still written because impute.py expects gene-names.txt
epochs="${EPOCHS:-400}"
n_states="${N_STATES:-5}"

if [ -z "${prefix}" ] || [ -z "${pathways}" ]; then
  echo "Usage: ./run_pathway.sh data/sample/ path/to/pathways.tsv"
  exit 1
fi

python validate_inputs.py "${prefix}" --pathways="${pathways}"

# Preprocess histology image.
echo "${pixel_size}" > "${prefix}pixel-size.txt"
python rescale.py "${prefix}" --image
python preprocess.py "${prefix}" --image

# Extract histology features and tissue mask.
python extract_features.py "${prefix}" --device="${device}"
python get_mask.py "${prefix}embeddings-hist.pickle" "${prefix}mask-small.png"

# Keep the original gene-selection step for backward compatibility with
# impute.py setup; pathway mode replaces the actual training target later.
python select_genes.py --n-top="${n_genes}" "${prefix}cnts.tsv" "${prefix}gene-names.txt"

# Rescale spot coordinates and spot radius.
python rescale.py "${prefix}" --locs --radius

# Train against pathway activity targets and save predictions to
# cnts-pathway-super/ without overwriting gene-level cnts-super/.
python impute.py "${prefix}" \
  --epochs="${epochs}" \
  --n-states="${n_states}" \
  --device="${device}" \
  --pathways="${pathways}"

# Optional visualization and pathway-feature clustering.
python plot_imputed.py "${prefix}" \
  --names=pathway-names.txt \
  --file-names=pathway-file-names.txt \
  --input-dir=cnts-pathway-super \
  --output-dir=cnts-pathway-super-plots

python cluster.py \
  --filter-size=8 \
  --min-cluster-size=20 \
  --n-clusters=10 \
  --mask="${prefix}mask-small.png" \
  "${prefix}embeddings-pathway.pickle" \
  "${prefix}clusters-pathway/"
