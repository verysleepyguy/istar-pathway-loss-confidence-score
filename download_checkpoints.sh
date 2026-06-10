#!/bin/bash
set -e

source_256="https://github.com/mahmoodlab/HIPT/raw/master/HIPT_4K/Checkpoints/vit256_small_dino.pth"
source_4k="https://github.com/mahmoodlab/HIPT/raw/master/HIPT_4K/Checkpoints/vit4k_xs_dino.pth"
target_256="checkpoints/vit256_small_dino.pth"
target_4k="checkpoints/vit4k_xs_dino.pth"

mkdir -p checkpoints
wget --tries=3 --timeout=60 ${source_256} -O ${target_256}
wget --tries=3 --timeout=60 ${source_4k} -O ${target_4k}

if [ ! -s "${target_256}" ] || [ ! -s "${target_4k}" ]; then
  echo "Checkpoint download failed or produced an empty file." >&2
  exit 1
fi
