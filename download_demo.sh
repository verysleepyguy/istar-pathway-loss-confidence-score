#!/bin/bash
set -e

prefix=$1

source_img="https://upenn.box.com/shared/static/yya0lvlur8aase29hvy630jd06r64tdn.jpg"
source_cnts="https://upenn.box.com/shared/static/kaoo8j31dx5lupyz8dctay7p5x3exqsa.tsv"
source_locs="https://upenn.box.com/shared/static/7nbnorlr2h6tkeyghjibqitztezkadwh.tsv"
source_radius="https://upenn.box.com/shared/static/a8655bmb02q9cqegndnwhcb5r0mqqphi.txt"
source_pixsize="https://upenn.box.com/shared/static/1stmq5ly6iqnljt0uq8rotlki5q8sjfs.txt"

target_img="${prefix}he-raw.jpg"
target_cnts="${prefix}cnts.tsv"
target_locs="${prefix}locs-raw.tsv"
target_radius="${prefix}radius-raw.txt"
target_pixsize="${prefix}pixel-size-raw.txt"

download_file () {
  source_url=$1
  target_file=$2
  mkdir -p "$(dirname "${target_file}")"
  wget --tries=3 --timeout=60 "${source_url}" -O "${target_file}"
  if [ ! -s "${target_file}" ]; then
    echo "Download failed or produced an empty file: ${target_file}" >&2
    exit 1
  fi
}

download_file "${source_img}" "${target_img}"
download_file "${source_cnts}" "${target_cnts}"
download_file "${source_locs}" "${target_locs}"
download_file "${source_radius}" "${target_radius}"
download_file "${source_pixsize}" "${target_pixsize}"

python validate_inputs.py "${prefix}"
