import argparse
import os
import sys

import pandas as pd
from PIL import Image


REQUIRED_RAW_FILES = [
    'he-raw.jpg',
    'cnts.tsv',
    'locs-raw.tsv',
    'pixel-size-raw.txt',
    'radius-raw.txt',
]


def check_nonempty(path, errors):
    if not os.path.exists(path):
        errors.append(f'Missing file: {path}')
        return False
    if os.path.getsize(path) == 0:
        errors.append(f'Empty file: {path}')
        return False
    return True


def validate_prefix(prefix):
    errors = []
    for name in REQUIRED_RAW_FILES:
        check_nonempty(os.path.join(prefix, name), errors)

    cnts_file = os.path.join(prefix, 'cnts.tsv')
    locs_file = os.path.join(prefix, 'locs-raw.tsv')
    image_file = os.path.join(prefix, 'he-raw.jpg')

    if os.path.exists(cnts_file) and os.path.getsize(cnts_file) > 0:
        try:
            cnts = pd.read_csv(cnts_file, sep='\t', index_col=0, nrows=5)
            if cnts.shape[1] == 0:
                errors.append(f'No gene columns found in {cnts_file}')
        except Exception as exc:
            errors.append(f'Could not read {cnts_file}: {exc}')

    if os.path.exists(locs_file) and os.path.getsize(locs_file) > 0:
        try:
            locs = pd.read_csv(locs_file, sep='\t', index_col=0, nrows=5)
            if locs.shape[1] < 2:
                errors.append(f'Expected at least two coordinate columns in {locs_file}')
        except Exception as exc:
            errors.append(f'Could not read {locs_file}: {exc}')

    if os.path.exists(image_file) and os.path.getsize(image_file) > 0:
        try:
            with Image.open(image_file) as img:
                img.verify()
        except Exception as exc:
            errors.append(f'Could not open {image_file} as an image: {exc}')

    return errors


def validate_pathways(pathways_file):
    errors = []
    if not check_nonempty(pathways_file, errors):
        return errors

    try:
        pathways = pd.read_csv(pathways_file, sep='\t', header=None, usecols=[0, 1])
    except Exception as exc:
        errors.append(f'Could not read {pathways_file}: {exc}')
        return errors

    if pathways.empty:
        errors.append(f'No pathway-gene rows found in {pathways_file}')
    if pathways.shape[1] < 2:
        errors.append(f'Expected long format: pathway_name<TAB>gene_name in {pathways_file}')
    return errors


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('prefix', type=str)
    parser.add_argument('--pathways', type=str, default=None)
    return parser.parse_args()


def main():
    args = get_args()
    prefix = args.prefix
    errors = validate_prefix(prefix)
    if args.pathways is not None:
        errors.extend(validate_pathways(args.pathways))

    if errors:
        print('Input validation failed:', file=sys.stderr)
        for error in errors:
            print(f'  - {error}', file=sys.stderr)
        print(
                '\nFor the demo, run: ./download_demo.sh data/demo/\n'
                'For GitHub/Colab, do not commit zero-byte placeholders; '
                'use the download scripts, Git LFS, or a release/Drive download.',
                file=sys.stderr)
        sys.exit(1)

    print(f'Input validation passed for {prefix}')


if __name__ == '__main__':
    main()
