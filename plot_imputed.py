import argparse
import sys

import numpy as np
import matplotlib.pyplot as plt

from utils import load_pickle, save_image, read_lines, load_image
# from visual import cmap_turbo_truncated


def plot_super(
        x, outfile, underground=None, truncate=None):

    x = x.copy()
    mask = np.isfinite(x)

    if truncate is not None:
        x -= np.nanmean(x)
        x /= np.nanstd(x) + 1e-12
        x = np.clip(x, truncate[0], truncate[1])

    x -= np.nanmin(x)
    x /= np.nanmax(x) + 1e-12

    cmap = plt.get_cmap('turbo')
    # cmap = cmap_turbo_truncated
    if underground is not None:
        under = underground.mean(-1, keepdims=True)
        under -= under.min()
        under /= under.max() + 1e-12

    img = cmap(x)[..., :3]
    if underground is not None:
        img = img * 0.5 + under * 0.5
    img[~mask] = 1.0
    img = (img * 255).astype(np.uint8)
    save_image(img, outfile)


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('prefix', type=str)
    parser.add_argument('--names', type=str, default='gene-names.txt')
    parser.add_argument('--input-dir', type=str, default='cnts-super')
    parser.add_argument('--output-dir', type=str, default='cnts-super-plots')
    parser.add_argument('--file-names', type=str, default=None)
    return parser.parse_args()


def main():

    args = get_args()
    prefix = args.prefix  # e.g. 'data/her2st/B1/'
    names = read_lines(f'{prefix}{args.names}')
    file_names = names
    if args.file_names is not None:
        file_names = read_lines(f'{prefix}{args.file_names}')
    mask = load_image(f'{prefix}mask-small.png') > 0

    for name, file_name in zip(names, file_names):
        cnts = load_pickle(f'{prefix}{args.input_dir}/{file_name}.pickle')
        cnts[~mask] = np.nan
        plot_super(cnts, f'{prefix}{args.output_dir}/{file_name}.png')


if __name__ == '__main__':
    main()
