import argparse
import os

import numpy as np
import pandas as pd
import squidpy as sq
from PIL import Image


def save_image(img, outfile):
    img = np.asarray(img)
    if img.max() <= 1:
        img = (img * 255).clip(0, 255)
    img = img.astype(np.uint8)
    Image.fromarray(img).save(outfile)


def get_counts(adata):
    x = adata.X
    if hasattr(x, 'toarray'):
        x = x.toarray()
    return pd.DataFrame(x, index=adata.obs_names, columns=adata.var_names)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--prefix', default='data/visium_mouse_brain/')
    parser.add_argument(
            '--pixel-size', type=float, default=0.5,
            help='Microns per pixel for the saved hires image.')
    args = parser.parse_args()

    prefix = args.prefix
    os.makedirs(prefix, exist_ok=True)

    adata = sq.datasets.visium_hne_adata()
    library_id = list(adata.uns['spatial'].keys())[0]
    spatial = adata.uns['spatial'][library_id]
    scalefactors = spatial['scalefactors']
    hires_scale = float(scalefactors['tissue_hires_scalef'])

    # iSTAR needs a histology image. Squidpy provides the 10x hires image.
    save_image(spatial['images']['hires'], os.path.join(prefix, 'he-raw.jpg'))

    # adata.obsm['spatial'] stores full-resolution pixel coordinates.
    # Because he-raw.jpg is the hires image, convert spot centers/radius into
    # that same hires pixel coordinate system before writing iSTAR inputs.
    locs = pd.DataFrame(
            adata.obsm['spatial'] * hires_scale,
            index=adata.obs_names,
            columns=['x', 'y']).round().astype(int)
    locs.to_csv(os.path.join(prefix, 'locs-raw.tsv'), sep='\t')

    cnts = get_counts(adata)
    cnts.to_csv(os.path.join(prefix, 'cnts.tsv'), sep='\t')

    radius = float(scalefactors['spot_diameter_fullres']) * hires_scale / 2
    with open(os.path.join(prefix, 'radius-raw.txt'), 'w') as file:
        file.write(str(int(round(radius))))
    with open(os.path.join(prefix, 'pixel-size-raw.txt'), 'w') as file:
        file.write(str(args.pixel_size))

    print('Created iSTAR inputs in', prefix)
    print('counts:', cnts.shape)
    print('hires_scale:', hires_scale)
    print('radius_raw_hires_pixels:', radius)


if __name__ == '__main__':
    main()
