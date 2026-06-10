import argparse
import os

import numpy as np
import pandas as pd


OUTPUT_COLUMNS = [
    'test_index',
    'pathway_name',
    'predicted_pathway_score',
    'raw_conformal_pval',
    'adjusted_pval',
    'confidence_score',
    'cluster_label',
]


def read_matrix(filename, label):
    df = pd.read_csv(filename, sep='\t', header=0, index_col=0)
    if df.empty:
        raise ValueError(f'{label} matrix is empty: {filename}')
    try:
        df = df.astype(float)
    except ValueError as exc:
        raise ValueError(f'{label} matrix must contain numeric values: {filename}') from exc
    return df


def read_pathway_names(filename):
    with open(filename, 'r') as file:
        names = [line.strip() for line in file if line.strip() != '']
    if len(names) == 0:
        raise ValueError(f'No pathway names found in {filename}')
    return names


def read_clusters(filename, expected_index, label):
    df = pd.read_csv(filename, sep='\t', header=0, index_col=0)
    if df.empty:
        raise ValueError(f'{label} cluster file is empty: {filename}')
    if df.shape[1] != 1:
        raise ValueError(
                f'{label} cluster file must have exactly one cluster-label column: '
                f'{filename}')
    clusters = df.iloc[:, 0]
    if len(clusters) != len(expected_index):
        raise ValueError(
                f'{label} cluster count ({len(clusters)}) does not match matrix '
                f'row count ({len(expected_index)})')
    if set(expected_index).issubset(set(clusters.index)):
        clusters = clusters.loc[expected_index]
    elif not clusters.index.equals(expected_index):
        clusters.index = expected_index
    return clusters.astype(str)


def validate_inputs(cal_true, cal_pred, test_pred, pathway_names):
    if cal_true.shape != cal_pred.shape:
        raise ValueError(
                'Calibration true and calibration predicted matrices must have '
                f'the same shape, got {cal_true.shape} and {cal_pred.shape}')
    if cal_true.shape[1] != len(pathway_names):
        raise ValueError(
                'Number of pathway names must match calibration matrix columns, '
                f'got {len(pathway_names)} names and {cal_true.shape[1]} columns')
    if test_pred.shape[1] != len(pathway_names):
        raise ValueError(
                'Number of pathway names must match test prediction columns, '
                f'got {len(pathway_names)} names and {test_pred.shape[1]} columns')


def bh_adjust(pvals):
    pvals = np.asarray(pvals, dtype=float)
    adjusted = np.full(pvals.shape, np.nan, dtype=float)
    finite = np.isfinite(pvals)
    if not finite.any():
        return adjusted

    finite_pvals = pvals[finite]
    try:
        from statsmodels.stats.multitest import multipletests
        adjusted[finite] = multipletests(finite_pvals, method='fdr_bh')[1]
        return adjusted
    except ImportError:
        pass

    order = np.argsort(finite_pvals)
    ranked = finite_pvals[order]
    n = len(ranked)
    ranked_adj = ranked * n / np.arange(1, n + 1)
    ranked_adj = np.minimum.accumulate(ranked_adj[::-1])[::-1]
    ranked_adj = np.clip(ranked_adj, 0, 1)
    out = np.empty_like(ranked_adj)
    out[order] = ranked_adj
    adjusted[finite] = out
    return adjusted


def confidence_from_adjusted_pvals(adjusted_pvals):
    adjusted_pvals = np.asarray(adjusted_pvals, dtype=float)
    confidence = np.full(adjusted_pvals.shape, np.nan, dtype=float)
    finite = np.isfinite(adjusted_pvals) & (adjusted_pvals > 0)
    if not finite.any():
        return confidence

    # UTOPIA-style conversion from adjusted p-value to a bounded confidence.
    numerator = -np.log10(adjusted_pvals[finite]) + np.log10(0.5)
    denominator = -np.log10(0.01) + np.log10(0.5)
    confidence[finite] = np.clip(numerator / denominator, 0, 1)
    return confidence


def compute_group_rows(
        pathway_name, pathway_idx, cluster_label, cal_true, cal_pred,
        test_pred, threshold, min_null):
    cal_true_values = cal_true[:, pathway_idx]
    cal_pred_values = cal_pred[:, pathway_idx]
    test_pred_values = test_pred[:, pathway_idx]

    null_mask = (
            np.isfinite(cal_true_values)
            & np.isfinite(cal_pred_values)
            & (cal_true_values <= threshold))
    null_pred = cal_pred_values[null_mask]

    raw_pvals = np.full(test_pred_values.shape, np.nan, dtype=float)
    if len(null_pred) >= min_null:
        for i, value in enumerate(test_pred_values):
            if np.isfinite(value):
                raw_pvals[i] = (
                        1 + np.sum(null_pred >= value)) / (1 + len(null_pred))

    return pd.DataFrame({
        'pathway_name': pathway_name,
        'predicted_pathway_score': test_pred_values,
        'raw_conformal_pval': raw_pvals,
        'cluster_label': cluster_label,
    })


def compute_confidence(
        cal_true, cal_pred, test_pred, pathway_names, threshold=0.05,
        min_null=10, cal_clusters=None, test_clusters=None):
    # UTOPIA-inspired conformal calibration for pathway activity predictions.
    # The null calibration set is calibration spots whose true pathway activity
    # is at or below the activity threshold.
    cal_true = cal_true.copy()
    cal_pred = cal_pred.copy()
    test_pred = test_pred.copy()
    cal_true.columns = pathway_names
    cal_pred.columns = pathway_names
    test_pred.columns = pathway_names

    rows = []
    if cal_clusters is None and test_clusters is None:
        cal_groups = {'global': np.ones(len(cal_true), dtype=bool)}
        test_groups = {'global': np.ones(len(test_pred), dtype=bool)}
    elif cal_clusters is not None and test_clusters is not None:
        labels = sorted(set(cal_clusters).intersection(set(test_clusters)))
        cal_groups = {
            label: (cal_clusters.to_numpy() == label)
            for label in labels}
        test_groups = {
            label: (test_clusters.to_numpy() == label)
            for label in labels}
    else:
        raise ValueError(
                'Provide both --cal-clusters and --test-clusters, or neither.')

    for pathway_idx, pathway_name in enumerate(pathway_names):
        for cluster_label in cal_groups:
            test_mask = test_groups[cluster_label]
            if not test_mask.any():
                continue
            group = compute_group_rows(
                    pathway_name=pathway_name,
                    pathway_idx=pathway_idx,
                    cluster_label=cluster_label,
                    cal_true=cal_true.to_numpy()[cal_groups[cluster_label]],
                    cal_pred=cal_pred.to_numpy()[cal_groups[cluster_label]],
                    test_pred=test_pred.to_numpy()[test_mask],
                    threshold=threshold,
                    min_null=min_null)
            group.insert(0, 'test_index', test_pred.index[test_mask].to_numpy())
            rows.append(group)

    if len(rows) == 0:
        raise ValueError('No test rows were available for confidence scoring.')

    result = pd.concat(rows, ignore_index=True)
    if cal_clusters is None:
        result['cluster_label'] = np.nan

    group_cols = ['pathway_name']
    if cal_clusters is not None:
        group_cols.append('cluster_label')

    result['adjusted_pval'] = np.nan
    groupby_key = group_cols[0] if len(group_cols) == 1 else group_cols
    for _, idx in result.groupby(groupby_key, dropna=False).groups.items():
        idx = list(idx)
        result.loc[idx, 'adjusted_pval'] = bh_adjust(
                result.loc[idx, 'raw_conformal_pval'].to_numpy())

    result['confidence_score'] = confidence_from_adjusted_pvals(
            result['adjusted_pval'].to_numpy())
    return result[OUTPUT_COLUMNS]


def write_outputs(result, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    result.to_csv(
            os.path.join(out_dir, 'pathway_confidence_results.tsv'),
            sep='\t', index=False)


def get_args():
    parser = argparse.ArgumentParser(
            description=(
                'UTOPIA-inspired post-hoc conformal confidence scoring for '
                'pathway-iSTAR predictions. Run after pathway prediction using '
                'wide TSV matrices with spots as rows and pathways as columns.'))
    parser.add_argument('--cal-true', required=True, type=str)
    parser.add_argument('--cal-pred', required=True, type=str)
    parser.add_argument('--test-pred', required=True, type=str)
    parser.add_argument('--pathway-names', required=True, type=str)
    parser.add_argument('--out-dir', required=True, type=str)
    parser.add_argument('--threshold', default=0.05, type=float)
    parser.add_argument('--cal-clusters', default=None, type=str)
    parser.add_argument('--test-clusters', default=None, type=str)
    parser.add_argument('--min-null', default=10, type=int)
    return parser.parse_args()


def main():
    args = get_args()
    pathway_names = read_pathway_names(args.pathway_names)
    cal_true = read_matrix(args.cal_true, 'Calibration true')
    cal_pred = read_matrix(args.cal_pred, 'Calibration predicted')
    test_pred = read_matrix(args.test_pred, 'Test predicted')

    validate_inputs(cal_true, cal_pred, test_pred, pathway_names)

    cal_clusters = None
    test_clusters = None
    if args.cal_clusters is not None or args.test_clusters is not None:
        if args.cal_clusters is None or args.test_clusters is None:
            raise ValueError(
                    'Cluster-stratified confidence requires both '
                    '--cal-clusters and --test-clusters.')
        cal_clusters = read_clusters(
                args.cal_clusters, cal_true.index, 'Calibration')
        test_clusters = read_clusters(
                args.test_clusters, test_pred.index, 'Test')

    result = compute_confidence(
            cal_true=cal_true,
            cal_pred=cal_pred,
            test_pred=test_pred,
            pathway_names=pathway_names,
            threshold=args.threshold,
            min_null=args.min_null,
            cal_clusters=cal_clusters,
            test_clusters=test_clusters)
    write_outputs(result, args.out_dir)
    print(f'Wrote pathway confidence outputs to {args.out_dir}')


if __name__ == '__main__':
    main()
