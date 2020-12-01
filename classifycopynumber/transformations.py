import numpy as np
import pandas as pd

def calculate_log_change(data, ploidy):
    normalize = (
        data.groupby('gene_name')['overlap_width']
        .sum().rename('sum_overlap_width').reset_index())

    data['total_raw_weighted'] = data['copy'] * data['overlap_width']

    data = data.groupby(['gene_name'])['total_raw_weighted'].sum().reset_index()
    data = data.merge(normalize)
    data['total_raw_mean'] = data['total_raw_weighted'] / data['sum_overlap_width']
    data['log_change'] = np.log2(data['total_raw_mean']/ploidy)
    return data[["gene_name", "total_raw_mean", "total_raw_weighted", "log_change"]]


def aggregate_adjacent(cnv, value_cols=(), stable_cols=(), length_normalized_cols=(), summed_cols=()):
    """ Aggregate adjacent segments with similar copy number state.

    see: https://github.com/amcpherson/remixt/blob/master/remixt/segalg.py

    Args:
        cnv (pandas.DataFrame): copy number table

    KwArgs:
        value_cols (list): list of columns to compare for equivalent copy number state
        stable_cols (list): columns for which values are the same between equivalent states
        length_normalized_cols (list): columns that are width normalized for equivalent states

    Returns:
        pandas.DataFrame: copy number with adjacent segments aggregated
    """

    # Group segments with same state
    cnv = cnv.sort_values(['chr', 'start'])
    cnv['chromosome_index'] = np.searchsorted(np.unique(cnv['chr']), cnv['chr'])
    cnv['diff'] = cnv[['chromosome_index'] + value_cols].diff().abs().sum(axis=1)
    cnv['is_diff'] = (cnv['diff'] != 0)
    cnv['cn_group'] = cnv['is_diff'].cumsum()

    def agg_segments(df):
        a = df[stable_cols].iloc[0]
        a['chr'] = df['chr'].iloc[0] # can't be min chr is str
        a['start'] = df['start'].min()
        a['end'] = df['end'].max()
        a['width'] = df['width'].sum()

        for col in length_normalized_cols:
            a[col] = (df[col] * df['width']).sum() / (df['width'].sum() + 1e-16)

        for col in summed_cols:
            a[col] = df[col].sum()
        return pd.Series(a)

    aggregated = cnv.groupby('cn_group').apply(agg_segments)
    for col in aggregated:
        aggregated[col] = aggregated[col].astype(cnv[col].dtype)

    return aggregated


