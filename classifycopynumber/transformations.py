import numpy as np
import pandas as pd
import wgs_analysis.algorithms.cnv


def generate_segmental_cn(filename, cn, stats_data):
    aggregated_cn_data = wgs_analysis.algorithms.cnv.aggregate_adjacent(
        cn,
        value_cols=['major_0', 'minor_0', 'major_1', 'minor_1', 'major_2', 'minor_2'],
        stable_cols=['major_0', 'minor_0', 'major_1', 'minor_1', 'major_2', 'minor_2', 'sample'],
        length_normalized_cols=['major_raw', 'minor_raw'],
    )
    # Clean up segs and write to disk
    aggregated_cn_data['ploidy'] = stats_data['ploidy']
    aggregated_cn_data['total_raw'] = aggregated_cn_data['major_raw'] + aggregated_cn_data['minor_raw']
    aggregated_cn_data['seg.mean'] = np.log2(aggregated_cn_data['total_raw'] / aggregated_cn_data['ploidy'])
    aggregated_cn_data['num.mark'] = (aggregated_cn_data['length'] / 500000).astype(int)
    aggregated_cn_data = aggregated_cn_data.rename(columns={'sample': 'ID', 'chromosome': 'chrom', 'start': 'loc.start', 'end': 'loc.end'})
    aggregated_cn_data = aggregated_cn_data[['ID', 'chrom', 'loc.start', 'loc.end', 'num.mark', 'seg.mean']]
    aggregated_cn_data['seg.mean'] = aggregated_cn_data['seg.mean'].fillna(np.exp(-8))
    aggregated_cn_data.loc[aggregated_cn_data['seg.mean'] == np.NINF, 'seg.mean'] = np.exp(-8)
    aggregated_cn_data.to_csv(filename, index=None, sep='\t')


def calculate_log_change(data, ploidy):
    normalize = (
        data.groupby('gene_name')['overlap_width']
        .sum().rename('sum_overlap_width').reset_index())

    data['total_raw_weighted'] = data['copy'] * data['overlap_width']

    data = data.groupby(['gene_name'])['total_raw_weighted'].sum().reset_index()
    data = data.merge(normalize)

    data['total_raw_mean'] = data['total_raw_weighted'] / data['sum_overlap_width']
    data['log_change'] = np.log2(data['total_raw_mean']/ploidy)
    data["ploidy"] = [ploidy] * len(data)
    return data


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


