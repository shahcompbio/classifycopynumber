import numpy as np
import pandas as pd
import wgs_analysis.algorithms.cnv


def generate_segmental_cn(filename, aggregated_cn_data, ploidy,  cn_col="copy", length_col="length"):

    aggregated_cn_data['ploidy'] = ploidy 
    aggregated_cn_data['seg.mean'] = np.log2(aggregated_cn_data[cn_col] / aggregated_cn_data['ploidy'])
    aggregated_cn_data['num.mark'] = (aggregated_cn_data[length_col] / 500000).astype(int)
    aggregated_cn_data = aggregated_cn_data.rename(columns={'sample': 'ID', 'chromosome': 'chrom', 'start': 'loc.start', 'end': 'loc.end'})
    aggregated_cn_data = aggregated_cn_data[['ID', 'chrom', 'loc.start', 'loc.end', 'num.mark', 'seg.mean']]
    aggregated_cn_data['seg.mean'] = aggregated_cn_data['seg.mean'].fillna(np.exp(-8))
    aggregated_cn_data.loc[aggregated_cn_data['seg.mean'] == np.NINF, 'seg.mean'] = np.exp(-8)
    aggregated_cn_data = _correct_seg_bin_ends(aggregated_cn_data)
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



def _correct_seg_bin_ends(data):
    data.loc[(data['chrom'] == '1') & (data['loc.end'] == 249500000), 'loc.end'] = 249250621
    data.loc[(data['chrom'] == '2') & (data['loc.end'] == 243500000), 'loc.end'] = 243199373
    data.loc[(data['chrom'] == '3') & (data['loc.end'] == 198500000), 'loc.end'] = 198022430
    data.loc[(data['chrom'] == '4') & (data['loc.end'] == 191500000), 'loc.end'] = 191154276
    data.loc[(data['chrom'] == '5') & (data['loc.end'] == 181000000), 'loc.end'] = 180915260
    data.loc[(data['chrom'] == '6') & (data['loc.end'] == 171500000), 'loc.end'] = 171115067
    data.loc[(data['chrom'] == '7') & (data['loc.end'] == 159500000), 'loc.end'] = 159138663
    data.loc[(data['chrom'] == '8') & (data['loc.end'] == 146500000), 'loc.end'] = 146364022
    data.loc[(data['chrom'] == '9') & (data['loc.end'] == 141500000), 'loc.end'] = 141213431
    data.loc[(data['chrom'] == '10') & (data['loc.end'] == 136000000), 'loc.end'] = 135534747
    data.loc[(data['chrom'] == '11') & (data['loc.end'] == 135500000), 'loc.end'] = 135006516
    data.loc[(data['chrom'] == '12') & (data['loc.end'] == 134000000), 'loc.end'] = 133851895
    data.loc[(data['chrom'] == '13') & (data['loc.end'] == 115500000), 'loc.end'] = 115169878
    data.loc[(data['chrom'] == '14') & (data['loc.end'] == 107500000), 'loc.end'] = 107349540
    data.loc[(data['chrom'] == '15') & (data['loc.end'] == 103000000), 'loc.end'] = 102531392
    data.loc[(data['chrom'] == '16') & (data['loc.end'] == 90500000), 'loc.end'] = 90354753
    data.loc[(data['chrom'] == '17') & (data['loc.end'] == 81500000), 'loc.end'] = 81195210
    data.loc[(data['chrom'] == '18') & (data['loc.end'] == 78500000), 'loc.end'] = 78077248
    data.loc[(data['chrom'] == '19') & (data['loc.end'] == 59500000), 'loc.end'] = 59128983
    data.loc[(data['chrom'] == '20') & (data['loc.end'] == 63500000), 'loc.end'] = 63025520
    data.loc[(data['chrom'] == '21') & (data['loc.end'] == 48500000), 'loc.end'] = 48129895
    data.loc[(data['chrom'] == '22') & (data['loc.end'] == 51500000), 'loc.end'] = 51304566
    data.loc[(data['chrom'] == 'X') & (data['loc.end'] == 155500000), 'loc.end'] = 155270560
    data.loc[(data['chrom'] == 'Y') & (data['loc.end'] == 59500000), 'loc.end'] = 59373566
    return data