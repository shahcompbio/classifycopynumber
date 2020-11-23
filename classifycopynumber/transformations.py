import numpy as np

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

        a['chr'] = df['chr'].min()
        a['start'] = df['start'].min()
        a['end'] = df['end'].max()
        a['width'] = df['width'].sum()

        for col in length_normalized_cols:
            a[col] = (df[col] * df['width']).sum() / (df['width'].sum() + 1e-16)

        for col in summed_cols:
            a[col] = df[col].sum()

        return a
    print(cnv)
    aggregated = cnv.groupby('cn_group').apply(agg_segments)

    for col in aggregated:
        aggregated[col] = aggregated[col].astype(cnv[col].dtype)

    return aggregated


