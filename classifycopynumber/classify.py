import wgs_analysis.algorithms.cnv
import numpy as np


def calculate_amp_percentile(cn, gene_cn):
    """ Calculate the percentile of each gene wrt the cumulative distribution
    of total copy number across the genome.

    Args:
        cn (pandas.DataFrame): segment copy number
        gene_cn (pandas.DataFrame): gene copy number

    Returns:
        pandas.DataFrame: per gene amp percentile
    """

    # Calculate cumulative distribution of total copy number
    cn = cn.copy()
    cn['segment_length'] = cn['end'] - cn['start']
    cn = cn[['total_raw', 'segment_length']].dropna().sort_values('total_raw')
    cn['fraction_genome'] = cn['segment_length'] / cn['segment_length'].sum()
    cn['cum_fraction_genome'] = 1 - cn['fraction_genome'].cumsum()

    # Search cumulative distribution for percentile of each genes total copy number
    gene_cn = gene_cn[['gene_id', 'total_raw_mean']].dropna().drop_duplicates()
    gene_raw_mean_ix = cn['total_raw'].searchsorted(gene_cn['total_raw_mean']) - 1 # [1, ..., n] - 1
    gene_cn['amp_percentile'] = cn['cum_fraction_genome'].values[gene_raw_mean_ix]

    return gene_cn[['gene_id', 'amp_percentile']]


def calculate_mean_cn(gene_cn, cn_cols):
    """ Calculate log2 change over ploidy in mean raw copy number

    Args:
        gene_cn (pandas.DataFrame): gene copy number
        cn_cols (list): list of cn columns to calculate mean for

    Returns:
        pandas.DataFrame: per gene log change over ploidy
    """

    gene_cn = gene_cn.copy()
    normalize = (
        gene_cn.groupby('gene_id')['overlap_width']
        .sum().rename('sum_overlap_width').reset_index())

    weighted_cols = []
    for col in cn_cols:
        gene_cn[f'{col}_weighted'] = gene_cn[col] * gene_cn['overlap_width']
        weighted_cols.append(f'{col}_weighted')

    gene_cn = gene_cn.groupby(['gene_id'])[weighted_cols].sum().reset_index()
    gene_cn = gene_cn.merge(normalize)

    mean_cols = []
    for col in cn_cols:
        gene_cn[f'{col}_mean'] = gene_cn[f'{col}_weighted'] / gene_cn['sum_overlap_width']
        mean_cols.append(f'{col}_mean')

    gene_cn = gene_cn[['gene_id', 'sum_overlap_width'] + mean_cols]

    return gene_cn


def calculate_hdel_width(gene_cn, hdel_cn_threshold=0.5):
    """ Calculate total length of overlapping hdel segments

    Args:
        gene_cn (pandas.DataFrame): gene copy number

    KwArgs:
        hdel_threshold (float): threshold on total copy number for homozygous deletion

    Returns:
        pandas.DataFrame: per gene log change over ploidy
    """

    hdel_data = gene_cn[gene_cn['total_raw'] < hdel_cn_threshold].copy()
    hdel_data = hdel_data.groupby(['gene_id'])['overlap_width'].sum().rename('hdel_width').reset_index()
    hdel_data['hdel_width'] = hdel_data['hdel_width'].fillna(0).astype(int)

    return hdel_data


def classify_cn_change(
        cn, genes,
        hlamp_percentile_threshold=0.02,
        amp_log_change_threshold=1,
        hdel_overlap_threshold=10000,
        loh_mean_cn_threshold=0.5,
        del_log_change_threshold=-1,
        hdel_cn_threshold=0.5,
    ):
    """ Classify CN changes and calculate gistic values

    Args:
        cn (pandas.DataFrame): table of copy number values
        genes (pandas.DataFrame): gene regions of interest

    KwArgs: filtering thresholds

    """
    cn = cn[np.isfinite(cn['total_raw'])]
    cn = cn[cn['length'] > 1e5] #
    cn = cn[cn['minor_readcount'] > 100] #
    cn['chromosome'] = cn['chromosome'].str.replace('chr', '') #

    # Calculate ploidy
    segment_length = cn['end'] - cn['start']
    ploidy = (cn['total_raw'] * segment_length).sum() / segment_length.sum()

    # Calculate gene copy number overlaps
    cn_cols = ['total_raw', 'minor_raw']
    if 'minor_raw' not in cn:
        cn_cols = ['total_raw']
    gene_cn = wgs_analysis.algorithms.cnv.calculate_gene_copy(cn, genes, cn_cols)

    # Calculate mean copy number across each gene
    mean_cn = calculate_mean_cn(gene_cn, cn_cols)
    mean_cn['log_change'] = np.log2(mean_cn['total_raw_mean'] / ploidy)

    # Call LOH from minor CN if possible
    if 'minor_raw_mean' in mean_cn:
        mean_cn['is_loh'] = mean_cn['minor_raw_mean'] < loh_mean_cn_threshold

    # Call HL amps as amplified type and percentile threshold
    amp_percentile = calculate_amp_percentile(cn, mean_cn)
    amp_percentile = amp_percentile.merge(genes[['gene_id', 'amplification_type']], how='right')
    amp_percentile['is_hlamp'] = (
        #amp_percentile['amplification_type'] &
        (amp_percentile['amp_percentile'] < hlamp_percentile_threshold))
    amp_percentile = amp_percentile.drop(['amplification_type'], axis=1)

    # Call HDels as deletion type and hdel overlap threshold
    hdel_width = calculate_hdel_width(gene_cn, hdel_cn_threshold=hdel_cn_threshold)
    hdel_width = hdel_width.merge(genes[['gene_id', 'deletion_type']], how='right')
    hdel_width['hdel_width'] = hdel_width['hdel_width'].fillna(0)
    hdel_width['is_hdel'] = (
        #hdel_width['deletion_type'] &
        (hdel_width['hdel_width'] > hdel_overlap_threshold))
    hdel_width = hdel_width.drop(['deletion_type'], axis=1)

    # Compile list of CN change
    cn_change = (
        genes.merge(mean_cn, on='gene_id', how='left')
        .merge(amp_percentile, on='gene_id', how='left')
        .merge(hdel_width, on='gene_id', how='left'))

    # Default state
    cn_change['gistic_value'] = 0

    # Amplified as log change greater than 1
    cn_change.loc[cn_change['log_change'] > amp_log_change_threshold, 'gistic_value'] = 1

    # Deleted state as LOH if possible otherwise log change threshold
    if 'is_loh' in cn_change:
        cn_change.loc[cn_change['is_loh'].fillna(False), 'gistic_value'] = -1
        cn_change['has_loh'] = True
    else:
        cn_change.loc[cn_change['log_change'] < del_log_change_threshold, 'gistic_value'] = -1
        cn_change['has_loh'] = False

    # High level amplified state from amp call
    cn_change.loc[cn_change['is_hlamp'].fillna(False), 'gistic_value'] = 2

    # Homozygous deleted state from amp call
    cn_change.loc[cn_change['is_hdel'].fillna(False), 'gistic_value'] = -2

    return cn_change

