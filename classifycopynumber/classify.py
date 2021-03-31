import wgs_analysis.algorithms.cnv
import numpy as np


gene_cols = [
    'gene_name',
    'chromosome',
    'gene_id',
    'gene_start',
    'gene_end',
]


def label_amplifications(data, ploidy, min_log_change=1):
    data = data[data['cn_type'] == 'amplification']

    amp_data = data.copy()
    normalize = (
        amp_data.groupby('gene_name')['overlap_width']
        .sum().rename('sum_overlap_width').reset_index())
    amp_data['total_raw_weighted'] = amp_data['total_raw'] * amp_data['overlap_width']
    amp_data = amp_data.groupby(['gene_name'])['total_raw_weighted'].sum().reset_index()
    amp_data = amp_data.merge(normalize)
    amp_data['total_raw_mean'] = amp_data['total_raw_weighted'] / amp_data['sum_overlap_width']
    amp_data['log_change'] = np.log2(amp_data['total_raw_mean'] / ploidy)
    amp_data['pass_filter'] =  amp_data['log_change'] > min_log_change

    cols = [
        'total_raw_mean',
        'log_change',
        'pass_filter',
    ]

    amp_data = amp_data.merge(data[gene_cols].drop_duplicates())[gene_cols + cols]

    return amp_data


def label_deletions(data, ploidy, min_overlap=10000):
    data = data[data['cn_type'] == 'deletion']

    hdel_data = data[data['total_raw'] < 0.5]
    hdel_data = hdel_data.groupby(['gene_name'])['overlap_width'].sum().rename('hdel_width').reset_index()
    hdel_data = hdel_data.merge(data[gene_cols].drop_duplicates(), how='right')
    hdel_data['hdel_width'] = hdel_data['hdel_width'].fillna(0).astype(int)

    hdel_data['pass_filter'] =  hdel_data['hdel_width'] > min_overlap

    cols = [
        'hdel_width',
        'pass_filter',
    ]

    hdel_data = hdel_data[gene_cols + cols]

    return hdel_data

