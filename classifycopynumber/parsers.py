import pandas as pd
import numpy as np

import classifycopynumber.transformations


def read_remixt(filename, max_ploidy=None, min_ploidy=None, max_divergence=0.5):
    with pd.HDFStore(filename) as store:
        stats = store['stats']
        stats = stats[stats['proportion_divergent'] < max_divergence]

        if max_ploidy is not None:
            stats = stats[stats['ploidy'] < max_ploidy]

        if min_ploidy is not None:
            stats = stats[stats['ploidy'] > min_ploidy]

        if len(stats.index) == 0:
            raise ValueError(f'failed to select correct ploidy for {filename}')

        stats = stats.sort_values('elbo').iloc[-1]
        
        init_id = stats['init_id']

        cn = store[f'/solutions/solution_{init_id}/cn']
        cn['segment_length'] = cn['end'] - cn['start'] + 1
        cn['length_ratio'] = cn['length'] / cn['segment_length']
        cn['total_raw'] = cn['major_raw'] + cn['minor_raw']
        cn['total_raw_e'] = cn['major_raw_e'] + cn['minor_raw_e']
        
        filtered_cn = cn.query('length > 100000')
        raw_mean_sq_err = ((filtered_cn['total_raw_e'] - filtered_cn['total_raw']) ** 2).mean()
        readcount_mean_sq_err = ((filtered_cn['readcount'] - filtered_cn['total_e']) ** 2).mean()
        raw_integer_mean_sq_err = ((filtered_cn['total_raw'] - filtered_cn['total_raw'].round()) ** 2).mean()

        mix = store[f'/solutions/solution_{init_id}/mix']

        stats['normal_proportion'] = mix[0]
        stats['raw_mean_sq_err'] = raw_mean_sq_err
        stats['readcount_mean_sq_err'] = readcount_mean_sq_err
        stats['raw_integer_mean_sq_err'] = raw_integer_mean_sq_err

        return cn, stats


def read_hmmcopy(filename, filter_normal=False):
    """ Read hmmcopy data, filter normal cells and aggregate into segments
    """
    data = pd.read_csv(
        filename,
        usecols=['chr', 'start', 'end', 'width', 'state', 'copy', 'reads', 'cell_id'],
        dtype={'chr': 'category', 'cell_id': 'category'})

    # Filter normal cells that are approximately diploid
    if filter_normal:
        cell_stats = (
            data[data['chr'].isin(autosomes)]
            .groupby('cell_id')['state']
            .agg(['mean', 'std'])
            .reset_index())

        cell_stats['is_normal'] = (
            (cell_stats['mean'] > 1.95) &
            (cell_stats['mean'] < 2.05) &
            (cell_stats['std'] < 0.01))

        data = data.merge(cell_stats[['cell_id', 'is_normal']], how='left')

        data = data[~data['is_normal']]

    # Aggregate cell copy number
    data = (
        data
        .groupby(['chr', 'start', 'end', 'width'])
        .agg({'state': 'median', 'copy': np.nanmean, 'reads': 'sum'})
        .reset_index())

    assert not data.duplicated(['chr', 'start', 'end']).any()

    # Aggregate cell copy number
    data = classifycopynumber.transformations.aggregate_adjacent(
        data,
        value_cols=['state'],
        stable_cols=['state'],
        length_normalized_cols=['copy'],
        summed_cols=['reads'],
    )

    data['chr'] = data['chr'].astype(str)

    return data


def read_gene_data(gtf):
    data = pd.read_csv(
        gtf,
        delimiter='\t',
        names=['chromosome', 'gene_start', 'gene_end', 'info'],
        usecols=[0,3,4,8],
        converters={'chromosome': str},
    )

    def extract_info(info):
        info_dict = {}
        for a in info.split('; '):
            k, v = a.split(' ')
            info_dict[k] = v.strip(';').strip('"')
        return info_dict
    
    data['info'] = data['info'].apply(extract_info)
    data['gene_id'] = data['info'].apply(lambda a: a['gene_id'])
    data['gene_name'] = data['info'].apply(lambda a: a['gene_name'])

    data = data.groupby(['chromosome', 'gene_id', 'gene_name']).agg({'gene_start':'min', 'gene_end':'max'}).reset_index()

    return data

