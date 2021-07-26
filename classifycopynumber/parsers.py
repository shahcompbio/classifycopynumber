import pandas as pd
import numpy as np
import os
import classifycopynumber.transformations
import pkg_resources
import yaml
from scgenome.utils import concat_with_categories


autosomes = [str(a) for a in range(1, 23)]


def read_remixt_parsed_csv(filename):
    cn = pd.read_csv(filename, sep='\t', converters={'chromosome': str}, low_memory=False)
    cn['total_raw'] = cn['major_raw'] + cn['minor_raw']

    metadata_filename = os.path.join(os.path.dirname(filename), 'meta.yaml')
    with open(metadata_filename, 'r') as f:
        metadata = yaml.load(f, Loader=yaml.FullLoader)
    ploidy = metadata['ploidy']

    return cn, {'ploidy':ploidy}


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

        cn['total_raw'] = cn['major_raw'] + cn['minor_raw']

        return cn, stats


def read_hmmcopy_files(filenames, filter_normal=False, group_label_col='cell_id', sample_ids=None):
    """ Read hmmcopy data, filter normal cells and aggregate into segments
    """
    dfs=[]
    for filename in filenames:
        data = pd.read_csv(
            filename,
            usecols=['chr', 'start', 'end', 'width', 'state', 'copy', 'reads', 'cell_id'],
            dtype={'chr': 'category', 'cell_id': 'category'})
        dfs.append(data)

    if len(dfs) > 1:
        data = concat_with_categories(dfs)

    data['sample_id'] = data['cell_id'].str.split('-', expand=True)[0]

    if sample_ids is not None:
        data = data[data['sample_id'].isin(sample_ids)]
        if data.empty:
            raise ValueError('no data after filtering samples')

    # HACK: remove the last segment of each chromosome
    chrom_ends = data.groupby('chr')['end'].max().rename('chrom_end').reset_index()
    data = data.merge(chrom_ends)
    data = data.query('end != chrom_end')
    data = data.drop(['chrom_end'], axis=1)

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
    aggregated_data = {}

    data = (
        data
        .groupby(['chr', 'start', 'end', 'width'], observed=True)
        .agg({'state': np.nanmedian, 'copy': np.nanmean, 'reads': np.nansum})
        .reset_index())

    assert not data.duplicated(['chr', 'start', 'end']).any()

    data = classifycopynumber.transformations.aggregate_adjacent(
        data,
        value_cols=['state'],
        stable_cols=['state'],
        length_normalized_cols=['copy'],
        summed_cols=['reads'],
    )

    ploidy = (data['copy'] * data['width']).sum() / data['width'].sum()

    data['chr'] = data['chr'].astype(str)
    data = data.rename(columns={'chr': 'chromosome'})

    data['total_raw'] = data['copy']

    return data.reset_index(), ploidy.mean()


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
        for a in info.strip(' ').split('; '):
            k, v = a.split(' ')
            info_dict[k] = v.strip(';').strip('"')
        return info_dict

    data['info'] = data['info'].apply(extract_info)
    data['gene_id'] = data['info'].apply(lambda a: a['gene_id'])
    data['gene_name'] = data['info'].apply(lambda a: a['gene_name'])

    data = data.groupby(['chromosome', 'gene_id', 'gene_name']).agg({'gene_start':'min', 'gene_end':'max'}).reset_index()

    return data


def _get_gene_lists():
    meta = os.path.join(os.path.dirname(__file__), 'metadata')

    return {
        'amp': pkg_resources.resource_stream(__name__, 'metadata/census_amps.csv'),  
        'del': pkg_resources.resource_stream(__name__, 'metadata/census_dels.csv'),
        'cancer_gene_census': pkg_resources.resource_stream(__name__, 'metadata/cancer_gene_census.csv'),
        'additional_genes': pkg_resources.resource_stream(__name__, 'metadata/additional_genes.csv'),
        'antigen_genes': pkg_resources.resource_stream(__name__, 'metadata/antigen_presenting_genes.csv'),
        'hr_genes': pkg_resources.resource_stream(__name__, 'metadata/hr_genes.txt'),
    }


default_additional_gene_lists = ('additional_genes', 'antigen_genes', 'hr_genes')


def compile_genes_of_interest(additional_gene_lists=default_additional_gene_lists):
    """ Compile a list of genes of interest.

    KwArgs:
        additional (tuple): list of additional gene sets.

    Returns:
        pandas.DataFrame: table of genes
    """

    gene_lists = _get_gene_lists()

    genes = []

    amp_genes = pd.read_csv(gene_lists['amp'], usecols=['Gene Symbol']).rename(columns={'Gene Symbol': 'gene_name'})
    amp_genes['cn_type'] = 'amplification'
    genes.append(amp_genes)

    del_genes = pd.read_csv(gene_lists['del'], usecols=['Gene Symbol']).rename(columns={'Gene Symbol': 'gene_name'})
    del_genes['cn_type'] = 'deletion'
    genes.append(del_genes)

    cgc_genes = pd.read_csv(gene_lists['cancer_gene_census'], usecols=['Gene Symbol']).rename(columns={'Gene Symbol': 'gene_name'})
    cgc_genes['cn_type'] = 'unspecified'
    genes.append(cgc_genes)

    for additional in additional_gene_lists:
        genes.append(pd.read_csv(gene_lists[additional]))

    genes = pd.concat(genes, ignore_index=True)

    genes['cn_type'] = genes['cn_type'].fillna('unspecified')
    genes = genes.drop_duplicates()

    # HACK: remove IG, TCR genes
    ig_genes = [
        'IGH',
        'IGK',
        'IGL',
        'TRA',
        'TRB',
        'TRD',
    ]
    genes = genes[~genes['gene_name'].isin(ig_genes)]

    gene_rename_file = pkg_resources.resource_stream(__name__, 'metadata/renames.csv')
    gene_rename = pd.read_csv(gene_rename_file).set_index('gene_name')['gtf_gene_name'].to_dict()
    genes['gene_name'] = genes['gene_name'].apply(lambda a: gene_rename.get(a, a))

    # Reshape to one gene per row
    genes = genes.set_index(['gene_name', 'cn_type']).assign(indicator=1)['indicator'].unstack(fill_value=0).reset_index()
    genes['amplification_type'] = genes['amplification'] == 1
    genes['deletion_type'] = genes['deletion'] == 1
    genes = genes[['gene_name', 'amplification_type', 'deletion_type']]

    return genes
