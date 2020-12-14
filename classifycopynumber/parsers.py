import pandas as pd
import numpy as np
import os
import classifycopynumber.transformations
import pkg_resources
import yaml
from scgenome.utils import concat_with_categories

def read_remixt_parsed_csv(filename):
    metadata = os.path.join(os.path.dirname(filename), "meta.yaml")
    cn = pd.read_csv(filename, sep="\t")
    cn['total_raw'] = cn['major_raw'] + cn['minor_raw']
    cn=cn.rename(columns={"total_raw":"copy"})
    with open(metadata, 'r') as f:
        yam = yaml.load(f, Loader=yaml.FullLoader)
    ploidy = yam["ploidy"]
    return cn, {"ploidy":ploidy}


def read_remixt(filename, max_ploidy=None, min_ploidy=None, max_divergence=0.5):
    with pd.HDFStore(filename) as store:
        stats = store['stats']
        brk_cn = store["brk_cn"]

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
        cn=cn.rename(columns={"total_raw":"copy"})

        return cn, stats


def read_hmmcopy_files(filenames, filter_normal=False, group_label_col='cell_id'):
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
        .groupby([ 'chr', 'start', 'end', 'width'])
        .agg({'state': 'median', 'copy': np.nanmean, 'reads': 'sum'})
        .reset_index())

    assert not data.duplicated(['chr', 'start', 'end']).any()

    data = classifycopynumber.transformations.aggregate_adjacent(
            data,
            value_cols=['state'],
            stable_cols=['state'],
            length_normalized_cols=['copy'],
            summed_cols=['reads'],
    )

    data['chr'] = data['chr'].astype(str)
    data = data.rename(columns={"chr":"chromosome"})
    ploidy = (data["copy"] * data["width"]).sum() / data["width"].sum()

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
        for a in info.split('; '):
            k, v = a.split(' ')
            info_dict[k] = v.strip(';').strip('"')
        return info_dict
    
    data['info'] = data['info'].apply(extract_info)
    data['gene_id'] = data['info'].apply(lambda a: a['gene_id'])
    data['gene_name'] = data['info'].apply(lambda a: a['gene_name'])

    data = data.groupby(['chromosome', 'gene_id', 'gene_name']).agg({'gene_start':'min', 'gene_end':'max'}).reset_index()

    return data

def _get_default_genes():

    meta = os.path.join(os.path.dirname(__file__), "metadata")
    return {"amp":pkg_resources.resource_stream(__name__, 'metadata/census_amps.csv'),  
    "del": pkg_resources.resource_stream(__name__, 'metadata/census_dels.csv'),
    "additional_genes": pkg_resources.resource_stream(__name__, 'metadata/additional_genes.csv'),
    "antigen_genes": pkg_resources.resource_stream(__name__, 'metadata/antigen_presenting_genes.csv'),
    "hr_genes": pkg_resources.resource_stream(__name__, 'metadata/hr_genes.txt'),
    }


def compile_genes_of_interest(gene_regions, amp_genes='default', 
    del_genes='default', 
    additional_genes='default', 
    antigen_genes='default',
    hr_genes='default'):

    if amp_genes == "default":
       amp_genes =  _get_default_genes()["amp"]

    if del_genes == "default":
       del_genes =  _get_default_genes()["del"]

    if additional_genes == "default":
       additional_genes =  _get_default_genes()["additional_genes"]

    if antigen_genes == "default":
       antigen_genes =  _get_default_genes()["antigen_genes"]
       
    if hr_genes == "default":
       hr_genes =  _get_default_genes()["hr_genes"]

    amp_genes = amp_genes = pd.read_csv(amp_genes, usecols=['Gene Symbol'])
    del_genes = pd.read_csv(del_genes, usecols=['Gene Symbol'])

    amp_genes['cn_type'] = 'amplification'
    del_genes['cn_type'] = 'deletion'

    cgc_genes = pd.concat([amp_genes, del_genes], ignore_index=False)
    cgc_genes = cgc_genes.rename(columns={'Gene Symbol': 'gene_name'})

    if additional_genes:
        cgc_genes = pd.concat([cgc_genes, pd.read_csv(additional_genes)])
    if hr_genes:
        cgc_genes = pd.concat([cgc_genes, pd.read_csv(hr_genes)])
    if antigen_genes:
        cgc_genes = pd.concat([cgc_genes, pd.read_csv(antigen_genes)])

    cgc_genes["gene_name"] = cgc_genes.gene_name.replace({'NSD3': 'WHSC1L1','MRE11': 'MRE11A','SEM1': 'SHFM1' })

    gene_rename = {
        'NSD3': 'WHSC1L1',
        'MRE11': 'MRE11A',
        'SEM1': 'SHFM1',
    }
    cgc_genes['gene_name'] = cgc_genes['gene_name'].apply(lambda a: gene_rename.get(a, a))


    cgc_genes = cgc_genes.merge(gene_regions, how='left')

    assert not cgc_genes['gene_start'].isnull().any()

    return cgc_genes
