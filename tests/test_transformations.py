from classifycopynumber.transformations import calculate_log_change,aggregate_adjacent
import pkg_resources
import pandas as pd
import numpy as np

def get_testdata():
    return {"log_change_testdata": pkg_resources.resource_stream(__name__, 'testdata/hmmcopy_A90554A_gene_cn.csv'), 
    "aggregate_adjacent_testdata":pkg_resources.resource_stream(__name__, 'testdata/A90554A_reads.csv')}

def test_calculate_log_change(gene_cn="default"):
    if gene_cn=="default":
        gene_cn = get_testdata()["log_change_testdata"]
    gene_cn = pd.read_csv(gene_cn)
    ploidy = 10

    log_gene_cn = calculate_log_change(gene_cn, ploidy)

    assert "log_change" in log_gene_cn.columns

def setup_test_aggregate_adjacent(filename):
    data = pd.read_csv(
    filename,
    usecols=['chr', 'start', 'end', 'width', 'state', 'copy', 'reads', 'cell_id'],
    dtype={'chr': 'category', 'cell_id': 'category'})

    data = (
        data
        .groupby([ 'chr', 'start', 'end', 'width'])
        .agg({'state': 'median', 'copy': np.nanmean, 'reads': 'sum'})
        .reset_index())
    return data


def test_aggregate_adjacent(hmmcopy_file="default"):
    if hmmcopy_file == "default":
        hmmcopy_file = get_testdata()["aggregate_adjacent_testdata"]

    unaggregated = setup_test_aggregate_adjacent(hmmcopy_file)
    aggregated = aggregate_adjacent(
        unaggregated, 
        value_cols=['state'],
        stable_cols=['state'],
        length_normalized_cols=['copy'],
        summed_cols=['reads'],
    )
    def check_states(cnv_chr):
        assert 0 not in cnv_chr.state.diff().tolist()
        
    assert all(aggregated.groupby("chr").apply(check_states))

test_calculate_log_change()