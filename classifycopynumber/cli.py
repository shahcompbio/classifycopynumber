import click
import classifycopynumber.classify
import classifycopynumber.parsers
import wgs_analysis.algorithms.cnv
import os
import pandas as pd
import yaml


@click.command()
@click.argument('genes_gtf')
@click.argument('cn_change_filename')
@click.option('--remixt_h5_filename', help='Remixt H5')
@click.option('--remixt_parsed_csv', help='Remixt H5')
@click.option('--hmmcopy_csv_filenames', help='hmmcopy csv',  multiple=True)
@click.option('--sample_ids', help='sample ids',  multiple=True)
def main(genes_gtf, cn_change_filename, remixt_h5_filename=None, remixt_parsed_csv=None, hmmcopy_csv_filenames=(), sample_ids=()):
    if remixt_h5_filename is None and len(hmmcopy_csv_filenames) == 0 and remixt_parsed_csv is None:
        raise click.ClickException('One of remixt_h5_filename, hmmcopy_csv_filename required')

    if remixt_parsed_csv is not None:
        cn, stats = classifycopynumber.parsers.read_remixt_parsed_csv(remixt_parsed_csv)

    if remixt_h5_filename is not None:
        cn, stats = classifycopynumber.parsers.read_remixt(remixt_h5_filename)

    elif len(hmmcopy_csv_filenames) > 0:
        if len(sample_ids) == 0:
            raise click.ClickException('sample_ids required if providing hmmcopy')
        cn = classifycopynumber.parsers.read_hmmcopy_files(hmmcopy_csv_filenames, sample_ids=sample_ids)

    genes = classifycopynumber.parsers.read_gene_data(genes_gtf)
    genes_of_interest = classifycopynumber.parsers.compile_genes_of_interest()

    genes = genes.merge(genes_of_interest, how='right')

    if genes['gene_start'].isnull().any():
        raise Exception(genes[genes['gene_start'].isnull()])

    cn_change = classifycopynumber.classify.classify_cn_change(cn, genes)

    cn_change.to_csv(cn_change_filename, index=False)
