import click
import classifycopynumber.classify
import classifycopynumber.parsers
import classifycopynumber.plots
import os

@click.command()
@click.argument('genes_gtf')
@click.argument('results_dir')

@click.option('--remixt_h5_filename', help='Remixt H5')
@click.option('--hmmcopy_csv_filename', help='hmmcopy csv')
# @click.option('--results_dir', help='dir for results')
@click.option('--plot', help='bool to generate plots')

def main(genes_gtf, results_dir, remixt_h5_filename=None, hmmcopy_csv_filename=None, plot=False):
    genes = classifycopynumber.parsers.read_gene_data(genes_gtf)

    if remixt_h5_filename is None and hmmcopy_csv_filename is None:
        raise click.ClickException('One of remixt_h5_filename, hmmcopy_csv_filename required')

    if remixt_h5_filename is not None:
        cn, _ = classifycopynumber.parsers.read_remixt(remixt_h5_filename)

    elif hmmcopy_csv_filename is not None:
        cn, _ = classifycopynumber.parsers.read_hmmcopy(hmmcopy_csv_filename)
    print("elkblkfb\n\n\n")
    cn = classifycopynumber.classify.aggregate_cn_data(cn_data)
    genes_of_interest = classifycopynumber.parsers.compile_genes_of_interest(genes)
    gene_cn = classifycopynumber.classify.calculate_gene_cn_data(cn, genes_of_interest)

    amp_genes = classifycopynumber.classify.label_amplifications(gene_cn)
    hdel_genes = classifycopynumber.classify.label_hdels(gene_cn)

    amp_genes.to_csv(os.path.join(results_dir, "amplified_genes.tsv"))
    hdel_genes.to_csv(os.path.join(results_dir, "deleted_genes.tsv"))

    if plots:
        plots_dir = os.path.join(results_dir, "plots")
        
        log_change_plot = os.path.join(plots_dir, "log_change.png")
        classifycopynumber.plots.plot_log_change(amp_genes, log_change_plot)

        log_change_plot = os.path.join(plots_dir, "amp_matrix.png")
        classifycopynumber.plots.plot_amp_matrix(amp_genes, log_change_plot)

        log_change_plot = os.path.join(plots_dir, "hdels.png")
        classifycopynumber.plots.plot_log_change(hdel_genes, log_change_plot)

    from IPython import embed; embed(); raise


# main("/home/mcphersa1/work/apolloanalysis/metadata/Homo_sapiens.GRCh37.73.gtf.gz",  "/juno/work/shah/abramsd/CODE/classif_copynumber/classifycopynumber/tests/results", hmmcopy_csv_filename="/work/shah/tantalus/SC-2640/results/hmmcopy/A96213A_reads.csv.gz", plot=True)
