import scanpy as sc
import re
from loguru import logger

from core.utils import remove_ids

class Plotting(object):

    """
    A class to deal with the matlib plotting for notebooks and as images. Instead of plotting the graphs we should
    send them to the web page to view.
    """

    def __init__(self, adata):
        """
        Initializes the Plotting class
            adata (AnnData): The AnnData object with single-cell data.
        """


        self.adata = adata

    def remove_missing_genes(self, gene_list, error_message):

        # Get the genes that are not matching from the error message to remove them.
        gene_pattern = r'\[([^]]*)\]'

        genes_to_remove = []

        # Find all words beginning with 'LOC' in the sentence
        matches = re.findall(gene_pattern, error_message)
        if matches:
            genes_to_remove = [item.strip(" '") for item in matches[0].split(",")]

        logger.info(f'removing genes {genes_to_remove}')
        ortho_list = remove_ids(gene_list, genes_to_remove)

        return ortho_list

    def run_dot_plots(self, gene_list, groupby='clusters', use_raw=None):


        print (gene_list)

        if gene_list:
            try:
                # Recalculate dendrogram for the current groupby
                sc.tl.dendrogram(self.adata, groupby=groupby)
                # Produce the dot and matrix plots
                sc.pl.matrixplot(self.adata, gene_list, groupby, dendrogram=True, cmap='Blues', standard_scale='var',
                                 colorbar_title='column scaled\nexpression', use_raw=use_raw)
                sc.pl.dotplot(self.adata, gene_list, groupby, dendrogram=True, use_raw=use_raw)

            except KeyError as e:

                # Extracting keys from the error message
                error_message = str(e)
                gene_list = self.remove_missing_genes(gene_list, error_message)

                if gene_list:
                    sc.pl.matrixplot(self.adata, gene_list, groupby, dendrogram=True, cmap='Blues',
                                     standard_scale='var', colorbar_title='column scaled\nexpression', use_raw=use_raw)
                    sc.pl.dotplot(self.adata, gene_list, groupby, dendrogram=True, use_raw=use_raw)
                else:
                    print("gene_list is ", gene_list, " empty! sorry")

        else:
            print("gene_list is ", gene_list, " empty! sorry")

        return gene_list

    def run_violin_plots(self, gene_list, groupby='clusters', use_raw=None):


        if gene_list:
            try:

                gene_groups = [gene_list[i:i + 3] for i in range(0, len(gene_list), 3)]

                # Plot each group of genes separately
                for gene_group in gene_groups:
                    sc.pl.violin(self.adata, keys=gene_group, groupby=groupby, use_raw=use_raw)

            except KeyError as e:

                # Extracting keys from the error message
                error_message = str(e)
                gene_list = self.remove_missing_genes(gene_list, error_message)

                if gene_list:
                    gene_groups = [gene_list[i:i + 3] for i in range(0, len(gene_list), 3)]
                    for gene_group in gene_groups:
                        sc.pl.violin(self.adata, keys=gene_group, groupby=groupby, use_raw=use_raw)
                else:
                    print("gene_list is ", gene_list, " empty! sorry")

            return gene_list

    def run_umap_plots(self, gene_list, use_raw=None):

        plot_group_size = 3  # Set the group size here

        if gene_list:
            try:
                gene_groups = [gene_list[i:i + plot_group_size] for i in range(0, len(gene_list), plot_group_size)]

                # Plot each group of genes separately
                for gene_group in gene_groups:
                    sc.pl.umap(self.adata, color=gene_group,  use_raw=use_raw)

            except KeyError as e:

                # Extracting keys from the error message
                error_message = str(e)
                # Remove missing genes from error message and return the list
                gene_list = self.remove_missing_genes(self, error_message)

                if gene_list:
                    gene_groups = [gene_list[i:i + plot_group_size] for i in range(0, len(gene_list), plot_group_size)]
                    for gene_group in gene_groups:
                        sc.pl.umap(self.adata, color=gene_group, use_raw=use_raw)
                else:
                    logger.info(f'gene_list is {gene_list} empty! sorry')

        return gene_list