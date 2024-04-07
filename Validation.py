import scanpy as sc

import pickle

import pandas as pd

import gseapy as gp
from gseapy import Biomart, barplot


adata = sc.datasets.ebi_expression_atlas("E-MTAB-10553")

# Load clustering of cells
file = open('color_nou.pkl', 'rb')
color = pickle.load(file)
file.close()

# Adding our new cluster labels to the anotated data
adata.obs["Class_own"]= pd.Categorical(color)

# Differential expression analysis with the new annotation
sc.tl.rank_genes_groups(adata, groupby='Class_own', method='wilcoxon')

# Access the results of differential expression analysis
result = adata.uns['rank_genes_groups']

# Extract marker genes for a specific cluster
marker_genes ={}
for i in set(color):
    marker_genes[i] = result['names'][i]

gene_sets=['GO_Biological_Process_2023','GO_Cellular_Component_2023','GO_Molecular_Function_2023','KEGG_2021_Human']

for cluster in set(color):
    bm = Biomart()
    queries ={'ensembl_gene_id': list(marker_genes[cluster]) } # need to be a dict object
    results = bm.query(dataset='hsapiens_gene_ensembl',
                    attributes=['ensembl_gene_id', 'external_gene_name', 'entrezgene_id', 'go_id'],
                    filters=queries)

    gene_list = list(set(results["external_gene_name"].dropna()))

    enr = gp.enrichr(gene_list=gene_list,
                    gene_sets=gene_sets,
                    organism='human',
                    outdir=None,
                    )

    ax = barplot(enr.results,
                column="Adjusted P-value",
                group='Gene_set',
                size=10,
                top_term=5,
                figsize=(3,5),
                color = {'GO_Cellular_Component_2023': 'salmon', 'GO_Biological_Process_2023':'orange','GO_Molecular_Function_2023':'darkblue','KEGG_2021_Human':'lightblue'},
                title = 'Over-representation analysis of cluster ' + str(cluster),
                ofname = 'ORA_'+ str(cluster) + '.png'
                )