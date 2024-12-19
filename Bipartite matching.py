from networkx.algorithms import bipartite
import networkx as nx
import pickle
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


file = open('TOT-color.pkl', 'rb')
colors = pickle.load(file)
file.close()

adata  = sc.datasets.ebi_expression_atlas("E-MTAB-10553")
color = [str(i) for i in colors]
adata.obs["Class_own"]= pd.Categorical(color)
cells_type = sc.get.obs_df(adata, keys='Factor Value[inferred cell type - ontology labels]')
cell_color = sc.get.obs_df(adata, keys='Class_own')

G = nx.Graph()

# Add nodes with the node attribute "bipartite"
G.add_nodes_from(list(set(cells_type.dropna().values)), bipartite=0)
G.add_nodes_from(list(set((color))), bipartite=1)

for c in set(color):
    idx = np.where(np.array(color) == c)[0]
    for ct in set(cells_type.dropna().values):
        G.add_edge(c,ct,weight = (cells_type[idx] == ct).sum()/len(color))
        #a = np.where(cells_type[:]==ct)
        #b=np.where(np.array(color)==c)
        #G.add_edge(c,ct,weight = (cells_type[idx] == ct).sum()/np.union1d(a,b).shape[0]) #Jaccard distance (intersection/union)
        #G.add_edge(c,ct,weight = (cells_type[idx] == ct).sum()/np.min(a.shape[0],b.shape[0]))

#matching = nx.bipartite.maximum_matching(G)
matching = nx.algorithms.matching.max_weight_matching(G, maxcardinality=True) 
print(matching)


ground_truth = list(set(cells_type.dropna().values))
predicted = list(set((color)))

# Initialize a matrix with NaN values
heatmap_data = pd.DataFrame(
    np.nan, index=ground_truth, columns=predicted
)

#for u, v in matching:
    #if u in ground_truth and v in predicted:
        #heatmap_data.loc[u, v] = G[u][v]['weight']
    #elif v in ground_truth and u in predicted:
        #heatmap_data.loc[v, u] = G[u][v]['weight']

for u in ground_truth:
    for v in predicted:
        heatmap_data.loc[u, v] = G[u][v]['weight']

print(heatmap_data)

plt.figure(figsize=(8, 6))
sns.heatmap(
    heatmap_data,
    annot=True,  
    cmap="YlGnBu",  
    linewidths=0.5,  
    linecolor='gray',  
    cbar_kws={'label': 'Weight'}  
)
plt.title("Heatmap of Matching Weights")
plt.xlabel("Predicted Clusters")
plt.ylabel("Ground Truth Clusters")
plt.show()

