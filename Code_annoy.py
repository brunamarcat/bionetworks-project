import numpy as np
from annoy import AnnoyIndex

from matplotlib import pyplot as plt   
import seaborn as sns

from scipy.io import mmread
from scipy.stats import norm
from scipy.sparse import csr_matrix, lil_matrix

# print("Before import")
# import scanpy as sc
# print("After import")

from sklearn.cluster import SpectralClustering, AgglomerativeClustering
from sklearn import manifold

from statsmodels.stats.multitest import multipletests

import pickle
import json 
from JSONEncoder import NumpyArrayEncoder
from Tree import Node


def calculate_split_gens(count_matrix, cel, old_phi, fwer):
    """
    Statistical function to calculate the genes for splitting current cluster of cells based 
    on the gene count differences. Based on paper "Demystifying 'drop-outs' in single-cell UMI data" 
    https://doi.org/10.1186/s13059-020-02096-y
    
    :param count_matrix: Sparse matrix that contains gene expression counts for each cell. 
    :param cel: Set of indexes of cells in current cluster.
    :param old_phi: List of indices representing genes used for previous split.
    :param fwer: Family-wise error rate.

    :return: List of genes indexes for the new split.
    """

    num_cel = len(cel)

    zero_proportion = np.array([1-count_matrix[gen, cel].getnnz(axis=1)/num_cel for gen in old_phi]).flatten()
    gene_mean = np.mean(count_matrix[old_phi, :],axis=1).mean()

    ex_pi = np.array(np.minimum(np.exp(-gene_mean),1-1e-10)).flatten()
    se = np.sqrt(ex_pi*(1-ex_pi))
    ex_pi = np.maximum(1e-10,ex_pi)
    se = se/np.sqrt(num_cel)
    prop_diff = zero_proportion - ex_pi

    zvalue = prop_diff/se
    p_values = 1 - norm.cdf(zvalue, loc=0, scale=1) 

    # Bonferroni correction for multiple testing
    new_phi_ind = np.where(multipletests(p_values, alpha=fwer, method='bonferroni')[0])[0]

    new_phi = np.array(old_phi)[new_phi_ind].tolist()

    return new_phi


def split_cells(count_matrix, phi, cells, method='SpectralClustering', number_neighbors=1000):
    """
    Function to split cells into two clusters using the specified clustering method and the normalized 
    Hamming distance as adjacency matrix.
    
    :param count_matrix: Sparse matrix that contains gene expression counts for each cell. 
    :param phi:  List of indices representing genes to be used for splitting.
    :param cells:  Set of indexes of cells in current cluster.
    :param method: Specifies clustering algorithm to be used. The default
    value is set to 'SpectralClustering'

    :return: Two sets of cell indexes representing the two clusters after splitting.
    """

    # Filter the matrix to keep only the genes in phi
    filt_matrix = count_matrix[phi, :][:,cells].copy()
    # Binarize the matrix
    filt_matrix[filt_matrix != 0] = 1
    
    # Initialize the Annoy index to store the normalized Hamming distance
    t = AnnoyIndex(len(phi), 'hamming')
    for i in range(len(cells)):
        t.add_item(i, filt_matrix[:,i].toarray().T[0])

    t.build(10) # 10 trees
    #t.save('test.ann')
    lim = 100 if len(cells) > 100 else len(cells)                   ## TODO: set minimum number of neighbors
    number_neighbors = np.maximum(lim, len(cells)/10).astype(int)   ##TODO: Adjust number of neighbors
    HS = lil_matrix((len(cells), len(cells)))    
    for i in range(len(cells)):
        row_index = t.get_nns_by_item(i, number_neighbors)
        HS[i, row_index] = 1                                ## TODO: Should be the normalized hamming similitude instead of 1
                                                            ##   (len(cells) - [t.get_distance(i, x) for x in row_index])/n

    # Symmetric matrix to store the normalized Hamming distance
    HS =  (HS + HS.T)/2
    ## TODO: Further normalization of the matrix?

    # Split the cells into two clusters using the specified clustering method
    if method == 'SpectralClustering':
        sc = SpectralClustering(2, affinity='precomputed', assign_labels='discretize', eigen_solver='amg')
        sc.fit_predict(HS)  
        C1 = np.array(cells)[sc.labels_ == 0]
        C2 = np.array(cells)[sc.labels_ == 1]

    elif method == 'AglomerativeClustering':
        model = AgglomerativeClustering(n_clusters=2, linkage="average", metric='precomputed')
        model.fit_predict(1-HS)
        C1 = np.array(cells)[model.labels_ == 0]
        C2 = np.array(cells)[model.labels_ == 1]

    return C1,C2


def split_cells_old(count_matrix, phi, cells, method='SpectralClustering'):
    """
    Function to split cells into two clusters using the specified clustering method and the normalized 
    Hamming distance as adjacency matrix.
    
    :param count_matrix: Sparse matrix that contains gene expression counts for each cell. 
    :param phi:  List of indices representing genes to be used for splitting.
    :param cells:  Set of indexes of cells in current cluster.
    :param method: Specifies clustering algorithm to be used. The default
    value is set to 'SpectralClustering'

    :return: Two sets of cell indexes representing the two clusters after splitting.
    """

    # Filter the matrix to keep only the genes in phi
    filt_matrix = count_matrix[phi, :][:,cells].copy()

    # Binarize the matrix
    filt_matrix[filt_matrix != 0] = 1

    # Make the matrix dense for efficient calculation of Hamming distance
    filt_matrix = filt_matrix.toarray()
    n = filt_matrix.shape[0]

    # Efficient method to calculate the normalized Hamming distance in the binarized matrix
    HS = np.matmul(filt_matrix.T,filt_matrix) 
    filt_matrix = -filt_matrix+1
    HS = (HS + np.matmul(filt_matrix.T,filt_matrix))/n 

    # Check range of values of the normalized Hamming distance matrix
    assert HS.max() <= 1
    assert HS.min() >= 0

    # Split the cells into two clusters using the specified clustering method
    if method == 'SpectralClustering':
        sc = SpectralClustering(2, affinity='precomputed', assign_labels='discretize')
        sc.fit_predict(HS)  
        C1 = np.array(cells)[sc.labels_ == 0]
        C2 = np.array(cells)[sc.labels_ == 1]

    elif method == 'AglomerativeClustering':
        model = AgglomerativeClustering(n_clusters=2, linkage="average", metric='precomputed')
        model.fit_predict(1-HS)
        C1 = np.array(cells)[model.labels_ == 0]
        C2 = np.array(cells)[model.labels_ == 1]

    return C1,C2



def tree_generator(tree, count_matrix, cells, old_phi, fwer,  method='SpectralClustering'):
    """
    Recursive function to split the current set of cells using specified clustering method.
 
    :param tree: It is a Node object representing current node in the clustering tree. 
    :param count_matrix: Sparse matrix that contains gene expression counts for each cell. 
    :param cells:  Set of indexes of cells in current cluster.
    :param old_phi:  List of indices representing genes used for previous split.
    :param fwer: Family-wise error rate.
    :param method: Specifies clustering algorithm to be used. The default
    value is set to 'SpectralClustering'
    """

    global c 

    # Calculate the genes for splitting the current cluster
    new_phi = calculate_split_gens(count_matrix, cells, old_phi, fwer)

    # If no genes are found for splitting, then the current cluster is a leaf node
    if new_phi == []:
        tree.gens_split = old_phi
        tree.cells = cells.tolist()
        tree.type = c
        colors[cells] = c
        c = c + 1
        tree.left = None
        tree.right = None

        print("*** Classe: ", c-1)
        print("  Celules: ", cells.tolist())

        return None
    
    # Split the set of cells into two clusters using the specified clustering method.
    C1,C2 = split_cells(count_matrix, new_phi, cells, method)

    # Expand the tree with the new split
    tree.gens_split = new_phi
    tree.cells = cells
    tree.type = None
    tree.left = Node()
    tree.right = Node()

    # Recursive call to split the left and right clusters
    tree_generator(tree.left, count_matrix, C1, new_phi, fwer, method)
    tree_generator(tree.right, count_matrix,C2, new_phi, fwer, method)

    return tree
    

if __name__ == "__main__":

    fwer = 0.05
    #dataset_name = "E-MTAB-10553.aggregated_filtered_counts.mtx"
    dataset_name = "matrix.mtx"

    ## Load the dataset
    count_matrix = mmread(dataset_name)
    count_matrix = csr_matrix(count_matrix)

    # print("Before load")
    # dataset_name = "E-MTAB-10553"
    # count_matrix = sc.datasets.ebi_expression_atlas(dataset_name)
    # print("After load")
    # count_matrix = count_matrix.X.T
    # print("After transpose")
    # count_matrix = csr_matrix(count_matrix)

    num_genes = count_matrix.shape[0]
    num_cells = count_matrix.shape[1]
    print("Num genes: ", num_genes)
    print("Num cells: ", num_cells)

    # Initialize the cluster array to store the cluster number of each cell
    c = 1
    colors = np.array([0]*num_cells)

    # Initialize the tree with a root node
    tree = Node()

    # Generate the tree using the Spectral Clustering method
    tree = tree_generator(tree, count_matrix, list(range(0,num_cells)), list(range(0,num_genes)), fwer, method='SpectralClustering')

    # Print the number of clusters generated
    print("Num clusters: ", c-1)

    # Save the tree generated in a JSON file
    json_data = json.dumps(tree.to_json(), cls=NumpyArrayEncoder)
    with open('tree.json', 'w') as f:
        f.write(json_data)

    # Save each cell cluster number in a pickle file
    file = open('color.pkl', 'wb')
    pickle.dump(colors, file)
    file.close()

    # Plot the t-SNE visualization of the cells colored by the clusters
    tsne = manifold.TSNE(n_components=2, random_state=0,perplexity=30)
    trans_data = tsne.fit_transform(count_matrix.transpose().toarray()).T
    d = {'x':trans_data[0], 'y':trans_data[1], 'class':[str(i) for i in colors]}
    sns.scatterplot(d, x='x',y='y', legend = 'full', hue='class',palette='Paired',alpha=0.5)
    plt.title('t-SNE visualization of cells colored by the clusters')
    plt.xlabel('t-SNE 1')
    plt.ylabel('t-SNE 2')
    plt.savefig('t-SNE.png')
    plt.show()