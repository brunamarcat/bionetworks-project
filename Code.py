from scipy.io import mmread
import numpy as np
import time
from scipy.sparse import csr_matrix, csc_matrix, coo_matrix, random

a = mmread('matrix.mtx')
a = csc_matrix(a)

num_gen = a.shape[1]
num_cel = a.shape[0]
fwer = 0.05

def calculate_phi(a, num_cel, num_gen, fwer):
    phi = []
    xg = [a[:, gen].sum()/num_cel for gen in range(num_gen)]

    x_b = np.mean(xg)
    ex_b = np.exp(-x_b)

    for gen in range(num_gen):
        pg = a[:, gen].getnnz(axis=0)/ num_cel
        zg = (pg - ex_b) / (pg * (1 - pg) / num_cel)
        if pg * num_gen < fwer:
            phi.append(gen)
    
    return phi

### Funció antiga. Descomentar per fer la comparació de temps.
# def calculate_phi(a, num_cel, num_gen, fwer):
#     phi = []
#     xg = [a.data[a.col == gen].sum() / num_cel for gen in range(num_gen)]
#
#     x_b = np.mean(xg)
#     ex_b = np.exp(-x_b)
#
#     for gen in range(num_gen):
#         pg = a.shape[0] - len(a.data[a.col == gen]) / num_cel
#         zg = (pg - ex_b) / (pg * (1 - pg) / num_cel)
#         if pg * num_gen < fwer:
#             phi.append(gen)
#
#     return phi


start = time.time()
phi = calculate_phi(a, num_cel, num_gen, fwer)
end = time.time()
print(end - start)

print(phi)


pass