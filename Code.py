from scipy.io import mmread
import numpy as np

a = mmread('matrix.mtx')

num_gen = a.shape[1]
#num_gen = 100
num_cel = a.shape[0]
fwer = 0.05

phi = list()
xg = [a.data[a.col ==gen].sum()/num_cel for gen in range(num_gen)]

x_b= np.mean(xg)
ex_b= np.exp(-x_b)


for gen in range(num_gen):
    #col_gen = a.getcol(gen).toarray().flatten()
    #pg = len(np.where(col_gen==0)[0])/num_cel
    pg = a.shape[0] - len(a.data[a.col == gen])/num_cel

    zg = (pg - ex_b)/(pg*(1-pg)/num_cel)

    if pg * num_gen < fwer:
        phi.append(gen)

print(phi)
pass