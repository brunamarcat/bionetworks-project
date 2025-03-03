import json
from Tree import Node
import numpy as np

with open('tree.json', 'r') as f:
    data_tree = json.load(f)

tree = Node()
tree.from_json(data_tree)

for i in range(0, tree.get_tree_depth() + 1):
    print('Level:', i)
    colors = np.zeros(len(tree.cells))-1
    c = 0
    for k in tree.get_nodes_at_level(i+1):
        colors[k.cells] = c
        c=c+1
        print('Node:', k.type)
    assert np.all(colors!=-1)
    # Analyze the colors at this depth



print('enter the node you want to search for:')