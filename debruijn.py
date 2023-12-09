import matplotlib.pyplot as plt
import networkx as nx
from pyvis.network import Network

def debruijnize(reads):
    nodes = {r[:-1] for r in reads}.union({r[1:] for r in reads})
    edges = [(r[:-1], r[1:]) for r in reads]
    return nodes, edges, list(nodes - {r[1:] for r in reads})

def build_k_mer(string, k):
    return [string[i:k + i] for i in range(len(string) - k + 1)]

def make_node_edge_map(edges):
    node_edge_map = {}
    for start, end in edges:
        node_edge_map.setdefault(start, []).append(end)
    return node_edge_map

def eulerian_trail(node_edge_map, start_vertex):
    nemap = node_edge_map.copy()
    result_trail = [start_vertex]

    # Loop to find the Eulerian trail
    while True:
        trail = []
        previous = start_vertex
        while True:
            if previous not in nemap:
                break
            next_node = nemap[previous].pop()
            if len(nemap[previous]) == 0:
                nemap.pop(previous, None)
            trail.append(next_node)
            if next_node == start_vertex:
                break
            previous = next_node

        # Inserting the found trail into the result trail
        index = result_trail.index(start_vertex)
        result_trail = result_trail[:index + 1] + trail + result_trail[index + 1:]
        
        # Finding new start vertex if there are remaining edges
        if not nemap:
            break
        found_new_start = False
        for n in result_trail:
            if n in nemap:
                start_vertex = n
                found_new_start = True
                break
        if not found_new_start:
            break
    return result_trail


def assemble_trail(trail):
    return trail[0][:-1] + ''.join(node[-1] for node in trail)




# Example usage of the functions
test_string = "HAIIIIIIIIIIIIIIIIIIIILMFDK:AFNLDKNGFKPDALNGDJSGNKJPFNGFSKJFGNLSFKJGNSFKLG"
k = 9
reads = build_k_mer(test_string, k)
# reads = list(set(reads))
print(reads)
G = debruijnize(reads)
node_edge_map = make_node_edge_map(G[1])
start = G[2][0] if G[2] else list(G[0])[0]
trail = eulerian_trail(node_edge_map, start)
assembled_sequence = assemble_trail(trail)

# Visualizing the graph and displaying the assembled sequence
print(assembled_sequence)
