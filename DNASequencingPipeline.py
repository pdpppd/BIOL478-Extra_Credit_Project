import random
import matplotlib.pyplot as plt
import networkx as nx
import math
from pyvis.network import Network
import matplotlib.pyplot as plt
from itertools import product
from graphviz import Digraph

def generate_dna_sequence(length, freq=None, repeat_length=None, num_repeats=None, repeat_type='tandem'):
    """
    Generates a random DNA sequence with the option to add repetitive sequences.
    
    Parameters:
    - length (int): The length of the DNA sequence to generate.
    - freq (dict): A dictionary with the frequency of each nucleotide (e.g., {'A': 0.3, 'T': 0.3, 'G': 0.2, 'C': 0.2}).
    - repeat_length (int): The length of the repetitive sequence to add.
    - num_repeats (int): The number of times the repetitive sequence should be added.
    - repeat_type (str): The type of repeats to add ('tandem', 'interspersed', or 'both').
    
    Returns:
    - str: A random DNA sequence.
    """
    
    # If no frequency is provided, assume equal distribution of nucleotides
    if freq is None:
        freq = {'A': 0.25, 'T': 0.25, 'G': 0.25, 'C': 0.25}
    
    # Calculate the total length of the repeat sequence to be inserted
    total_repeat_length = repeat_length * num_repeats if repeat_length and num_repeats else 0
    
    # Generate the basic sequence without the repeats
    base_length = length - total_repeat_length
    nucleotides = ''.join([k * int(v * base_length) for k, v in freq.items()])
    base_sequence = ''.join(random.sample(nucleotides, len(nucleotides)))

    # If we have a repeat sequence to add, generate it
    if repeat_length and num_repeats:
        # Create the repeat sequence
        repeat_sequence = ''.join(random.choices('ATGC', k=repeat_length)) * num_repeats

        if repeat_type == 'tandem':
            # Insert the tandem repeat sequence at a random position within the base sequence
            insert_pos = random.randint(0, len(base_sequence))
            sequence = base_sequence[:insert_pos] + repeat_sequence + base_sequence[insert_pos:]
        
        elif repeat_type == 'interspersed':
            # Insert each repeat sequence at a random position within the base sequence
            sequence = base_sequence
            for _ in range(num_repeats):
                single_repeat_sequence = ''.join(random.choices('ATGC', k=repeat_length))
                insert_pos = random.randint(0, len(sequence))
                sequence = sequence[:insert_pos] + single_repeat_sequence + sequence[insert_pos:]
        
        elif repeat_type == 'both':
            # Insert the tandem repeat sequence at a random position within the base sequence
            insert_pos_tandem = random.randint(0, len(base_sequence))
            sequence = base_sequence[:insert_pos_tandem] + repeat_sequence + base_sequence[insert_pos_tandem:]
            # Additionally, intersperse more repeats within the entire sequence
            for _ in range(num_repeats):
                single_repeat_sequence = ''.join(random.choices('ATGC', k=repeat_length))
                insert_pos_interspersed = random.randint(0, len(sequence))
                sequence = sequence[:insert_pos_interspersed] + single_repeat_sequence + sequence[insert_pos_interspersed:]
    else:
        sequence = base_sequence

    # Return the sequence trimmed to the desired length
    return sequence[:length]

def generate_reads(dna_sequence, read_length, coverage, error_rate):
    sequence_length = len(dna_sequence)
    total_reads = int(coverage * sequence_length / read_length)
    reads = []

    for _ in range(total_reads):
        start = random.randint(0, sequence_length - read_length)
        read = dna_sequence[start:start + read_length]
        read_with_error = introduce_errors(read, error_rate)
        reads.append(read_with_error)

    return reads

def introduce_errors(read, error_rate):
    nucleotides = ['A', 'C', 'G', 'T']
    read_with_error = ""

    for base in read:
        if random.random() < error_rate:
            read_with_error += random.choice([nuc for nuc in nucleotides if nuc != base])
        else:
            read_with_error += base

    return read_with_error

def generate_kmers(reads, kmer_length):
    kmers = []

    for read in reads:
        for i in range(len(read) - kmer_length + 1):
            kmer = read[i:i + kmer_length]
            kmers.append(kmer)

    return kmers

# Function to create a de Bruijn graph from a list of reads (or kmers)
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

# Function to find an Eulerian trail in the graph
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

# Function to assemble the sequence from the Eulerian trail
def assemble_trail(trail):
    if not trail:
        return ""
    result = trail[0][:-1]
    for node in trail:
        result += node[-1]
    return result

# Visualization Functions
def visualize_debruijn_matplotlib(nodes, edges):
    G = nx.DiGraph()
    G.add_nodes_from(nodes)
    G.add_edges_from(edges)
    nx.draw(G, with_labels=True, node_color='lightblue', node_size=2000, font_size=10, font_weight='bold')
    plt.show()

def visualize_debruijn_pyvis(nodes, edges, file_name='debruijn_graph.html'):
    net = Network(bgcolor="#222222", font_color="white", notebook=True)  # Set notebook=True if you are using a Jupyter notebook
    net.add_nodes(nodes)  # Add nodes to the network
    net.add_edges(edges)
    net.show_buttons(filter_=['physics'])  # Add edges to the network
    net.show(file_name)


def visualize_debruijn_graphviz(nodes, edges, file_path='debruijn_graph'):
    dot = Digraph(comment='De Bruijn Graph')
    for node in nodes:
        dot.node(node)
    for edge in edges:
        dot.edge(edge[0], edge[1])
    
    dot.render(file_path, format='png', cleanup=True)
    return file_path + '.png'

# Compute Levenshtein Distance of the Generated and Reconstructed DNA 
def levenshtein_distance(s1, s2):
    # Initialize matrix of zeros
    rows = len(s1)+1
    cols = len(s2)+1
    distance = [[0 for _ in range(cols)] for _ in range(rows)]

    # Populate matrix of distances
    for i in range(1, rows):
        distance[i][0] = i
    for j in range(1, cols):
        distance[0][j] = j

    # Compute Levenshtein distance
    for col in range(1, cols):
        for row in range(1, rows):
            if s1[row-1] == s2[col-1]:
                cost = 0
            else:
                cost = 1
            distance[row][col] = min(distance[row-1][col] + 1,      # deletion
                                     distance[row][col-1] + 1,      # insertion
                                     distance[row-1][col-1] + cost) # substitution

    # Return Levenshtein distance
    return distance[rows-1][cols-1]

def predict_number_of_contigs(N, L, G):
    coverage = N * L / G
    return N * math.exp(-coverage)

#------------------------------------------#
# The Pipeline
dna_sequence = generate_dna_sequence(1000, repeat_length=0, num_repeats=0, repeat_type='tandem')
print(dna_sequence)
print('___')

simulated_reads = generate_reads(dna_sequence, read_length=20, coverage = 10, error_rate = 0.05)
print(f"Generated {len(simulated_reads)} reads")
print('___')
print(simulated_reads)
print('___')

kmers = generate_kmers(simulated_reads, kmer_length=10)
kmers = set(kmers)
print(kmers)
print('___')

G = debruijnize(kmers)
node_edge_map = make_node_edge_map(G[1])
start = G[2][0] if G[2] else list(G[0])[0]
trail = eulerian_trail(node_edge_map, start)
assembled_sequence = assemble_trail(trail)

# Visualizing the graph and displaying the assembled sequence
visualize_debruijn_graphviz(G[0], G[1])
print(assembled_sequence)
print(dna_sequence)
print('___')

print("Distance: ", levenshtein_distance(dna_sequence, assembled_sequence))