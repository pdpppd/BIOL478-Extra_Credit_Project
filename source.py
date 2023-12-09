# Importing necessary libraries
import random
import matplotlib.pyplot as plt
import networkx as nx
import math
from pyvis.network import Network
from itertools import product
from graphviz import Digraph

# Function to generate a random DNA sequence with optional repeated sequences
def generate_dna_sequence(length, freq=None, repeat_length=None, num_repeats=None, repeat_type='tandem'):
    # Set default nucleotide frequencies to 25% each if none are provided
    if freq is None:
        freq = {'A': 0.25, 'T': 0.25, 'G': 0.25, 'C': 0.25}

    # Determine the total length of all repeat sequences combined
    total_repeat_length = repeat_length * num_repeats if repeat_length and num_repeats else 0

    # Generate the base DNA sequence without the repeats, based on the specified frequencies
    base_length = length - total_repeat_length
    nucleotides = ''.join([k * int(v * base_length) for k, v in freq.items()])
    base_sequence = ''.join(random.sample(nucleotides, len(nucleotides)))

    # Generate and insert repeat sequences, if specified, into the base DNA sequence
    if repeat_length and num_repeats:
        repeat_sequence = ''.join(random.choices('ATGC', k=repeat_length)) * num_repeats
        if repeat_type == 'tandem':
            # Insert the tandem repeat sequence at a random location in the base sequence
            insert_pos = random.randint(0, len(base_sequence))
            sequence = base_sequence[:insert_pos] + repeat_sequence + base_sequence[insert_pos:]
        elif repeat_type == 'interspersed':
            # Insert each repeat sequence at different random positions within the base sequence
            sequence = base_sequence
            for _ in range(num_repeats):
                single_repeat_sequence = ''.join(random.choices('ATGC', k=repeat_length))
                insert_pos = random.randint(0, len(sequence))
                sequence = sequence[:insert_pos] + single_repeat_sequence + sequence[insert_pos:]
        elif repeat_type == 'both':
            # Insert both tandem and interspersed repeats into the base sequence
            insert_pos_tandem = random.randint(0, len(base_sequence))
            sequence = base_sequence[:insert_pos_tandem] + repeat_sequence + base_sequence[insert_pos_tandem:]
            for _ in range(num_repeats):
                single_repeat_sequence = ''.join(random.choices('ATGC', k=repeat_length))
                insert_pos_interspersed = random.randint(0, len(sequence))
                sequence = sequence[:insert_pos_interspersed] + single_repeat_sequence + sequence[insert_pos_interspersed:]
    else:
        sequence = base_sequence

    # Return the final DNA sequence, trimmed to the desired length
    return sequence[:length]

# Function to simulate sequencing reads from the given DNA sequence
def generate_reads(dna_sequence, read_length, coverage, error_rate):
    sequence_length = len(dna_sequence)
    # Calculate the total number of reads needed to achieve the desired coverage
    total_reads = int(coverage * sequence_length / read_length)
    reads = []

    # Generate reads by selecting random subsequences of the DNA sequence
    for _ in range(total_reads):
        start = random.randint(0, sequence_length - read_length)
        read = dna_sequence[start:start + read_length]
        # Introduce random sequencing errors into each read
        read_with_error = introduce_errors(read, error_rate)
        reads.append(read_with_error)

    return reads

# Function to introduce random errors into a DNA read based on a specified error rate
def introduce_errors(read, error_rate):
    nucleotides = ['A', 'C', 'G', 'T']
    read_with_error = ""

    # Iterate over each nucleotide in the read
    for base in read:
        # Randomly introduce an error based on the error rate
        if random.random() < error_rate:
            # Replace the current nucleotide with a different random nucleotide
            read_with_error += random.choice([nuc for nuc in nucleotides if nuc != base])
        else:
            # Keep the original nucleotide if no error is introduced
            read_with_error += base

    return read_with_error

# Function to generate kmers (subsequences of a fixed length) from a set of reads
def generate_kmers(reads, kmer_length):
    kmers = []

    # Extract all possible kmers of the specified length from each read
    for read in reads:
        for i in range(len(read) - kmer_length + 1):
            kmer = read[i:i + kmer_length]
            kmers.append(kmer)

    return kmers

# Function to create a de Bruijn graph from a list of kmers
def debruijnize(reads):
    # Create nodes for the graph, each representing a unique k-1 mer
    nodes = {r[:-1] for r in reads}.union({r[1:] for r in reads})
    # Create edges for the graph, each representing a kmer
    edges = [(r[:-1], r[1:]) for r in reads]
    # Return the nodes, edges, and starting nodes for the graph
    return nodes, edges, list(nodes - {r[1:] for r in reads})

def build_k_mer(string, k):
    # Generate all kmers for a given string and k value
    return [string[i:k + i] for i in range(len(string) - k + 1)]

def make_node_edge_map(edges):
    # Create a mapping from each node to its outgoing edges
    node_edge_map = {}
    for start, end in edges:
        node_edge_map.setdefault(start, []).append(end)
    return node_edge_map

# Function to find an Eulerian trail in a graph
def eulerian_trail(node_edge_map, start_vertex):
    nemap = node_edge_map.copy()
    result_trail = [start_vertex]

    # Construct the Eulerian trail using Fleury's algorithm
    while True:
        trail = []
        previous = start_vertex
        while True:
            # Break the loop if there are no more edges for the current node
            if previous not in nemap:
                break
            # Select an edge to traverse and remove it from the map
            next_node = nemap[previous].pop()
            if len(nemap[previous]) == 0:
                nemap.pop(previous, None)
            trail.append(next_node)
            # Complete a loop if the start vertex is reached
            if next_node == start_vertex:
                break
            previous = next_node

        # Merge the new trail into the main result trail
        index = result_trail.index(start_vertex)
        result_trail = result_trail[:index + 1] + trail + result_trail[index + 1:]
        
        # Find a new starting point for the next trail segment
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

# Function to assemble a DNA sequence from an Eulerian trail of kmers
def assemble_trail(trail):
    if not trail:
        return ""
    # Start with the prefix of the first kmer
    result = trail[0][:-1]
    # Append the last nucleotide of each subsequent kmer
    for node in trail:
        result += node[-1]
    return result

# Visualization Functions
# Function to visualize a de Bruijn graph using matplotlib
def visualize_debruijn_matplotlib(nodes, edges):
    G = nx.DiGraph()
    G.add_nodes_from(nodes)
    G.add_edges_from(edges)
    # Draw the graph with node labels and specified visual properties
    nx.draw(G, with_labels=True, node_color='lightblue', node_size=2000, font_size=10, font_weight='bold')
    plt.show()

# Function to visualize a de Bruijn graph using pyvis (for Jupyter notebooks)
def visualize_debruijn_pyvis(nodes, edges, file_name='debruijn_graph.html'):
    net = Network(bgcolor="#222222", font_color="white", notebook=True)
    net.add_nodes(nodes)
    net.add_edges(edges)
    net.show_buttons(filter_=['physics'])
    net.show(file_name)

# Function to visualize a de Bruijn graph using graphviz
def visualize_debruijn_graphviz(nodes, edges, file_path='debruijn_graph'):
    dot = Digraph(comment='De Bruijn Graph')
    for node in nodes:
        dot.node(node)
    for edge in edges:
        dot.edge(edge[0], edge[1])
    # Render the graph to a file and return the file path
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

