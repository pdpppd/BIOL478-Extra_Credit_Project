from DNASequencingPipeline import *
import matplotlib.pyplot as plt
from itertools import product
from tqdm import tqdm

# Assuming you have the necessary functions and libraries imported

def run_pipeline(kmer_length, error_rate, coverage, num_repeats, repeat_length):
    # Generate the DNA sequence
    dna_sequence = generate_dna_sequence(300, repeat_length=repeat_length, num_repeats=num_repeats, repeat_type='tandem')

    # Generate reads
    simulated_reads = generate_reads(dna_sequence, read_length=20, coverage=coverage, error_rate=error_rate)

    # Generate kmers
    kmers = generate_kmers(simulated_reads, kmer_length=kmer_length)
    kmers = set(kmers)

    # Construct the De Bruijn graph and assemble the sequence
    G = debruijnize(kmers)
    node_edge_map = make_node_edge_map(G[1])
    if not G[0]:
        # Return a special value or handle the serror appropriately
        return float('inf')
    start = G[2][0] if G[2] else list(G[0])[0]
    trail = eulerian_trail(node_edge_map, start)
    assembled_sequence = assemble_trail(trail)

    # Calculate the Levenshtein distance
    distance = levenshtein_distance(dna_sequence, assembled_sequence)
    return distance


# Base case parameters
base_kmer_length = 10
base_error_rate = 0.01
base_coverage = 300
base_num_repeats = 2
base_repeat_length = 10

# Parameter ranges
kmer_lengths = range(5, 36)  # 5 to 35
error_rates = [i * 0.005 for i in range(31)]  # 0% to 15% in 0.5% increments
coverages = range(10, 1011, 33)  # 10x to 1000x
num_repeats_list = range(1, 31)  # 1 to 30
repeat_lengths = range(5, 61, 2)  # 5 to 60

# Function to run experiments and plot results
def experiment_and_plot(variable_param, variable_values):
    distances = []
    for value in tqdm(variable_values, desc=f"Testing {variable_param}"):
        # Update the variable parameter
        params = {
            'kmer_length': base_kmer_length,
            'error_rate': base_error_rate,
            'coverage': base_coverage,
            'num_repeats': base_num_repeats,
            'repeat_length': base_repeat_length
        }
        params[variable_param] = value

        # Run the pipeline
        distance = run_pipeline(**params)
        distances.append(distance)

    # Plotting
    plt.plot(variable_values, distances, marker='o')
    plt.xlabel(variable_param)
    plt.ylabel('Levenshtein Distance')
    plt.title(f'Levenshtein Distance vs {variable_param}')
    plt.show()

# Running experiments for each parameter
experiment_and_plot('kmer_length', kmer_lengths)
experiment_and_plot('error_rate', error_rates)
experiment_and_plot('coverage', coverages)
experiment_and_plot('num_repeats', num_repeats_list)
experiment_and_plot('repeat_length', repeat_lengths)