from DNASequencingPipeline import *
# PIPELINE FOR DNA SEQUENCE GENERATION, READ SIMULATION, AND ASSEMBLY

# Step 1: Generate a DNA sequence
# - Length: 100 bases (customize as needed)
# - Repeat_length and num_repeats: Set to 0 for no specific repeats
# - Repeat_type: 'tandem' (can be 'tandem', 'interspersed', or 'both')
# - Optional: freq (dictionary) for nucleotide frequencies (A, T, G, C)
# Example freq: {'A': 0.3, 'T': 0.2, 'G': 0.2, 'C': 0.3}
# Uncomment and modify the freq line below to use custom nucleotide frequencies.
# freq = {'A': 0.3, 'T': 0.2, 'G': 0.2, 'C': 0.3}
dna_sequence = generate_dna_sequence(100, repeat_length=0, num_repeats=0, repeat_type='tandem') #, freq=freq)
print("Generated DNA Sequence:")
print(dna_sequence)
print('___')

# Step 2: Generate simulated reads from the DNA sequence
# - Read_length: Length of each read (e.g., 50 bases)
# - Coverage: Average number of times a base is read (e.g., 100x coverage)
# - Error_rate: Probability of a base being read incorrectly (e.g., 0 for no errors)
simulated_reads = generate_reads(dna_sequence, read_length=50, coverage=100, error_rate=0)
print(f"Generated {len(simulated_reads)} reads.")
print('Simulated Reads:')
print(simulated_reads)
print('___')

# Step 3: Generate k-mers from the simulated reads
# - Kmer_length: Length of each k-mer (e.g., 10 bases)
kmers = generate_kmers(simulated_reads, kmer_length=10)
kmers = set(kmers)  # Remove duplicates
print("Unique K-mers generated:")
print(kmers)
print('___')

# Step 4: Construct a de Bruijn graph from the k-mers
G = debruijnize(kmers)
node_edge_map = make_node_edge_map(G[1])
start = G[2][0] if G[2] else list(G[0])[0]

# Step 5: Find an Eulerian trail in the graph and assemble the sequence
trail = eulerian_trail(node_edge_map, start)
assembled_sequence = assemble_trail(trail)

# Visualize the de Bruijn graph (saved to file and optionally displayed)
# There are 3 options for visualization, graphviz will display the clearest graph, pyvis makes a cool interactive graph that you can play around with, matplotlib is just there if the other two do not work
visualize_debruijn_graphviz(G[0], G[1])
#visualize_debruijn_pyvis(G[0], G[1])
#visualize_debruijn_matplotlib(G[0], G[1])
print("Assembled DNA Sequence:")
print(assembled_sequence)
print("Original DNA Sequence:")
print(dna_sequence)
print('___')

# Step 6: Compute the Levenshtein distance between the original and assembled sequences
print("Levenshtein Distance between Original and Assembled Sequences:", levenshtein_distance(dna_sequence, assembled_sequence))
