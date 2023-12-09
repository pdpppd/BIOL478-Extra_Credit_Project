# DNA Sequence Analysis Pipeline

This pipeline simulates the process of DNA sequencing and assembly. It includes generating a DNA sequence, simulating the generation of reads from this sequence, constructing a de Bruijn graph from the reads, and assembling the sequence back from the graph.

## Features

- DNA sequence generation with customizable length and nucleotide frequency.
- Simulating read generation with adjustable read length, coverage, and error rate.
- K-mer generation from simulated reads.
- Construction of a de Bruijn graph from the k-mers.
- Finding an Eulerian trail in the graph for sequence assembly.
- Visualizing the de Bruijn graph.
- Calculating the Levenshtein distance between the original and assembled sequences.

## Installation

Before using the pipeline, you need to install the necessary Python libraries. 

1. Clone the repository or download the source code.
3. Run the following command to install the required libraries:

   ```bash
   pip install -r requirements.txt


## Usage 

- Open pipeline.py in an IDE of choice and see the comments
- You can modify the parameters as you like and then run the program to construct the graph and compare the initial seq and assembled seq using leveschtein distance. 

## Output

- The pipeline outputs the original DNA sequence, the reads, the unique k-mers, the assembled sequence, and the Levenshtein distance between the original and assembled sequences.
- It also outputs an image of the debruijn graph. 

## Visualization

- The graph visualized with graphviz is saved as an image
- The graph visualized with pyvis will be saved as a html page
- The graph visualized with matplotlib will automatically appear

## Evaluation

- testing.py is used to run experiments with this code. It tries different combinations of parameters and evaluates it's effect on accuracy. 
- A PDF document is included for submission purposes. 