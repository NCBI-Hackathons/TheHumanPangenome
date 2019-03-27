from collections import defaultdict
import logging
import statistics
import sys
import re
import pandas as pd

pack_file_path = sys.argv[1]
gGFF_file_path = sys.argv[2]

node_to_gene_map = {}

node_pat = re.compile(r"([0-9]+)\[([0-9]+):([0-9]+)\]([\+\-\?])")
gene_name_pat = re.compile('gene_name "(.+?)"; ')
pack_df = pd.read_csv(pack_file_path, sep = "\t")

with open(gGFF_file_path) as gGFF_file:
    for line in gGFF_file:
        subgraph, source, f_type, score, phase, attributes = line.rstrip().split("\t")
        nodes = subgraph.split(',')
        gene_name_search = re.search(gene_name_pat, attributes)
        if gene_name_search:
            gene_name = gene_name_search.group(1)
        else:
            logging.warning("Could not find gene name in line: {}".format(line.rstrip()))
            continue
        for node in nodes:
            result = re.fullmatch(node_pat, node)
            if result:
                node_id = int(result.group(1))
                start = int(result.group(2))
                stop = int(result.group(3))
                strand = result.group(4)
                if node_id in node_to_gene_map:
                    logging.warning("Node {} already in map for gene {}, overwriting with {}".format(node_id, node_to_gene_map[node_id], gene_name))
                node_to_gene_map[node_id] = gene_name
            else:
                logging.warning("Could not parse node {}".format(node))

# sum coverage over node IDs for each gene
node_summed_coverage = pack_df.groupby('node.id').sum()['coverage']
node_lengths = pack_df.groupby('node.id').size()
node_avg_cov = node_summed_coverage / node_lengths
# filter to nodes with coverage
node_avg_cov = node_avg_cov.loc[node_avg_cov > 0]

gene_coverage = defaultdict(list)
for node_id in node_avg_cov.index:
    if node_id in node_to_gene_map:
        gene_coverage[node_to_gene_map[node_id]].append(node_avg_cov[node_id])

for gene in gene_coverage:
    print("{}\t{}".format(gene, statistics.mean(gene_coverage[gene])))
