import argparse
from gfa import GFAGraph

parser = argparse.ArgumentParser(prog='generate-paths.py', description=__doc__)
parser.add_argument('gfa', metavar='GFA', help='Input graph (gfa)')
parser.add_argument('name', metavar='NAME', help='Name of the primary sequence')
prgs = parser.parse_args()

# read the graph
graph = GFAGraph(varying_overlaps=True, filename=args.gfa)

# generate the paths
graph.generate_paths(primary_name=args.name)

