from heapq import *
from collections import defaultdict

class Node:
	def __init__(self, node_id, orientation):
		self.node_id = node_id
		self.orientation = orientation
	def __eq__(self, other):
		return (self.node_id == other.node_id) and (self.orientation == other.orientation)
	def __ne__(self, other):
		return not self.__eq__(other)
	def __str__(self):
		return ''.join(['Node(', str(self.node_id), ',', str(self.orientation), ')'])
	def __hash__(self):
		return hash((self.node_id, self.orientation))

class GFAGraph:
	def __init__(self, varying_overlaps=False, filename=None):
		"""
		Construct an empty GFA Graph (assume a de-Bruijn graph).
		varying_overlaps -- node overlaps can vary
		filename -- read graph from given file
		"""
		self._graph = defaultdict(set)
		self._nodes = {}
		self._varying_overlaps = varying_overlaps
		self._const_overlap = -1
		self._edge_overlaps = {}
		self._header = None
		self._name_to_id = {}
		self._id_to_name = {}		
		self._next_id = 0
		if filename is not None:
			# read graph from file
			self._read_from_file(filename)
		self._paths = defaultdict(list)

	def _get_id(self, node_name):
		"""
		Given a node name, return its id in the graph.
		node_name -- name of the node
		"""
		if node_name in self._name_to_id:
			node_id = self._name_to_id[node_name]
		else:
			self._name_to_id[node_name] = self._next_id
			self._id_to_name[self._next_id] = node_name
			node_id = self._next_id
			self._next_id += 1
		return node_id

	def _get_name(self, node_id):
		"""
		Given a node id, return its name.
		node_id -- id of the node
		"""
		try:
			return self._id_to_name[node_id]
		except KeyError:
			raise ValueError('No node with ID: ' + str(node_id))

	def get_overlap(self, start_name, start_orientation, end_name, end_orientation):
		"""
		Get node overlap.
		"""
		if not self.has_link(start_name, start_orientation, end_name, end_orientation):
			raise RuntimeError('Link does not exist.')
		if self._varying_overlaps:
			start_node = Node(self._get_id(start_name), start_orientation)
			end_node = Node(self._get_id(end_name), end_orientation)
			return self._edge_overlaps[(start_node, end_node)]
		else:
			return self._const_overlap

	def get_sequence(self, node_name):
		"""
		Get the sequence of a node.
		node_name -- name of the node.
		"""
		if not node_name in self:
			raise RuntimeError('node ' + node_name + ' is not in graph.')
		node_id = self._name_to_id[node_name]
		return self._nodes[node_id]

	def set_header(self, header_line):
		"""
		Add a header.
		"""
		self._header = header_line

	def get_header(self):
		"""
		Get header line.
		"""
		return self._header

	def size(self):
		"""
		Return the number of nodes.
		"""
		return len(self._nodes)

	def add_node(self, node_name, sequence):
		"""
		Add a node to the graph.
		node_id -- node identifier
		sequence -- sequence represented by this node
		"""
		node_id = self._get_id(node_name)
		self._nodes[node_id] = sequence

	def outdegree_of(self, node_name, orientation):
		"""
		Get the outdegree.
		node_name -- name of the node to check.
		orientation -- orientation of the node
		"""
		if node_name not in self:
			raise RuntimeError('Node ' + node_name + ' is not in graph.')
		node_id = self._get_id(node_name)
		node = Node(node_id, orientation)
		return len(self._graph[node])

	def add_link(self, from_name, from_orientation, to_name, to_orientation, overlap):
		"""
		Add a link from "from_id" to "end_id".
		from_id -- start node
		from_orientation -- whether the sequence is reverse complemented
		to_id -- end node
		to_orientation -- whether the end node is reverse complemented
		"""
		if not self._varying_overlaps:
			if (self._const_overlap != -1) and (self._const_overlap != overlap):
				raise RuntimeError('Edge overlaps are not allowed to vary.')
		from_id = self._get_id(from_name)
		to_id = self._get_id(to_name)
		f_n = Node(from_id, from_orientation)
		t_n = Node(to_id, to_orientation)
		self._graph[f_n].add(t_n)
		if self._varying_overlaps:
			self._edge_overlaps[(f_n, t_n)] = overlap
		else:
			self._const_overlap = overlap


	def has_link(self, from_name, from_orientation, to_name, to_orientation):
		"""
		Check whether graph contains the given link.
		from_name -- name of start node
		from_orientation -- orientation of start node
		to_name -- name of end node
		to_orientation -- orientation of end node
		"""
		if not from_name in self:
			return False
		start_id = self._name_to_id[from_name]
		start_node = Node(start_id, from_orientation)
		for end_node in self._graph[start_node]:
			end_name = self._get_name(end_node.node_id)
			end_orientation = end_node.orientation
			if (end_name == to_name) and (end_orientation == to_orientation):
				return True
		return False

	def _read_from_file(self, filename):
		"""
		Read graph from given file.
		filename -- gfa-file to read from
		"""
		reverse = {'+':False, '-':True}
		l = 0
		print('Read graph from file: ' + filename + ' ...')
		for line in open(filename, 'r'):
			l += 1
			if line.startswith('H'):
				fields = line.split()
				if len(fields) > 1:
					self._header = fields[1]
			if line.startswith('S'):
				fields = line.split()
				node_name = fields[1]
				sequence = fields[2]
				self.add_node(node_name, sequence)
				continue
			if line.startswith('L'):
				fields = line.split()
				overlap = int(fields[5].split('M')[0])
				from_orientation = reverse[fields[2]]
				to_orientation = reverse[fields[4]]
				self.add_link(fields[1], from_orientation, fields[3], to_orientation, overlap)
				continue
			if line.startswith('P'):
				fields = line.split()
				path_name = fields[1]
				path_nodes = fields[2].split(',')
				for path_node in path_nodes:
					node_name = path_node[:-1]
					node_orientation = path_node[-1] == '-'
					node_id = self._name_to_id[node_name]
					self._paths[path_name].append(Node(node_id, node_orientation))
		if not self._varying_overlaps:
			assert len(self._edge_overlaps) == 0
		print('Verify ...')
		if not self.valid():
			raise RuntimeError('Read an invalid graph.')
	
	def write_to_file(self, filename):
		"""
		write the graph to a file.
		filename -- name of file to write data to
		"""
#		if not self.valid():
#			raise RuntimeError('Graph format is not valid.')
		reverse = {False:'+', True:'-'}
		with open(filename, 'w') as outfile:
			if self._header is not None:
				header_line = 'H\t' + self._header + '\n'
			else:
				header_line = 'H\n'
			outfile.write(header_line)
			for node_id in self._nodes:
				node_name = self._id_to_name[node_id]
				line = '\t'.join(['S', node_name, self._nodes[node_id]]) + '\n'
				outfile.write(line)
				# write the forward edges
				forward_node = Node(node_id, False)
				for end_node in self._graph[forward_node]:
					end_name = self._id_to_name[end_node.node_id]
					end_orientation = end_node.orientation
					if not self._varying_overlaps:
						cigar = str(self._const_overlap) + 'M'
					else:
						cigar = str(self._edge_overlaps[(forward_node, end_node)]) + 'M'
					line = '\t'.join(['L', node_name, reverse[False], end_name, reverse[end_orientation], cigar]) + '\n'
					outfile.write(line)
				# write reverse edges
				reverse_node = Node(node_id, True)
				for end_node in self._graph[reverse_node]:
					end_name = self._id_to_name[end_node.node_id]
					end_orientation = end_node.orientation
					if not self._varying_overlaps:
						cigar = str(self._const_overlap) + 'M'
					else:
						cigar = str(self._edge_overlaps[(reverse_node, end_node)]) + 'M'
					line = '\t'.join(['L', node_name, reverse[True], end_name, reverse[end_orientation], cigar]) + '\n'
					outfile.write(line)


	def double_edges(self):
		"""
		Also include all reverse edges in the graph (doubles all edges)
		Neccessary to properly construct subgraphs
		""" 
		for node in self._nodes:
			forward_node = Node(node, False)
			reverse_node = Node(node, True)
			for end_node in self._graph[forward_node]:
				rev_end_node = Node(end_node.node_id, not end_node.orientation)
				self._graph[rev_end_node].add(reverse_node)
				if self._varying_overlaps:
					self._edge_overlaps[(rev_end_node, reverse_node)] = self._edge_overlaps[(forward_node, end_node)]
			for end_node in self._graph[reverse_node]:
				rev_end_node = Node(end_node.node_id, not end_node.orientation)
				self._graph[rev_end_node].add(forward_node)
				if self._varying_overlaps:
					self._edge_overlaps[(rev_end_node, forward_node)] = self._edge_overlaps[(reverse_node, end_node)]

	def contains_all(self, node_names):
		"""
		Check whether graph contains all nodes.
		"""
		for node_name in node_names:
			if not node_name in self:
				return False
		return True

	def __contains__(self, node_name):
		"""
		Check whether graph contains a node.
		"""
		if node_name in self._name_to_id:
			node_id = self._name_to_id[node_name]
			if node_id in self._nodes:
				return True
			else:
				return False
		else:
			return False

	def valid(self):
		"""
		Verify the graph. All nodes defined in edges should
		have been defined.
		"""
		for node, links in self._graph.items():
			if node.node_id not in self._nodes.keys():
				return False
			for to_node in links:
				if to_node.node_id not in self._nodes.keys():
					return False
		return True

	# TODO
	def generate_paths(self, primary_name):
		"""
		Traverse the graph and generate paths.
		"""
		start_id = self._paths[primary_name][0].node_id
		end_id = self._paths[primary_name][-1].node_id
		start_name = self._id_to_name[start_id]
		# perform DFS
		nodes_to_process = [start_id, [primary_name, start_name]]
		paths = []
		while nodes_to_process:
			(node, path) = nodes_to_process.pop()
			for end_node in self._graph[node]:
				end_node_id = self._name_to_id[end_node.node_id]
				if end_id == end_node_id:
					paths.append(path)
				else:
					nodes_to_process.append((end_id, [end_node]))
		print(paths)
