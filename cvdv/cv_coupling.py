# This code is part of Qiskit.
#
# (C) Copyright IBM 2017, 2019.
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.

"""
Directed graph object for representing coupling between physical qubits.

The nodes of the graph correspond to physical qubits (represented as integers) and the
directed edges indicate which physical qubits are coupled and the permitted direction of
CNOT gates. The object has a distance function that can be used to map quantum circuits
onto a device with this coupling.
"""

import math
from typing import List

import rustworkx as rx
from rustworkx.visualization import graphviz_draw

from qiskit.transpiler.exceptions import CouplingError

from graphviz import Graph
from PIL import Image


class CVCouplingMap:
    """
    Directed graph specifying fixed coupling.

    Nodes correspond to physical qubits (integers) and directed edges correspond
    to permitted CNOT gates, with source and destination corresponding to control
    and target qubits, respectively.
    """

    __slots__ = (
        "description",
        "graph",
        "_dist_matrix",
        "_qubit_list",
        "_qumode_list",
        "_size",
        "_is_symmetric",
    )

    def __init__(self, couplinglist=None, description=None):
        """
        Create coupling graph. By default, the generated coupling has no nodes.

        Args:
            couplinglist (list or None): An initial coupling graph, specified as
                an adjacency list containing couplings, e.g. [[0,1], [0,2], [1,2]].
                It is required that nodes are contiguously indexed starting at 0.
                Missed nodes will be added as isolated nodes in the coupling map.
            description (str): A string to describe the coupling map.
        """
        self.description = description
        # the coupling map graph
        self.graph = rx.PyGraph()
        # a dict of dicts from node pairs to distances
        self._dist_matrix = None
        # a sorted list of physical qubits (integers) in this coupling map
        self._qubit_list = None
        # a sorted list of physical qumodes (integers) in this coupling map
        self._qumode_list = None
        # number of qubits in the graph
        self._size = None
        self._is_symmetric = None

        if couplinglist is not None:
            self.graph.extend_from_edge_list([tuple(x) for x in couplinglist])

    def size(self):
        """Return the number of physical qubits in this graph."""
        if self._size is None:
            self._size = len(self.graph)
        return self._size

    def get_edges(self):
        """
        Gets the list of edges in the coupling graph.

        Returns:
            Tuple(int,int): Each edge is a pair of physical qubits.
        """
        return self.graph.edge_list()

    def __iter__(self):
        return iter(self.graph.edge_list())

    # def add_physical_qubit(self, physical_qubit):
    #     """Add a physical qubit to the coupling graph as a node.

    #     physical_qubit (int): An integer representing a physical qubit.

    #     Raises:
    #         CouplingError: if trying to add duplicate qubit
    #     """
    #     if not isinstance(physical_qubit, int):
    #         raise CouplingError("Physical qubits should be integers.")
    #     if physical_qubit in self.physical_qubits:
    #         raise CouplingError(
    #             f"The physical qubit {physical_qubit} is already in the coupling graph"
    #         )
    #     self.graph.add_node(physical_qubit)
    #     self._dist_matrix = None  # invalidate
    #     self._qubit_list = None  # invalidate
    #     self._size = None  # invalidate

    # def add_physical_qumode(self, physical_qumode):
    #     """Add a physical qumode to the coupling graph as a node.

    #     physical_qumode (int): An integer representing a physical qumode.

    #     Raises:
    #         CouplingError: if trying to add duplicate qumode
    #     """
    #     if not isinstance(physical_qumode, int):
    #         raise CouplingError("Physical qumodes should be integers.")
    #     if physical_qumode in self.physical_qumodes:
    #         raise CouplingError(
    #             f"The physical qumode {physical_qumode} is already in the coupling graph"
    #         )
    #     self.graph.add_node(physical_qumode)
    #     self._dist_matrix = None  # invalidate
    #     self._qumode_list = None  # invalidate
    #     self._size = None  # invalidate


    def add_physical_qubit(self, physical_qubit):
        """Add a physical qubit to the coupling graph as a node.

        physical_qubit (int): An integer representing a physical qubit.

        Raises:
            CouplingError: if trying to add duplicate qubit
        """
        if not isinstance(physical_qubit, int):
            raise CouplingError("Physical qubits should be integers.")
        if physical_qubit in self.physical_qubits:
            raise CouplingError(
                f"The physical qubit {physical_qubit} is already in the coupling graph"
            )
        # Add the node with a type label 'qubit'
        self.graph.add_node({"index": physical_qubit, "type": "qubit"})
        self._dist_matrix = None  # invalidate
        self._qubit_list = None  # invalidate
        self._size = None  # invalidate

    def add_physical_qumode(self, physical_qumode):
        """Add a physical qumode to the coupling graph as a node.

        physical_qumode (int): An integer representing a physical qumode.

        Raises:
            CouplingError: if trying to add duplicate qumode
        """
        if not isinstance(physical_qumode, int):
            raise CouplingError("Physical qumodes should be integers.")
        if physical_qumode in self.physical_qumodes:
            raise CouplingError(
                f"The physical qumode {physical_qumode} is already in the coupling graph"
            )
        # Add the node with a type label 'qumode'
        self.graph.add_node({"index": physical_qumode, "type": "qumode"})
        self._dist_matrix = None  # invalidate
        self._qumode_list = None  # invalidate
        self._size = None  # invalidate



    def add_edge(self, src, dst, edge_type="qumode-qumode", cost=1, weight=1, **attrs):
        """
        Add directed edge to coupling graph.

        Args:
            src (int): source node (qubit or qumode)
            dst (int): destination node (qubit or qumode)
        """
        # Get index of src and dst
        for index, node in enumerate(self.graph.nodes()):
            if node["index"] == src:
                src = index
            elif node["index"] == dst:
                dst = index

        edge_attrs = {"type": edge_type, "weight": weight, "cost": cost}
        edge_attrs.update(attrs)  # Add any additional attributes

        # Add edge to the graph
        self.graph.add_edge(src, dst, edge_attrs)
        self._dist_matrix = None  # invalidate
        self._is_symmetric = None  # invalidate


    # @property
    # def physical_qubits(self):
    #     """Returns a sorted list of physical_qubits"""
    #     if self._qubit_list is None:
    #         self._qubit_list = self.graph.node_indexes()
    #     return self._qubit_list

    @property
    def physical_qubits(self):
        """Returns a sorted list of physical_qubits"""
        if self._qubit_list is None:
            self._qubit_list = sorted(
                node["index"] for node in self.graph.nodes() if node["type"] == "qubit"
            )
        return self._qubit_list

    @property
    def physical_qumodes(self):
        """Returns a sorted list of physical_qumodes"""
        if self._qumode_list is None:
            self._qumode_list = sorted(
                node["index"] for node in self.graph.nodes() if node["type"] == "qumode"
            )
        return self._qumode_list


    def is_connected(self):
        """
        Test if the graph is connected.

        Return True if connected, False otherwise
        """
        try:
            return rx.is_weakly_connected(self.graph)
        except rx.NullGraph:
            return False

    def neighbors(self, physical_qubit):
        """Return the nearest neighbors of a physical qubit.

        Directionality matters, i.e. a neighbor must be reachable
        by going one hop in the direction of an edge.
        """
        return self.graph.neighbors(physical_qubit)

    @property
    def distance_matrix(self):
        """Return the distance matrix for the coupling map.

        For any qubits where there isn't a path available between them the value
        in this position of the distance matrix will be ``math.inf``.
        """
        self.compute_distance_matrix()
        return self._dist_matrix

    def compute_distance_matrix(self):
        """Compute the full distance matrix on pairs of nodes.

        The distance map self._dist_matrix is computed from the graph using
        all_pairs_shortest_path_length. This is normally handled internally
        by the :attr:`~qiskit.transpiler.CouplingMap.distance_matrix`
        attribute or the :meth:`~qiskit.transpiler.CouplingMap.distance` method
        but can be called if you're accessing the distance matrix outside of
        those or want to pre-generate it.
        """
        if self._dist_matrix is None:
            self._dist_matrix = rx.graph_distance_matrix(
                self.graph, null_value=math.inf
            )

    # def distance(self, physical_qubit1, physical_qubit2):
    #     """Returns the undirected distance between physical_qubit1 and physical_qubit2.

    #     Args:
    #         physical_qubit1 (int): A physical qubit
    #         physical_qubit2 (int): Another physical qubit

    #     Returns:
    #         int: The undirected distance

    #     Raises:
    #         CouplingError: if the qubits do not exist in the CouplingMap
    #     """
    #     if physical_qubit1 >= self.size():
    #         raise CouplingError(f"{physical_qubit1} not in coupling graph")
    #     if physical_qubit2 >= self.size():
    #         raise CouplingError(f"{physical_qubit2} not in coupling graph")
    #     self.compute_distance_matrix()
    #     res = self._dist_matrix[physical_qubit1, physical_qubit2]
    #     if res == math.inf:
    #         raise CouplingError(f"No path from {physical_qubit1} to {physical_qubit2}")
    #     return int(res)


    def distance(self, physical_node1, physical_node2):
        """Returns the distance between two physical nodes as the sum of the cost of the edges in the path.
        
        Args:
            physical_node1 (int): A physical qubit or qumode.
            physical_node2 (int): Another physical qubit or qumode.
        
        Returns:
            float: The total cost of the path between the two nodes.
        
        Raises:
            CouplingError: If the nodes do not exist in the coupling map or no path is found.
        """
        # Validate if the nodes exist in the graph
        if physical_node1 >= self.size() or physical_node2 >= self.size():
            raise CouplingError(f"One or both nodes ({physical_node1}, {physical_node2}) are not in the coupling graph.")
        
        # Find the shortest path using `all_shortest_paths` (returns a generator)
        paths_generator = rx.all_shortest_paths(self.graph, physical_node1, physical_node2)

        # Convert the generator to a list of paths
        paths = list(paths_generator)

        if not paths:
            raise CouplingError(f"No path found between {physical_node1} and {physical_node2}")
        
        # Pick the first path (since all_shortest_paths returns multiple paths if there are multiple shortest ones)
        path = paths[0]

        # Now sum the costs along the edges of the path
        total_cost = 0
        for i in range(len(path) - 1):
            src = path[i]
            dst = path[i + 1]
            edge_cost = self.get_edge_cost((src, dst))
            total_cost += edge_cost
        
        return total_cost

    def get_edge_cost(self, edge):
        """Helper function to safely retrieve the cost of an edge, with a default value."""
        edge_data = self.graph.get_edge_data(edge[0], edge[1])
        
        # Ensure edge data exists and has a 'cost' field, defaulting to 1 if not.
        if edge_data is not None:
            return edge_data.get('cost', 1)  # Default cost = 1 if 'cost' is not set
        else:
            return 1  # Default to 1 if no edge data exists


    def shortest_undirected_path(self, physical_qubit1, physical_qubit2):
        """Returns the shortest undirected path between physical_qubit1 and physical_qubit2.

        Args:
            physical_qubit1 (int): A physical qubit
            physical_qubit2 (int): Another physical qubit
        Returns:
            List: The shortest undirected path
        Raises:
            CouplingError: When there is no path between physical_qubit1, physical_qubit2.
        """
        paths = rx.digraph_dijkstra_shortest_paths(
            self.graph, source=physical_qubit1, target=physical_qubit2, as_undirected=True
        )
        if not paths:
            raise CouplingError(
                f"Nodes {str(physical_qubit1)} and {str(physical_qubit2)} are not connected"
            )
        return paths[physical_qubit2]

    @property
    def is_symmetric(self):
        """
        Test if the graph is symmetric.

        Return True if symmetric, False otherwise
        """
        if self._is_symmetric is None:
            self._is_symmetric = self._check_symmetry()
        return self._is_symmetric

    def make_symmetric(self):
        """
        Convert uni-directional edges into bi-directional.
        """

        edges = self.get_edges()
        edge_set = set(edges)
        for src, dest in edges:
            if (dest, src) not in edge_set:
                self.graph.add_edge(dest, src, None)
        self._dist_matrix = None  # invalidate
        self._is_symmetric = None  # invalidate

    def _check_symmetry(self):
        """
        Calculates symmetry

        Returns:
            Bool: True if symmetric, False otherwise
        """
        return self.graph.is_symmetric()

    def reduce(self, mapping, check_if_connected=True):
        """Returns a reduced coupling map that
        corresponds to the subgraph of qubits
        selected in the mapping.

        Args:
            mapping (list): A mapping of reduced qubits to device
                qubits.
            check_if_connected (bool): if True, checks that the reduced
                coupling map is connected.

        Returns:
            CouplingMap: A reduced coupling_map for the selected qubits.

        Raises:
            CouplingError: Reduced coupling map must be connected.
        """

        inv_map = [None] * (max(mapping) + 1)
        for idx, val in enumerate(mapping):
            inv_map[val] = idx

        reduced_cmap = []

        for edge in self.get_edges():
            if edge[0] in mapping and edge[1] in mapping:
                reduced_cmap.append([inv_map[edge[0]], inv_map[edge[1]]])

        # Note: using reduced_coupling_map.graph is significantly faster
        # than calling add_physical_qubit / add_edge.
        reduced_coupling_map = CouplingMap()
        for node in range(len(mapping)):
            reduced_coupling_map.graph.add_node(node)
        reduced_coupling_map.graph.extend_from_edge_list([tuple(x) for x in reduced_cmap])

        if check_if_connected and not reduced_coupling_map.is_connected():
            raise CouplingError("coupling_map must be connected.")

        return reduced_coupling_map

    @classmethod
    def from_full(cls, num_qubits, bidirectional=True) -> "CouplingMap":
        """Return a fully connected coupling map on n qubits."""
        cmap = cls(description="full")
        if bidirectional:
            cmap.graph = rx.generators.directed_mesh_graph(num_qubits)
        else:
            edge_list = []
            for i in range(num_qubits):
                for j in range(i):
                    edge_list.append((j, i))
            cmap.graph.extend_from_edge_list(edge_list)
        return cmap

    @classmethod
    def from_line(cls, num_qubits, bidirectional=True) -> "CouplingMap":
        """Return a coupling map of n qubits connected in a line."""
        cmap = cls(description="line")
        cmap.graph = rx.generators.directed_path_graph(num_qubits, bidirectional=bidirectional)
        return cmap

    @classmethod
    def from_ring(cls, num_qubits, bidirectional=True) -> "CouplingMap":
        """Return a coupling map of n qubits connected to each of their neighbors in a ring."""
        cmap = cls(description="ring")
        cmap.graph = rx.generators.directed_cycle_graph(num_qubits, bidirectional=bidirectional)
        return cmap

    @classmethod
    def from_grid(cls, num_rows, num_columns, bidirectional=True) -> "CouplingMap":
        """Return a coupling map of qubits connected on a grid of num_rows x num_columns."""
        cmap = cls(description="grid")
        cmap.graph = rx.generators.directed_grid_graph(
            num_rows, num_columns, bidirectional=bidirectional
        )
        return cmap

    @classmethod
    def from_heavy_hex(cls, distance, bidirectional=True) -> "CouplingMap":
        """Return a heavy hexagon graph coupling map.

        A heavy hexagon graph is described in:

        https://journals.aps.org/prx/abstract/10.1103/PhysRevX.10.011022

        Args:
            distance (int): The code distance for the generated heavy hex
                graph. The value for distance can be any odd positive integer.
                The distance relates to the number of qubits by:
                :math:`n = \\frac{5d^2 - 2d - 1}{2}` where :math:`n` is the
                number of qubits and :math:`d` is the ``distance`` parameter.
            bidirectional (bool): Whether the edges in the output coupling
                graph are bidirectional or not. By default this is set to
                ``True``
        Returns:
            CouplingMap: A heavy hex coupling graph
        """
        cmap = cls(description="heavy-hex")
        cmap.graph = rx.generators.directed_heavy_hex_graph(distance, bidirectional=bidirectional)
        return cmap

    @classmethod
    def from_heavy_square(cls, distance, bidirectional=True) -> "CouplingMap":
        """Return a heavy square graph coupling map.

        A heavy square graph is described in:

        https://journals.aps.org/prx/abstract/10.1103/PhysRevX.10.011022

        Args:
            distance (int): The code distance for the generated heavy square
                graph. The value for distance can be any odd positive integer.
                The distance relates to the number of qubits by:
                :math:`n = 3d^2 - 2d` where :math:`n` is the
                number of qubits and :math:`d` is the ``distance`` parameter.
            bidirectional (bool): Whether the edges in the output coupling
                graph are bidirectional or not. By default this is set to
                ``True``
        Returns:
            CouplingMap: A heavy square coupling graph
        """
        cmap = cls(description="heavy-square")
        cmap.graph = rx.generators.directed_heavy_square_graph(
            distance, bidirectional=bidirectional
        )
        return cmap

    @classmethod
    def from_hexagonal_lattice(cls, rows, cols, bidirectional=True) -> "CouplingMap":
        """Return a hexagonal lattice graph coupling map.

        Args:
            rows (int): The number of rows to generate the graph with.
            cols (int): The number of columns to generate the graph with.
            bidirectional (bool): Whether the edges in the output coupling
                graph are bidirectional or not. By default this is set to
                ``True``
        Returns:
            CouplingMap: A hexagonal lattice coupling graph
        """
        cmap = cls(description="hexagonal-lattice")
        cmap.graph = rx.generators.directed_hexagonal_lattice_graph(
            rows, cols, bidirectional=bidirectional
        )
        return cmap

    def largest_connected_component(self):
        """Return a set of qubits in the largest connected component."""
        return max(rx.weakly_connected_components(self.graph), key=len)

    def connected_components(self) -> List["CouplingMap"]:
        """Separate a :Class:`~.CouplingMap` into subgraph :class:`~.CouplingMap`
        for each connected component.

        The connected components of a :class:`~.CouplingMap` are the subgraphs
        that are not part of any larger subgraph. For example, if you had a
        coupling map that looked like::

            0 --> 1   4 --> 5 ---> 6 --> 7
            |     |
            |     |
            V     V
            2 --> 3

        then the connected components of that graph are the subgraphs::

            0 --> 1
            |     |
            |     |
            V     V
            2 --> 3

        and::

            4 --> 5 ---> 6 --> 7

        For a connected :class:`~.CouplingMap` object there is only a single connected
        component, the entire :class:`~.CouplingMap`.

        This method will return a list of :class:`~.CouplingMap` objects, one for each connected
        component in this :class:`~.CouplingMap`. The data payload of each node in the
        :attr:`~.CouplingMap.graph` attribute will contain the qubit number in the original
        graph. This will enables mapping the qubit index in a component subgraph to
        the original qubit in the combined :class:`~.CouplingMap`. For example::

            from qiskit.transpiler import CouplingMap

            cmap = CouplingMap([[0, 1], [1, 2], [2, 0], [3, 4], [4, 5], [5, 3]])
            component_cmaps = cmap.connected_components()
            print(component_cmaps[1].graph[0])

        will print ``3`` as index ``0`` in the second component is qubit 3 in the original cmap.

        Returns:
            list: A list of :class:`~.CouplingMap` objects for each connected
                components. The order of this list is deterministic but
                implementation specific and shouldn't be relied upon as
                part of the API.
        """
        # Set payload to index
        for node in self.graph.node_indices():
            self.graph[node] = node
        components = rx.weakly_connected_components(self.graph)
        output_list = []
        for component in components:
            new_cmap = CouplingMap()
            new_cmap.graph = self.graph.subgraph(sorted(component))
            output_list.append(new_cmap)
        return output_list

    def __str__(self):
        """Return a string representation of the coupling graph."""
        string = ""
        if self.get_edges():
            string += "["
            string += ", ".join([f"[{src}, {dst}]" for (src, dst) in self.get_edges()])
            string += "]"
        return string

    def __eq__(self, other):
        """Check if the graph in ``other`` has the same node labels and edges as the graph in
        ``self``.

        This function assumes that the graphs in :class:`.CouplingMap` instances are connected.

        Args:
            other (CouplingMap): The other coupling map.

        Returns:
            bool: Whether or not other is isomorphic to self.
        """
        if not isinstance(other, CouplingMap):
            return False
        return set(self.graph.edge_list()) == set(other.graph.edge_list())

    # def draw(self):
    #     """Draws the coupling map.

    #     This function calls the :func:`~rustworkx.visualization.graphviz_draw` function from the
    #     ``rustworkx`` package to draw the :class:`CouplingMap` object.

    #     Returns:
    #         PIL.Image: Drawn coupling map.

    #     """

    #     # Distinguish node attributes based on their type
    #     node_attr = {}
    #     for i, node in enumerate(self.graph.nodes()):
    #         if node["type"] == "qubit":
    #             node_attr[i] = {"color": "blue", "shape": "ellipse", "label": f"Qubit {node['index']}"}
    #         elif node["type"] == "qumode":
    #             node_attr[i] = {"color": "green", "shape": "rectangle", "label": f"Qumode {node['index']}"}


    #     return graphviz_draw(self.graph, node_attr=node_attr, method="neato")


    def draw(self, label_color="red"):
        """Draws the coupling map with node and edge attributes, with customizable label colors."""

        dot = Graph(engine='neato')

        # Distinguish node attributes based on their type
        for i, node in enumerate(self.graph.nodes()):
            if node["type"] == "qubit":
                dot.node(str(i), label=f"Qubit {node['index']}", color="blue", shape="ellipse")
            elif node["type"] == "qumode":
                dot.node(str(i), label=f"Qumode {node['index']}", color="green", shape="rectangle")

        # Add edges with attributes like weight and cost, and set label color
        for src, dst in self.graph.edge_list():  # Only source and destination are returned
            attrs = self.graph.get_edge_data(src, dst)  # Retrieve edge attributes
            if attrs:
                # label = f"type={attrs.get('type', 'N/A')}, weight={attrs.get('weight', 1)}, cost={attrs.get('cost', 0)}"
                label = f"{attrs.get('cost', 0)}"
            else:
                label = ""  # If there are no attributes, leave label empty

            # Add the edge with a custom label color
            dot.edge(str(src), str(dst), label=label, fontcolor=label_color)

        # Render the graph to a file and then open it using PIL
        dot.render(filename='cv_map', format='png', cleanup=True)
        image = Image.open('cv_map.png')
        return image




"""
# Test program

cv_map = CVCouplingMap(description="simple test map")

# # Add qubits (using IDs distinct from qumodes)
# cv_map.add_physical_qubit(0)
# cv_map.add_physical_qubit(1)
# cv_map.add_physical_qubit(2)
# cv_map.add_physical_qubit(3)

# # Add qumodes (using IDs distinct from qubits)
# cv_map.add_physical_qumode(100)
# cv_map.add_physical_qumode(101)
# cv_map.add_physical_qumode(102)
# cv_map.add_physical_qumode(103)

# # Add edges between qubits
# cv_map.add_edge(0, 1)
# cv_map.add_edge(1, 2)
# cv_map.add_edge(2, 3)

# # Add edges between qumodes
# cv_map.add_edge(100, 101)
# cv_map.add_edge(101, 102)
# cv_map.add_edge(102, 103)

# # Add edges between qubit and qumode
# cv_map.add_edge(1, 101)

# Construct Superconducting coupling map
## Add qubits and qumode
m = 3 # N = m * m
N = m * m # number of qubits = number of qumodes
for i in range(0, N):
    cv_map.add_physical_qumode(i)
for i in range(N, N * 2):
    cv_map.add_physical_qubit(i)

## Connnect qumodes
for i in range(0, m): # connect rows
    pos = i * m
    cv_map.add_edge(pos, pos + 1)
    cv_map.add_edge(pos + 1, pos + 2)
for j in range(0, m): # connect columns
    pos = j
    cv_map.add_edge(pos, pos + m)
    cv_map.add_edge(pos + m, pos + m * 2)

## Connect qumodes and qubits
for i in range(0, N):
    cv_map.add_edge(i, i + N)


# print(f"Distance bwtween qubit 0 and qubit 3: {cv_map.distance(0, 3)}")

# Visualize the coupling map
image = cv_map.draw()
# image.show()  # Display the image
image.save('cv_map.png')
"""