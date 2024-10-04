import networkx as nx
import dag
import sys
import os
sys.path.append(os.path.abspath('/common/home/mg1998/CVDV/bosonic-qiskit'))
import c2qa
from cvdv.cv_coupling import CVCouplingMap
import parser_cvdvqasm
import pdb
import heapq
import math





r""" Cost Function
'''
H_{basic} = \sum_{gate \in F} D[\pi(gate.q_1)][\pi(gate.q2)]

H_{decay}=\frac{1}{\left|{F}\right|}\sum_{gate \in F} D[\pi(gate.q_1)][\pi(gate.q2)]
                    + W*\frac{1}{\left|{E}\right|} \sum_{gate \in E} D[\pi(gate.q_1)][\pi(gate.q2)]

H_{decay} = max(decay(SWAP.q_1), decay(SWAP.q_2)) {
                    \frac{1}{\left|{F}\right|} \sum_{gate \in F} D[\pi(gate.q_1)][\pi(gate.q2)]\\
                    + W *\frac{1}{\left|{E}\right|} \sum_{gate \in E} D[\pi(gate.q_1)][\pi(gate.q2)]
                    }
'''


# Input
'''
## DAG
## Coupling Map
## Latency of Each Gate
'''
coupling_graph = CVCouplingMap()

# Output
'''
## Initial Maping
## Physical Circuit
    - Original gates
    - SWAP gates
## Evaluation Circuit
    - Gate count
    - Circuit duration
'''

"""



class CVDVQuantumCircuitMapper:
    def __init__(self, qc_type="super_conducting", num_qubits=3, num_qumodes=3):
        self.qc_type = qc_type
        self.num_qubits = num_qubits
        self.num_qumodes = num_qumodes
        self.coupling_graph = self.gen_coupling_graph()
        self.distance_matrix = self.gen_distance_matrix()

    def gen_coupling_graph(self):
        cv_map = CVCouplingMap(description=f"{self.num_qubits}*{self.num_qumodes} {self.qc_type} coupling map")
        if self.qc_type == "super_conducting":
            if self.num_qubits != self.num_qumodes:
                raise ValueError("Number of qubits should be equal to the number of qumodes in superconducting architecture")
            self.setup_superconducting(cv_map)
        return cv_map    


def gen_coupling_graph(qc_type = "super_conducting", num_qubits = 3, num_qumodes = 3):

    # Construct Superconducting coupling map
    cv_map = CVCouplingMap(description = f"3*3 superconducting coupling map")

    if qc_type == "super_conducting":
        if num_qubits != num_qumodes:
            raise ValueError("In super conducting architecture, the number of qubits should be equale to the number of qumodes")

        ## Add qubits and qumode
        m = num_qubits # N = m * m
        N = m * m # number of qubits = number of qumodes
        for i in range(0, N):
            cv_map.add_physical_qumode(i)
        for i in range(N, N * 2):
            cv_map.add_physical_qubit(i)

        ## Connnect qumodes
        for i in range(0, m): # connect rows
            pos = i * m
            cv_map.add_edge(pos, pos + 1, edge_type="qumode-qumode", weight=1, cost=1)
            cv_map.add_edge(pos + N, pos + 1 + N, edge_type="qubit-qubit", weight=1, cost=9)        # Virtual Link
            cv_map.add_edge(pos + 1, pos + 2, edge_type="qumode-qumode", weight=1, cost=1)
            cv_map.add_edge(pos + 1 + N, pos + 2 + N, edge_type="qubit-qubit", weight=1, cost=9)    # Virtual Link
        for j in range(0, m): # connect columns
            pos = j
            cv_map.add_edge(pos, pos + m, edge_type="qumode-qumode", weight=1, cost=1)
            cv_map.add_edge(pos + N, pos + m + N, edge_type="qubit-qubit", weight=1, cost=9)
            cv_map.add_edge(pos + m, pos + m * 2, edge_type="qumode-qumode", weight=1, cost=1)      
            cv_map.add_edge(pos + m + N, pos + m * 2 + N, edge_type="qubit-qubit", weight=1, cost=9)

        ## Connect qumodes and qubits
        for i in range(0, N):
            cv_map.add_edge(i, i + N, edge_type="qumode-qubit", weight=1, cost=2)


    # print(f"Distance bwtween qubit 0 and qubit 3: {cv_map.distance(0, 3)}")

    # Visualize the coupling map
    image = cv_map.draw()
    print(cv_map.distance(8, 9))
    # image.show()  # Display the image
    image.save('cv_map_superconducting.png')

    return cv_map



def heuristic_search(front_layer, ini_mapping, dist_matrix, cir_dag, coupling_graph):
    return 0


def cal_distance(coupling_graph: CVCouplingMap, node1: int, node2: int):
    return coupling_graph.distance(node1, node2)
# Distance matrix   get_dis()


# Get the distance between two nodes from the matrix distance
# def get_distance(node1: int, node2: int):
#     return 0


def gen_distance_matrix(coupling_graph: CVCouplingMap):
    # Assuming CVCouplingMap has an attribute nodes_count that gives the number of nodes
    num_nodes = coupling_graph.size()

    # Initialize the distance matrix with default values (0s, infs, or appropriate values)
    distance_matrix = [[float('inf')] * num_nodes for _ in range(num_nodes)]

    # Populate the distance matrix with the distances between each pair of nodes
    for i in range(num_nodes):
        for j in range(num_nodes):
            if i == j:
                distance_matrix[i][j] = 0  # Distance from a node to itself is 0
            else:
                distance_matrix[i][j] = cal_distance(coupling_graph, i, j)

    return distance_matrix


def gen_ini_map(qumodes_count, qubits_count, qc_type = "super_conducting", map_alg = "index"):
    if map_alg == "index":
        log_to_phy_mapping = {} # logical -> physical
        phy_to_log_mapping = {} # physical -> logical
        nodes_count = qumodes_count + qubits_count
        for i in range(nodes_count):
            log_to_phy_mapping[f"qm{i}"] = i
        for i in range(nodes_count):
            phy_to_log_mapping[i] = f"qm{i}"

    return log_to_phy_mapping, phy_to_log_mapping



def swap_alg_baseline(dag_input: nx.DiGraph, coupling_graph: CVCouplingMap, mapping_log_to_phy, mapping_phy_to_log):
    # Distance_matrix
    distance_matrix = gen_distance_matrix(coupling_graph)

    # Assuming G is a typo and should be dag
    node_depths = {node: dag_input.nodes[node]['depth'] for node in dag_input.nodes}
    # Sort nodes by depth
    sorted_nodes = sorted(node_depths.items(), key=lambda x: x[1])  # nodes in DAG

    # Determine the maximum depth
    longest_depth = max(node_depths.values())

    # Initialize data structure to hold layers
    layers = [[] for _ in range(longest_depth + 1)]

    # Distribute nodes into layers based on their depth
    for node, depth in sorted_nodes:
        if depth > 0:   # exclude qumodes node
            layers[depth].append(node)

    # Process each layer to evaluate and perform swaps if necessary
    for i, layer in enumerate(layers):  # nodes in DAG
        front_layer = layer
        swap_candidate_list = []

        gate_qumodes_list = []   # List of qumodes pair which need to be swapped to be connected in current front_layer
        # Example process for identifying swap candidates
        for gate in front_layer:    # nodes in front_layer of DAG
            # node means Gate in DAG, so extract qumodes/qubits applied in this gate
            gate_qumodes = []   # A simple list to record the gate name and corresponding qumodes
            gate_qumodes.append(gate)
            qumodes = dag_input.nodes[gate]["qumodes"]
            gate_qumodes.append(qumodes)
            # print(qumodes)
            gate_qumodes_list.append(gate_qumodes)

        print(gate_qumodes_list)   
        for gate_qumodes in gate_qumodes_list:
            gate = gate_qumodes[0]
            qumodes = gate_qumodes[1]
            while distance_matrix[mapping_log_to_phy[qumodes[0]]][mapping_log_to_phy[qumodes[1]]] > 1:
                # print(f"DISTANCE: {distance_matrix[mapping_log_to_phy[qumodes[0]]][mapping_log_to_phy[qumodes[1]]]}")
                # Apply SWAP on these qumodes
                qumode_to_swap_0 = qumodes[0]
                qumode_to_swap_1 = qumodes[1]
                neighbors = coupling_graph.neighbors(mapping_log_to_phy[qumode_to_swap_0])    # neighbors of qumodes[0]
                for neighbor in neighbors:
                    print(neighbor)
                    if distance_matrix[neighbor][mapping_log_to_phy[qumode_to_swap_1]] < distance_matrix[mapping_log_to_phy[qumode_to_swap_0]][mapping_log_to_phy[qumode_to_swap_1]]:
                        # Do SWAP
                        print(f"SWAP Physical {neighbor} and {mapping_log_to_phy[qumode_to_swap_0]}")
                        ## Update Mapping
                        temp_phy_qumode = mapping_log_to_phy[qumode_to_swap_0]
                        temp_log_neighbor = mapping_phy_to_log[neighbor]
                        mapping_log_to_phy[qumode_to_swap_0] = neighbor
                        mapping_log_to_phy[temp_log_neighbor] = temp_phy_qumode
                        mapping_phy_to_log[neighbor] = qumode_to_swap_0
                        mapping_phy_to_log[temp_phy_qumode] = temp_log_neighbor
                        print(f"mapping_log_to_phy[{qumode_to_swap_0}] = {neighbor} mapping_log_to_phy[{temp_log_neighbor}] = {temp_phy_qumode}  mapping_phy_to_log[{neighbor}] = {qumode_to_swap_0} mapping_phy_to_log[{temp_phy_qumode}] = {temp_log_neighbor}")
                        
                        ## Insert SWAP gates into DAG
                        SWAP = [temp_log_neighbor, qumode_to_swap_0]
                        dag_input.nodes[gate]['SWAPs_list'].append(SWAP)
                        print("uuuuuuuu")
                        print(dag_input.nodes[gate])

                        break


    return dag_input


def if_gate_can_execute_on_device(gate, dag_input, distance_matrix, mapping_log_to_phy, mapping_phy_to_log):    # It means the gate can execute directly without SWAP
    qumodes = dag_input.nodes[gate]["qumodes"]
    print(f"qumodes: {qumodes[0]} and {qumodes[1]}")
    distance = distance_matrix[mapping_log_to_phy[qumodes[0]]][mapping_log_to_phy[qumodes[1]]]
    # pdb.set_trace()
    if distance > 1:
        print("Cannot Execute")
        return False
    else:
        print("Executable")
        return True

def get_successor_gate(dag_input, gate):
    successors = list(dag_input.successors(gate))
    return successors

def if_gate_dependence_solve(dag_input, gate):
    return dag_input.in_degree(gate) == 0



# Obtain swap candidate list (in physical qumode format) and total distance need to be swapped currently
def obtain_SWAPs_and_totalDist(front_layer, dag_input, coupling_graph, distance_matrix, mapping_log_to_phy, mapping_phy_to_log):
    swap_candidate_list = []
    total_dist = 0

    for gate in front_layer:
        qumodes = dag_input.nodes[gate]["qumodes"]  # Find all(2) qumodes of this gate
        total_dist += distance_matrix[mapping_log_to_phy[qumodes[0]]][mapping_log_to_phy[qumodes[1]]]
        for qumode in qumodes:
            phy_qumode = mapping_log_to_phy[qumode]
            neighbors = coupling_graph.neighbors(phy_qumode)    # Find all neighbors of this qumode
            for neighbor in neighbors: 
                swap = []
                swap.append(phy_qumode)
                swap.append(neighbor)
                swap.sort()
                if swap not in swap_candidate_list:
                    swap_candidate_list.append(swap)

    print("swap_candidate_list:")
    print(swap_candidate_list)
    return swap_candidate_list, total_dist

def heuristic_cost(front_layer, dag_input, swap, coupling_graph, distance_matrix, mapping_log_to_phy, mapping_phy_to_log, total_dist, decay, heuristic="basic", weight_extended_set=0.5):
    cost = 0
    cost_basic = 0
    cost_lookahead = 0
    phy_q_swap_0 = swap[0]
    phy_q_swap_1 = swap[1]

    # Calculate basic cost
    for gate in front_layer:
        qumodes = dag_input.nodes[gate]["qumodes"]  # Find all(2) qumodes of this gate
        phy_qumode_0 = mapping_log_to_phy[qumodes[0]]
        phy_qumode_1 = mapping_log_to_phy[qumodes[1]]
        original_distance = distance_matrix[phy_qumode_0][phy_qumode_1]
        
        if phy_qumode_0 == phy_q_swap_0:
            change = distance_matrix[phy_q_swap_1][phy_qumode_1] - original_distance
            cost_basic += change
        elif phy_qumode_0 == phy_q_swap_1:
            change = distance_matrix[phy_q_swap_0][phy_qumode_1] - original_distance
            cost_basic += change
        elif phy_qumode_1 == phy_q_swap_0:
            change = distance_matrix[phy_qumode_0][phy_q_swap_1] - original_distance
            cost_basic += change
        elif phy_qumode_1 == phy_q_swap_1: 
            change = distance_matrix[phy_qumode_0][phy_q_swap_0] - original_distance
            cost_basic += change

    # Calculate lookahead cost
    ## construct extended set
    extended_set = []
    for gate in front_layer:
        for successor_gate in get_successor_gate(dag_input, gate):
            extended_set.append(successor_gate)
    # pdb.set_trace()


    if heuristic == "basic":
        cost = cost_basic
    elif heuristic == "lookahead": 
        for gate in extended_set:
            qumodes = dag_input.nodes[gate]["qumodes"]  # Find all(2) qumodes of this gate
            phy_qumode_0 = mapping_log_to_phy[qumodes[0]]
            phy_qumode_1 = mapping_log_to_phy[qumodes[1]]
            original_distance = distance_matrix[phy_qumode_0][phy_qumode_1]
            
            if phy_qumode_0 == phy_q_swap_0:
                change = distance_matrix[phy_q_swap_1][phy_qumode_1] - original_distance
                cost_lookahead += change * weight_extended_set
            elif phy_qumode_0 == phy_q_swap_1:
                change = distance_matrix[phy_q_swap_0][phy_qumode_1] - original_distance
                cost_lookahead += change * weight_extended_set
            elif phy_qumode_1 == phy_q_swap_0:
                change = distance_matrix[phy_qumode_0][phy_q_swap_1] - original_distance
                cost_lookahead += change * weight_extended_set
            elif phy_qumode_1 == phy_q_swap_1: 
                change = distance_matrix[phy_qumode_0][phy_q_swap_0] - original_distance
                cost_lookahead += change * weight_extended_set

        cost = cost_basic / len(front_layer) + cost_lookahead / max(len(extended_set), 1)
    elif heuristic == "decay":
        cost = max(decay[phy_q_swap_0], decay[phy_q_swap_1]) * cost_basic / len(front_layer) + cost_lookahead / max(len(extended_set), 1)

    # pdb.set_trace()

    print(f"cost: {cost}")
    return cost

def gen_cvdvqasm(final_gate_list, filename: str):
    with open(filename, 'w') as fp:
        for gate in final_gate_list:
            if gate[0] == 'BS':
                num_qm1 = gate[1]
                num_qm2 = gate[2]
                fp.write(f"bsCV qmr[{num_qm1}], qmr[{num_qm2}]; // bsCV {gate[3]}, {gate[4]}\n")
            elif gate[0] == 'SWAP':
                fp.write(f"SWAP {gate[1]}, {gate[2]}, qmr[{gate[3]}], qmr[{gate[4]}]; // SWAP {gate[5]}, {gate[6]}\n")
            else:
                raise ValueError(f"Unrecognized instruction: {gate[0]}")


def cal_cir_duration(num_phy_nodes, final_gate_list, gate_latency):
    cir_duration = 0
    # Initialize time stamp for each qumode/qubit
    time_nodes = []
    for i in range(num_phy_nodes):
        time_nodes.append(0)

    # Schedule
    for gate in final_gate_list:
        if gate[0] == 'BS':
            num_qm1 = int(gate[1])
            num_qm2 = int(gate[2])
            start_time = max(time_nodes[num_qm1], time_nodes[num_qm2])
            end_time = start_time + gate_latency['BS']
            time_nodes[num_qm1] = end_time
            time_nodes[num_qm2] = end_time
        elif gate[0] == 'SWAP':
            num_qm1 = int(gate[3])
            num_qm2 = int(gate[4])
            start_time = max(time_nodes[num_qm1], time_nodes[num_qm2])
            end_time = start_time + gate_latency['SWAP']
            time_nodes[num_qm1] = end_time
            time_nodes[num_qm2] = end_time
        else:
            raise ValueError(f"Unrecognized instruction: {gate[0]}")

    for time in time_nodes:
        cir_duration = max(cir_duration, time)
    return cir_duration


def swap_sabre(dag_input: nx.DiGraph, coupling_graph: CVCouplingMap, mapping_log_to_phy, mapping_phy_to_log, gate_latency, heuristic="basic", weight_extended_set=0.5, decay_increment=0.01, decay_reduction=0.005):
    front_layer =[]
    final_gate_list = []    # contain original gates and SWAP gates, with the execution order
    num_nodes = coupling_graph.size()

    # Distance_matrix
    distance_matrix = gen_distance_matrix(coupling_graph)

    # Decay value of each qumodes
    decay = []
    ## Initialization
    for i in range(num_nodes):
        decay.append(1)

    # Delete qumode nodes in DAG
    qumode_nodes = []
    for node in dag_input.nodes:
        if node[0] == 'q' and node[1] == 'm':
            qumode_nodes.append(node)
    for node in qumode_nodes:
        dag_input.remove_node(node)
    
    # Put nodes whose dependence solved into the front layer
    for node in dag_input.nodes:
        if dag_input.in_degree(node) == 0:
            print(f"Add node into front layer:{node}")
            front_layer.append(node)

    while len(front_layer) != 0:
        execute_gate_list = []
        # Find ready gates
        for gate in front_layer:
            if if_gate_can_execute_on_device(gate, dag_input, distance_matrix, mapping_log_to_phy, mapping_phy_to_log):
                print(f"Can execute: {gate}")
                execute_gate_list.append(gate)
        
        if len(execute_gate_list) != 0:
            for gate in execute_gate_list:
                front_layer.remove(gate)
                # Construct instruction
                gate_num = str(gate)
                gate_type = dag_input.nodes[gate]["gate_type"]
                gate_qumodes = dag_input.nodes[gate]["qumodes"]
                gate_instr = []
                # gate_instr.append(gate_num)
                gate_instr.append(gate_type)
                gate_instr.append(mapping_log_to_phy[gate_qumodes[0]])  # physical qumodes
                gate_instr.append(mapping_log_to_phy[gate_qumodes[1]])
                gate_instr.append(gate_qumodes[0])                      # logical qumodes
                gate_instr.append(gate_qumodes[1])
                final_gate_list.append(gate_instr)

                successor_gate_list = get_successor_gate(dag_input, gate)
                dag_input.remove_node(gate)
                print(f"DAG remove gate {gate}")
                print(f"successor_gate_list: {successor_gate_list}")
                for successor_gate in successor_gate_list:
                    if if_gate_dependence_solve(dag_input, successor_gate):
                        front_layer.append(successor_gate)
                    else:
                        # predecessors = dag_input.predecessors(successor_gate)
                        print(f"{successor_gate} cannot execute")
        else:
            scores = {}
            # Update the reduction of decay value
            for i in range(len(decay)):
                decay[i] = max(1, decay[i] - decay_reduction)

            swap_candidate_list, total_dist = obtain_SWAPs_and_totalDist(front_layer, dag_input, coupling_graph, distance_matrix, mapping_log_to_phy, mapping_phy_to_log)
            print(f"total distance: {total_dist}")
            for swap in swap_candidate_list:
                print(f"Swap {swap[0]} and {swap[1]}")
                scores[str(swap)] = heuristic_cost(front_layer, dag_input, swap, coupling_graph, distance_matrix, mapping_log_to_phy, mapping_phy_to_log, total_dist, decay, heuristic, weight_extended_set)
            
            # Find the SWAP with minimal score
            # pdb.set_trace()
            priority_queue = [(value, key) for key, value in scores.items()]
            heapq.heapify(priority_queue)
            sorted_items = [heapq.heappop(priority_queue) for _ in range(len(priority_queue))]

            # Select the swap with minimum cost
            _, swap = sorted_items[0]
            swap_list = swap.strip("[]").split(",")
            phy_q_swap_0 = int(swap_list[0])
            phy_q_swap_1 = int(swap_list[1])
            # Update Decay Value
            decay[phy_q_swap_0] += decay_increment
            decay[phy_q_swap_1] += decay_increment
            # Perform SWAP
            print(f"SWAP {phy_q_swap_0} and {phy_q_swap_1}")
            swap_instr = []
            swap_instr.append('SWAP')
            swap_instr.append('2pi')
            swap_instr.append('-0.5pi')
            swap_instr.append(phy_q_swap_0)     # physical qumodes
            swap_instr.append(phy_q_swap_1)
            swap_instr.append(mapping_phy_to_log[phy_q_swap_0])     # logical qumodes
            swap_instr.append(mapping_phy_to_log[phy_q_swap_1])
            final_gate_list.append(swap_instr)
            ## Update Mapping
            temp_log_q_swap_0 = mapping_phy_to_log[phy_q_swap_0]
            temp_log_q_swap_1 = mapping_phy_to_log[phy_q_swap_1]
            mapping_phy_to_log[phy_q_swap_0] = temp_log_q_swap_1
            mapping_phy_to_log[phy_q_swap_1] = temp_log_q_swap_0
            mapping_log_to_phy[temp_log_q_swap_0] = phy_q_swap_1
            mapping_log_to_phy[temp_log_q_swap_1] = phy_q_swap_0


    print(final_gate_list)
    gen_cvdvqasm(final_gate_list, "sabre_example.cvdvqasm")
    cir_duration = cal_cir_duration(num_nodes, final_gate_list, gate_latency)
    return final_gate_list, cir_duration

def perform_swaps(dag, coupling_graph, swap_candidate_list):
    # Placeholder function for swap logic
    # This could involve finding the shortest path to swap, applying the swap, and updating the DAG
    pass



    



def dag_to_qasm():


    return 0


def run():
    return 0, 0



def gen_test_dag():
    gates = [{"gate": "G1", "gate_type": "BS", "qumodes": ["qm1", "qm2"], "SWAPs_list":[]}, 
        {"gate": "G2", "gate_type": "BS", "qumodes": ["qm1", "qm3"], "SWAPs_list":[]}, 
        {"gate": "G3", "gate_type": "BS", "qumodes": ["qm1", "qm4"], "SWAPs_list":[]}, 
        {"gate": "G4", "gate_type": "BS", "qumodes": ["qm5", "qm6"], "SWAPs_list":[]}, 
        {"gate": "G5", "gate_type": "BS", "qumodes": ["qm5", "qm4"], "SWAPs_list":[]}]

    G, qumode_depth, gate_depth = dag.gen_dag(7, gates)
    dag.draw_dag(G, qumode_depth, gate_depth, "DAG_before_swaps")
    return G, qumode_depth, gate_depth


# Initialization
## INPUT
coupling_graph = gen_coupling_graph()                       # Coupling Graph
distance_matrix = gen_distance_matrix(coupling_graph)       # Distance Matrix
print(f"dis:{distance_matrix[6][5]}")
mapping_log_to_phy, mapping_phy_to_log = gen_ini_map(9, 9)    # Initial Mapping
cir_dag_input, qumode_depth, gate_depth = gen_test_dag()
gate_latency = {"BS": 1, "SWAP": 1}
# cir_dag_output = swap_alg_baseline(cir_dag_input, coupling_graph, mapping_log_to_phy, mapping_phy_to_log)
# dag.draw_dag(cir_dag_output, qumode_depth, gate_depth, "DAG_after_swaps")



# cir_dag_output, ini_mapping = run()


_, cir_duration = swap_sabre(cir_dag_input, coupling_graph, mapping_log_to_phy, mapping_phy_to_log, gate_latency, "decay")
print(f"Circuit Duration: {cir_duration}")