import networkx as nx
import matplotlib.pyplot as plt


# To calculate the depth of each gate(node)
def longest_path_to_node(G, node, memo={}):
    if node in memo:
        return memo[node]
    if G.in_degree(node) == 0:
        # Base case: No incoming edges, return 0
        return 0
    # Recursive case: 1 + maximum of longest paths to predecessors
    max_length = 0
    for predecessor in G.predecessors(node):
        max_length = max(max_length, longest_path_to_node(G, predecessor, memo))
    memo[node] = 1 + max_length
    return memo[node]


# To generate graph by gates
'''
def gen_dag(num_qumodes, gates, fig_name = "grid_dag.png"):
    # Create a directed graph
    G = nx.DiGraph()

    # Add nodes for each qmode
    num_qumodes = num_qumodes
    qumode_depth = {}
    gate_depth = {}
    gate_positions = []

    # Add nodes for qumodes
    for i in range(num_qumodes):
        G.add_node(f"qm{i}", cur_G=f"qm{i}")    # cur_G means the current gate such qumode is applied by
        qumode_depth[f"qm{i}"] = (i, 0)  # Position qumodes in a single row

    # Add new gate
    gates = gates

    # Initialize Y position for gates
    y_pos = 1

    # Add gate nodes and their connections
    for gate in gates:
        gate_name = gate["gate"]
        print(gate_name)
        qumodes = gate["qumodes"]
        SWAPs_list = gate["SWAPs_list"]
        
        # Compute X position as average of the connected qumodes
        x_pos = sum(qumode_depth[qm][0] for qm in qumodes) / len(qumodes)
        gate_depth[gate_name] = (x_pos, y_pos)
        gate_positions.append(x_pos)  # Store positions for determining next row's y position

        G.add_node(gate_name, qumodes=qumodes, SWAPs_list=SWAPs_list)
        
        # Add edges based on hypothetical dependencies
        for qm in qumodes:
            G.add_edge(G.nodes[qm]["cur_G"], gate_name)
            G.nodes[qm]["cur_G"] = gate_name

        # Increase y_pos for next gate row
        y_pos += 1


    G = calculate_depth(G)

    # Draw the Graph

    # Consolidate all positions
    pos = {**qumode_depth, **gate_depth}

    # Drawing the graph
    nx.draw(G, pos, with_labels=False, node_color='skyblue', node_size=3000, font_weight='bold')

    # Check if the graph is a DAG
    print("Is the graph a DAG?", nx.is_directed_acyclic_graph(G))

    # Custom labels with qumodes and gates
    labels = {node: f"{node}\n {', '.join(data.get('qumodes', []))}\n depth={data.get('depth')}\n {', '.join(data.get('SWAPs_list', []))}" for node, data in G.nodes(data=True)}
    nx.draw_networkx_labels(G, pos, labels=labels)

    # Display the graph
    plt.show()

    # Save the graph to a file
    plt.savefig(fig_name)



    # print(longest_path_to_node(G, 'G5'))


    return G

'''

def gen_dag(num_qumodes, gates):
    # Create a directed graph
    G = nx.DiGraph()

    # Add nodes for each qumode
    num_qumodes = num_qumodes
    qumode_depth = {}
    gate_depth = {}
    gate_positions = []

    # Add nodes for qumodes
    for i in range(num_qumodes):
        G.add_node(f"qm{i}", cur_G=f"qm{i}")  # cur_G means the current gate such qumode is applied by
        qumode_depth[f"qm{i}"] = (i, 0)  # Position qumodes in a single row

    # Add new gate
    gates = gates

    # Initialize Y position for gates
    y_pos = 1

    # Add gate nodes and their connections
    for gate in gates:
        gate_name = gate["gate"]
        print(gate_name)
        gate_type = gate["gate_type"]
        qumodes = gate["qumodes"]
        SWAPs_list = gate["SWAPs_list"]
        
        # Compute X position as average of the connected qumodes
        x_pos = sum(qumode_depth[qm][0] for qm in qumodes) / len(qumodes)
        gate_depth[gate_name] = (x_pos, y_pos)
        gate_positions.append(x_pos)  # Store positions for determining next row's y position

        G.add_node(gate_name, gate_type=gate_type, qumodes=qumodes, SWAPs_list=SWAPs_list)
        
        # Add edges based on hypothetical dependencies
        for qm in qumodes:
            G.add_edge(G.nodes[qm]["cur_G"], gate_name)
            G.nodes[qm]["cur_G"] = gate_name

        # Increase y_pos for next gate row
        y_pos += 1

    G = calculate_depth(G)  # Assuming calculate_depth is defined elsewhere

    # Return the graph and positions
    return G, qumode_depth, gate_depth


def draw_dag(G, qumode_depth, gate_depth, fig_name="grid_dag.png"):
    # Consolidate all positions
    pos = {**qumode_depth, **gate_depth}

    # Draw the graph
    nx.draw(G, pos, with_labels=False, node_color='skyblue', node_size=3000, font_weight='bold')

    # Check if the graph is a DAG
    print("Is the graph a DAG?", nx.is_directed_acyclic_graph(G))

    # Custom labels with qumodes and gates
    labels = {node: f"{node}\n {', '.join(data.get('qumodes', []))}\n depth={data.get('depth')}\n {', '.join(map(str, data.get('SWAPs_list', [])))}" for node, data in G.nodes(data=True)}
    nx.draw_networkx_labels(G, pos, labels=labels, font_size=8)

    # Display the graph
    plt.show()

    # Save the graph to a file
    plt.savefig(fig_name)



def calculate_depth(G):
    for node in G.nodes:
        depth = longest_path_to_node(G, node)
        G.nodes[node]["depth"] = depth
    
    return G


def insert_node_into_edge(G, src, dst, new_node):
    # Check if the edge exists
    if graph.has_edge(src, dst):
        # Remove the original edge
        graph.remove_edge(src, dst)
        
        # Add the new node
        graph.add_node(new_node)
        
        # Add new edges
        graph.add_edge(src, new_node)
        graph.add_edge(new_node, dst)
    else:
        raise ValueError("The specified edge does not exist in the graph.")




# gates = [{"gate": "G1", "gate_type": "BS", "qumodes": ["qm1", "qm2"]}, 
#         {"gate": "G2", "gate_type": "BS", "qumodes": ["qm1", "qm3"]}, 
#         {"gate": "G3", "gate_type": "BS", "qumodes": ["qm1", "qm4"]}, 
#         {"gate": "G4", "gate_type": "BS", "qumodes": ["qm5", "qm6"]}, 
#         {"gate": "G5", "gate_type": "BS", "qumodes": ["qm5", "qm4"]}]

# G = gen_dag(7, gates)
