import json
import pickle
import sys
import os
script_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(script_dir)
sys.path.append(parent_dir)
from cvdv.cv_coupling import CVCouplingMap


def gen_gate_latency():
    # Dictionary to store
    gate_latency = {"BS": 1, "SWAP": 1}

    # Serialize and write to a file
    with open('gate_latency.json', 'w') as file:
        json.dump(gate_latency, file)


# Generate Coupling Graph
## linked_qubits means whether virtual links between qubits are needed (for superconducting)
def gen_coupling_graph(qc_type = "superconducting", m_qubits = 3, m_qumodes = 3, linked_qubits = True):  # N = m * m

    # Construct Superconducting coupling map
    cv_map = CVCouplingMap(description = f"3*3 superconducting coupling map")

    if qc_type == "superconducting":
        if m_qubits != m_qumodes:
            raise ValueError("In superconducting architecture, the number of qubits should be equale to the number of qumodes")

        ## Add qubits and qumode
        m = m_qubits # N = m * m
        N = m * m # number of qubits = number of qumodes
        for i in range(0, N):
            cv_map.add_physical_qumode(i)
        for i in range(N, N * 2):
            cv_map.add_physical_qubit(i)

        ## Connnect qumodes/qubits
        for i in range(0, m): # connect rows
            pos = i * m
            cv_map.add_edge(pos, pos + 1, edge_type="qumode-qumode", weight=1, cost=1)
            cv_map.add_edge(pos + 1, pos + 2, edge_type="qumode-qumode", weight=1, cost=1)
            if linked_qubits:
                cv_map.add_edge(pos + N, pos + 1 + N, edge_type="qubit-qubit", weight=1, cost=9)        # Virtual Link
                cv_map.add_edge(pos + 1 + N, pos + 2 + N, edge_type="qubit-qubit", weight=1, cost=9)    # Virtual Link
        for j in range(0, m): # connect columns
            pos = j
            cv_map.add_edge(pos, pos + m, edge_type="qumode-qumode", weight=1, cost=1)
            cv_map.add_edge(pos + m, pos + m * 2, edge_type="qumode-qumode", weight=1, cost=1) 
            if linked_qubits:   
                cv_map.add_edge(pos + N, pos + m + N, edge_type="qubit-qubit", weight=1, cost=9)  
                cv_map.add_edge(pos + m + N, pos + m * 2 + N, edge_type="qubit-qubit", weight=1, cost=9)

        ## Connect qumodes and qubits
        for i in range(0, N):
            cv_map.add_edge(i, i + N, edge_type="qumode-qubit", weight=1, cost=2)


    # print(f"Distance bwtween qubit 0 and qubit 3: {cv_map.distance(0, 3)}")

    # Visualize the coupling map
    image = cv_map.draw()
    # print(cv_map.distance(8, 9))
    # image.show()  # Display the image
    image.save(f"cv_map_{qc_type}_{m_qubits}*{m_qumodes}_linked-qubits-{linked_qubits}.png")

    with open(f"cv_map_{qc_type}_{m_qubits}*{m_qumodes}_linked-qubits-{linked_qubits}.pkl", 'wb') as f:
        pickle.dump(cv_map, f)








gen_gate_latency()
gen_coupling_graph()
