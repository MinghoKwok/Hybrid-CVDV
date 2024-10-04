from typing_extensions import Self, Any, Literal #For annotations only
import math
from collections import deque, defaultdict
import numpy as np
import depq
import queue
import functools

#class for graph nodes. Used to represent single gates and phase polynomials
#notes for use in comparisons: 
# * Node.nodeType compares types of nodes, ignoring parameters and additional data
# * Node.is_equivalent_to() compares if 2 nodes have equivalent data, including or excluding ins/outs based on a given parameter
# * == operator behavior is python default

class Node:
    #note that ins and outs need to be initialized before proper use, see connect_to() and connect_from()
    sig = 0 #<--- debug purposes only
    
    def __init__(this, qubits:list[int], ins:(dict[int,Self]|None) = None, outs:(dict[int,Self]|None) = None):
        this.qubits = qubits.copy()
        this.nodeType = "generic"
        this.sig = Node.sig
        Node.sig += 1
        if ins != None:
            this.ins = ins.copy()
        else:
            this.ins:dict[int,Node] = dict()
        
        if outs != None:
            this.outs = outs.copy()
        else:
            this.outs:dict[int,Node] = dict()

    #makes edge directing this node to 1 other node
    def connect_to(this, target_node:Self, qubit:int):
        if this.outs.get(qubit) != None:
            this.outs[qubit].ins.pop(qubit)
        if target_node.ins.get(qubit) != None:
            target_node.ins[qubit].outs.pop(qubit)
        this.outs[qubit] = target_node
        target_node.ins[qubit] = this

    def connect_from(this, target_node:Self, qubit:int):
        target_node.connect_to(this, qubit)

    #find previous node applying to a qubit
    def prev(this, qubit:int):
        return this.ins.get(qubit)
    
    #find next node applying to a qubit
    def next(this, qubit:int):
        return this.outs.get(qubit)
    
    #attach this node after a set of previous nodes, while remaking the graph edges correctly
    def attach_after(this, prev_nodes:dict[int,Self]):
        for qubit, node in prev_nodes.items():
            this.connect_to(node.outs[qubit], qubit)
            this.connect_from(node, qubit)

    #remove this node from its current circuit, while remaking the graph edges correctly
    def remove_from_circuit(this):
        for q in set(this.ins.keys()) & set(this.outs.keys()):
            this.prev(q).connect_to(this.next(q))
        for q, node in this.ins.items():
            del node.outs[q]
        for q, node in this.outs.items():
            del node.ins[q]
        this.ins = dict()
        this.outs = dict()

    def to_instr(this):
        return "<?>"
    
    def __repr__(this) -> str:
        return this.to_instr() + f" #{this.sig}"
    
    #checks if this and the other gate is equivalent, ignoring input/outputs nodes. 
    #should not be used to replace ==/__eq__, as that should be reserved for default object compare behavior
    def is_equivalent_to(this, other_node:Self, ignore_qubits:bool = False):
        if this.__class__ != other_node.__class__:
            return False
        return this.nodeType == other_node.nodeType and (ignore_qubits or this.qubits == other_node.qubits)
    
    def copy_disconnected(this):
        return Node(this.qubits)
    
#subclass that represents gates
class Gate(Node):
    def __init__(this, gateType:str, qubits:list[int], params:list[float] = []):
        super().__init__(qubits)
        this.nodeType = gateType
        this.params = params.copy()

    def to_instr(this):
        out = this.nodeType
        if this.params != None and len(this.params) > 0:
            out += f"({','.join([str(p) for p in this.params])})"
        out += " "+",".join(f"q[{q}]" for q in this.qubits)
        return out + ';'
    
    def copy_disconnected(this):
        return Gate(this.nodeType, this.qubits, this.params)
    
    def is_equivalent_to(this, other_node:Self, ignore_qubits:bool = False):
        if not super().is_equivalent_to(other_node, ignore_qubits=ignore_qubits):
            return False
        return this.params == other_node.params

#local use, mainly
def guassian_elim(matrix:np.ndarray, n:int):
    id_matrix = np.identity(n, dtype=int)
    cnot_list = []
    #print("lcs")
    checked_cols = set()
    unchecked_cols = set(range(n))
    while len(unchecked_cols) > 0:
        col_min = min(np.sum((matrix ^ id_matrix)[:, a]) for a in unchecked_cols)
        valid_cols = np.nonzero(np.sum(matrix ^ id_matrix, 0) == col_min)[0].reshape(-1)
        for c in valid_cols:
            if c in unchecked_cols:
                col = c
                break

        unchecked_cols.remove(col)

        diag_one = False if matrix[col, col] == 0 else True
        for row in unchecked_cols:
            #print(matrix)
            if matrix[row, col] == 0:
                continue
            if diag_one == False:
                matrix[col,:] ^= matrix[row,:]
                cnot_list.append((row, col))
                #print((row,col))
                diag_one = True
            matrix[row,:] ^= matrix[col,:]
            cnot_list.append((col, row))
            #print((col,row))
            #print(matrix)
        for row in checked_cols:
            #print(matrix)
            if matrix[row, col] == 0:
                continue
            if not diag_one:
                raise Exception("Non-singular matrix")
            matrix[row,:] ^= matrix[col,:]
            cnot_list.append((col, row))
            #print(matrix)
        checked_cols.add(col)
    #print("-")
    return cnot_list

#subclass that represents a phase polynomial
class PhasePoly(Node):
    def __init__(this, qubits:list[int], rotations:list[tuple[float,set[int]]], affineOut:dict[int,set[int]]):
        super().__init__(qubits)
        this.rotations = rotations.copy()
        this.affineOut = affineOut.copy()
        this.nodeType = "phasePoly"
        this.possible_circuits:set[Circuit] = set()
        pass
    
    @functools.total_ordering
    class State():
        def __init__(this, r_angles, parity_table, output_table, active_cols):
            this.parity_table:np.ndarray = parity_table.copy()
            this.output_table:np.ndarray = output_table.copy()
            this.prev_state:Self = None
            this.prev_ops = []
            this.ops_list_len:int = 0
            this.active_col:set[int]
            this.r_angles:list[float] = r_angles.copy()
            this._cost_val = None

            '''print(">",r_angles)
            print(parity_table)'''

            #TODO eliminate single row parities
            ready_parities = np.nonzero(np.sum(this.parity_table, 0) == 1)[0]
            offset = 0
            for col in sorted(ready_parities):
                row = np.nonzero(this.parity_table[:, col-offset])[0][0]
                this.prev_ops.append(("rz",[row], [this.r_angles[col-offset]]))
                this.ops_list_len += 1
                this.parity_table = np.delete(this.parity_table, col-offset, 1)
                this.r_angles = np.delete(this.r_angles, col-offset)
                offset += 1
            this.active_col = set(c - len([None for p in ready_parities if p < c]) for c in active_cols if c not in ready_parities)

            '''print(this.r_angles, this.prev_ops)
            print(this.parity_table)'''

        def cost(this):
            if this._cost_val == None:
                ops_len = this.ops_list_len
                parity_sum = np.sum(this.parity_table)
                guassian_len = len(guassian_elim(this.output_table.copy(), this.output_table.shape[0]))

                this._cost_val = (ops_len + parity_sum + guassian_len, parity_sum, guassian_len, -ops_len)
            return this._cost_val
        
        def full_ops_list(this):
            out_seq = []
            ptr = this
            while ptr is not None:
                out_seq = ptr.prev_ops + out_seq
                ptr = ptr.prev_state
            return out_seq

        def __eq__(this, other:Self):
            return this.cost() == other.cost()
        def __lt__(this, other:Self):
            return this.cost() < other.cost()
    
    def synthesize_row_search(this, buffer_size:int = -1, ends_checked:int = 10):
        
        parity_table = np.array([])
        output_table = np.array([])
        for _,s in this.rotations:
            col = np.array([(1 if q in s else 0) for q in this.qubits]).reshape(-1,1)
            if len(parity_table) > 0:
                parity_table = np.append(parity_table, col, 1)
            else:
                parity_table = col
        for w in this.qubits:
            col = np.array([(1 if q in this.affineOut[w] else 0) for q in this.qubits]).reshape(-1,1)
            if len(output_table) > 0:
                output_table = np.append(output_table, col, 1)
            else:
                output_table = col
        
        initial_state = PhasePoly.State([r[0] for r in this.rotations], parity_table, output_table, set())

        buffer:depq.DEPQ[PhasePoly.State] = depq.DEPQ()
        end_states:list[PhasePoly.State] = []
        buffer.insert(initial_state, initial_state.cost())
        #TODO

        while not buffer.is_empty():
            #print("bufsize:", buffer.size())

            state:PhasePoly.State = buffer.poplast()[0]
            #print("current state cost:", state.cost())

            active_cols = state.active_col
            if len(state.active_col) == 0 and len(state.parity_table) > 0:
                active_cols = set(range(state.parity_table.shape[1]))

            if len(state.parity_table) <= 0 or state.parity_table.shape[1] <= 0:
                end_states.append(state)
                if len(end_states) >= ends_checked:
                    break
                continue

            active_rows = np.nonzero(np.sum(state.parity_table[:,list(sorted(active_cols))], 1) > 0)[0]
            for i in active_rows:
                for j in active_rows:
                    if i == j:
                        continue
                    helped_cols = set(c for c in active_cols if state.parity_table[i,c] == 1 and state.parity_table[j,c] == 1)
                    if len(helped_cols) == 0:
                        continue
                    new_pt = state.parity_table.copy()
                    new_pt[i] ^= new_pt[j] 
                    new_ot = state.output_table.copy()
                    new_ot[i] ^= new_ot[j] 
                    new_state = PhasePoly.State(state.r_angles, new_pt, new_ot, helped_cols)
                    new_state.prev_state = state
                    new_state.prev_ops.insert(0, ("cx", [i,j], []))
                    new_state.ops_list_len += state.ops_list_len + 1
                    if buffer.size() < buffer_size or buffer_size == -1:
                        buffer.insert(new_state, new_state.cost())
                    elif buffer.high() > new_state.cost(): 
                        buffer.popfirst()
                        buffer.insert(new_state, new_state.cost())

        '''print(this)
        for s in end_states:
            print(s.cost())
            out_seq = []
            ptr = s
            while ptr is not None:
                out_seq = ptr.prev_ops + out_seq
                ptr = ptr.prev_state
            print([t[1] for t in out_seq] + guassian_elim(s.output_table.copy(), len(this.qubits)))'''

        end_states.sort()
        out = Circuit(qubits=this.qubits)
        final_state = end_states[0]
        for g in final_state.full_ops_list():
            out.append_node(Gate(g[0], [this.qubits[i] for i in g[1]], g[2]))
        for i,j in guassian_elim(final_state.output_table.copy(), len(this.qubits)):
            out.append_node(Gate("cx", [this.qubits[j], this.qubits[i]], []))

        '''print(out)
        print(this)
        print(out.convert_to_phasePoly())'''

        return out
            

    @functools.total_ordering
    class QueueObject():
        def __init__(this):
            this.parity_table:np.ndarray
            this.output_table:np.ndarray
            this.ops_list:list[tuple[int,int]]
            this._cost_val = None

        def cost(this):
            if this._cost_val == None:
                this._cost_val = len(this.ops_list) + len(guassian_elim(this.output_table.copy(), this.output_table.shape[0]))
                #this._cost_val = len(this.ops_list) + np.sum(this.parity_table) + len(guassian_elim(this.output_table.copy(), this.output_table.shape[0]))
            return this._cost_val

        def __eq__(this, other:Self):
            return this.cost() == other.cost()
        def __lt__(this, other:Self):
            return this.cost() < other.cost()

    #outputs a circuit of gates equivalent to the phase polynomial
    def synthesize_by_col(this):
        out = Circuit(qubits=this.qubits)
        r = np.array([[3,4],[0,1]])
        np.append(r, np.array([7,9]).reshape(-1,1), 1)
        num_qubits = len(this.qubits)

        #qubit_order_map = dict([(v,k) for k,v in enumerate(this.qubits)])
        
        r_angles = [r[0] for r in this.rotations]
        parity_table = np.array([])
        output_table = np.array([])
        for _,s in this.rotations:
            col = np.array([(1 if q in s else 0) for q in this.qubits]).reshape(-1,1)
            if len(parity_table) > 0:
                parity_table = np.append(parity_table, col, 1)
            else:
                parity_table = col
        for w in this.qubits:
            col = np.array([(1 if q in this.affineOut[w] else 0) for q in this.qubits]).reshape(-1,1)
            if len(output_table) > 0:
                output_table = np.append(output_table, col, 1)
            else:
                output_table = col

        def insert_cnot(i,j):
            out.append_node(Gate("cx", [this.qubits[i],this.qubits[j]]))
            parity_table[i] ^= parity_table[j]
            output_table[i] ^= output_table[j]

        weights = np.identity(num_qubits)

        def implement_rotation(col):
            nonlocal parity_table, r_angles, out, weights
            if sum(parity_table[:,col]) != 1:
                buffer = queue.PriorityQueue()
                item = PhasePoly.QueueObject()
                item.parity_table = parity_table.copy()
                item.output_table = output_table.copy()
                item.ops_list = []
                buffer.put(item)

                item:PhasePoly.QueueObject
                while True:
                    #print(buffer.qsize())
                    item = buffer.get()
                    lines = np.where(item.parity_table[:, col])[0]
                    #print("new", item.ops_list)
                    #print(item.parity_table)
                    #print("lines", lines)
                    if len(lines) == 1:
                        break
                    for i in lines:
                        for j in lines:
                            if i == j:
                                continue
                            new_item = PhasePoly.QueueObject()
                            new_item.parity_table = item.parity_table.copy()
                            #print(">testy", i, j)
                            new_item.parity_table[i] ^= new_item.parity_table[j]
                            new_item.output_table = item.output_table.copy()
                            new_item.output_table[i] ^= new_item.output_table[j]
                            new_item.ops_list = item.ops_list.copy()
                            new_item.ops_list.append((i,j))
                            #print(new_item.parity_table)
                            buffer.put(new_item)

                for (i, j) in item.ops_list:
                    insert_cnot(i, j)
                    #print(parity_table)
                    #print(output_table)
                    #print("~")
            index = np.nonzero(parity_table[:,col])[0][0]
            out.append_node(Gate("rz", [this.qubits[index]], [r_angles[col]])) 
            parity_table = np.delete(parity_table, col, 1)
            r_angles = np.delete(r_angles, col)

        #print(r_angles)
        '''print(parity_table)
        print(output_table)
        print("-")'''
        val = np.sum(parity_table, 0) == np.amin(np.sum(parity_table, 0))
        cols = np.argwhere(val)
        #print("testy", val, cols)
        while len(parity_table) > 0 and parity_table.shape[1]:
            cols = np.argwhere(sum(parity_table) == np.amin(sum(parity_table)))
            for i,_ in enumerate(this.qubits):
                ones = cols[np.where(parity_table[i, cols])[0]]
                if len(ones):
                    cols = ones
                    if len(cols) == 1:
                        break
            implement_rotation(cols[0][0])
            '''print(parity_table)
            print(output_table)
            print(";")'''

        circuit_l = guassian_elim(output_table, num_qubits)
        
        #print("LU", circuit_l, circuit_u)
        #print("")
        for i,j in circuit_l:
            out.append_node(Gate("cx", [this.qubits[j],this.qubits[i]]))
        '''print(out)
        print(this)
        print(out.convert_to_phasePoly())'''

        return out

    #synthesize, but the results replace the current circuit in place
    def systhesize_in_place(this, synthesis_type:Literal["pick_column_row_heap", "row_heap"], buffer_size:int = -1, ends_checked:int = 10, only_replace_on_improve:bool = True):
        result:Circuit
        if synthesis_type == "pick_column_row_heap":
            result = this.synthesize_by_col()
        elif synthesis_type == "row_heap":
            result = this.synthesize_row_search(buffer_size, ends_checked)
        
        if only_replace_on_improve:
            for circ in this.possible_circuits:
                if len(result.get_sequence()) > len(circ.get_sequence()):
                    result = circ

        result.replace_at(Circuit(ins=this.ins, outs=this.outs))

    def to_instr(this):
        out = "|" + ",".join([f"[{i}]" for i in this.qubits])+"> -> e^i("

        rotations = sorted(this.rotations, key=lambda x:list(sorted(x[1])))
        out += "+".join([str(t)+"*".join([f"[{i}]" for i in q]) for t,q in rotations]) + ")|"
        out += ",".join(["*".join([f"[{i}]" for i in this.affineOut[q]]) for q in this.qubits])
        return out + '>;'
    
    def copy_disconnected(this):
        out = PhasePoly(this.qubits, this.rotations, this.affineOut)
        out.possible_circuits = this.possible_circuits.copy()
        return out
    
    def is_equivalent_to(this, other_node:Self, ignore_qubits:bool = False):
        if not super().is_equivalent_to(other_node, ignore_qubits = False):
            return False
        return all(this.affineOut[q] == other_node.affineOut[q] for q in this.qubits) and set(this.rotations) == set(other_node.rotations)

#actual circuit class. Should work as a graph.
class Circuit:
    # For input params, its either (size:int) or (qubits:list[int]) for a new circuit, or (ins:dict[int,Node], outs,dict[int,Node]) for a subcircuit. 
    
    # Note: It is suggested against making multiple subcircuits on the same circuit due to how the ins/outs system works.
    # Any alternate solutions to define circuits in ways other than exclusive ins/outs are welcome, but remember that this should still support empty circuits.
    def __init__(this, size:(int|None) = None, qubits:(list[int]|None) = None, ins:(dict[int,Node]|None) = None, outs:(dict[int,Node]|None) = None):
        if ins != None and outs != None:
            if set(ins.keys()) != set(outs.keys()):
                raise Exception("Invalid")
            this.qubits = sorted(list(ins.keys()))
            this.ins = ins.copy()
            this.outs = outs.copy()
            this.is_subcircuit = True #temp measure
            #this.get_sequence()
            return
        if (size == None) == (qubits == None):
            raise Exception("Invalid")
        if qubits != None:
            this.qubits = qubits
        else:
            this.qubits = list(range(size))
        this.ins:dict[int,Node] = dict()
        this.outs:dict[int,Node] = dict()
        this.is_subcircuit = False
        for q in this.qubits:
            this.ins[q] = Node([q])
            this.outs[q] = Node([q])
            this.ins[q].connect_to(this.outs[q], q)

    #returns all nodes in the circuit as a list of nodes, ordered with dependency in mind.
    def get_sequence(this):
        sequence:list[Node] = []
        tracked = set()
        buffer:deque[Node] = deque()
        for q, in_node in this.ins.items():
            tracked.add(in_node)
            node = in_node.next(q)
            if in_node not in node.ins.values():
                raise Exception(f"Improper Graph Structure: {this.qubits}, {this.ins} -> {this.outs}")
            if node in this.ins.values():
                raise Exception(f"Improper Graph Structure: {this.qubits}, {this.ins} -> {this.outs}")
            if node not in buffer and node not in tracked:
                buffer.append(node)

        skip_streak = 0
        #print("\nget_seq buffer:")
        while len(buffer) > skip_streak:
            #print(buffer)
            target = buffer.popleft()
            if not all(prev_node in tracked for prev_node in target.ins.values()) and target not in this.outs.values():
                skip_streak += 1
                buffer.append(target)
                continue
            skip_streak = 0
            if target not in this.outs.values():
                for node in target.outs.values():
                    if target not in node.ins.values():
                        raise Exception(f"Improper Graph Structure: {this.qubits}, {this.ins} -> {this.outs}")
                    if node not in buffer and node not in tracked:
                        buffer.append(node)
                sequence.append(target)
            tracked.add(target)

        if len(buffer) != 0:
            raise Exception(f"Improper Graph Structure: {this.qubits}, {this.ins} -> {this.outs}\n buffer:{buffer}")
        return sequence

    def append_node(this, new_node:Node):
        if any(q not in this.qubits for q in new_node.qubits):
            raise Exception(f"Invalid Qubit input")
        
        prev_nodes:dict[int,Node] = dict()
        for q in new_node.qubits:
            prev_nodes[q] = this.outs[q].prev(q)
        new_node.attach_after(prev_nodes)

    #TODO: Untested
    def append_circuit(this, new_circuit:Self):
        if any(q not in this.qubits for q in new_circuit.qubits):
            raise Exception("Invalid Qubit input")
        
        for q,n in new_circuit.ins.items():
            this.outs[q].prev(q).connect_to(n.next(q), q)
        for q,n in new_circuit.outs.items():
            n.prev(q).connect_to(this.outs[q], q)

    
    #adds gate to end of circuit
    #deprecated in favor of append_node(Gate(...)) 
    '''
    def append_new_gate(this, gateType:str, qubits:list[int] = None, params:list[float] = None):
        if any(q not in this.qubits for q in qubits):
            raise Exception("Invalid Qubit input")
        
        new_gate:Gate = Gate(gateType, qubits, params)
        prev_nodes:dict[int,Node] = dict()
        for q in qubits:
            prev_nodes[q] = this.outs[q].prev(q)
        new_gate.attach_after(prev_nodes)
    '''
    #remaps the circuits qubits (very jank)
    #In particlur, the function input parameter is very awkward to use. See remap_demo.py
    def remap(this, mapping):
        if this.is_subcircuit:
            raise Exception("Illegal Function") #yes, this is horrible design. I'll refactor this eventually.
        
        def remap_list(in_list:list[int]):
            out_list:list[int] = []
            for num in in_list: 
                out_list.append(mapping(num))
            return out_list
        
        def remap_dict(in_dict:dict[int]):
            out_dict:dict[int] = dict()
            for i, node in in_dict.items():
                new_i = mapping(i)
                if new_i in out_dict.keys():
                    raise Exception("Mapping error: QuBit Overlap")
                out_dict[new_i] = node
            return out_dict

        this.qubits = remap_list(this.qubits)

        for in_node in this.ins.values():
            in_node.qubits = remap_list(in_node.qubits)
            in_node.outs = remap_dict(in_node.outs)
        this.ins = remap_dict(this.ins)
        
        sequence = this.get_sequence()
        for node in sequence:
            node.qubits = remap_list(node.qubits)
            node.ins = remap_dict(node.ins)
            node.outs = remap_dict(node.outs)
        
        for out_node in this.outs.values():
            out_node.qubits = remap_list(out_node.qubits)
            out_node.ins = remap_dict(out_node.ins)
        this.outs = remap_dict(this.outs)
        
        this.get_sequence() # this is to verify the circuit structure, for debug purposes

    #creates a copy, with different node objects
    def copy(this):
        out_circuit = Circuit(qubits=this.qubits)
        for node in this.get_sequence():
            out_circuit.append_node(node.copy_disconnected())
        return out_circuit
    
    #copy the target circuit to the given inputs/outputs. 
    #Note that the parameter circuit is the circuit to be replaced and must match the given circuit's qubits.
    #the circuit object the function is called at will be unchanged.
    #Below function is completely untested as of note, use at your own risk (for now)
    def replace_at(this, subcircuit_to_replace:Self):
        if not (set(this.qubits) == set(subcircuit_to_replace.qubits)):
            raise Exception("Invalid Qubits for Parameters")
        sequence = this.get_sequence()
        for node in sequence:
            if not set(node.qubits).issubset(set(this.qubits)):
                raise Exception("Circuit to be placed is not independent")
        
        for q in subcircuit_to_replace.qubits:
            subcircuit_to_replace.ins[q].connect_to(subcircuit_to_replace.outs[q], q)
        for node in sequence:
            if not set(node.qubits).issubset(set(subcircuit_to_replace.qubits)):
                raise Exception("Circuit to replace is not independent")
            subcircuit_to_replace.append_node(node.copy_disconnected())

    #TODO
    def find_match(this, inputs:dict[int,Node]):
        subcirc:Circuit; mapping:list[tuple[int,int]]
        return subcirc, mapping

    #Checks to see if this and another circuit/subcircuit match gates. remap is an optional function to remap other_circuit
    def structural_match(this, other_circuit:Self, remap = None):
        if remap == None:
            remap = lambda x : x
        for q in this.qubits:
            c1ptr = this.ins[q].next(q)
            c2ptr = other_circuit.ins[remap(q)].next(remap(q))
            while True:
                if c1ptr == this.outs[q] and c2ptr == other_circuit.outs[remap(q)]:
                    break
                if c1ptr != this.outs[q] and c2ptr != other_circuit.outs[remap(q)]:
                    if not c1ptr.is_equivalent_to(c2ptr):
                        return False
                else:
                    return False
                c1ptr = c1ptr.next(q)
                c2ptr = c2ptr.next(remap(q))
        return True

    #Untested
    def partition_to_phasePoly(this):
        tracked = set()
        legal_gates = {"cx", "cnot", "x", "not", "rz"}

        for node in this.get_sequence():
            if node in tracked:
                continue
            tracked.add(node)
            if node.nodeType not in legal_gates:
                continue

            #print(f"$ {node}")

            seq_nums = dict((b,a) for a,b in enumerate(this.get_sequence()))

            current_ins:dict[int,Node] = node.ins.copy()
            current_outs:dict[int,Node] = node.outs.copy()
            pending:set[tuple[Node,bool]] = set()
            buffer = deque()
            for q in node.qubits:
                buffer.extend([(q, node.next(q), True), (q, node.prev(q), False)])

            def parse_anchor(qubit:int, node:Node, direction:bool):
                ptr = node
                if direction:
                    current_bound = current_outs
                else:
                    current_bound = current_ins
                while ptr != None and ptr.nodeType in legal_gates:
                    if ptr.nodeType in {"cnot", "cx"}:
                        other_q = ptr.qubits[0] if ptr.qubits[0] != qubit else ptr.qubits[1] 
                        if other_q not in current_ins.keys():
                            min_out_num = min(seq_nums.get(n,math.inf) for n in current_outs.values())
                            max_in_num = max(seq_nums.get(n,-1) for n in current_ins.values())

                            #print("min_out_num", min_out_num, "max_in_num", max_in_num, ptr, other_q)

                            #anti self-dependency measures

                            stop_flag = False
                            backbuffer = deque([ptr.prev(other_q)])
                            tracked_b = set()
                            while len(backbuffer) > 0:
                                #print("backbuffer", [f"({q},{seq_nums.get(q)})" for q in backbuffer])
                                item = backbuffer.popleft()
                                tracked_b.add(item)
                                if item == None:
                                    raise Exception("Out of bounds (this shouldn't happen)")
                                if item in current_outs.values():
                                    #print(f"back stopped at {item}")
                                    stop_flag = True
                                    break
                                if item in seq_nums.keys() and seq_nums[item] > min_out_num:
                                    for o in item.ins.values():
                                        if o not in tracked_b and o not in backbuffer:
                                            backbuffer.append(o)
                            if stop_flag:
                                pending.add((ptr,direction))
                                break
                            forwardbuffer = deque([ptr.next(other_q)])
                            tracked_f = set()
                            while len(forwardbuffer) > 0:
                                #print("forwardbuffer", [f"({q},{seq_nums.get(q)})" for q in forwardbuffer])
                                item = forwardbuffer.popleft()
                                tracked_f.add(item)
                                if item == None:
                                    raise Exception("Out of bounds (this shouldn't happen)")
                                if item in current_ins.values():
                                    #print(f"forward stopped at {item}")
                                    stop_flag = True
                                    break
                                if item in seq_nums.keys() and seq_nums[item] < max_in_num:
                                    for o in item.outs.values():
                                        if o not in tracked_f and (o not in forwardbuffer):
                                            forwardbuffer.append(o)
                            if stop_flag:
                                pending.add((ptr,direction))
                                break

                            #print("added")
                            current_outs[other_q] = ptr.next(other_q)
                            buffer.append((other_q, ptr.next(other_q), True))
                            current_ins[other_q] = ptr.prev(other_q)
                            buffer.append((other_q, ptr.prev(other_q), False))

                        elif ptr in [i[0] for i in pending]:
                            _, dir = [p for p in pending if p[0] == ptr][0]
                            pending.remove((ptr,dir))
                            #print("pended",other_q,ptr,dir)
                            if dir:
                                current_outs[other_q] = ptr.next(other_q)
                                buffer.append((other_q, current_outs[other_q], dir))
                            else:
                                current_ins[other_q] = ptr.prev(other_q)
                                buffer.append((other_q, current_ins[other_q], dir))
                        else:
                            pending.add((ptr,direction))
                            #print("pending",ptr,direction)
                            current_bound[qubit] = ptr
                            break
                    
                    tracked.add(ptr)
                    if direction:
                        ptr = ptr.next(qubit)
                    else:
                        ptr = ptr.prev(qubit)
                    current_bound[qubit] = ptr
                    '''print(buffer)
                    print("cin",current_ins)
                    print("cout",current_outs)
                    print("pending",pending)
                    print("@", Circuit(ins=current_ins, outs=current_outs))'''
                if ptr == None:
                    raise Exception("Out of bounds (this shouldn't happen)")

            while len(buffer) > 0:
                '''print()
                print("keys:", current_ins.keys())
                print("buffer:",buffer)
                print("pending:",pending)
                print("current_ins:",current_ins)
                print("current_outs:",current_outs)
                
                sub_seq = Circuit(ins=current_ins, outs=current_outs).get_sequence()
                for node in this.get_sequence():
                    if node in sub_seq:
                        print(">"+str(node))
                    else:
                        print(node)'''
                
                q,n,d = buffer.popleft()
                parse_anchor(q,n,d)
            
            '''
            print("\nFinal:")
            print("buffer:",buffer)
            print("pending:",pending)
            print("current_ins:",current_ins)
            print("current_outs:",current_outs)
            print(Circuit(ins=current_ins, outs=current_outs))
            sub_seq = Circuit(ins=current_ins, outs=current_outs).get_sequence()
            for node in this.get_sequence():
                if node in sub_seq:
                    print(">"+str(node))
                else:
                    print(node)
            '''
            Circuit(ins=current_ins, outs=current_outs).replace_with_phasePoly()
            '''print(Circuit(ins=current_ins, outs=current_outs))
            this.get_sequence()
            print("~~~~~~~~~")'''


    #Note: I tried to run the below on quartz/circuit/nam_rm_circs/gf2^32_mult.qasm
    #It did not terminate despite runnning for more than a few days.

    #Call only for circuit consisting only of Rz, CNOT/CX, and NOT/X gates
    #Untested, use at your own risk
    def convert_to_phasePoly(this):
        subcircuit_copy = Circuit(qubits=this.qubits)

        hanging_nots = defaultdict(lambda:False)
        rotations:list[tuple[float, set[int]]] = []
        affineOut:dict[int, set[int]] = dict()
        for q in this.qubits:
            affineOut[q] = {q}

        out = Circuit(qubits = this.qubits)
        for node in this.get_sequence():
            node:Gate
            if node.nodeType == "cx" or node.nodeType == "cnot":
                q0,q1 = node.qubits[0:2]
                affineOut[q1] = affineOut[q0] ^ affineOut[q1]
                hanging_nots[q1] = hanging_nots[q0] ^ hanging_nots[q1]
                subcircuit_copy.append_node(Gate('cx', [q0,q1]))
            elif node.nodeType == "x" or node.nodeType == "cnot":
                hanging_nots[node.qubits[0]] = not hanging_nots[node.qubits[0]]
            elif node.nodeType == "rz":
                q = node.qubits[0]
                deg = node.params[0]
                if hanging_nots[q]:
                    deg = -deg
                subcircuit_copy.append_node(Gate('rz', [q],[deg]))
                parity = affineOut[q].copy()
                for i in range(len(rotations)):
                    if rotations[i][1] == parity:
                        rotations[i] = (rotations[i][0]+deg, parity)
                        break
                else:
                    rotations.append((deg, parity))
            else:
                raise Exception("Illegal gate in subcircuit for phase polynomial")

        rotations = [t for t in rotations if t[0] != 0.0]

        poly_qubits = this.qubits.copy()
        
        for q in this.qubits:
            if any((q in r) for _,r in rotations):
                continue
            if any((q != k and q in v) for k,v in affineOut.items()):
                continue
            if affineOut[q] != {q}:
                continue
            poly_qubits.remove(q)
            affineOut.pop(q)
        
        phasepoly_gate = PhasePoly(poly_qubits,rotations,affineOut)
        phasepoly_gate.possible_circuits.add(subcircuit_copy)
        out.append_node(phasepoly_gate)

        for k,v in hanging_nots.items():
            if v:
                out.append_node(Gate("x",[k]))
        return out
    
    def replace_with_phasePoly(this):
        pp = this.convert_to_phasePoly()
        '''print("+=+")
        print(this)
        print("PP: ",pp)'''
        pp.replace_at(this)

    def s_and_t_to_rz(this):
        for g in this.get_sequence():
            g:Gate
            if g.nodeType == "t":
                g.nodeType = "rz"
                g.params = [math.pi/4]
            elif g.nodeType == "tdg":
                g.nodeType = "rz"
                g.params = [-math.pi/4]
            elif g.nodeType == "s":
                g.nodeType = "rz"
                g.params = [math.pi/2]
            elif g.nodeType == "sdg":
                g.nodeType = "rz"
                g.params = [-math.pi/2]
    
    def __repr__(this) -> str:
        out = "\n".join([n.__repr__() for n in this.get_sequence()])
        return "\n>"+out
        