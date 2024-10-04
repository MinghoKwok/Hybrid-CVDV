from typing import Any
from collections import deque
import math
import qiskit
import sys
import os

sys.path.append(os.path.abspath('/common/home/mg1998/CVDV/bosonic-qiskit'))

import c2qa



def read_qasm(filename: str):
    out: c2qa.CVCircuit
    out = None
    # Register preparation
    qbr = None
    qmr = None
    cr = None

    with open(filename, 'r') as fp:
        for line in fp:
            comment_index = line.find(r"//")
            if comment_index != -1:
                line = line[:comment_index]
            
            if line.strip() == "":
                continue

            command, args_str = line.strip().strip(";").split(" ", 1)

            print("\ncommand: " + command)
            print("args_str: " + args_str)

            if command in ['qbreg', 'qmreg', 'creg']:
                if command == "qbreg":
                    arg = int(args_str.strip('qbr[]'))
                    print(f"qb_arg: {arg}")
                    qbr = qiskit.QuantumRegister(arg)
                    continue
                elif command == "creg":
                    arg = int(args_str.strip('cr[]'))
                    print(f"c_arg: {arg}")
                    cr = qiskit.ClassicalRegister(arg)
                    continue
                elif command == "qmreg":
                    arg = int(args_str.strip('qmr[]'))
                    print(f"qm_arg: {arg}")
                    qmr = c2qa.QumodeRegister(num_qumodes = 8)
                    continue
            else:
                if out == None:
                    out = c2qa.CVCircuit(qmr, qbr, cr)
                    print(f"==========={out.qregs[0]}")
                
                args = [arg.strip() for arg in args_str.split(',')]

                print([arg for arg in args])

                if command == "bsCV":
                    cmd = f"out.cv_bs(math.pi * 2, {args[2]}, {args[3]})"
                    exec(cmd)
                    continue
                elif command == "rCV":
                    cmd = f"out.cv_r(theta = - math.pi/2, qumode={args[1]})"
                    exec(cmd)
                    continue
                elif command == "cdHybrid":
                    cmd = f"out.cv_c_d(alpha=1.0, qumode={args[0]}, qubit={args[1]})"
                    continue
                elif command == "measure":
                    # cmd = "out.cv_measure(" + "qmr[0]" + "," + "cr[0]" + ")"
                    # exec(cmd)
                    out.measure(qbr[0], cr[0])

                



                    


            # if command in ["bsCV", "rCV"]:
            #     out.append_node(circuits.Gate(command, qbit_vals, params))
            # elif command == "measure" and cbits is not None:
            #     out.append_node(circuits.Measure(qbit_vals[0], cbit_vals[0]))
            # elif command == "cdHybrid":
            #     out.append_node(circuits.HybridGate(command, qbit_vals))
            # else:
            #     raise ValueError(f"Unrecognized command: {command}")
                
    return out


def write_qasm(circuit: c2qa.CVCircuit, filename: str):
    print(circuit.data)
    print(f"\nqmregs: {len(circuit.qmregs)}")
    print(f"\n_qubit_regs: {len(circuit._qubit_regs)}")


    with open(filename, 'w') as fp:
        qbr_size = len(circuit.qregs[0]) if circuit.qregs else 0
        cr_size = len(circuit.cregs[0]) if circuit.cregs else 0
        qmreg_size = len(circuit.qmregs[0]) if circuit.qmregs else 0

        if qbr_size > 0:
            fp.write(f"qbreg qbr[{qbr_size}];\n")
        if cr_size > 0:
            fp.write(f"creg cr[{cr_size}];\n")
        if qmreg_size > 0:
            fp.write(f"qmreg qmr[{qmreg_size}];\n")

        for inst, qargs, cargs in circuit.data:
            # Check the instruction name and convert it to the appropriate QASM command
            if inst.name == 'BS':
                theta = inst.params[0]  # Assuming `inst.params[0]` holds the theta value for `bsCV`
                fp.write(f"bsCV {theta}, 0, {qargs[0]._register.name}[{qargs[0]._index}], {qargs[1]._register.name}[{qargs[1]._index}];\n")
            elif inst.name == 'R':
                theta = inst.params[0]  # Assuming `inst.params[0]` holds the theta value for `rCV`
                fp.write(f"rCV {theta}, {qargs[0]._register.name}[{qargs[0]._index}];\n")
            elif inst.name == 'cv_c_d':
                alpha = inst.params[0]  # Assuming `inst.params[0]` holds the alpha value for `cdHybrid`
                fp.write(f"cdHybrid {qargs[0]._register.name}[{qargs[0]._index}], {qargs[1]._register.name}[{qargs[1]._index}];\n")
            elif inst.name == 'measure':
                fp.write(f"measure {qargs[0]._register.name}[{qargs[0]._index}], {cargs[0]._register.name}[{cargs[0]._index}];\n")
            else:
                raise ValueError(f"Unrecognized instruction: {inst.name}")





'''
def math_str_to_float(math_str: str):
    try:
        out = float(math_str)
        return out
    except ValueError:
        pass
    
    tokens: deque[tuple[str, Any]] = deque()
    ptr = 0
    while ptr < len(math_str):
        c: str = math_str[ptr]
        if c in "-*/":
            tokens.append((c, None))
        elif c.isspace():
            ptr += 1
            continue
        else:
            if c == "p" and len(math_str) - ptr >= 2 and math_str[ptr:ptr+2] == "pi":
                tokens.append(("num", math.pi))
                ptr += 1
            elif c in "1234567890.":
                start_ptr = ptr
                while ptr < len(math_str) - 1 and math_str[ptr + 1] in "1234567890.":
                    ptr += 1
                tokens.append(("num", float(math_str[start_ptr:ptr + 1])))
            else:
                raise Exception("Parsing Error")
            if len(tokens) > 1 and tokens[-2][0] == "-":
                _, num = tokens.pop()
                tokens.pop()
                tokens.append(("num", -num))
        ptr += 1

    while len(tokens) > 1:
        if len(tokens) == 2 or len(tokens) <= 0:
            raise Exception("Parsing Error")
        if tokens[0][0] == tokens[2][0] == "num":
            _, num1 = tokens.popleft()
            op, _ = tokens.popleft()
            _, num2 = tokens.popleft()
            if op == "*":
                tokens.appendleft(("num", num1 * num2))
            elif op == "/":
                tokens.appendleft(("num", num1 / num2))
            else:
                raise Exception("Parsing Error")
            continue
        raise Exception("Parsing Error")
    
    return tokens[0][1]
'''

# read_qasm("swap_example.cvdvqasm")
write_qasm(read_qasm("swap_example.cvdvqasm"), "convert.cvdvqasm")
