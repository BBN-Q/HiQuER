#!/usr/bin/env python3

#
#Copyright 2022 Raytheon BBN
#
#Licensed under the Apache License, Version 2.0 (the "License");
#you may not use this file except in compliance with the License.
#
#You may obtain a copy of the License at
#   http://www.apache.org/licenses/LICENSE-2.0
#
#Unless required by applicable law or agreed to in writing, software
#distributed under the License is distributed on an "AS IS" BASIS,
#WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#See the License for the specific language governing permissions and
#limitations under the License.
#

#Script for rotation synthesis to the Clifford+T set using pyLIQTR's
#implementation of grid search.

import cirq
from cirq.contrib.qasm_import import circuit_from_qasm
from pyLIQTR.gate_decomp.cirq_transforms import *
from decimal import Decimal
import numpy as np 
import texttable
import argparse, sys, json

#Gates in the Clifford+T set

cliffordT = set([
    cirq.H,
    cirq.T,
    cirq.X,
    cirq.Y,
    cirq.Z,
    cirq.S,
    cirq.T,
    cirq.CNOT])
cliffordT |= set([cirq.inverse(g) for g in cliffordT])

def moment_nonCT_gates(moment):
    '''All non-CliffordT gates in a cirq.moment.
    '''
    nonCT = [op.gate for op in moment.operations if op.gate not in cliffordT]
    return any(nonCT)

def rebuild_moment_rzs(moment, precision=10):
    subcircuit = cirq.Circuit()
    #for op in tqdm(moment.operations, desc="Operations", leave=False):
    for op in moment.operations:
        if op.gate in cliffordT:
            subcircuit.append(op)
        else:
            #print(op.gate)
            subcircuit.append(decompose_diagonal_cirq(Decimal(op.gate._rads), precision, op.qubits[0]))
    return subcircuit

def cliffordT_circuit(circuit, precision=10):
    new_circuit = cirq.Circuit()
    #for jj in tqdm(range(len(circuit.moments)), desc="Moments"):
    for jj in range(len(circuit.moments)):
        if moment_nonCT_gates(circuit.moments[jj]):
            subs = rebuild_moment_rzs(circuit.moments[jj], precision=precision)
            new_circuit.append(subs)
        else:
            new_circuit.append(circuit.moments[jj])
    return new_circuit

def circuit_stats(circuit):
    
    gates = [op.gate for op in circuit.all_operations()]
    
    nQB = len(circuit.all_qubits())
    depth = len(cirq.Circuit(circuit.all_operations()))
    
    nGates = len(gates)
    nCNOT  = sum([g == cirq.CNOT for g in gates])
    nT     = sum([g == cirq.T or g == cirq.inverse(cirq.T) for g in gates])
    nNC    = sum([g not in cliffordT for g in gates])

    stats = {"Number of qubits":    nQB,
             "Depth":               depth,
             "Number of gates":     nGates,
             "Number of CNOTs":     nCNOT,
             "T gate count":        nT,
             "Non-clifford count":  nNC}
    
    return stats

def print_stats(stats, names=None):
    if not names:
        names = [f"Circuit {n}" for n in len(stats)]
    if len(names) < len(stats):
        N = len(names)
        for j in range(N, len(stats)):
            names.append(f"Circuit {j}")

    ncols = len(stats)
    table = Texttable()
    table.set_cols_align(["l"] + ncols*['c'])
    table.set_cols_dtype(["t"] + ncols*['i'])
    table.add_rows([[k] + [s[k] for s in stats] for k in stats[0].keys()])
    print(table.draw())

def read_qasm(file):
    #f = open(file, "r")
    lines = [line for line in file.readlines() if "reset" not in line] #nasty hack
    qasm = "".join(lines)
    file.close()
    return circuit_from_qasm(qasm)

parser = argparse.ArgumentParser(description='Rotation Synthesis')

parser.add_argument("-f", '--file', type=argparse.FileType('r'), default='-')
parser.add_argument("-p", '--precision')
parser.add_argument("-e", '--error', default=1e-4)
parser.add_argument("-s", '--stats', action='store_true')

args = parser.parse_args()

if __name__ == "__main__":

    qasm = read_qasm(args.file)


    initial_stats = circuit_stats(qasm)

    if args.precision:
        prec = int(args.precision )
    else:
        #Guess that there are on average 200 T gates per non-Clifford
        prec = int(np.abs(np.log10(0.5*args.error/(initial_stats["Non-clifford count"]*200))))

    circ_ct = cliffordT_circuit(qasm, precision=prec)

    new_stats = circuit_stats(circ_ct)

    sys.stdout.write(cirq.qasm(circ_ct))

    if args.stats:

        stats = {"Initial Circuit": initial_stats, "Final Circuit": new_stats}

        with open("stats.json", "w") as f:
            json.dump(stats, f, indent=4)









