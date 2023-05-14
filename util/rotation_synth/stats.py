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
from rotation_synth import read_qasm, circuit_stats
import argparse, sys, json

parser = argparse.ArgumentParser(description='Circuit Stats')

parser.add_argument("-f", '--file', type=argparse.FileType('r'), default='-')
parser.add_argument("-j", '--json', action='store_true', default=False)

args = parser.parse_args()

if __name__ == "__main__":

    qasm = read_qasm(args.file)
    stats = circuit_stats(qasm)
    
    if True:
        json.dump(stats, sys.stdout, indent=4)
        exit(0)
        
    print("---- Circuit Statistics ----")
    L = max([len(x) for x in stats.keys()])
    for k,v in stats.items():
        print(k + ":" + " "*(L+4-len(k)) + str(v))
        
    
    




    