#!/usr/bin/env python3
import argparse
import json
import numpy as np
import sys

# # MSD factories
# 
# from Table 1 in Litinski's Magic State Distillation: Not as Costly as You Think (https://arxiv.org/pdf/1905.06903.pdf)

factory_1 = {
    "p_phys": 1e-4,
    "p_out": 4.4 * 1e-8,
    "qubits": 810,
    "cycles": 18.1
}

factory_2 = {
    "p_phys": 1e-4,
    "p_out": 9.3 * 1e-10,
    "qubits": 1150,
    "cycles": 18.1
}

factory_3 = {
    "p_phys": 1e-4,
    "p_out": 1.9 * 1e-11,
    "qubits": 2070,
    "cycles": 30.0
}

factory_4 = {
    "p_phys": 1e-4,
    "p_out": 2.4 * 1e-15,
    "qubits": 16400,
    "cycles": 90.3
}

factory_5 = {
    "p_phys": 1e-4,
    "p_out": 6.3 * 1e-25,
    "qubits": 18600,
    "cycles": 67.8
}

factory_6 = {
    "p_phys": 1e-3,
    "p_out": 4.5 * 1e-8,
    "qubits": 4620,
    "cycles": 42.6
}

factory_7 = {
    "p_phys": 1e-3,
    "p_out": 1.4 * 1e-10,
    "qubits": 43300,
    "cycles": 130
}

factory_8 = {
    "p_phys": 1e-3,
    "p_out": 2.6 * 1e-11,
    "qubits": 46800,
    "cycles": 157
}

factory_9 = {
    "p_phys": 1e-3,
    "p_out": 2.7 * 1e-12,
    "qubits": 30700,
    "cycles": 82.5
}

factory_10 = {
    "p_phys": 1e-3,
    "p_out": 3.3 * 1e-14,
    "qubits": 39100,
    "cycles": 97.5
}

factory_11 = {
    "p_phys": 1e-3,
    "p_out": 4.5 * 1e-20,
    "qubits": 73400,
    "cycles": 128
}

factory_12 = {
    "p_phys": 1e-4,

    "p_out": 1.5 * 1e-9,
    "qubits": 762,
    "cycles": 36.2
}

factory_13 = {
    "p_phys": 1e-3,
    "p_out": 6.1 * 1e-10,
    "qubits": 7780,
    "cycles": 469
}

factories = [factory_1, factory_2, factory_3, factory_4, factory_5, factory_6,
             factory_7, factory_8, factory_9,
             factory_10, factory_11, factory_12, factory_13]


# Functions to calculate required distance and data block parameters

def overhead_d_func(scheme, num_logical_qubits, T_count):
    if scheme == "fast":
        tiles = 2 * num_logical_qubits + np.sqrt(8 * num_logical_qubits) + 1
        timesteps = T_count

    if scheme == "compact":
        tiles = 1.5 * num_logical_qubits + 3
        timesteps = 9 * T_count

    overhead_d = {"space (d^2)": tiles, "time (d)": timesteps,
                  "spacetime (d^3)": tiles * timesteps}

    return overhead_d


def overhead_phys_func(overhead_d, p_phys, t_cycle):
    tiles = overhead_d["space (d^2)"]
    timesteps = overhead_d["time (d)"]

    d = 1
    while tiles * timesteps * 0.1 * (100 * p_phys) ** ((d + 1) / 2) > 0.01:
        d = d + 2

    qubits = tiles * 2 * d ** 2
    time = timesteps * t_cycle * d

    overhead_phys = {"space (physical qubits)": qubits, "time (sec)": time,
                     "spacetime (qubitsecs)": qubits * time,
                     "distance_used": d}

    return overhead_phys


# Function to determine best MSD factory

def best_factory_func(p_phys, T_count, d, scheme):
    if scheme == "compact":
        rate = 9

    if scheme == "fast":
        rate = 1

    # this is hacky
    best_factory = {"qubits": 1e20}

    num_facs = 1
    while best_factory["qubits"] == 1e20:
        for factory in factories:
            if factory["p_phys"] == p_phys:
                if factory['p_out'] * T_count < 0.01:
                    if (factory['cycles'] / d) / num_facs < rate:
                        if factory['qubits'] < best_factory['qubits']:
                            best_factory = factory

        num_facs = num_facs + 1

    return best_factory, num_facs


# Calculating the overall overhead

def overheads_func(scheme, num_logical_qubits, T_count, p_phys, t_cycle):
    overhead_d = overhead_d_func(scheme, num_logical_qubits, T_count)
    overhead_phys = overhead_phys_func(overhead_d, p_phys, t_cycle)
    best_factory, num_facs = best_factory_func(p_phys, T_count, overhead_phys["distance_used"], scheme)

    overheads = {
        "physical qubits": int(
            overhead_phys['space (physical qubits)'] + best_factory[
                'qubits'] * num_facs),
        "time (sec)": overhead_phys['time (sec)'],
        "time (min)": overhead_phys['time (sec)'] / 60,
        "time (hours)": overhead_phys['time (sec)'] / 3600,
        "time (days)": overhead_phys['time (sec)'] / (24 * 3600),
        "time (months)": overhead_phys['time (sec)'] / (24 * 3600 * 30),
        "time (years)": overhead_phys['time (sec)'] / (24 * 3600 * 30 * 12),
        "distance used": overhead_phys["distance_used"],
        "fraction factories": np.round((best_factory['qubits'] * num_facs) / (
                overhead_phys['space (physical qubits)'] + best_factory['qubits'] * num_facs), 2),
        "number of factories": num_facs,
        "factory used": best_factory
    }

    return overheads


def compute_t_layers(lattice_surgery_output):
    circuit = lattice_surgery_output['3. Circuit after T depth reduction']
    if circuit is None:
        print("Error: Circuit after T depth reduction: null, exiting...", file=sys.stderr)
        sys.exit(-1)

    n_qubits = int(lattice_surgery_output['3. Circuit after T depth reduction']['n'])
    t_layers = lattice_surgery_output['3. Circuit after T depth reduction']['T layers']
    t_counts = []
    if t_layers is not None:
        for layer in t_layers:
            t_counts.append(len(layer))

    return {"n_qubits": n_qubits, "Commuting T-layers": len(t_counts),
            "T-count": sum(t_counts), "T-count/layer": t_counts}


def qre(lattice_surgery_output, config):
    t_layers = compute_t_layers(lattice_surgery_output)
    n_qubits = t_layers["n_qubits"]
    t_counts = t_layers["T-count/layer"]
    t_count = sum(t_counts)
    commuting_t_layers_no = len(t_counts)
    logical_circuit = {"n_qubits": n_qubits, "T-count": t_count, "Commuting T-layers": commuting_t_layers_no}
    # logical_circuit["T-count/layer"] = t_layers["T-count/layer"]
    overheads = overheads_func(config["scheme"], n_qubits, t_count, config["p_g"], config["cycle_time"])

    return logical_circuit, overheads


# Initialize parser
parser = argparse.ArgumentParser(description='Quantum resource estimator',
                                 epilog='(c) 2023 softwareQ Inc. All rights reserved.')

# Adding required arguments
parser.add_argument("-f", "--file", help="Reads from file (by default reads from stdin)")
parser.add_argument("-c", "--config", help="QECC physical parameters configuration (by default, uses a fast scheme)")

# Read arguments from command line
args = parser.parse_args()

if __name__ == '__main__':
    lattice_surgery_json = {}
    if args.file:
        with open(args.file) as f:
            lattice_surgery_json = json.loads(f.read())
    else:
        lattice_surgery_json = json.load(sys.stdin)

    default_scheme = {
        "scheme": "fast",
        "p_g": 1e-3,
        "cycle_time": 1e-6
    }
    if args.config:
        with open(args.config) as f_config:
            config_json = json.loads(f_config.read())
    else:
        config_json = default_scheme

    logical, physical = qre(lattice_surgery_json, config_json)

    final_result = {"logical": logical, "qecc": config_json, "physical": physical}
    print(json.dumps(final_result, indent=4))
