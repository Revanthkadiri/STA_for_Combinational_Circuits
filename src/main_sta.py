import argparse  # parsing command-line arguments
from collections import Counter  # Counting occurrences
import time

tic = time.perf_counter()


def read_ckt(input_path):
    """Reads Circuit details from Bench file"""
    inputs = {}  # Creating Dictionary for Inputs
    outputs = {}  # Creating Dictionary for Outputs
    gates = {}  # Creating Dictionary for Gates
    output_ids_counter = Counter()  # To count OUTPUT IDs
    input_ids_counter = Counter()  # To count INPUT IDs
    primary_inputs = set()  # Initialize primary inputs set
    direct_connections = []  # Initialize direct connections list

    # Reads and analyses the file
    with open(input_path, 'r') as file:
        for line in file:
            line = line.strip()
            if not line or line.startswith("#"):  # Skip comments starting with "#" and empty lines
                continue

            # If the line starts with the keyword INPUT
            if line.startswith("INPUT"):
                input_id = line.split("(")[1].rstrip(")")
                inputs[input_id] = {"outputs": []}
                input_ids_counter[input_id] += 1
                primary_inputs.add(input_id)

            # If the line starts with the keyword OUTPUT
            elif line.startswith("OUTPUT"):
                output_id = line.split("(")[1].rstrip(")")
                outputs[output_id] = {"inputs": []}  # Ensure outputs have an "inputs" key initialized as a list
                output_ids_counter[output_id] += 1

            # Circuit Connections
            else:
                parts = line.split(" = ")
                gate_output, operation_with_inputs = parts[0], parts[1]
                operation, gate_inputs_str = operation_with_inputs.split("(")
                gate_inputs = gate_inputs_str.rstrip(")").split(", ")
                gates[gate_output] = {"type": operation, "inputs": gate_inputs, "outputs": []}

    # Updates each gate's outputs according to its inputs
    for gate_id, gate_info in gates.items():
        for input_id in gate_info["inputs"]:
            if input_id in gates:
                gates[input_id]["outputs"].append(gate_id)
            elif input_id in inputs:
                inputs[input_id]["outputs"].append(gate_id)

    # Connects primary outputs to the gate(s) generating them
    for gate_id, gate_info in gates.items():
        if gate_id in outputs:  # Correctly connect outputs to their generating gate
            outputs[gate_id]["inputs"].append(gate_id)

    # Identify direct connections between inputs and outputs
    for output_id in outputs:
        if output_id in inputs:
            direct_connections.append((f"INPUT_{output_id}", f"OUTPUT_{output_id}"))


    return inputs, outputs, gates, output_ids_counter, input_ids_counter, primary_inputs, direct_connections

def write_ckt(inputs, outputs, gates, input_ids_counter, output_ids_counter, output_file_path):
    """Write circuit details"""
    with open(output_file_path, 'w') as file:
        file.write(f"{sum(input_ids_counter.values())} primary inputs\n")
        file.write(f"{sum(output_ids_counter.values())} primary outputs\n")

        gate_types = {}
        for gate_id, gate_info in gates.items():
            gate_type = gate_info["type"]
            gate_types[gate_type] = gate_types.get(gate_type, 0) + 1  # Count and update repetition of each gate type.
        for gate_type, count in gate_types.items():
            if count >> 1:
                file.write(f"{count} {gate_type} gates\n")
            else:
                file.write(f"{count} {gate_type} gate\n")  # Write the count of each gate type to the output file

        # Write fanout information
        file.write("\nFanout...\n")
        for gate_id, gate_info in gates.items():
            output_parts = []  # Creating a list to obtain output_parts

            # Check and collect gate outputs
            if gate_info["outputs"]:
                unique_outputs = list(dict.fromkeys(gate_info["outputs"]))
                outputs_list = ', '.join([f"{gates[out_id]['type']}-{out_id}" for out_id in unique_outputs])
                output_parts.append(outputs_list)

            # Verify and store primary output without any fanout
            if gate_id in outputs:
                output_parts.append(f"OUTPUT-{gate_id}")

            # Writing to file if anything is available to print
            if output_parts:
                file.write(f"{gate_info['type']}-{gate_id}: {', '.join(output_parts)}\n")

            else:
                # Print only Gates if it has no output(s)
                file.write(f"{gate_info['type']}-{gate_id}\n")

        # Write fan-in information
        file.write("\nFanin...\n")
        for gate_id, gate_info in gates.items():
            unique_inputs = list(dict.fromkeys(gate_info["inputs"]))
            inputs_list = ', '.join(
                [f"INPUT-{inp_id}" if inp_id in inputs else f"{gates[inp_id]['type']}-{inp_id}" for inp_id in
                 unique_inputs])
            file.write(f"{gate_info['type']}-{gate_id}: {inputs_list}\n")


# Reading the Liberty File
def open_NLDM(input_path):
    with open(input_path, 'r') as file:
        return file.read()

def read_delays(NLDM_file):
    """Read delay delatils from liberty file"""
    lines = NLDM_file.split('\n')
    NLDM_details = []  # emtpy list to store input slew, load capacitance and delay
    cell = None  # checker to check if it is inside cell
    if_delay = False  # checker to check if it is inside delay
    multi_line = False  # For iterating through multiple line in delay
    delays = ""  # empty string to extract delay values

    for line in lines:
        line = line.strip()
        if line.startswith("cell ("):
            if cell:
                NLDM_details.append(cell)
            cell_name = line.split("(")[1].split(")")[0].strip()
            cell = {"cell": cell_name, "input_slews": "", "load_cap": "", "delay": []}
        elif "cell_delay" in line and cell is not None:
            if_delay = True
        elif line.startswith("index_1") and if_delay:  # condition to get input slew
            cell["input_slews"] = line.split('"')[1]
        elif line.startswith("index_2") and if_delay:  # condition to get load capacitance
            cell["load_cap"] = line.split('"')[1]
        elif if_delay and line.startswith("values"):  # condition to get delays
            delays = " " + line.split("values")[1].lstrip().replace('(', '').replace('"', '').replace(', \\', '\n')
            if not line.endswith(";"):
                multi_line = True
            else:
                cell["delay"].append(delays)
                if_delay = False
        elif multi_line:
            delays += " " + line.replace(', \\', '\n').replace('"', '').replace('(', '').replace(')', '').replace(';',
                                                                                                                  '')
            if line.endswith(";"):
                multi_line = False
                if_delay = False
                cell["delay"].append(delays)
        if line == "}":
            if_delay = False
            multi_line = False
    if cell:
        NLDM_details.append(cell)
    return NLDM_details  # return NLDM details

def write_delays(NLDM_details, output_path):
    """Writes delay LUT"""
    with open(output_path, 'w') as file:
        for cell in NLDM_details:
            file.write(f"cell: {cell['cell']}\n")
            file.write(f"input slews: {cell['input_slews']}\n")
            file.write(f"load cap: {cell['load_cap']}\n")
            file.write("delays:")
            file.write("\n")
            for delay_line in cell['delay']:
                write_delay = [line + ";" for line in delay_line.split('\n') if line.strip()]
                file.write('\n'.join(write_delay) + "\n")
            file.write("\n")

def read_slews(NLDM_file):
    """Read Slews and capacitance's of cell from liberty file"""
    lines = NLDM_file.split('\n')
    NLDM_details = []  # emtpy list to store input slew, load capacitance and slew
    cell = None  # checker to check if it is inside cell
    if_slew = False  # checker to check if it is inside slew
    multi_line = False  # For iterating through multiple line in delay
    output_slews = ""  # empty string to extract slew values
    cell_capacitance_load = {}
    current_cell = None

    for line in lines:
        line = line.strip()
        if line.startswith("cell ("):
            current_cell = line.split(" ")[1].strip("()")
            if cell:
                NLDM_details.append(cell)
            cell_name = line.split("(")[1].split(")")[0].strip()
            cell = {"cell": cell_name, "input_slews": "", "load_cap": "", "slews": []}
        elif "output_slew" in line and cell is not None:
            if_slew = True
        elif line.startswith("index_1") and if_slew:  # condition to get input slews
            cell["input_slews"] = line.split('"')[1]
        elif line.startswith("index_2") and if_slew:  # condition to get load capacitance
            cell["load_cap"] = line.split('"')[1]
        elif line.startswith("capacitance") and current_cell:
            capacitance_value = float(line.split(":")[1].strip(';'))
            cell_capacitance_load[current_cell.split('2')[0]] = capacitance_value
            # cell_capacitance
            current_cell = None  # Reset after capturing capacitance
        elif if_slew and line.startswith("values"):  # condition to get output slews
            output_slews = " " + line.split("values")[1].lstrip().replace('(', '').replace('"', '').replace(', \\',
                                                                                                            '\n')
            if not line.endswith(";"):
                multi_line = True
            else:
                cell["slews"].append(output_slews)
                if_slew = False
        elif multi_line:
            output_slews += " " + line.replace(', \\', '\n').replace('"', '').replace('(', '').replace(')', '').replace(
                ';', '')
            if line.endswith(";"):
                multi_line = False
                if_slew = False
                cell["slews"].append(output_slews)
        if line == "}":
            if_slew = False
            multi_line = False
    if cell:  # Append last cell value if any
        NLDM_details.append(cell)
    # print(cell_capacitance_load)
    return NLDM_details, cell_capacitance_load

def write_slews(NLDM_details, output_path):
    """Writes Slew LUT"""
    with open(output_path, 'w') as file:
        for cell in NLDM_details:
            file.write(f"cell: {cell['cell']}\n")
            file.write(f"input slews: {cell['input_slews']}\n")
            file.write(f"load cap: {cell['load_cap']}\n")
            file.write("slews:")
            file.write("\n")
            for slew_line in cell['slews']:
                slew_lines_with_semicolon = [line + ";" for line in slew_line.split('\n') if line.strip()]
                file.write('\n'.join(slew_lines_with_semicolon) + "\n")
            file.write("\n")

def gate_dict(output_ids_counter):
    """ Feteching Parameters for individual gates: fanin, fanouts, input & output connection details """
    gate_dict = {}  # Initialize dictionary for gates
    fain_dict = {}  # Initialize dictionary for fanins
    output_connected_gates = {}  # Initialize dictionary for output connections
    input_connected_gates = {}  # Initialize dictionary for input connections
    fanout_list = []  # Initialize fanout_list

    for gate_id, gate_info in gates.items():
        key = f"{gate_info['type']}-{gate_id}"
        output_parts_gates = [] # Stores Fanouts of gates
        output_parts_outputs = []  # stores nodes which has OUTPUT node as Fanout
        output_parts = []  # Stores all Fanouts of gates
        input_parts = []  # Stores Fanins of gates

        # Process outputs
        if gate_info["outputs"]:
            # print(gate_info["outputs"])
            unique_outputs = [gates[out_id]['type'] + '-' + out_id for out_id in gate_info["outputs"]]
            # print(unique_outputs)
            outputs_list = ', '.join(unique_outputs)
            output_parts_gates.append(outputs_list)
            # print(output_parts)
            fanout_list.append(unique_outputs)  # Add to fanout_list as a set

        if gate_id in outputs:
            while output_ids_counter[gate_id]:
                output_parts_outputs.append(f"OUTPUT-{gate_id}")
                fanout_list.append({f"OUTPUT-{gate_id}"})
                output_ids_counter[gate_id] = output_ids_counter[gate_id] - 1
        if output_parts_gates:
            combined_outputs = output_parts_gates[0]
            # If there are output parts in output_parts_outputs, append them
            if output_parts_outputs:
                combined_outputs += ', ' + ', '.join(output_parts_outputs)
            # Add the combined string to output_parts
            output_parts.append(combined_outputs)
        else:
            # If output_parts_gates is empty but there are items in output_parts_outputs, use them instead
            if output_parts_outputs:
                output_parts.append(', '.join(output_parts_outputs))

        gate_dict[key] = output_parts
        # print(fanout_list)

        # Process inputs
        if gate_info["inputs"]:
            unique_inputs = [f"INPUT-{inp_id}" if inp_id in inputs else gates[inp_id]['type'] + '-' + inp_id
                             for inp_id in list(dict.fromkeys(gate_info["inputs"]))]
            inputs_list = ', '.join(unique_inputs)
            input_parts.append(inputs_list)

            for inp in unique_inputs:
                if inp.startswith("INPUT"):
                    if inp in input_connected_gates:
                        input_connected_gates[inp].append(key)
                    else:
                        input_connected_gates[inp] = [key]

        fain_dict[key] = input_parts
    for gate, gate_outputs in gate_dict.items():
        if gate_outputs:  # Check if gate_outputs is not empty
            input_list = gate_outputs[0].split(', ')
            for input_source in input_list:
                if "OUTPUT" in input_source:
                    if input_source in output_connected_gates:
                        output_connected_gates[input_source].append(gate)
                    else:
                        output_connected_gates[input_source] = [gate]

    fanout_list = list(gate_dict.values())
    return gate_dict, fanout_list, fain_dict, input_connected_gates, output_connected_gates


def load_cap(cell_capacitance_load, gate_dict, fanout_list, input_connected_gates):
    """ Computing Load Capacitance for all gates """

    cap_list = []  

    def is_connected_to_output(fanout_group):
        """ Checking for gates connected to the OUTPUT node """
        return any('OUTPUT' in fanout for fanout in fanout_group)

    for fo in fanout_list:
        fanout_cap_values = []
        for fanout_group in fo:
            fanout_group = fanout_group.split(', ')  # Split the fanout string into individual fanouts
            connected_to_output = is_connected_to_output(fanout_group)
            for single_fanout in fanout_group:
                gate_type = single_fanout.split('-')[0]
                if gate_type in ["BUFF", "NOT", "OUTPUT"]:
                    single_fanout = "INV_X1" if gate_type != "BUFF" else "BUF_X1"  # modify names as to match with cell names in liberty file
                gate_type = single_fanout.split('-')[0]
                capacitance_value = cell_capacitance_load.get(gate_type,
                                                              0)  # Directly use the value from the dictionary
                if connected_to_output and gate_type == "INV_X1":
                    capacitance_value *= 4

                fanout_cap_values.append(capacitance_value)
        cap_list.append(sum(fanout_cap_values))  # Add the nodes total capacitance to the list

    final_dict = dict(zip(gate_dict.keys(), cap_list))  # map gate and capcitance

    for input_gate, connections in input_connected_gates.items():  # Setting up capcitances for INPUT nodes
        cap_values = []  # To store capacitance values for this input
        for connection in connections:
            gate_type = connection.split('-')[0]
            if gate_type in ["BUFF", "NOT", "OUTPUT"]:
                connection = "INV_X1" if gate_type != "BUFF" else "BUF_X1"
            gate_type = connection.split('-')[0]
            capacitance_value = cell_capacitance_load.get(gate_type, 0)
            cap_values.append(capacitance_value)
            final_dict[input_gate] = sum(cap_values)  # Update the final dict with total capacitance for the input
    return final_dict


def bilinear_interpolation(C, t, C1, C2, t1, t2, v11, v12, v21, v22):
    """Performs Bilinear Interpolation"""

    dC = (C - C1) / (C2 - C1)
    dt = (t - t1) / (t2 - t1)
    return (v11 * (1 - dC) * (1 - dt) + v21 * dC * (1 - dt) + v12 * (1 - dC) * dt + v22 * dC * dt)


def get_delay(cell_data, input_slew, load_cap):
    """
    Finds the delay for a given input_slew and load_cap using bilinear interpolation, considering edge cases.
    """
    # Data preparation for interpolation
    input_slews = [float(slew) for slew in cell_data['input_slews'].split(',')]
    load_caps = [float(cap) for cap in cell_data['load_cap'].split(',')]

    delay_values = []  # list to store delays in a linear format
    for delay_line in cell_data['delay']:
        slew_values = delay_line.replace('\n', ',').split(',')
        delay_values.extend([float(value.strip()) for value in slew_values if value.strip()])

    slew_L = max([i for i in range(len(input_slews)) if input_slews[i] <= input_slew] + [
        0])  # Finding left most index for input slew
    slew_R = min([i for i in range(len(input_slews)) if input_slews[i] >= input_slew] + [
        len(input_slews) - 1])  # Finding right most index for input slew
    cap_L = max(
        [j for j in range(len(load_caps)) if load_caps[j] <= load_cap] + [0])  # Finding left most index for load cap
    cap_R = min([j for j in range(len(load_caps)) if load_caps[j] >= load_cap] + [
        len(load_caps) - 1])  # Finding right most index for load cap

    # Handle edge cases by adjusting indices
    if input_slew < input_slews[0]: slew_L, slew_R = 0, 1  # if input slew is less than minimum value in lUT
    if input_slew > input_slews[-1]: slew_L, slew_R = len(input_slews) - 2, len(
        input_slews) - 1  # if input slew is more than maximum value in lUT
    if load_cap < load_caps[0]: cap_L, cap_R = 0, 1  # if load cap is less than minimum value in lUT
    if load_cap > load_caps[-1]: cap_L, cap_R = len(load_caps) - 2, len(
        load_caps) - 1  # if load cap is more than maximum value in lUT

    # Retrieve values for interpolation
    v11 = delay_values[(slew_L * len(load_caps)) + cap_L]  # Get delay values wrt slew_L and cap_L
    v12 = delay_values[(slew_L * len(load_caps)) + cap_R]  # Get delay values wrt slew_L and cap_R
    v21 = delay_values[(slew_R * len(load_caps)) + cap_L]  # Get delay values wrt slew_R and cap_L
    v22 = delay_values[(slew_R * len(load_caps)) + cap_R]  # Get delay values wrt slew_R and cap_R

    delay = bilinear_interpolation(input_slew, load_cap, input_slews[slew_L], input_slews[slew_R], load_caps[cap_L],
                                   load_caps[cap_R], v11, v12, v21, v22)  # Perform Bilinear Interpolation
    return delay


def get_slew(cell_data, input_slew, load_cap):
    # Data preparation for interpolation
    input_slews = ([float(slew) for slew in cell_data['input_slews'].split(',')])
    load_caps = ([float(cap) for cap in cell_data['load_cap'].split(',')])

    slews = []  # list to store slews in a linear format
    for slew_line in cell_data['slews']:
        slew_values = slew_line.replace('\n', ',').split(',')
        slews.extend([float(slew.strip()) for slew in slew_values if slew.strip()])

    slew_L = max([i for i in range(len(input_slews)) if input_slews[i] <= input_slew] + [
        0])  # Finding left most index for input slew
    slew_R = min([i for i in range(len(input_slews)) if input_slews[i] >= input_slew] + [
        len(input_slews) - 1])  # Finding right most index for input slew
    cap_L = max(
        [j for j in range(len(load_caps)) if load_caps[j] <= load_cap] + [0])  # Finding left most index for load cap
    cap_R = min([j for j in range(len(load_caps)) if load_caps[j] >= load_cap] + [
        len(load_caps) - 1])  # Finding right most index for load cap

    # Handle edge cases by adjusting indices
    if input_slew < input_slews[0]: slew_L, slew_R = 0, 1  # if input slew is less than minimum value in lUT
    if input_slew > input_slews[-1]: slew_L, slew_R = len(input_slews) - 2, len(
        input_slews) - 1  # if input slew is more than maximum value in lUT
    if load_cap < load_caps[0]: cap_L, cap_R = 0, 1  # if load cap is less than minimum value in lUT
    if load_cap > load_caps[-1]: cap_L, cap_R = len(load_caps) - 2, len(
        load_caps) - 1  # if load cap is more than maximum value in lUT

    v11 = slews[(slew_L * len(load_caps)) + cap_L]  # Get delay values wrt slew_L and cap_L
    v12 = slews[(slew_L * len(load_caps)) + cap_R]  # Get delay values wrt slew_L and cap_R
    v21 = slews[(slew_R * len(load_caps)) + cap_L]  # Get delay values wrt slew_R and cap_L
    v22 = slews[(slew_R * len(load_caps)) + cap_R]  # Get delay values wrt slew_R and cap_R

    output_slew = bilinear_interpolation(input_slew, load_cap, input_slews[slew_L], input_slews[slew_R],
                                         load_caps[cap_L],
                                         load_caps[cap_R], v11, v12, v21, v22)  # Perform Bilinear Interpolation
    return output_slew


def infer_cell_name(gate):
    """
    TO match Bench file cell name with liberty file cell name
    """
    cell_type = gate.split('-')[0]
    cell_name_map = {
        'NAND': 'NAND2_X1',
        'NOR': 'NOR2_X1',
        'AND': 'AND2_X1',
        'OR': 'OR2_X1',
        'XOR': 'XOR2_X1',
        'INV': 'INV_X1',
        'NOT': 'INV_X1',
        'BUFF': 'BUF_X1',
    }

    return cell_name_map.get(cell_type, None)


def forward_traversal(gates_fanins, gate_capacitances, slew_LUT, delay_LUT, input_connected_gates,
                      direct_connections):
    all_slews = {gate: 0.002 for gate in
                 input_connected_gates}  # Dictionary to store calculated slews INPUT nodes are set to 0.002 ns
    all_delays = {gate: 0 for gate in
                  input_connected_gates}  # Dictionary to store calculated delays of each fanin to a gate INPUT nodes are set to 0
    all_max_delays = {gate: 0 for gate in
                      input_connected_gates}  # Dictionary to store calculated delays INPUT nodes are set to 0
    all_arrival_times = {gate: 0 for gate in
                         input_connected_gates}  # Dictionary to store calculated arrival times INPUT nodes are set to 0

    def find_cell_data(lut, cell_name):
        """ Returns data for the given cell name with liberty file name"""
        for cell_data in lut:
            if cell_data['cell'] == cell_name:
                return cell_data
        return None

    gates_processed = set(input_connected_gates.keys())  # keeps track of gates processes
    gates_to_process = set(gates_fanins.keys()) - gates_processed  # keeps track of gates to be processed

    while gates_to_process:
        for gate in gates_to_process.copy():
            fanins = gates_fanins[gate][0].split(', ')
            if all(fanin in gates_processed for fanin in fanins):  # Check if all fanins have been processed
                cell_name = infer_cell_name(gate)
                cell_data_delay = find_cell_data(delay_LUT, cell_name)
                cell_data_slew = find_cell_data(slew_LUT, cell_name)
                load_cap = gate_capacitances[gate]

                adjustment_factor = len(fanins) / 2 if len(fanins) > 2 else 1  # adjustment factor condtion

                gate_delays = {fanin: get_delay(cell_data_delay, all_slews[fanin], load_cap) * adjustment_factor for
                               fanin in fanins}  # calculated gate delays
                critical_fanin = max(fanins, key=lambda fanin: all_arrival_times[fanin] + gate_delays[
                    fanin])  # Finds critial fanin: which has maximum arrival + delay times

                max_delay = gate_delays[critical_fanin]

                all_delays[gate] = gate_delays
                all_max_delays[gate] = max_delay
                all_arrival_times[gate] = all_arrival_times[critical_fanin] + gate_delays[critical_fanin]

                gate_slew = get_slew(cell_data_slew, all_slews[critical_fanin], load_cap) * adjustment_factor
                all_slews[gate] = gate_slew

                gates_processed.add(gate)
                gates_to_process.remove(gate)
    for output, connections in output_connected_gates.items():
        max_arrival_time = max(all_arrival_times[gate] for gate in connections)
        all_arrival_times[output] = max_arrival_time
    for direct in direct_connections:  # Used to set arrival times of IN-OUT nodes
        for direct_arrival in direct:
            all_arrival_times[direct_arrival] = 0
            all_delays[direct_arrival] = 0

    return all_slews, all_delays, all_max_delays, all_arrival_times


def find_circuit_delay(arrival_time):
    """ Circuit maximum delay of a circuit"""
    output_arrival_times = [time for gate, time in arrival_time.items() if
                            gate.startswith('OUTPUT-')]  # Arrival times of all OUTPUT nodes is taken
    circuit_delay = max(output_arrival_times)  # Max arrival times of all OUTPUT nodes is taken

    return circuit_delay


def topological_sort(fanout_dict, input_connected_gates):
    """Performs DFS to sort gates from OUTPUT to INPUT"""

    extended_fanout_dict = fanout_dict.copy()
    for input_gate, outputs in input_connected_gates.items():
        extended_fanout_dict[input_gate] = outputs

    fanin_dict = {}
    for gate, outputs in extended_fanout_dict.items():
        for output in outputs:
            if output not in fanin_dict:
                fanin_dict[output] = []
            fanin_dict[output].append(gate)

    no_incoming = {gate for gate in extended_fanout_dict.keys() if gate not in fanin_dict}
    sorted_nodes = []

    while no_incoming:
        node = no_incoming.pop()
        sorted_nodes.append(node)
        if node in extended_fanout_dict:
            for successor in extended_fanout_dict[node]:
                fanin_dict[successor].remove(node)
                if not fanin_dict[successor]:
                    no_incoming.add(successor)
                    del fanin_dict[successor]

    return sorted_nodes


def backward_traversal(fanout_dict, gate_delays_list, circuit_delay, network, direct_connections):
    """ Performs backward traversal to find Required Arrival Time"""
    connections = {}
    for gate, inputs in network.items():
        inputs_list = inputs[0].split(', ')
        for input_node in inputs_list:
            if input_node not in connections:
                connections[input_node] = []
            connections[input_node].append(gate)
    for key, value in fanout_dict.items():
        if not isinstance(value, list):
            fanout_dict[key] = [value]

    extended_network = fanout_dict.copy()

    sorted_nodes = topological_sort(extended_network, input_connected_gates)

    rat_values = {node: None for node in sorted_nodes}  # Initialize RAT values for all nodes

    for node in sorted_nodes:  # Set RAT for OUTPUT nodes based on the circuit delay
        if node.startswith('OUTPUT'):
            rat_values[node] = circuit_delay

    for direct in direct_connections:  # Set RAT for IN-OUT nodes based on the circuit delay
        for direct_arrival in direct:
            rat_values[direct_arrival] = circuit_delay

    # Calculate RAT values in reverse topological order
    for node in reversed(sorted_nodes):
        if node in fanout_dict:
            min_rat_among_fanouts = float('inf')

            for fanout in fanout_dict[node]:
                if fanout in rat_values and rat_values[fanout] is not None:
                    delay_to_fanout = gate_delays_list.get(fanout, {}).get(node,
                                                                           0)  # Default delay is 0 if not specified
                    fanout_rat = rat_values[fanout] - delay_to_fanout
                    if fanout_rat < min_rat_among_fanouts:
                        min_rat_among_fanouts = fanout_rat
            rat_values[node] = min_rat_among_fanouts

        elif node in connections:
            min_rat_among_gates = float('inf')  # Initialize to infinity for comparison
            for gate in connections[node]:  # Iterate over gates connected to the current INPUT node
                if gate in rat_values:  # Ensure the gate has a calculated RAT value
                    delay_to_gate = gate_delays_list.get(gate, {}).get(node,
                                                                       0)  # Get delay from INPUT node to gate, defaulting to 0 if not specified
                    gate_rat = rat_values[gate] - delay_to_gate  # Calculate adjusted RAT for the gate
                    min_rat_among_gates = min(min_rat_among_gates,
                                              gate_rat)  # Update the minimum RAT among gates connected to this INPUT node
            rat_values[node] = min_rat_among_gates if min_rat_among_gates != float(
                'inf') else None  # Set the RAT for the INPUT node, ensuring it's not set to infinity

    return rat_values


def critical_path(fanin_dict, output_connected_gates, slack):
    """To find critical path of the circuit by considering min slews of each node"""
    min_slack_output = min(output_connected_gates.keys(), key=lambda x: slack[x])  # Finds min slack itterating through all connected nodes
    critical_path = []  # List to store nodes in Critical path

    # Function to recursively trace back from a gate to its fan-in gates or inputs
    def trace_back(gate_or_input):
        """ Traces back the path from OUTPUT to INPUT nodes"""
        if gate_or_input.startswith('INPUT'):  # If INPUT, return it
            return [gate_or_input]
        predecessors = fanin_dict[gate_or_input]
        min_slack_predecessor = min(predecessors, key=lambda x: slack[x] if x in slack else float('inf'))
        return trace_back(min_slack_predecessor) + [gate_or_input]

    if min_slack_output in output_connected_gates:
        connected_gate = output_connected_gates[min_slack_output][0]
        critical_path = trace_back(connected_gate)
        critical_path.append(min_slack_output)

    return critical_path


def sorting_output(node):
    """Sorts Nodes"""
    if node.startswith('INPUT'):
        return 0, node
    elif node.startswith('OUTPUT'):
        return 1, node
    else:
        return 2, node


if __name__ == "__main__":
    """Set parser to read input command"""
    parser = argparse.ArgumentParser(description="Parse and process circuit and NLDM files.")
    parser.add_argument('--read_ckt', help='Path to the circuit (.bench) file')
    parser.add_argument('--read_nldm', help='Path to the NLDM (.lib) file')
    args = parser.parse_args()

    inputs, outputs, gates, output_ids_counter, input_ids_counter, primary_inputs, direct_connections = read_ckt(
        args.read_ckt)  #  calls read_ckt fn
    # print("prited gates")
    output_path_ckt = f"../outputs/ckt_details.txt"  # Path to write circuit details
    write_ckt(inputs, outputs, gates, input_ids_counter, output_ids_counter, output_path_ckt)  # write circuit details
    primary_outputs = set(outputs.keys())
    NLDM_content = open_NLDM(args.read_nldm)  # Open liberty file
    NLDM_details_delays = read_delays(NLDM_content)  # Read delays from liberty file
    output_path_delays = f"../outputs/delay_LUT.txt"   #  Path to write delay LUT
    write_delays(NLDM_details_delays, output_path_delays)  # write delay LUT
    NLDM_details_slews, cell_capacitance = read_slews(NLDM_content)  # Read slews from liberty file
    output_path_slews = f"../outputs/slew_LUT.txt"  #  Path to write slew LUT
    write_slews(NLDM_details_slews, output_path_slews)  # write slew LUT
    fanout_dict, fanout_list, fain_dict, input_connected_gates, output_connected_gates = gate_dict(output_ids_counter)  # Finds basic circuit connection
    load_capacitance = load_cap(cell_capacitance, fanout_dict, fanout_list, input_connected_gates)  # find load capacitances for each node

    sorted_gates = topological_sort(fain_dict, input_connected_gates)  # Sort nodes in reverse topological order

    all_slews, all_delays, all_max_delays, all_arrival_times = forward_traversal(fain_dict, load_capacitance,
                                                                                 NLDM_details_slews,
                                                                                 NLDM_details_delays,
                                                                                 input_connected_gates,
                                                                                 direct_connections)  # Performs forward traversal
    circuit_delay = find_circuit_delay(all_arrival_times)
    
    fanout_dict = {gate: fanout[0].split(', ') for gate, fanout in fanout_dict.items() if
                   fanout}  # Splits Fanout dictionary in required format

    rat = backward_traversal(fanout_dict, all_delays, circuit_delay * 1.1, fain_dict, direct_connections)  # Performs forward traversal
    
    slack = {node: rat[node] - all_arrival_times[node] for node in rat if rat[node] is not None}  # Finds slack for each node

    sorted_slack = dict(sorted(slack.items(), key=lambda item: sorting_output(item[0])))  #  Sorts slack for writing Format
    fanin_dict = {gate: fanins[0].split(', ') for gate, fanins in fain_dict.items()}  # Splits Fanin dictionary in required format

    critical_path = critical_path(fanin_dict, output_connected_gates, slack)  # Finds critical path of the circuit

    output_path_traversal = f"../outputs/ckt_traversal.txt"  # Path to write Circuit Traversal Data
    with open(output_path_traversal, "w") as file:
        # Write circuit delay
        file.write(f"Circuit delay: {circuit_delay * 1e3:.4f} ps\n\n")

        # Write the gate slacks
        file.write("Gate slacks:\n")
        for gate, value in sorted_slack.items():
            file.write(f"{gate}: {value * 1e3:.4f} ps\n")

        # Write the critical path
        file.write("\nCritical path:\n")
        file.write(", ".join(critical_path))

toc = time.perf_counter()
print(toc - tic)
