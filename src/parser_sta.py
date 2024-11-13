# STA Parser Multi-Block

# Bash Terminal Commands Supported:

# [1] python3.7 parser.py --read_ckt c17.bench <-- Change as per bench file requirement
# Function: Provides a .txt file output for Primary Input, Outputs, Gates, Fan-In and Fanout fetched from the Bench File.

# [2] python3.7 parser.py --delays --read_nldm sample_NLDM.lib <-- Change as per Liberty file requirement
# Function: Provides .txt file output for the "delays" fetched from the Liberty File.

# [3] python3.7 parser.py --slews --read_nldm sample_NLDM.lib <-- Change as per Liberty file requirement
# Function: Provides .txt file output for the "slews" fetched from the Liberty File.

# Libraries Used:
import argparse  # parsing command-line arguments
from collections import Counter  # Counting occurrences


# [1] To read and process Bench file data


def read_ckt(input_path):
    inputs = {}  # Creating Dictionary for Inputs
    outputs = {}  # Creating Dictionary for Outputs
    gates = {}  # Creating Dictionary for Gates
    output_ids_counter = Counter()  # To count OUTPUT IDs
    input_ids_counter = Counter()  # To count INPUT IDs

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

            # If the line starts with the keyword OUTPUT
            elif line.startswith("OUTPUT"):
                output_id = line.split("(")[1].rstrip(")")
                outputs[output_id] = {"input": None}
                output_ids_counter[output_id] += 1
            # Circuit Connections
            else:
                parts = line.split(" = ")
                gate_output, operation_with_inputs = parts[0], parts[1]
                operation, gate_inputs_str = operation_with_inputs.split("(")
                gate_inputs = gate_inputs_str.rstrip(")").split(", ")
                gates[gate_output] = {"type": operation, "inputs": gate_inputs, "outputs": []}
    #  Updates each gate's outputs according to its inputs
    for gate_id, gate_info in gates.items():
        for input_id in gate_info["inputs"]:
            if input_id in gates:
                gates[input_id]["outputs"].append(gate_id)
            elif input_id in inputs:
                inputs[input_id]["outputs"].append(gate_id)
    # Connects primary outputs to the gate(s) generating them
    for output_id, output_info in outputs.items():
        for gate_id, gate_info in gates.items():
            if output_id in gate_info["outputs"]:
                output_info["input"] = gate_id

    return inputs, outputs, gates, output_ids_counter, input_ids_counter  # Return to the extracted values


# Function to write circuit data to a specific output file with counts of primary inputs & outputs
def write_ckt(inputs, outputs, gates, input_ids_counter, output_ids_counter, output_file_path):
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


# Reading the NLDM Liberty File
def open_NLDM(input_path):
    with open(input_path, 'r') as file:
        return file.read()


# [2] Read Delay Values from the liberty file
def read_delays(NLDM_file):
    lines = NLDM_file.split('\n')
    NLDM_details = []   # emtpy list to store input slew, load capacitance and delay
    cell = None     # checker to check if it is inside cell
    if_delay = False    # checker to check if it is inside delay
    multi_line = False  # For iterating through multiple line in delay
    delays = ""     # empty string to extract delay values

    for line in lines:
        line = line.strip()
        if line.startswith("cell ("):
            if cell:
                NLDM_details.append(cell)
            cell_name = line.split("(")[1].split(")")[0].strip()
            cell = {"cell": cell_name, "input_slews": "", "load_cap": "", "delay": []}
        elif "cell_delay" in line and cell is not None:
            if_delay = True
        elif line.startswith("index_1") and if_delay:       # condition to get input slew
            cell["input_slews"] = line.split('"')[1]
        elif line.startswith("index_2") and if_delay:       # condition to get load capacitance
            cell["load_cap"] = line.split('"')[1]
        elif if_delay and line.startswith("values"):        # condition to get delays
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
    return NLDM_details     # return NLDM details


# Write Delay Values to the delay_LUT.txt file
def write_delays(NLDM_details, output_path):
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


# [3] Read Slew Values from the liberty file
def read_slews(NLDM_file):
    lines = NLDM_file.split('\n')
    NLDM_details = []   # emtpy list to store input slew, load capacitance and slew
    cell = None     # checker to check if it is inside cell
    if_slew = False     # checker to check if it is inside slew
    multi_line = False  # For iterating through multiple line in delay
    output_slews = ""     # empty string to extract slew values

    for line in lines:
        line = line.strip()
        if line.startswith("cell ("):
            if cell:
                NLDM_details.append(cell)
            cell_name = line.split("(")[1].split(")")[0].strip()
            cell = {"cell": cell_name, "input_slews": "", "load_cap": "", "slews": []}
        elif "output_slew" in line and cell is not None:
            if_slew = True
        elif line.startswith("index_1") and if_slew:        # condition to get input slews
            cell["input_slews"] = line.split('"')[1]
        elif line.startswith("index_2") and if_slew:        # condition to get load capacitance
            cell["load_cap"] = line.split('"')[1]
        elif if_slew and line.startswith("values"):     # condition to get output slews
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
    return NLDM_details


# Write Slew Values to the slew_LUT.txt file


def write_slews(NLDM_details, output_path):
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

# Main function to parse the commands
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    # Command Declarations to run in terminal window
    parser.add_argument('--read_ckt')
    parser.add_argument('--read_nldm')
    parser.add_argument("--delays", action="store_true")
    parser.add_argument('--slews', action="store_true")

    args = parser.parse_args()
    # parser for reading and writing Bench file
    if args.read_ckt:
        inputs, outputs, gates, output_ids_counter, input_ids_counter = read_ckt(args.read_ckt)
        output_path = f"../outputs/ckt_details.txt"  # Path to ckt_details.txt file
        write_ckt(inputs, outputs, gates, input_ids_counter, output_ids_counter, output_path)
    # parser for reading and writing NLDM file for delays
    if args.read_nldm and args.delays:
        NLDM_file = open_NLDM(args.read_nldm)
        NLDM_details = read_delays(NLDM_file)
        output_path = f"../outputs/delay_LUT.txt"  # Path to delay_LUT.txt file
        write_delays(NLDM_details, output_path)
    # parser for reading and writing NLDM file for slews
    if args.read_nldm and args.slews:
        NLDM_file = open_NLDM(args.read_nldm)
        NLDM_details = read_slews(NLDM_file)
        output_path = f"../outputs/slew_LUT.txt"  # Path to slew_LUT.txt file
        write_slews(NLDM_details, output_path)
