**<h1>STA FOR COMBINATIONAL CIRCUITS**

**TEAM: IDEAL IDLE</h1>**


A Python-based tool called STA Combinational is intended to extract
and process important data from Non-Linear Delay Model (NLDM) liberty
files (**.lib**) and digital circuit benchmark files (**.bench**). It
makes the process of Static Timing Analysis (STA) for any combinational 
logic easier by running commands for:

1. **Extraction:** of primary inputs, outputs, gates, fan-ins, fan-outs, delays, and slews.

2. a) **Extraction:** of primary inputs, outputs, gates, fan-ins, fan-outs, delays, and slews.

   b) **Computation:** of Circuit Delay, Gate Slacks and Critical Path.

**<h3>FEATURES:</h3>**

-   **Circuit Analysis:** Extracts circuit details like number of
    primary inputs, number of primary outputs and related circuit gate
    details from '**.bench**' file.

-   **Delay Extraction:** Extracts input_slews, load_capacitance and
    delays from a '**.lib**'

-   **Slew Extraction:** Extracts input_slews, load_capacitance and
    output_slews from a '**.lib**'
	
-   **Circuit Delay Computation:** Computes Circuit Delay (in ps) of a netlist 
	'**.bench**' file with the cell parameters fetched from a'**.lib**' file.
	
-   **Gate Slacks Computation:** Computes slack for every gate(in ps) of a netlist 
	'**.bench**' file with the cell parameters fetched from a'**.lib**' file.
	
-   **Critical Path Computation** Computes Critical Path for a netlist 
	'**.bench**' file with the cell parameters fetched from a'**.lib**' file.
	
**<h3>DEPENDENCIES:</h3>**

-   Python 3.7 or higher

-   Libraries:

    -   '**argparse**' for parsing command-line arguments.

    -   '**collections**' for utilizing the Counter class.
	
**<h3>MODES SUPPORTED FOR SCRIPT(S):</h3>**

1. Extraction: parser.py

2. Computation/Extraction: main_sta.py 

- Running (2) will automatically generate outputs from (1).

**<h3>USAGE:</h3>**

- Bash Terminal Commands Supported for a given mode:

**(1) Instructions to run parser_sta.py**

**(a)To parse and analyze a .bench file:**

>> python3.7 parser_sta.py --read_ckt (path to .bench file)

[Change path as per Bench File (.bench) selection]

Function: Provides a .txt file output for Primary Input, Outputs, Gates,
Fan-In and Fanout fetched from the Bench File.

**(b)To extract delays from a .lib file:**

>> python3.7 parser_sta.py --delays --read_nldm (path to .lib file)

[Change path as per Liberty file (.lib) selection]

Function: Provides .txt file output for the \"delays\" fetched from the
Liberty File for each gate.

**(c)To extract slews from a .lib file:**

>> python3.7 parser_sta.py --slews --read_nldm (path to .lib file)

[Change path as per Liberty file (.lib) selection]

Function: Provides .txt file output for the \"slews\" fetched from the
Liberty File for each gate.

**(2): Instructions to run main_sta.py**

**To compute Circuit Delay, Gate Slacks and Critical Path and extract circuit details, slews and delays from a selected .bench and .lib file:**

>> python3.7 main_sta.py --read_ckt (path to .bench file) --read_nldm (path to .lib file)

[Change path as per Liberty file (.lib) and Bench File (.bench) selection]

Function: Provides .txt file output for the \"circuit delay\", \"gate slacks\" and \"critical path\"
processed from the selected bench file referenced with the parameters fetched from the Liberty File for each gate.

**<h3>EXAMPLE(S):</h3>**

**(a)To analyze a benchmark file named c17.bench and extract circuit details:**

>> python3.7 parser_sta.py --read_ckt c17.bench

**(b)To fetch slews information from a liberty file named sample_NLDM.lib:**

>> python3.7 parser_sta.py --slews --read_nldm sample_NLDM.lib

**(c)To fetch delay information from a liberty file named sample_NLDM.lib:**

>> python3.7 parser_sta.py --delays --read_nldm sample_NLDM.lib

**(d)To compute Circuit Delay, Gate Slacks and Critical Path for a selected c17.bench and sample_NLDM.lib file:**

>> python3.7 main_sta.py --read_ckt c17.bench --read_nldm sample_NLDM.lib

**<h3>OUTPUTS</h3>**

The generated output files after running any of the above Bash Commands are stored under the folder '**outputs**' in the parent directory.
