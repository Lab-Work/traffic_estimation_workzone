This is the equivalent OD matrices demand.

input files:
flows.txt, paras.txt, turns.txt

OD demand output file:
I57_matrices.txt

OD demand output file format: see file matrices.txt
Each demand matrix for a time interval is specified in a block. Blocks are separated by blank lines.

For each block,
demand_matrix_id demand_matrix_name
vehicle_type_id vehicle_type_string
start_time (in format 15:00:00)
duration (in format 00:05:00)
from_centroid_id to_centroid_id #_vehicles_in_duration
…
all from_id to all to_centroid_id #_vehicles_in_duration

