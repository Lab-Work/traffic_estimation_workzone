__author__ = 'Yanning Li'
"""
This script reads the demand data for I57 and convert it to OD matrix format. The reason for the conversion is that
AIMSUN only allows exporting the trajectory data if the traffic demand is specified by OD matrices instead of traffic states.

The input files are:
- para.txt: state_name, vehicle_type, start_time, duration
- flow.txt: the inflow to main entrance and onramps

The output is a file containing an OD matrix.

Notes:
- OD matrix is more informative than traffic states. To reconstruct the OD matrix from the less informative traffic
  states, we assume: all offramp flows origin from the main entrance; all onramp flows head to the main exit.
- Theoretically, we can use the turn.txt and flow.txt to compute the OD pairs. However, an easier and more accurate way
  is to simply simulate the traffic in AIMSUN using the traffic states demand. Measure the offramp flows and specify
  the offramp flow manually.
- Hence the procedure is:
    - Open AIMSUN and add all centroids to entrances and exits. Create a dictionary in setup_centroid():
        centroid[sec_id] = centroid_id
    - Simulate in AIMSUN and get the offramp flow. Use function generate_off_flow to manually set the offramp flow
    - Read flow.txt and para.txt and save in dict
    - write the OD matrix.
"""

# 1. The vehicle trajectory data can only be exported when the traffic demand is specified by OD matrices.
# 2. We have tested that we can simply compute convert the turn ratio to OD matrices.
# 3. The off ramp flow all comes from the flow entered from the mainstream, the flow is roughly set by checking the
#    simulated data in AIMSUN
# 4. On ramp flows all travel to the exit of the freeway.


import numpy as np
import sys


# set up the centroid
# centroid[section_id] = centroid_id
centroid  = dict()


main_entrance = 51413
main_exit = 473
onramp = 519
duration = 5    # each state is 5 min

# keep an ordered matrix names
matrix_names = []

# para[matrix_name] = ['car', 'start_time', 'duration']
para = dict()
# inflow[matrix_name] = [sec_id, flow]
inflow = dict()


# set up centroid dictionary
# centroid[sec_id] = centroid_id
def setup_centroid():

    # main entrance
    centroid[51413] = 51512

    # main exit
    centroid[473] = 51515

    # on ramp
    centroid[519] = 51517


# Read the parameter file to get the matrix names and the car type, start time, and duration of each state.
def read_para(para_file):

    f = open(para_file)

    for line in f:

        line = line.strip()

        items = line.split(',')

        # keep an ordered name list
        matrix_names.append(items[0])

        para[items[0]] = [items[1], items[2], items[3]]

    f.close()


# Read the flow file and save in inflow[matrix_name] = the amount of inflow
def read_flow(flow_file):

    f = open(flow_file)

    for line in f:

        line = line.strip()

        items = line.split(',')
        name = items[0]
        sec_id = int(items[1])
        flow = float(items[2])

        if name not in inflow.keys():
            inflow[name] = []

        inflow[name].append([sec_id, flow])

    f.close()


# write the OD matrix
# - if the flow is from freeway entrance, specify the destination to off ramp centroid.
# - if the flow is from on ramps, specify the destination as the freeway off ramp.
def write_OD_matrix():

    f = open('I57_matrices.txt', 'wb')

    # random set matrix id
    tmp_id = 1
    for name in matrix_names:

        # the first line of a OD matrix block
        # od_matrix_id od_matrix_name
        f.write('{0} {1}\n'.format( tmp_id ,name))
        tmp_id += 1

        # second line:
        # car_type_id car_type_string
        f.write('{0}\n'.format(para[name][0]))

        # third line:
        # start_time
        f.write('{0}\n'.format(para[name][1]))

        # forth line:
        # write the duration
        f.write('{0}\n'.format(para[name][2]))

        # write the OD values
        # from_centroid_id to_centroid_id number_of_vehicles
        for flow in inflow[name]:

            # all flow goes to the main exit.
            # Note, OD matrix specifies the number of cars sent during the interval, not the flow (veh/h)
            f.write('{0} {1} {2}\n'.format(centroid[flow[0]], centroid[main_exit], flow[1]*duration/60.0))

        # separate states by blank line
        f.write('\n')

    f.close()





def main(argv):

    # first set up the centroid ids
    setup_centroid()

    # read the in flow and state information
    read_para('paras.txt')
    read_flow('flows.txt')

    # reorganize the data and write in the standard OD matrix format.
    write_OD_matrix()





if __name__ == "__main__":
    sys.exit(main(sys.argv))


