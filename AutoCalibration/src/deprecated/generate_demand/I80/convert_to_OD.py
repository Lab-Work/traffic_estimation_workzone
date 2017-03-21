__author__ = 'Yanning Li'
"""
This script reads the demand data for I80 and convert it to OD matrix format. The reason for the conversion is that
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


import numpy as np
import sys


# set up the centroid
# centroid[section_id] = centroid_id
centroid  = dict()

# TODO: change the centroid ids once get the AIMSUN to work
main_entrance = 21216
main_exit = 13008
entrance_sec_ids = [30109, 30118, 21183, 21192, 21201, 1357]
exit_sec_ids = [23587, 15804, 343, 1039, 1501]
off_ramp_4 = 1501
duration = 5    # each state is 5 min

# keep an ordered matrix names
matrix_names = []

# para[matrix_name] = ['car', 'start_time', 'duration']
para = dict()
# inflow[matrix_name] = [sec_id, flow]
inflow = dict()

# off ramp flow for cars
# offramp[offramp_sec_id] = [flow] (veh/hr) 30 states list
off_flow = dict()
# trucks exist only at the last ramp; all others are 0
truck_off_flow_ramp_4 = None


# set up centroid dictionary
def setup_centroid():

    # main entrance
    centroid[21216] = 41142

    # onramp 1_1
    centroid[30109] = 41144
    # offramp 1
    centroid[23587] = 41146
    # onramp 1_2
    centroid[30118] = 41148

    # offramp 2
    centroid[15804] = 41150
    # onramp 2
    centroid[21183] = 41152

    # offramp 3_1
    centroid[343] = 41154
    # onramp 3_1
    centroid[21192] = 41156

    # offramp 3_2
    centroid[1039] = 41158
    # onramp 3_2
    centroid[21201] = 41160

    # offramp 4
    centroid[1501] = 41162
    # onramp 4
    centroid[1357] = 41164

    # main exit
    centroid[13008] = 41166


# manually generate the off ramp flows.
# The flow is approximated from the simulation using original traffic states.
def generate_off_flow():

    global truck_off_flow_ramp_4

    # offramp 1
    off_flow[23587] = np.concatenate(( np.ones((1,9),float)*30 ,
                                       np.ones((1,21),float)*50) , None )

    # offramp 2
    off_flow[15804] = np.concatenate(( np.ones((1,9),float)*20 ,
                                       np.ones((1,6),float)*40 ,
                                       np.ones((1,15),float)*50) , None )

    # offramp 3_1
    off_flow[343] = np.concatenate(( np.ones((1,10),float)*20 ,
                                     np.ones((1,5),float)*40 ,
                                     np.ones((1,15),float)*60) , None )

    # offramp 3_2
    off_flow[1039] = np.concatenate(( np.ones((1,15),float)*30 ,
                                      np.ones((1,15),float)*50) , None )

    # offramp 4
    off_flow[1501] = np.concatenate(( np.ones((1,10),float)*20 ,
                                      np.ones((1,20),float)*26) , None )

    truck_off_flow_ramp_4 = np.concatenate(( np.ones((1,10),float)*40 ,
                                             np.ones((1,20),float)*54) , None )


# read the parameter file
def read_para(para_file):

    f = open(para_file)

    for line in f:

        line = line.strip()

        items = line.split(',')

        # keep an ordered name list
        matrix_names.append(items[0])

        para[items[0]] = [items[1], items[2], items[3]]

    f.close()


# read the flow file
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

    f = open('I80_matrices.txt', 'wb')

    # random set matrix id
    tmp_id = 1
    for name in matrix_names:

        tup = name.split('_')
        vehicle_type = tup[2].strip()

        # print(tup)
        # find out the state id (which will also be used as the matrix id) from the name
        matrix_id = int(tup[1][5:])

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

            if flow[0] == main_entrance:
                # print vehicle_type
                # print len(vehicle_type)
                if vehicle_type == 'car':
                    # print 'flow: {0}'.format(flow)
                    flow_left = flow[1]

                    # first split flows to all offramp
                    for ramp_exit in exit_sec_ids:
                        f.write('{0} {1} {2}\n'.format( centroid[flow[0]], centroid[ramp_exit], off_flow[ramp_exit][matrix_id-1]*duration/60.0 ))
                        flow_left -= off_flow[ramp_exit][matrix_id-1]

                    # the flow left will be sent to the main exit
                    f.write('{0} {1} {2}\n'.format( centroid[flow[0]], centroid[main_exit], flow_left*duration/60.0 ))

                elif vehicle_type == 'truck':
                    flow_left = flow[1]

                    # first split flows to all offramp
                    # trucks exit on offramp4, and 0 on all other ramps
                    for ramp_exit in exit_sec_ids:
                        if ramp_exit != off_ramp_4:
                            f.write('{0} {1} {2}\n'.format( centroid[flow[0]], centroid[ramp_exit], 0 ))
                        else:
                            f.write('{0} {1} {2}\n'.format( centroid[flow[0]], centroid[ramp_exit], truck_off_flow_ramp_4[matrix_id-1]*duration/60.0 ))
                            flow_left -= truck_off_flow_ramp_4[matrix_id-1]

                    # flow left is sent to the main exit
                    f.write('{0} {1} {2}\n'.format( centroid[flow[0]], centroid[main_exit], flow_left*duration/60.0 ))

            else:
                # all vehicles that come from the onramp go to the main exit
                f.write('{0} {1} {2}\n'.format(centroid[flow[0]], centroid[main_exit], flow[1]*duration/60.0))

                # specify the od paris from the onramp to offramp all as 0.
                for ramp_exit in exit_sec_ids:
                    f.write('{0} {1} {2}\n'.format(centroid[flow[0]], centroid[ramp_exit], 0))

        # separate states by blank line
        f.write('\n')

    f.close()





def main(argv):

    # first set up the centroid ids
    setup_centroid()

    # generate off ramp flows
    generate_off_flow()

    # read the in flow and state information
    read_para('paras.txt')
    read_flow('flows.txt')

    # reorganize the data and write in the standard OD matrix format.
    write_OD_matrix()





if __name__ == "__main__":
    sys.exit(main(sys.argv))


