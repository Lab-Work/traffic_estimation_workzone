__author__ = 'Yanning Li'
"""
Yanning Li, updated June 19, 2016

This is a utility function, which converts the data into other formats.
- It converts the input_mainflow and input_onramp to the format used for specifying the traffic demand in AIMSUN.


"""


import sys
import csv
from collections import OrderedDict


# the main function
def main(argv):

    convert_pars('input_paras.txt')

    # convert_EB_to_flow('SB1.csv', 51413, 'cong_state')

    # you may put multiple onramp files
    main_flow_file = 'input_mainflow.txt'

    # truck ratio is for main flow
    truck_ratio = 0.32
    # truck_ratio = 0.0
    generate_flows_with_new_format(main_flow_file,
                   ['input_onramp.txt'], truck_ratio)

    # may use multiple offramp files
    # generate_turns(['input_offramp1.txt',
    #                 'input_offramp2.txt',
    #                'input_offramp3_1.txt',
    #               'input_offramp3_2.txt',
    #                 'input_offramp4.txt'])





# This function read the following files:
# main_flow.txt, onramp_1.txt, onramp_2.txt, and offramp.txt files
# Output a stitched file: flows.txt
def generate_flows(main_flow_file, onramp_file_list, main_flow_truck_ratio):

    # flow[new_statename] = [(sec_id, inflow(veh/hr)),(sec_id, inflow)]
    flow_car = {}
    flow_truck = {}

    # read the main flow file
    f_mainflow = open(main_flow_file, 'r')

    next(f_mainflow)

    for line in f_mainflow:

        line = line.strip()
        items = line.split(',')

        # create new states
        car_state_name = items[0] + '_car'
        truck_state_name = items[0] + '_truck'

        flow_car[car_state_name] = []
        flow_truck[truck_state_name] = []

        car_inflow = (int(items[1]), float(items[2])*12*(1-main_flow_truck_ratio))
        truck_inflow = (int(items[1]), float(items[2])*12*main_flow_truck_ratio )

        flow_car[car_state_name].append(car_inflow)
        flow_truck[truck_state_name].append(truck_inflow)

    f_mainflow.close()

    # now read all the ramps
    # we split because for the main onramp flow, we assume a truck ratio,
    # while for the other onramps, we assume only cars
    for ramp_file in onramp_file_list:

        f_onramp = open(ramp_file, 'r')

        next(f_onramp)

        for line in f_onramp:

            line = line.strip()
            items = line.split(',')

            car_state_name = items[0] + '_car'
            truck_state_name = items[0] + '_truck'

            # states created before
            flow_car[car_state_name].append( (int(items[1]), float(items[2])*12) )
            # just assume zero truck comes on
            flow_truck[truck_state_name].append( ( int(items[1]), 0 ))

        f_onramp.close()

    # now save the dict in to a file
    f_flow_all = open('flows.txt', 'wb')
    # write cars
    for state in flow_car.keys():

        for tup in flow_car[state]:
            # statename, section id, flow
            f_flow_all.write( state + ',' + str(tup[0]) + ',' + str(tup[1]) + '\n')

    # write trucks
    for state in flow_truck.keys():

        for tup in flow_truck[state]:
            # statename, section id, flow
            f_flow_all.write( state + ',' + str(tup[0]) + ',' + str(tup[1]) + '\n')


# This function read the following files:
# main_flow.txt, onramp_1.txt, onramp_2.txt, and offramp.txt files
# Two changes from previous function:
# 1. In onramps files, state name are complete: congflow_state1_car/truck
# 2. In onramp files, flow are veh/hr instead of count/5min
# Output a stitched file: flows.txt
def generate_flows_with_new_format(main_flow_file, onramp_file_list, main_flow_truck_ratio):

    # flow[new_statename] = [(sec_id, inflow(veh/hr)),(sec_id, inflow)]
    flow_car = OrderedDict()
    flow_truck = OrderedDict()

    # read the main flow file
    f_mainflow = open(main_flow_file, 'r')

    next(f_mainflow)

    for line in f_mainflow:

        line = line.strip()
        items = line.split(',')

        # create new states
        car_state_name = items[0] + '_car'
        truck_state_name = items[0] + '_truck'

        # key defined here
        flow_car[car_state_name] = []
        flow_truck[truck_state_name] = []

        car_inflow = (int(items[1]), float(items[2])*12*(1-main_flow_truck_ratio)*float(items[3]))
        truck_inflow = (int(items[1]), float(items[2])*12*main_flow_truck_ratio*float(items[3]) )

        flow_car[car_state_name].append(car_inflow)
        flow_truck[truck_state_name].append(truck_inflow)

    f_mainflow.close()

    # now read all the ramps
    # we split because for the main onramp flow, we assume a truck ratio,
    # while for the other onramps, we assume only cars
    for ramp_file in onramp_file_list:

        f_onramp = open(ramp_file, 'r')

        next(f_onramp)

        for line in f_onramp:

            line = line.strip()
            items = line.split(',')

            # find out which state is this
            name = items[0].split('_')[2]

            if name == 'car':
                car_state_name = items[0]
                # states created before
                flow_car[car_state_name].append( (int(items[1]), float(items[2])) )
            elif name == 'truck':
                truck_state_name = items[0]
                flow_truck[truck_state_name].append( ( int(items[1]), float(items[2]) ))

        f_onramp.close()

    # now save the dict in to a file
    f_flow_all = open('flows.txt', 'wb')
    # write cars
    for state in flow_car.keys():

        for tup in flow_car[state]:
            # statename, section id, flow
            f_flow_all.write( state + ',' + str(tup[0]) + ',' + str(tup[1]) + '\n')

    # write trucks
    for state in flow_truck.keys():

        for tup in flow_truck[state]:
            # statename, section id, flow
            f_flow_all.write( state + ',' + str(tup[0]) + ',' + str(tup[1]) + '\n')




# This function generates turns
# You may put multiple turns file in, and then stitch them together.
def generate_turns(offramp_file_list):

    f_towrite = open('turns.txt', 'wb')
    writer = csv.writer(f_towrite, delimiter=',')

    for offramp_file in offramp_file_list:

        f_raw = open(offramp_file, 'r')
        next(f_raw)

        for line in f_raw:

            line = line.strip()
            items = line.split(',')

            if len(items) == 3:
                # the section IDs
                from_sec_id = items[0]
                to_sec_id = items[1]
                off_sec_id = items[2]

            else:
                # cong_state#_car/truck, offratio
                state_name = items[0]
                offratio = float(items[1])

                # save into turns.txt file
                writer.writerow([state_name, from_sec_id, to_sec_id, str(100-offratio)])
                writer.writerow([state_name, from_sec_id, off_sec_id, str(offratio)])

        f_raw.close()

    f_towrite.close()





# convert the EB data into flows
def convert_EB_to_flow(EB_file, section_id, state_name):

    f_raw = open(EB_file, 'r')
    next(f_raw)

    f_towrite = open('input_mainflow.txt', 'wb')
    writer = csv.writer(f_towrite, delimiter=',')

    writer.writerow(['state', 'section_id', 'count/5min'])

    i = 1
    for line in f_raw:
        line = line.strip()
        items = line.split(',')

        writer.writerow([ state_name+str(int(i)),
                          str(section_id),
                          items[2]])
        i += 1

    f_raw.close()
    f_towrite.close()


# This function reads and just converts the paras
def convert_pars(pars_file_path):

     # Convert the parameter file
    f_par_in = open(pars_file_path, 'r')
    f_par_out = open('paras.txt', 'wb')

    orig_state_name = []

    while True:
        line = f_par_in.readline()
        if not (bool(line)):
            break

        line.strip()
        items = line.split(',')

        orig_state_name.append(items[0])
        items[3] = items[3].strip()

        f_par_out.write(items[0] + '_car' + ',' +
                        items[1] + ',' +
                        items[2] + ',' +
                        items[3] + '\n')
        f_par_out.write(items[0] + '_truck' + ',' +
                        '56 Truck' + ',' +
                        items[2] + ',' +
                        items[3] + '\n')

    f_par_in.close()
    f_par_out.close()







if __name__ == "__main__":
    sys.exit(main(sys.argv))

