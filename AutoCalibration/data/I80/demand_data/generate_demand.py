__author__ = 'Yanning'

# This is a utility function. It converts all kinds of data to other formats
# For now, the main purpose is that we need to calibrate the road network in segments.
# And this file provides the functions that can easily generate the inflow we need for the main
# entrance and the onramps, as well as the percentage of the offramps.


import sys
import csv
from collections import OrderedDict


# the main function
def main(argv):

    # convert_pars('input_paras.txt')

    # convert_EB_to_flow('EB3.csv', 21216, 'cong_state')

    # you may put multiple onramp files
    main_flow_file = 'input_mainflow.txt'
    truck_ratio = 0.27
    generate_flows_with_new_format(main_flow_file,
                   ['input_onramp1_1.txt',
                    'input_onramp1_2.txt',
                    'input_onramp2.txt',
                    'input_onramp3_1.txt',
                    'input_onramp3_2.txt',
                    'input_onramp4.txt'], truck_ratio)

    # may use multiple offramp files
    generate_turns(['input_offramp1.txt',
                    'input_offramp2.txt',
                    'input_offramp3_1.txt',
                    'input_offramp3_2.txt',
                    'input_offramp4.txt'])


# This function read the following files:
# main_flow.txt, onramp_1.txt, onramp_2.txt, and offramp.txt files
# Output a stitched file: flows.txt
def generate_flows(main_flow_file, onramp_file_list, main_flow_truck_ratio):

    # flow[new_statename] = [(sec_id, inflow(veh/hr)),(sec_id, inflow)]
    flow_car = OrderedDict()
    flow_truck = OrderedDict()

    # read the main flow file
    f_mainflow = open(main_flow_file, 'r')

    next(f_mainflow)

    # reduce the inflow after 15 states
    counter = 1
    for line in f_mainflow:

        line = line.strip()
        items = line.split(',')

        # create new states
        car_state_name = items[0] + '_car'
        truck_state_name = items[0] + '_truck'

        flow_car[car_state_name] = []
        flow_truck[truck_state_name] = []

        # if counter < 16:
        #     car_inflow = (int(items[1]), float(items[2])*12*(1-main_flow_truck_ratio))
        #     truck_inflow = (int(items[1]), float(items[2])*12*main_flow_truck_ratio )
        # else:
        # add 0.7 discount to clear the congestion
        car_inflow = (int(items[1]), float(items[2])*12*(1-main_flow_truck_ratio)*0.0 )
        truck_inflow = (int(items[1]), float(items[2])*12*main_flow_truck_ratio*0.0 )


        flow_car[car_state_name].append(car_inflow)
        flow_truck[truck_state_name].append(truck_inflow)

        counter += 1

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

        # tune the inflow
        # the last item is the percentage
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
# In each offramp file, first row is from_sec, to_sec, off_sec. Then each row is statename, off percentage.
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

