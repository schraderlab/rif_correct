from __future__ import print_function# inports python3's print function
if hasattr(__builtins__, 'raw_input'): #checks if python2
    input = raw_input # if python2 sets input function to python3's input function

import sys
import numpy as np
from scipy import stats
import csv



def load_predicted_start_sites(file_name):
    f=open(file_name)
    records=[]                                    #there is no need for creating a list before, the next step automatically does it and so even if this line is missing it will work fine
    records=f.readlines()
    f.close()

    fields=[]
    for line in records:
        field=line.split()
        fields.append(field)
    return fields

# ***** better way to read in CSV files, but no time to debug *****
# def load_predicted_start_sites(file_name):
#     with open(file_name) as csvfile:
#         csvreader = csv.reader(csvfile, delimiter="\t")
#
#         students = []
#         for row in csvreader:
#             students.append(row)
#     return students


#making program usable with 2 or 3 input types (with and without transcriptome file)
try:
    sys.argv[3]
except IndexError:
    if __name__=="__main__":
         start_sites=load_predicted_start_sites(sys.argv[1]) #cds_V5_20190801_2_working_correct.csv
         data=load_predicted_start_sites(sys.argv[2]) #dataset.csv
         arguments = 2
else:
    if __name__=="__main__":
          start_sites=load_predicted_start_sites(sys.argv[1]) #cds_V5_20190801_2_working_correct.csv
          data=load_predicted_start_sites(sys.argv[2]) #dataset.csv
          operon_info=load_predicted_start_sites(sys.argv[3]) #all_TSS_sites_hybrid_nospace.csv
          arguments = 3


internal_counter=0
leading_counter=0

#comment in the below script when PIPE > is figured out
# time_points_set = [0, 1, 2, 4, 8, 15, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30] #place holder for above to try with full 20 data points
# time_points = ['x','x','x','x','x','x','x','x','x','x','x','x','x','x','x','x','x','x','x','x' ]
# for x in range(0,len(time_points_set)): #comment in the below script when PIPE > is figured out
#     time_points[x] = time_points_set[x] #comment in the below script when PIPE > is figured out

filename_output = eval(input(" In Quotation Marks: Enter your filename. Filename must end in .csv (output is tab deliminted CSV file): "))

time_points = ['x','x','x','x','x','x','x','x','x','x','x','x','x','x','x','x','x','x','x','x' ] #list of 20 x's
time_points_set = eval(input("In ascending order enter your time points (in minutes) seperated by commas (no quotes): ")) #in put time_points
if len(time_points_set) < 3: # need a minimum of 3 time points
    print("")
    print("a minimum of three timpoints are needed, program stopped!")
    print("")
    exit()
if len(time_points_set) > 20: # need no more than 20 time points
    print("")
    print("a max of twenty timpoints are needed, program stopped!")
    print("")
    exit()

for x in range(0,len(time_points_set)): # combine inputted time points with time_points blank list
    time_points[x] = time_points_set[x]

elong_rate_input = eval(input("Enter Polymerase mRNA elongation rate in nucleotides per Second (E. coli has been measured at 25, C. Crescentus has been measured at 19.36) entering a value of 0 nt/sec will stop any lag correction: "))

# if elong_rate_input <= 0:
#     elong_rate = 1E-100
# if elong_rate_input > 0:
#     elong_rate = elong_rate_input*60

elong_rate = elong_rate_input*60


flag_r_value = eval(input("Enter at what r value would you like the transformed linear regression to be flagged above (r = -1 is perfect fit, r = 0 is a flat line, r = 1 is perfect INVERSE correlation) r = -0.65 is recommended: "))

cut_off_percent = eval(input("Enter at what percent of the 0 minute RPKM for each gene would you like to set as a cutoff (100% to 0%). For each individual gene RPKM datapoints below the cutoff will be excluded (2.71828 is recommended): "))

if cut_off_percent <= 0:
    cut_off = np.log(1E-100)
if cut_off_percent > 0:
    cut_off = np.log(cut_off_percent)



blank_x = ['x']
blank_simple = ['simple']
blank_operon = ['n']

new_list = []
new_list_temp = []

if arguments == 2: # this if when you do not have transcriptome file hence 2 arguments
    for x in start_sites:
        try:
            invalue = int(x[2])
        except ValueError:
            pass
        else:
            new_list_temp = []
            new_list_temp.append(str(x[0])) #operon
            new_list_temp.append(x[1]) #gene name/
            new_list_temp.append(int(x[2])) # start
            new_list_temp.append(int(x[3])) # / stop
            new_list_temp.append(abs(int(x[2])-int(x[3]))) # gene length
            new_list_temp.append(x[4]) # sense /
            new_list_temp.append(blank_operon[0]) # in operon
            new_list_temp.append(x[1]) # start gene
            new_list_temp.append(x[1]) # stop gene
            if x[4] == '+': #correct for reverse order of start pos and stop pos in transciptome file vs operon file
                new_list_temp.append(int(x[2])) # start position
                new_list_temp.append(int(x[3])) # stop position
            elif x[4] == '-': #correct for reverse order of start pos and stop pos in transciptome file vs operon file
                new_list_temp.append(int(x[3])) # start position
                new_list_temp.append(int(x[2])) # stop position
            else:
                break
            new_list_temp.append(blank_x[0]) # inserting 'X' blank for TSS start
            # if y[7] != 'x':
            #     new_list_temp.append(int(y[7])) # No Transcription file no TSS start position
            new_list_temp.append(blank_simple[0]) #simple/complex auto insert 'simple' because no transcription file
            new_list.append(new_list_temp)


if arguments == 3: # this if when you have transcriptome file hence 3 arguments
    for x in start_sites:
        counter = 0
        try:
            invalue = int(x[2])
        except ValueError:
            pass
        else:
            for y in operon_info:
                if x[0] == y[0]:
                    counter =+1
                    if counter == 1:
                        new_list_temp = []
                        new_list_temp.append(str(x[0])) #operon
                        new_list_temp.append(x[1]) #gene name/
                        new_list_temp.append(int(x[2])) # start
                        new_list_temp.append(int(x[3])) # / stop
                        new_list_temp.append(abs(int(x[2])-int(x[3]))) # gene length
                        new_list_temp.append(x[4]) # sense /
                        new_list_temp.append(y[6]) # in operon
                        new_list_temp.append(y[1]) # start gene
                        new_list_temp.append(y[2]) # stop gene
                        new_list_temp.append(int(y[3])) # start position
                        new_list_temp.append(int(y[4])) # stop position
                        if y[7] == 'x':
                            new_list_temp.append(y[7]) # start
                        if y[7] != 'x':
                            new_list_temp.append(int(y[7])) # start
                        new_list_temp.append(y[14]) # simple/complex
                        new_list.append(new_list_temp)
                    else:
                        print("Broken----------------------------->")

gene_operon = []
operon_pos = ['x']


for x in new_list:
    lag_dist = ['x']
    if x[1] == x[7]:
        operon_pos[0] = 1
    if x[1] != x[7]:
        operon_pos[0] += 1
    if x[5] == '+':
        if x[11] != 'x':
            lag_dist[0] = (x[3] - x[11])/2 #w/ TSS: changed from x[2] to x[3] to change 5' to 3'
            gene_operon.append(x + lag_dist + operon_pos)
        if x[11] == 'x':
            lag_dist[0] = (x[3] - x[9])/2 #w/o TSS:changed from x[2] to x[3] to change 5' to 3'
            gene_operon.append(x + lag_dist + operon_pos)
    if x[5] == '-':
        if x[11] != 'x':
            lag_dist[0] = (x[11] - x[2])/2 #w/ TSS: changed from x[3] to x[2] to change 5' to 3'
            gene_operon.append(x + lag_dist + operon_pos)
        if x[11] == 'x':
            lag_dist[0] = (x[9] - x[2])/2 #w/o TSS:changed from x[3] to x[2] to change 5' to 3'
            gene_operon.append(x + lag_dist + operon_pos)

# transcript length
# for x in new_list:
#     lag_dist = ['x']
#     if x[1] == x[7]:
#         operon_pos[0] = 1
#     if x[1] != x[7]:
#         operon_pos[0] += 1
#     if x[5] == '+':
#         if x[11] != 'x':
#             lag_dist[0] = (x[10] - x[11]) #TSS to stop
#             print  x[1],",",lag_dist[0]
#         if x[11] == 'x':
#             lag_dist[0] = (x[10] - x[9]) #start to stop
#             print  x[1],",",lag_dist[0]
#     if x[5] == '-':
#         if x[11] != 'x':
#             lag_dist[0] = (x[11] - x[10]) #TSS to stop
#             print  x[1],",",lag_dist[0]
#         if x[11] == 'x':
#             lag_dist[0] = (x[9] - x[10]) #start to stop
#             print x[1],",",lag_dist[0]
#
# exit()

log_counter = 0
data_log = []

for x in data: # log transforming the RPKM of the data, and getting rid of any genes that drop below detection then come back
    prime = ['x']
    second = ['x']
    third = ['x']
    forth = ['x']
    fifth = ['x']
    sixth = ['x']
    timepoint_7 = ['x']
    timepoint_8 = ['x']
    timepoint_9 = ['x']
    timepoint_10 = ['x']
    timepoint_11= ['x']
    timepoint_12 = ['x']
    timepoint_13 = ['x']
    timepoint_14 = ['x']
    timepoint_15 = ['x']
    timepoint_16 = ['x']
    timepoint_17 = ['x']
    timepoint_18 = ['x']
    timepoint_19 = ['x']
    timepoint_20 = ['x']
    if x[1] == 'x':
        add = ['x'] * 21
        add[0] = x[0]
        add[1] = 'x'
        dumb = ['x']
        dumb[0] = 'x'
        data_log.append(add+dumb)
        log_counter = log_counter + 1
    elif x[2] == 'x':
        add = ['x'] * 21
        add[0] = x[0]
        prime[0] = float(x[1])
        prime_adjust = np.log((prime[0] / prime[0]) * 100)
        add[1] = prime_adjust
        data_log.append(add+prime)
        log_counter = log_counter + 1
    elif x[3] == 'x':
        add = ['x'] * 21
        add[0] = x[0]
        prime[0] = float(x[1])
        prime_adjust = np.log((prime[0] / prime[0]) * 100)
        second[0] = float(x[2])
        second_adjust = np.log((second[0] / prime[0]) * 100)
        add[1] = prime_adjust
        add[2] = second_adjust
        data_log.append(add+prime)
        log_counter = log_counter + 1
    elif x[4] == 'x':
        add = ['x'] * 21
        add[0] = x[0]
        prime[0] = float(x[1])
        prime_adjust = np.log((prime[0] / prime[0]) * 100)
        second[0] = float(x[2])
        second_adjust = np.log((second[0] / prime[0]) * 100)
        third[0] = float(x[3])
        third_adjust = np.log((third[0] / prime[0]) * 100)
        add[1] = prime_adjust
        add[2] = second_adjust
        add[3] = third_adjust
        data_log.append(add+prime)
        log_counter = log_counter + 1
    elif x[5] == 'x':
        add = ['x'] * 21
        add[0] = x[0]
        prime[0] = float(x[1])
        prime_adjust = np.log((prime[0] / prime[0]) * 100)
        second[0] = float(x[2])
        second_adjust = np.log((second[0] / prime[0]) * 100)
        third[0] = float(x[3])
        third_adjust = np.log((third[0] / prime[0]) * 100)
        forth[0] = float(x[4])
        forth_adjust = np.log((forth[0] / prime[0]) * 100)
        add[1] = prime_adjust
        add[2] = second_adjust
        add[3] = third_adjust
        add[4] = forth_adjust
        data_log.append(add+prime)
        log_counter = log_counter + 1
    elif x[6] == 'x':
        add = ['x'] * 21
        add[0] = x[0]
        prime[0] = float(x[1])
        prime_adjust = np.log((prime[0] / prime[0]) * 100)
        second[0] = float(x[2])
        second_adjust = np.log((second[0] / prime[0]) * 100)
        third[0] = float(x[3])
        third_adjust = np.log((third[0] / prime[0]) * 100)
        forth[0] = float(x[4])
        forth_adjust = np.log((forth[0] / prime[0]) * 100)
        fifth[0] = float(x[5])
        fifth_adjust = np.log((fifth[0] / prime[0]) * 100)
        add[1] = prime_adjust
        add[2] = second_adjust
        add[3] = third_adjust
        add[4] = forth_adjust
        add[5] = fifth_adjust
        data_log.append(add+prime)
        log_counter = log_counter + 1
    elif x[7] == 'x':
        add = ['x'] * 21
        add[0] = x[0]
        prime[0] = float(x[1])
        prime_adjust = np.log((prime[0] / prime[0]) * 100)
        second[0] = float(x[2])
        second_adjust = np.log((second[0] / prime[0]) * 100)
        third[0] = float(x[3])
        third_adjust = np.log((third[0] / prime[0]) * 100)
        forth[0] = float(x[4])
        forth_adjust = np.log((forth[0] / prime[0]) * 100)
        fifth[0] = float(x[5])
        fifth_adjust = np.log((fifth[0] / prime[0]) * 100)
        sixth[0] = float(x[6])
        sixth_adjust = np.log((sixth[0] / prime[0]) * 100)
        add[1] = prime_adjust
        add[2] = second_adjust
        add[3] = third_adjust
        add[4] = forth_adjust
        add[5] = fifth_adjust
        add[6] = sixth_adjust
        data_log.append(add+prime)
        log_counter = log_counter + 1
    elif x[8] == 'x':
        add = ['x'] * 21
        add[0] = x[0]
        prime[0] = float(x[1])
        prime_adjust = np.log((prime[0] / prime[0]) * 100)
        second[0] = float(x[2])
        second_adjust = np.log((second[0] / prime[0]) * 100)
        third[0] = float(x[3])
        third_adjust = np.log((third[0] / prime[0]) * 100)
        forth[0] = float(x[4])
        forth_adjust = np.log((forth[0] / prime[0]) * 100)
        fifth[0] = float(x[5])
        fifth_adjust = np.log((fifth[0] / prime[0]) * 100)
        sixth[0] = float(x[6])
        sixth_adjust = np.log((sixth[0] / prime[0]) * 100)
        timepoint_7[0] = float(x[7])
        timepoint_7_adjust = np.log((timepoint_7[0] / prime[0]) * 100)
        add[1] = prime_adjust
        add[2] = second_adjust
        add[3] = third_adjust
        add[4] = forth_adjust
        add[5] = fifth_adjust
        add[6] = sixth_adjust
        add[7] = timepoint_7_adjust
        data_log.append(add+prime)
        log_counter = log_counter + 1
    elif x[9] == 'x':
        add = ['x'] * 21
        add[0] = x[0]
        prime[0] = float(x[1])
        prime_adjust = np.log((prime[0] / prime[0]) * 100)
        second[0] = float(x[2])
        second_adjust = np.log((second[0] / prime[0]) * 100)
        third[0] = float(x[3])
        third_adjust = np.log((third[0] / prime[0]) * 100)
        forth[0] = float(x[4])
        forth_adjust = np.log((forth[0] / prime[0]) * 100)
        fifth[0] = float(x[5])
        fifth_adjust = np.log((fifth[0] / prime[0]) * 100)
        sixth[0] = float(x[6])
        sixth_adjust = np.log((sixth[0] / prime[0]) * 100)
        timepoint_7[0] = float(x[7])
        timepoint_7_adjust = np.log((timepoint_7[0] / prime[0]) * 100)
        timepoint_8[0] = float(x[8])
        timepoint_8_adjust = np.log((timepoint_8[0] / prime[0]) * 100)
        add[1] = prime_adjust
        add[2] = second_adjust
        add[3] = third_adjust
        add[4] = forth_adjust
        add[5] = fifth_adjust
        add[6] = sixth_adjust
        add[7] = timepoint_7_adjust
        add[8] = timepoint_8_adjust
        data_log.append(add+prime)
        log_counter = log_counter + 1
    elif x[10] == 'x':
        add = ['x'] * 21
        add[0] = x[0]
        prime[0] = float(x[1])
        prime_adjust = np.log((prime[0] / prime[0]) * 100)
        second[0] = float(x[2])
        second_adjust = np.log((second[0] / prime[0]) * 100)
        third[0] = float(x[3])
        third_adjust = np.log((third[0] / prime[0]) * 100)
        forth[0] = float(x[4])
        forth_adjust = np.log((forth[0] / prime[0]) * 100)
        fifth[0] = float(x[5])
        fifth_adjust = np.log((fifth[0] / prime[0]) * 100)
        sixth[0] = float(x[6])
        sixth_adjust = np.log((sixth[0] / prime[0]) * 100)
        timepoint_7[0] = float(x[7])
        timepoint_7_adjust = np.log((timepoint_7[0] / prime[0]) * 100)
        timepoint_8[0] = float(x[8])
        timepoint_8_adjust = np.log((timepoint_8[0] / prime[0]) * 100)
        timepoint_9[0] = float(x[9])
        timepoint_9_adjust = np.log((timepoint_9[0] / prime[0]) * 100)
        add[1] = prime_adjust
        add[2] = second_adjust
        add[3] = third_adjust
        add[4] = forth_adjust
        add[5] = fifth_adjust
        add[6] = sixth_adjust
        add[7] = timepoint_7_adjust
        add[8] = timepoint_8_adjust
        add[9] = timepoint_9_adjust
        data_log.append(add+prime)
        log_counter = log_counter + 1
    elif x[11] == 'x':
        add = ['x'] * 21
        add[0] = x[0]
        prime[0] = float(x[1])
        prime_adjust = np.log((prime[0] / prime[0]) * 100)
        second[0] = float(x[2])
        second_adjust = np.log((second[0] / prime[0]) * 100)
        third[0] = float(x[3])
        third_adjust = np.log((third[0] / prime[0]) * 100)
        forth[0] = float(x[4])
        forth_adjust = np.log((forth[0] / prime[0]) * 100)
        fifth[0] = float(x[5])
        fifth_adjust = np.log((fifth[0] / prime[0]) * 100)
        sixth[0] = float(x[6])
        sixth_adjust = np.log((sixth[0] / prime[0]) * 100)
        timepoint_7[0] = float(x[7])
        timepoint_7_adjust = np.log((timepoint_7[0] / prime[0]) * 100)
        timepoint_8[0] = float(x[8])
        timepoint_8_adjust = np.log((timepoint_8[0] / prime[0]) * 100)
        timepoint_9[0] = float(x[9])
        timepoint_9_adjust = np.log((timepoint_9[0] / prime[0]) * 100)
        timepoint_10[0] = float(x[10])
        timepoint_10_adjust = np.log((timepoint_10[0] / prime[0]) * 100)
        add[1] = prime_adjust
        add[2] = second_adjust
        add[3] = third_adjust
        add[4] = forth_adjust
        add[5] = fifth_adjust
        add[6] = sixth_adjust
        add[7] = timepoint_7_adjust
        add[8] = timepoint_8_adjust
        add[9] = timepoint_9_adjust
        add[10] = timepoint_10_adjust
        data_log.append(add+prime)
        log_counter = log_counter + 1
    elif x[12] == 'x':
        add = ['x'] * 21
        add[0] = x[0]
        prime[0] = float(x[1])
        prime_adjust = np.log((prime[0] / prime[0]) * 100)
        second[0] = float(x[2])
        second_adjust = np.log((second[0] / prime[0]) * 100)
        third[0] = float(x[3])
        third_adjust = np.log((third[0] / prime[0]) * 100)
        forth[0] = float(x[4])
        forth_adjust = np.log((forth[0] / prime[0]) * 100)
        fifth[0] = float(x[5])
        fifth_adjust = np.log((fifth[0] / prime[0]) * 100)
        sixth[0] = float(x[6])
        sixth_adjust = np.log((sixth[0] / prime[0]) * 100)
        timepoint_7[0] = float(x[7])
        timepoint_7_adjust = np.log((timepoint_7[0] / prime[0]) * 100)
        timepoint_8[0] = float(x[8])
        timepoint_8_adjust = np.log((timepoint_8[0] / prime[0]) * 100)
        timepoint_9[0] = float(x[9])
        timepoint_9_adjust = np.log((timepoint_9[0] / prime[0]) * 100)
        timepoint_10[0] = float(x[10])
        timepoint_10_adjust = np.log((timepoint_10[0] / prime[0]) * 100)
        timepoint_11[0] = float(x[11])
        timepoint_11_adjust = np.log((timepoint_11[0] / prime[0]) * 100)
        add[1] = prime_adjust
        add[2] = second_adjust
        add[3] = third_adjust
        add[4] = forth_adjust
        add[5] = fifth_adjust
        add[6] = sixth_adjust
        add[7] = timepoint_7_adjust
        add[8] = timepoint_8_adjust
        add[9] = timepoint_9_adjust
        add[10] = timepoint_10_adjust
        add[11] = timepoint_11_adjust
        data_log.append(add+prime)
        log_counter = log_counter + 1
    elif x[13] == 'x':
        add = ['x'] * 21
        add[0] = x[0]
        prime[0] = float(x[1])
        prime_adjust = np.log((prime[0] / prime[0]) * 100)
        second[0] = float(x[2])
        second_adjust = np.log((second[0] / prime[0]) * 100)
        third[0] = float(x[3])
        third_adjust = np.log((third[0] / prime[0]) * 100)
        forth[0] = float(x[4])
        forth_adjust = np.log((forth[0] / prime[0]) * 100)
        fifth[0] = float(x[5])
        fifth_adjust = np.log((fifth[0] / prime[0]) * 100)
        sixth[0] = float(x[6])
        sixth_adjust = np.log((sixth[0] / prime[0]) * 100)
        timepoint_7[0] = float(x[7])
        timepoint_7_adjust = np.log((timepoint_7[0] / prime[0]) * 100)
        timepoint_8[0] = float(x[8])
        timepoint_8_adjust = np.log((timepoint_8[0] / prime[0]) * 100)
        timepoint_9[0] = float(x[9])
        timepoint_9_adjust = np.log((timepoint_9[0] / prime[0]) * 100)
        timepoint_10[0] = float(x[10])
        timepoint_10_adjust = np.log((timepoint_10[0] / prime[0]) * 100)
        timepoint_11[0] = float(x[11])
        timepoint_11_adjust = np.log((timepoint_11[0] / prime[0]) * 100)
        timepoint_12[0] = float(x[12])
        timepoint_12_adjust = np.log((timepoint_12[0] / prime[0]) * 100)
        add[1] = prime_adjust
        add[2] = second_adjust
        add[3] = third_adjust
        add[4] = forth_adjust
        add[5] = fifth_adjust
        add[6] = sixth_adjust
        add[7] = timepoint_7_adjust
        add[8] = timepoint_8_adjust
        add[9] = timepoint_9_adjust
        add[10] = timepoint_10_adjust
        add[11] = timepoint_11_adjust
        add[12] = timepoint_12_adjust
        data_log.append(add+prime)
        log_counter = log_counter + 1
    elif x[14] == 'x':
        add = ['x'] * 21
        add[0] = x[0]
        prime[0] = float(x[1])
        prime_adjust = np.log((prime[0] / prime[0]) * 100)
        second[0] = float(x[2])
        second_adjust = np.log((second[0] / prime[0]) * 100)
        third[0] = float(x[3])
        third_adjust = np.log((third[0] / prime[0]) * 100)
        forth[0] = float(x[4])
        forth_adjust = np.log((forth[0] / prime[0]) * 100)
        fifth[0] = float(x[5])
        fifth_adjust = np.log((fifth[0] / prime[0]) * 100)
        sixth[0] = float(x[6])
        sixth_adjust = np.log((sixth[0] / prime[0]) * 100)
        timepoint_7[0] = float(x[7])
        timepoint_7_adjust = np.log((timepoint_7[0] / prime[0]) * 100)
        timepoint_8[0] = float(x[8])
        timepoint_8_adjust = np.log((timepoint_8[0] / prime[0]) * 100)
        timepoint_9[0] = float(x[9])
        timepoint_9_adjust = np.log((timepoint_9[0] / prime[0]) * 100)
        timepoint_10[0] = float(x[10])
        timepoint_10_adjust = np.log((timepoint_10[0] / prime[0]) * 100)
        timepoint_11[0] = float(x[11])
        timepoint_11_adjust = np.log((timepoint_11[0] / prime[0]) * 100)
        timepoint_12[0] = float(x[12])
        timepoint_12_adjust = np.log((timepoint_12[0] / prime[0]) * 100)
        timepoint_13[0] = float(x[13])
        timepoint_13_adjust = np.log((timepoint_13[0] / prime[0]) * 100)
        add[1] = prime_adjust
        add[2] = second_adjust
        add[3] = third_adjust
        add[4] = forth_adjust
        add[5] = fifth_adjust
        add[6] = sixth_adjust
        add[7] = timepoint_7_adjust
        add[8] = timepoint_8_adjust
        add[9] = timepoint_9_adjust
        add[10] = timepoint_10_adjust
        add[11] = timepoint_11_adjust
        add[12] = timepoint_12_adjust
        add[13] = timepoint_13_adjust
        data_log.append(add+prime)
        log_counter = log_counter + 1
    elif x[15] == 'x':
        add = ['x'] * 21
        add[0] = x[0]
        prime[0] = float(x[1])
        prime_adjust = np.log((prime[0] / prime[0]) * 100)
        second[0] = float(x[2])
        second_adjust = np.log((second[0] / prime[0]) * 100)
        third[0] = float(x[3])
        third_adjust = np.log((third[0] / prime[0]) * 100)
        forth[0] = float(x[4])
        forth_adjust = np.log((forth[0] / prime[0]) * 100)
        fifth[0] = float(x[5])
        fifth_adjust = np.log((fifth[0] / prime[0]) * 100)
        sixth[0] = float(x[6])
        sixth_adjust = np.log((sixth[0] / prime[0]) * 100)
        timepoint_7[0] = float(x[7])
        timepoint_7_adjust = np.log((timepoint_7[0] / prime[0]) * 100)
        timepoint_8[0] = float(x[8])
        timepoint_8_adjust = np.log((timepoint_8[0] / prime[0]) * 100)
        timepoint_9[0] = float(x[9])
        timepoint_9_adjust = np.log((timepoint_9[0] / prime[0]) * 100)
        timepoint_10[0] = float(x[10])
        timepoint_10_adjust = np.log((timepoint_10[0] / prime[0]) * 100)
        timepoint_11[0] = float(x[11])
        timepoint_11_adjust = np.log((timepoint_11[0] / prime[0]) * 100)
        timepoint_12[0] = float(x[12])
        timepoint_12_adjust = np.log((timepoint_12[0] / prime[0]) * 100)
        timepoint_13[0] = float(x[13])
        timepoint_13_adjust = np.log((timepoint_13[0] / prime[0]) * 100)
        timepoint_14[0] = float(x[14])
        timepoint_14_adjust = np.log((timepoint_14[0] / prime[0]) * 100)
        add[1] = prime_adjust
        add[2] = second_adjust
        add[3] = third_adjust
        add[4] = forth_adjust
        add[5] = fifth_adjust
        add[6] = sixth_adjust
        add[7] = timepoint_7_adjust
        add[8] = timepoint_8_adjust
        add[9] = timepoint_9_adjust
        add[10] = timepoint_10_adjust
        add[11] = timepoint_11_adjust
        add[12] = timepoint_12_adjust
        add[13] = timepoint_13_adjust
        add[14] = timepoint_14_adjust
        data_log.append(add+prime)
        log_counter = log_counter + 1
    elif x[16] == 'x':
        add = ['x'] * 21
        add[0] = x[0]
        prime[0] = float(x[1])
        prime_adjust = np.log((prime[0] / prime[0]) * 100)
        second[0] = float(x[2])
        second_adjust = np.log((second[0] / prime[0]) * 100)
        third[0] = float(x[3])
        third_adjust = np.log((third[0] / prime[0]) * 100)
        forth[0] = float(x[4])
        forth_adjust = np.log((forth[0] / prime[0]) * 100)
        fifth[0] = float(x[5])
        fifth_adjust = np.log((fifth[0] / prime[0]) * 100)
        sixth[0] = float(x[6])
        sixth_adjust = np.log((sixth[0] / prime[0]) * 100)
        timepoint_7[0] = float(x[7])
        timepoint_7_adjust = np.log((timepoint_7[0] / prime[0]) * 100)
        timepoint_8[0] = float(x[8])
        timepoint_8_adjust = np.log((timepoint_8[0] / prime[0]) * 100)
        timepoint_9[0] = float(x[9])
        timepoint_9_adjust = np.log((timepoint_9[0] / prime[0]) * 100)
        timepoint_10[0] = float(x[10])
        timepoint_10_adjust = np.log((timepoint_10[0] / prime[0]) * 100)
        timepoint_11[0] = float(x[11])
        timepoint_11_adjust = np.log((timepoint_11[0] / prime[0]) * 100)
        timepoint_12[0] = float(x[12])
        timepoint_12_adjust = np.log((timepoint_12[0] / prime[0]) * 100)
        timepoint_13[0] = float(x[13])
        timepoint_13_adjust = np.log((timepoint_13[0] / prime[0]) * 100)
        timepoint_14[0] = float(x[14])
        timepoint_14_adjust = np.log((timepoint_14[0] / prime[0]) * 100)
        timepoint_15[0] = float(x[15])
        timepoint_15_adjust = np.log((timepoint_15[0] / prime[0]) * 100)
        add[1] = prime_adjust
        add[2] = second_adjust
        add[3] = third_adjust
        add[4] = forth_adjust
        add[5] = fifth_adjust
        add[6] = sixth_adjust
        add[7] = timepoint_7_adjust
        add[8] = timepoint_8_adjust
        add[9] = timepoint_9_adjust
        add[10] = timepoint_10_adjust
        add[11] = timepoint_11_adjust
        add[12] = timepoint_12_adjust
        add[13] = timepoint_13_adjust
        add[14] = timepoint_14_adjust
        add[15] = timepoint_15_adjust
        data_log.append(add+prime)
        log_counter = log_counter + 1
    elif x[17] == 'x':
        add = ['x'] * 21
        add[0] = x[0]
        prime[0] = float(x[1])
        prime_adjust = np.log((prime[0] / prime[0]) * 100)
        second[0] = float(x[2])
        second_adjust = np.log((second[0] / prime[0]) * 100)
        third[0] = float(x[3])
        third_adjust = np.log((third[0] / prime[0]) * 100)
        forth[0] = float(x[4])
        forth_adjust = np.log((forth[0] / prime[0]) * 100)
        fifth[0] = float(x[5])
        fifth_adjust = np.log((fifth[0] / prime[0]) * 100)
        sixth[0] = float(x[6])
        sixth_adjust = np.log((sixth[0] / prime[0]) * 100)
        timepoint_7[0] = float(x[7])
        timepoint_7_adjust = np.log((timepoint_7[0] / prime[0]) * 100)
        timepoint_8[0] = float(x[8])
        timepoint_8_adjust = np.log((timepoint_8[0] / prime[0]) * 100)
        timepoint_9[0] = float(x[9])
        timepoint_9_adjust = np.log((timepoint_9[0] / prime[0]) * 100)
        timepoint_10[0] = float(x[10])
        timepoint_10_adjust = np.log((timepoint_10[0] / prime[0]) * 100)
        timepoint_11[0] = float(x[11])
        timepoint_11_adjust = np.log((timepoint_11[0] / prime[0]) * 100)
        timepoint_12[0] = float(x[12])
        timepoint_12_adjust = np.log((timepoint_12[0] / prime[0]) * 100)
        timepoint_13[0] = float(x[13])
        timepoint_13_adjust = np.log((timepoint_13[0] / prime[0]) * 100)
        timepoint_14[0] = float(x[14])
        timepoint_14_adjust = np.log((timepoint_14[0] / prime[0]) * 100)
        timepoint_15[0] = float(x[15])
        timepoint_15_adjust = np.log((timepoint_15[0] / prime[0]) * 100)
        timepoint_16[0] = float(x[16])
        timepoint_16_adjust = np.log((timepoint_16[0] / prime[0]) * 100)
        add[1] = prime_adjust
        add[2] = second_adjust
        add[3] = third_adjust
        add[4] = forth_adjust
        add[5] = fifth_adjust
        add[6] = sixth_adjust
        add[7] = timepoint_7_adjust
        add[8] = timepoint_8_adjust
        add[9] = timepoint_9_adjust
        add[10] = timepoint_10_adjust
        add[11] = timepoint_11_adjust
        add[12] = timepoint_12_adjust
        add[13] = timepoint_13_adjust
        add[14] = timepoint_14_adjust
        add[15] = timepoint_15_adjust
        add[16] = timepoint_16_adjust
        data_log.append(add+prime)
        log_counter = log_counter + 1
    elif x[18] == 'x':
        add = ['x'] * 21
        add[0] = x[0]
        prime[0] = float(x[1])
        prime_adjust = np.log((prime[0] / prime[0]) * 100)
        second[0] = float(x[2])
        second_adjust = np.log((second[0] / prime[0]) * 100)
        third[0] = float(x[3])
        third_adjust = np.log((third[0] / prime[0]) * 100)
        forth[0] = float(x[4])
        forth_adjust = np.log((forth[0] / prime[0]) * 100)
        fifth[0] = float(x[5])
        fifth_adjust = np.log((fifth[0] / prime[0]) * 100)
        sixth[0] = float(x[6])
        sixth_adjust = np.log((sixth[0] / prime[0]) * 100)
        timepoint_7[0] = float(x[7])
        timepoint_7_adjust = np.log((timepoint_7[0] / prime[0]) * 100)
        timepoint_8[0] = float(x[8])
        timepoint_8_adjust = np.log((timepoint_8[0] / prime[0]) * 100)
        timepoint_9[0] = float(x[9])
        timepoint_9_adjust = np.log((timepoint_9[0] / prime[0]) * 100)
        timepoint_10[0] = float(x[10])
        timepoint_10_adjust = np.log((timepoint_10[0] / prime[0]) * 100)
        timepoint_11[0] = float(x[11])
        timepoint_11_adjust = np.log((timepoint_11[0] / prime[0]) * 100)
        timepoint_12[0] = float(x[12])
        timepoint_12_adjust = np.log((timepoint_12[0] / prime[0]) * 100)
        timepoint_13[0] = float(x[13])
        timepoint_13_adjust = np.log((timepoint_13[0] / prime[0]) * 100)
        timepoint_14[0] = float(x[14])
        timepoint_14_adjust = np.log((timepoint_14[0] / prime[0]) * 100)
        timepoint_15[0] = float(x[15])
        timepoint_15_adjust = np.log((timepoint_15[0] / prime[0]) * 100)
        timepoint_16[0] = float(x[16])
        timepoint_16_adjust = np.log((timepoint_16[0] / prime[0]) * 100)
        timepoint_17[0] = float(x[17])
        timepoint_17_adjust = np.log((timepoint_17[0] / prime[0]) * 100)
        add[1] = prime_adjust
        add[2] = second_adjust
        add[3] = third_adjust
        add[4] = forth_adjust
        add[5] = fifth_adjust
        add[6] = sixth_adjust
        add[7] = timepoint_7_adjust
        add[8] = timepoint_8_adjust
        add[9] = timepoint_9_adjust
        add[10] = timepoint_10_adjust
        add[11] = timepoint_11_adjust
        add[12] = timepoint_12_adjust
        add[13] = timepoint_13_adjust
        add[14] = timepoint_14_adjust
        add[15] = timepoint_15_adjust
        add[16] = timepoint_16_adjust
        add[17] = timepoint_17_adjust
        data_log.append(add+prime)
        log_counter = log_counter + 1
    elif x[19] == 'x':
        add = ['x'] * 21
        add[0] = x[0]
        prime[0] = float(x[1])
        prime_adjust = np.log((prime[0] / prime[0]) * 100)
        second[0] = float(x[2])
        second_adjust = np.log((second[0] / prime[0]) * 100)
        third[0] = float(x[3])
        third_adjust = np.log((third[0] / prime[0]) * 100)
        forth[0] = float(x[4])
        forth_adjust = np.log((forth[0] / prime[0]) * 100)
        fifth[0] = float(x[5])
        fifth_adjust = np.log((fifth[0] / prime[0]) * 100)
        sixth[0] = float(x[6])
        sixth_adjust = np.log((sixth[0] / prime[0]) * 100)
        timepoint_7[0] = float(x[7])
        timepoint_7_adjust = np.log((timepoint_7[0] / prime[0]) * 100)
        timepoint_8[0] = float(x[8])
        timepoint_8_adjust = np.log((timepoint_8[0] / prime[0]) * 100)
        timepoint_9[0] = float(x[9])
        timepoint_9_adjust = np.log((timepoint_9[0] / prime[0]) * 100)
        timepoint_10[0] = float(x[10])
        timepoint_10_adjust = np.log((timepoint_10[0] / prime[0]) * 100)
        timepoint_11[0] = float(x[11])
        timepoint_11_adjust = np.log((timepoint_11[0] / prime[0]) * 100)
        timepoint_12[0] = float(x[12])
        timepoint_12_adjust = np.log((timepoint_12[0] / prime[0]) * 100)
        timepoint_13[0] = float(x[13])
        timepoint_13_adjust = np.log((timepoint_13[0] / prime[0]) * 100)
        timepoint_14[0] = float(x[14])
        timepoint_14_adjust = np.log((timepoint_14[0] / prime[0]) * 100)
        timepoint_15[0] = float(x[15])
        timepoint_15_adjust = np.log((timepoint_15[0] / prime[0]) * 100)
        timepoint_16[0] = float(x[16])
        timepoint_16_adjust = np.log((timepoint_16[0] / prime[0]) * 100)
        timepoint_17[0] = float(x[17])
        timepoint_17_adjust = np.log((timepoint_17[0] / prime[0]) * 100)
        timepoint_18[0] = float(x[18])
        timepoint_18_adjust = np.log((timepoint_18[0] / prime[0]) * 100)
        add[1] = prime_adjust
        add[2] = second_adjust
        add[3] = third_adjust
        add[4] = forth_adjust
        add[5] = fifth_adjust
        add[6] = sixth_adjust
        add[7] = timepoint_7_adjust
        add[8] = timepoint_8_adjust
        add[9] = timepoint_9_adjust
        add[10] = timepoint_10_adjust
        add[11] = timepoint_11_adjust
        add[12] = timepoint_12_adjust
        add[13] = timepoint_13_adjust
        add[14] = timepoint_14_adjust
        add[15] = timepoint_15_adjust
        add[16] = timepoint_16_adjust
        add[17] = timepoint_17_adjust
        add[18] = timepoint_18_adjust
        data_log.append(add+prime)
        log_counter = log_counter + 1
    elif x[20] == 'x':
        add = ['x'] * 21
        add[0] = x[0]
        prime[0] = float(x[1])
        prime_adjust = np.log((prime[0] / prime[0]) * 100)
        second[0] = float(x[2])
        second_adjust = np.log((second[0] / prime[0]) * 100)
        third[0] = float(x[3])
        third_adjust = np.log((third[0] / prime[0]) * 100)
        forth[0] = float(x[4])
        forth_adjust = np.log((forth[0] / prime[0]) * 100)
        fifth[0] = float(x[5])
        fifth_adjust = np.log((fifth[0] / prime[0]) * 100)
        sixth[0] = float(x[6])
        sixth_adjust = np.log((sixth[0] / prime[0]) * 100)
        timepoint_7[0] = float(x[7])
        timepoint_7_adjust = np.log((timepoint_7[0] / prime[0]) * 100)
        timepoint_8[0] = float(x[8])
        timepoint_8_adjust = np.log((timepoint_8[0] / prime[0]) * 100)
        timepoint_9[0] = float(x[9])
        timepoint_9_adjust = np.log((timepoint_9[0] / prime[0]) * 100)
        timepoint_10[0] = float(x[10])
        timepoint_10_adjust = np.log((timepoint_10[0] / prime[0]) * 100)
        timepoint_11[0] = float(x[11])
        timepoint_11_adjust = np.log((timepoint_11[0] / prime[0]) * 100)
        timepoint_12[0] = float(x[12])
        timepoint_12_adjust = np.log((timepoint_12[0] / prime[0]) * 100)
        timepoint_13[0] = float(x[13])
        timepoint_13_adjust = np.log((timepoint_13[0] / prime[0]) * 100)
        timepoint_14[0] = float(x[14])
        timepoint_14_adjust = np.log((timepoint_14[0] / prime[0]) * 100)
        timepoint_15[0] = float(x[15])
        timepoint_15_adjust = np.log((timepoint_15[0] / prime[0]) * 100)
        timepoint_16[0] = float(x[16])
        timepoint_16_adjust = np.log((timepoint_16[0] / prime[0]) * 100)
        timepoint_17[0] = float(x[17])
        timepoint_17_adjust = np.log((timepoint_17[0] / prime[0]) * 100)
        timepoint_18[0] = float(x[18])
        timepoint_18_adjust = np.log((timepoint_18[0] / prime[0]) * 100)
        timepoint_19[0] = float(x[19])
        timepoint_19_adjust = np.log((timepoint_19[0] / prime[0]) * 100)
        add[1] = prime_adjust
        add[2] = second_adjust
        add[3] = third_adjust
        add[4] = forth_adjust
        add[5] = fifth_adjust
        add[6] = sixth_adjust
        add[7] = timepoint_7_adjust
        add[8] = timepoint_8_adjust
        add[9] = timepoint_9_adjust
        add[10] = timepoint_10_adjust
        add[11] = timepoint_11_adjust
        add[12] = timepoint_12_adjust
        add[13] = timepoint_13_adjust
        add[14] = timepoint_14_adjust
        add[15] = timepoint_15_adjust
        add[16] = timepoint_16_adjust
        add[17] = timepoint_17_adjust
        add[18] = timepoint_18_adjust
        add[19] = timepoint_19_adjust
        data_log.append(add+prime)
        log_counter = log_counter + 1
    elif x[21] == 'x':
        add = ['x'] * 21
        add[0] = x[0]
        prime[0] = float(x[1])
        prime_adjust = np.log((prime[0] / prime[0]) * 100)
        second[0] = float(x[2])
        second_adjust = np.log((second[0] / prime[0]) * 100)
        third[0] = float(x[3])
        third_adjust = np.log((third[0] / prime[0]) * 100)
        forth[0] = float(x[4])
        forth_adjust = np.log((forth[0] / prime[0]) * 100)
        fifth[0] = float(x[5])
        fifth_adjust = np.log((fifth[0] / prime[0]) * 100)
        sixth[0] = float(x[6])
        sixth_adjust = np.log((sixth[0] / prime[0]) * 100)
        timepoint_7[0] = float(x[7])
        timepoint_7_adjust = np.log((timepoint_7[0] / prime[0]) * 100)
        timepoint_8[0] = float(x[8])
        timepoint_8_adjust = np.log((timepoint_8[0] / prime[0]) * 100)
        timepoint_9[0] = float(x[9])
        timepoint_9_adjust = np.log((timepoint_9[0] / prime[0]) * 100)
        timepoint_10[0] = float(x[10])
        timepoint_10_adjust = np.log((timepoint_10[0] / prime[0]) * 100)
        timepoint_11[0] = float(x[11])
        timepoint_11_adjust = np.log((timepoint_11[0] / prime[0]) * 100)
        timepoint_12[0] = float(x[12])
        timepoint_12_adjust = np.log((timepoint_12[0] / prime[0]) * 100)
        timepoint_13[0] = float(x[13])
        timepoint_13_adjust = np.log((timepoint_13[0] / prime[0]) * 100)
        timepoint_14[0] = float(x[14])
        timepoint_14_adjust = np.log((timepoint_14[0] / prime[0]) * 100)
        timepoint_15[0] = float(x[15])
        timepoint_15_adjust = np.log((timepoint_15[0] / prime[0]) * 100)
        timepoint_16[0] = float(x[16])
        timepoint_16_adjust = np.log((timepoint_16[0] / prime[0]) * 100)
        timepoint_17[0] = float(x[17])
        timepoint_17_adjust = np.log((timepoint_17[0] / prime[0]) * 100)
        timepoint_18[0] = float(x[18])
        timepoint_18_adjust = np.log((timepoint_18[0] / prime[0]) * 100)
        timepoint_19[0] = float(x[19])
        timepoint_19_adjust = np.log((timepoint_19[0] / prime[0]) * 100)
        timepoint_20[0] = float(x[20])
        timepoint_20_adjust = np.log((timepoint_20[0] / prime[0]) * 100)
        add[1] = prime_adjust
        add[2] = second_adjust
        add[3] = third_adjust
        add[4] = forth_adjust
        add[5] = fifth_adjust
        add[6] = sixth_adjust
        add[7] = timepoint_7_adjust
        add[8] = timepoint_8_adjust
        add[9] = timepoint_9_adjust
        add[10] = timepoint_10_adjust
        add[11] = timepoint_11_adjust
        add[12] = timepoint_12_adjust
        add[13] = timepoint_13_adjust
        add[14] = timepoint_14_adjust
        add[15] = timepoint_15_adjust
        add[16] = timepoint_16_adjust
        add[17] = timepoint_17_adjust
        add[18] = timepoint_18_adjust
        add[19] = timepoint_19_adjust
        add[20] = timepoint_20_adjust
        data_log.append(add+prime)
        log_counter = log_counter + 1
    else:
        print('broke!!!!!!!!!!!!!!!!!!!!!!!')
        print(x)

combined_genes = []
combined_counter = 0

for x in gene_operon: # combine RPKM data with the gene data
    for y in data_log:
        if x[1] == y[0]:
            combined_genes.append(x + y[1:])
            combined_counter += 1

#next 2 lines set the header of the output file with titles, and times
csv_results_table = [["gene","regres. slope","regres. Y Int","regres. R value","Half-Life Calc (minutes)","# of data pts. After cuttoff","# of leading data points affected by lag calc","ln adjusted value of last data point used in regres.","In Operon (y or n)","Operon #","complex or simple","Total # of data points","R value flag","Half-life Flag"]]
csv_results_table[-1].extend(time_points_set)
data_log_corrected = []
time_corrected = []
operon_count_total = 0
#elong_rate = 19.36*60 # = nt/min.
# elong_rate = 22.143*60


#determining time points that are below cutoff and makes sure they are sequential (non-continous readings)
for x in range(0,len(combined_genes)):

    z = 'x'
    length = 'x'
    new_time = ['x']
    if combined_genes[x][15] == 'x' and combined_genes[x][16] == 'x': # no possible time points before cutoff
        i = 0
        length = 0
    if combined_genes[x][15] != 'x' and combined_genes[x][16] == 'x': # one possible time point before cutoff
        i = 1
        if combined_genes[x][15] >= cut_off:
            length = 1
        elif combined_genes[x][15] < cut_off:
            length = 0
        else:
            "Broken ------------------> 1"
    if combined_genes[x][16] != 'x' and combined_genes[x][17] == 'x': # 2 possible time points before cutoff
        i = 2
        if combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off: # 2 earliest points above cutoff
            length = 2
        elif combined_genes[x][16] < cut_off and combined_genes[x][15] >= cut_off: # first point above cutoff
            length = 1
        elif combined_genes[x][15] < cut_off: # first point isn't above cutoff
            length = 0
        else:
            "Broken ------------------> 2"
    if combined_genes[x][17] != 'x' and combined_genes[x][18] == 'x': # 3 possible time points before cutoff
        i = 3
        if combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off: # all 3 points above cutoff
            length = 3
        elif combined_genes[x][17] < cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off: # 2 earliest points above cutoff
            length = 2
        elif combined_genes[x][16] < cut_off and combined_genes[x][15] >= cut_off: # first point above cutoff
            length = 1
        elif combined_genes[x][15] < cut_off: # first point isn't above cutoff
            length = 0
        else:
            "Broken ------------------> 3"
    if combined_genes[x][18] != 'x' and combined_genes[x][19] == 'x':
        i = 4
        if combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 4
        elif combined_genes[x][18] < cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 3
        elif combined_genes[x][17] < cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 2
        elif combined_genes[x][16] < cut_off and combined_genes[x][15] >= cut_off:
            length = 1
        elif combined_genes[x][15] < cut_off:
            length = 0
        else:
            "Broken ------------------> 4"
    if combined_genes[x][19] != 'x' and combined_genes[x][20] == 'x':
        i = 5
        if combined_genes[x][19] >= cut_off and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 5
        elif combined_genes[x][19] < cut_off and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 4
        elif combined_genes[x][18] < cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 3
        elif combined_genes[x][17] < cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 2
        elif combined_genes[x][16] < cut_off and combined_genes[x][15] >= cut_off:
            length = 1
        elif combined_genes[x][15] < cut_off:
            length = 0
        else:
            "Broken ------------------> 5"
    if combined_genes[x][20] != 'x' and combined_genes[x][21] == 'x':
        i = 6
        if combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 6
        elif combined_genes[x][20] < cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 5
        elif combined_genes[x][19] < cut_off and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 4
        elif combined_genes[x][18] < cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 3
        elif combined_genes[x][17] < cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 2
        elif combined_genes[x][16] < cut_off and combined_genes[x][15] >= cut_off:
            length = 1
        elif combined_genes[x][15] < cut_off:
            length = 0
        else:
            "Broken ------------------> 6"
    if combined_genes[x][21] != 'x' and combined_genes[x][22] == 'x':
        i = 7
        if combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 7
        elif combined_genes[x][21] < cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 6
        elif combined_genes[x][20] < cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 5
        elif combined_genes[x][19] < cut_off and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 4
        elif combined_genes[x][18] < cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 3
        elif combined_genes[x][17] < cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 2
        elif combined_genes[x][16] < cut_off and combined_genes[x][15] >= cut_off:
            length = 1
        elif combined_genes[x][15] < cut_off:
            length = 0
        else:
            "Broken ------------------> 7"
    if combined_genes[x][22] != 'x' and combined_genes[x][23] == 'x':
        i = 8
        if combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 8
        elif combined_genes[x][22] < cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 7
        elif combined_genes[x][21] < cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 6
        elif combined_genes[x][20] < cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 5
        elif combined_genes[x][19] < cut_off and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 4
        elif combined_genes[x][18] < cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 3
        elif combined_genes[x][17] < cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 2
        elif combined_genes[x][16] < cut_off and combined_genes[x][15] >= cut_off:
            length = 1
        elif combined_genes[x][15] < cut_off:
            length = 0
        else:
            "Broken ------------------> 8"
    if combined_genes[x][23] != 'x' and combined_genes[x][24] == 'x':
        i = 9
        if combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 9
        elif combined_genes[x][23] < cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 8
        elif combined_genes[x][22] < cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 7
        elif combined_genes[x][21] < cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 6
        elif combined_genes[x][20] < cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 5
        elif combined_genes[x][19] < cut_off and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 4
        elif combined_genes[x][18] < cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 3
        elif combined_genes[x][17] < cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 2
        elif combined_genes[x][16] < cut_off and combined_genes[x][15] >= cut_off:
            length = 1
        elif combined_genes[x][15] < cut_off:
            length = 0
        else:
            "Broken ------------------> 9"
    if combined_genes[x][24] != 'x' and combined_genes[x][25] == 'x':
        i = 10
        if combined_genes[x][24] >= cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 10
        elif combined_genes[x][24] < cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 9
        elif combined_genes[x][23] < cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 8
        elif combined_genes[x][22] < cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 7
        elif combined_genes[x][21] < cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 6
        elif combined_genes[x][20] < cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 5
        elif combined_genes[x][19] < cut_off and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 4
        elif combined_genes[x][18] < cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 3
        elif combined_genes[x][17] < cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 2
        elif combined_genes[x][16] < cut_off and combined_genes[x][15] >= cut_off:
            length = 1
        elif combined_genes[x][15] < cut_off:
            length = 0
        else:
            "Broken ------------------> 10"
    if combined_genes[x][25] != 'x' and combined_genes[x][26] == 'x':
        i = 11
        if combined_genes[x][25] >= cut_off and combined_genes[x][24] >= cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 11
        elif combined_genes[x][25] < cut_off and combined_genes[x][24] >= cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 10
        elif combined_genes[x][24] < cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 9
        elif combined_genes[x][23] < cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 8
        elif combined_genes[x][22] < cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 7
        elif combined_genes[x][21] < cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 6
        elif combined_genes[x][20] < cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 5
        elif combined_genes[x][19] < cut_off and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 4
        elif combined_genes[x][18] < cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 3
        elif combined_genes[x][17] < cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 2
        elif combined_genes[x][16] < cut_off and combined_genes[x][15] >= cut_off:
            length = 1
        elif combined_genes[x][15] < cut_off:
            length = 0
        else:
            "Broken ------------------> 11"
    if combined_genes[x][26] != 'x' and combined_genes[x][27] == 'x':
        i = 12
        if combined_genes[x][26] >= cut_off and combined_genes[x][25] >= cut_off and combined_genes[x][24] >= cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 12
        elif combined_genes[x][26] < cut_off and combined_genes[x][25] >= cut_off and combined_genes[x][24] >= cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 11
        elif combined_genes[x][25] < cut_off and combined_genes[x][24] >= cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 10
        elif combined_genes[x][24] < cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 9
        elif combined_genes[x][23] < cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 8
        elif combined_genes[x][22] < cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 7
        elif combined_genes[x][21] < cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 6
        elif combined_genes[x][20] < cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 5
        elif combined_genes[x][19] < cut_off and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 4
        elif combined_genes[x][18] < cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 3
        elif combined_genes[x][17] < cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 2
        elif combined_genes[x][16] < cut_off and combined_genes[x][15] >= cut_off:
            length = 1
        elif combined_genes[x][15] < cut_off:
            length = 0
        else:
            "Broken ------------------> 12"
    if combined_genes[x][27] != 'x' and combined_genes[x][28] == 'x':
        i = 13
        if combined_genes[x][27] >= cut_off and combined_genes[x][26] >= cut_off and combined_genes[x][25] >= cut_off and combined_genes[x][24] >= cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 13
        elif combined_genes[x][27] < cut_off and combined_genes[x][26] >= cut_off and combined_genes[x][25] >= cut_off and combined_genes[x][24] >= cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 12
        elif combined_genes[x][26] < cut_off and combined_genes[x][25] >= cut_off and combined_genes[x][24] >= cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 11
        elif combined_genes[x][25] < cut_off and combined_genes[x][24] >= cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 10
        elif combined_genes[x][24] < cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 9
        elif combined_genes[x][23] < cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 8
        elif combined_genes[x][22] < cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 7
        elif combined_genes[x][21] < cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 6
        elif combined_genes[x][20] < cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 5
        elif combined_genes[x][19] < cut_off and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 4
        elif combined_genes[x][18] < cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 3
        elif combined_genes[x][17] < cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 2
        elif combined_genes[x][16] < cut_off and combined_genes[x][15] >= cut_off:
            length = 1
        elif combined_genes[x][15] < cut_off:
            length = 0
        else:
            "Broken ------------------> 13"
    if combined_genes[x][28] != 'x' and combined_genes[x][29] == 'x':
        i = 14
        if combined_genes[x][28] >= cut_off and combined_genes[x][27] >= cut_off and combined_genes[x][26] >= cut_off and combined_genes[x][25] >= cut_off and combined_genes[x][24] >= cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 14
        elif combined_genes[x][28] < cut_off and combined_genes[x][27] >= cut_off and combined_genes[x][26] >= cut_off and combined_genes[x][25] >= cut_off and combined_genes[x][24] >= cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 13
        elif combined_genes[x][27] < cut_off and combined_genes[x][26] >= cut_off and combined_genes[x][25] >= cut_off and combined_genes[x][24] >= cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 12
        elif combined_genes[x][26] < cut_off and combined_genes[x][25] >= cut_off and combined_genes[x][24] >= cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 11
        elif combined_genes[x][25] < cut_off and combined_genes[x][24] >= cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 10
        elif combined_genes[x][24] < cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 9
        elif combined_genes[x][23] < cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 8
        elif combined_genes[x][22] < cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 7
        elif combined_genes[x][21] < cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 6
        elif combined_genes[x][20] < cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 5
        elif combined_genes[x][19] < cut_off and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 4
        elif combined_genes[x][18] < cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 3
        elif combined_genes[x][17] < cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 2
        elif combined_genes[x][16] < cut_off and combined_genes[x][15] >= cut_off:
            length = 1
        elif combined_genes[x][15] < cut_off:
            length = 0
        else:
            "Broken ------------------> 14"
    if combined_genes[x][29] != 'x' and combined_genes[x][30] == 'x':
        i = 15
        if combined_genes[x][29] >= cut_off and combined_genes[x][28] >= cut_off and combined_genes[x][27] >= cut_off and combined_genes[x][26] >= cut_off and combined_genes[x][25] >= cut_off and combined_genes[x][24] >= cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 15
        elif combined_genes[x][29] < cut_off and combined_genes[x][28] >= cut_off and combined_genes[x][27] >= cut_off and combined_genes[x][26] >= cut_off and combined_genes[x][25] >= cut_off and combined_genes[x][24] >= cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 14
        elif combined_genes[x][28] < cut_off and combined_genes[x][27] >= cut_off and combined_genes[x][26] >= cut_off and combined_genes[x][25] >= cut_off and combined_genes[x][24] >= cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 13
        elif combined_genes[x][27] < cut_off and combined_genes[x][26] >= cut_off and combined_genes[x][25] >= cut_off and combined_genes[x][24] >= cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 12
        elif combined_genes[x][26] < cut_off and combined_genes[x][25] >= cut_off and combined_genes[x][24] >= cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 11
        elif combined_genes[x][25] < cut_off and combined_genes[x][24] >= cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 10
        elif combined_genes[x][24] < cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 9
        elif combined_genes[x][23] < cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 8
        elif combined_genes[x][22] < cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 7
        elif combined_genes[x][21] < cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 6
        elif combined_genes[x][20] < cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 5
        elif combined_genes[x][19] < cut_off and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 4
        elif combined_genes[x][18] < cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 3
        elif combined_genes[x][17] < cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 2
        elif combined_genes[x][16] < cut_off and combined_genes[x][15] >= cut_off:
            length = 1
        elif combined_genes[x][15] < cut_off:
            length = 0
        else:
            "Broken ------------------> 15"
    if combined_genes[x][30] != 'x' and combined_genes[x][31] == 'x':
        i = 16
        if combined_genes[x][30] >= cut_off and combined_genes[x][29] >= cut_off and combined_genes[x][28] >= cut_off and combined_genes[x][27] >= cut_off and combined_genes[x][26] >= cut_off and combined_genes[x][25] >= cut_off and combined_genes[x][24] >= cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 16
        elif combined_genes[x][30] < cut_off and combined_genes[x][29] >= cut_off and combined_genes[x][28] >= cut_off and combined_genes[x][27] >= cut_off and combined_genes[x][26] >= cut_off and combined_genes[x][25] >= cut_off and combined_genes[x][24] >= cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 15
        elif combined_genes[x][29] < cut_off and combined_genes[x][28] >= cut_off and combined_genes[x][27] >= cut_off and combined_genes[x][26] >= cut_off and combined_genes[x][25] >= cut_off and combined_genes[x][24] >= cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 14
        elif combined_genes[x][28] < cut_off and combined_genes[x][27] >= cut_off and combined_genes[x][26] >= cut_off and combined_genes[x][25] >= cut_off and combined_genes[x][24] >= cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 13
        elif combined_genes[x][27] < cut_off and combined_genes[x][26] >= cut_off and combined_genes[x][25] >= cut_off and combined_genes[x][24] >= cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 12
        elif combined_genes[x][26] < cut_off and combined_genes[x][25] >= cut_off and combined_genes[x][24] >= cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 11
        elif combined_genes[x][25] < cut_off and combined_genes[x][24] >= cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 10
        elif combined_genes[x][24] < cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 9
        elif combined_genes[x][23] < cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 8
        elif combined_genes[x][22] < cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 7
        elif combined_genes[x][21] < cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 6
        elif combined_genes[x][20] < cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 5
        elif combined_genes[x][19] < cut_off and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 4
        elif combined_genes[x][18] < cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 3
        elif combined_genes[x][17] < cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 2
        elif combined_genes[x][16] < cut_off and combined_genes[x][15] >= cut_off:
            length = 1
        elif combined_genes[x][15] < cut_off:
            length = 0
        else:
            "Broken ------------------> 16"
    if combined_genes[x][31] != 'x' and combined_genes[x][32] == 'x':
        i = 17
        if combined_genes[x][31] >= cut_off and combined_genes[x][30] >= cut_off and combined_genes[x][29] >= cut_off and combined_genes[x][28] >= cut_off and combined_genes[x][27] >= cut_off and combined_genes[x][26] >= cut_off and combined_genes[x][25] >= cut_off and combined_genes[x][24] >= cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 17
        elif combined_genes[x][31] < cut_off and combined_genes[x][30] >= cut_off and combined_genes[x][29] >= cut_off and combined_genes[x][28] >= cut_off and combined_genes[x][27] >= cut_off and combined_genes[x][26] >= cut_off and combined_genes[x][25] >= cut_off and combined_genes[x][24] >= cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 16
        elif combined_genes[x][30] < cut_off and  combined_genes[x][29] >= cut_off and combined_genes[x][28] >= cut_off and combined_genes[x][27] >= cut_off and combined_genes[x][26] >= cut_off and combined_genes[x][25] >= cut_off and combined_genes[x][24] >= cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 15
        elif combined_genes[x][29] < cut_off and combined_genes[x][28] >= cut_off and combined_genes[x][27] >= cut_off and combined_genes[x][26] >= cut_off and combined_genes[x][25] >= cut_off and combined_genes[x][24] >= cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 14
        elif combined_genes[x][28] < cut_off and combined_genes[x][27] >= cut_off and combined_genes[x][26] >= cut_off and combined_genes[x][25] >= cut_off and combined_genes[x][24] >= cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 13
        elif combined_genes[x][27] < cut_off and combined_genes[x][26] >= cut_off and combined_genes[x][25] >= cut_off and combined_genes[x][24] >= cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 12
        elif combined_genes[x][26] < cut_off and combined_genes[x][25] >= cut_off and combined_genes[x][24] >= cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 11
        elif combined_genes[x][25] < cut_off and combined_genes[x][24] >= cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 10
        elif combined_genes[x][24] < cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 9
        elif combined_genes[x][23] < cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 8
        elif combined_genes[x][22] < cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 7
        elif combined_genes[x][21] < cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 6
        elif combined_genes[x][20] < cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 5
        elif combined_genes[x][19] < cut_off and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 4
        elif combined_genes[x][18] < cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 3
        elif combined_genes[x][17] < cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 2
        elif combined_genes[x][16] < cut_off and combined_genes[x][15] >= cut_off:
            length = 1
        elif combined_genes[x][15] < cut_off:
            length = 0
        else:
            "Broken ------------------> 17"
    if combined_genes[x][32] != 'x' and combined_genes[x][33] == 'x':
        i = 18
        if combined_genes[x][32] >= cut_off and combined_genes[x][31] >= cut_off and combined_genes[x][30] >= cut_off and combined_genes[x][29] >= cut_off and combined_genes[x][28] >= cut_off and combined_genes[x][27] >= cut_off and combined_genes[x][26] >= cut_off and combined_genes[x][25] >= cut_off and combined_genes[x][24] >= cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 18
        elif combined_genes[x][32] < cut_off and combined_genes[x][31] >= cut_off and combined_genes[x][30] >= cut_off and combined_genes[x][29] >= cut_off and combined_genes[x][28] >= cut_off and combined_genes[x][27] >= cut_off and combined_genes[x][26] >= cut_off and combined_genes[x][25] >= cut_off and combined_genes[x][24] >= cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 17
        elif combined_genes[x][31] < cut_off and combined_genes[x][30] >= cut_off and combined_genes[x][29] >= cut_off and combined_genes[x][28] >= cut_off and combined_genes[x][27] >= cut_off and combined_genes[x][26] >= cut_off and combined_genes[x][25] >= cut_off and combined_genes[x][24] >= cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 16
        elif combined_genes[x][30] < cut_off and  combined_genes[x][29] >= cut_off and combined_genes[x][28] >= cut_off and combined_genes[x][27] >= cut_off and combined_genes[x][26] >= cut_off and combined_genes[x][25] >= cut_off and combined_genes[x][24] >= cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 15
        elif combined_genes[x][29] < cut_off and combined_genes[x][28] >= cut_off and combined_genes[x][27] >= cut_off and combined_genes[x][26] >= cut_off and combined_genes[x][25] >= cut_off and combined_genes[x][24] >= cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 14
        elif combined_genes[x][28] < cut_off and combined_genes[x][27] >= cut_off and combined_genes[x][26] >= cut_off and combined_genes[x][25] >= cut_off and combined_genes[x][24] >= cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 13
        elif combined_genes[x][27] < cut_off and combined_genes[x][26] >= cut_off and combined_genes[x][25] >= cut_off and combined_genes[x][24] >= cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 12
        elif combined_genes[x][26] < cut_off and combined_genes[x][25] >= cut_off and combined_genes[x][24] >= cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 11
        elif combined_genes[x][25] < cut_off and combined_genes[x][24] >= cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 10
        elif combined_genes[x][24] < cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 9
        elif combined_genes[x][23] < cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 8
        elif combined_genes[x][22] < cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 7
        elif combined_genes[x][21] < cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 6
        elif combined_genes[x][20] < cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 5
        elif combined_genes[x][19] < cut_off and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 4
        elif combined_genes[x][18] < cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 3
        elif combined_genes[x][17] < cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 2
        elif combined_genes[x][16] < cut_off and combined_genes[x][15] >= cut_off:
            length = 1
        elif combined_genes[x][15] < cut_off:
            length = 0
        else:
            "Broken ------------------> 18"
    if combined_genes[x][33] != 'x' and combined_genes[x][34] == 'x':
        i = 19
        if combined_genes[x][33] >= cut_off and combined_genes[x][32] >= cut_off and combined_genes[x][31] >= cut_off and combined_genes[x][30] >= cut_off and combined_genes[x][29] >= cut_off and combined_genes[x][28] >= cut_off and combined_genes[x][27] >= cut_off and combined_genes[x][26] >= cut_off and combined_genes[x][25] >= cut_off and combined_genes[x][24] >= cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 19
        elif combined_genes[x][33] < cut_off and combined_genes[x][32] >= cut_off and combined_genes[x][31] >= cut_off and combined_genes[x][30] >= cut_off and combined_genes[x][29] >= cut_off and combined_genes[x][28] >= cut_off and combined_genes[x][27] >= cut_off and combined_genes[x][26] >= cut_off and combined_genes[x][25] >= cut_off and combined_genes[x][24] >= cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 18
        elif combined_genes[x][32] < cut_off and combined_genes[x][31] >= cut_off and combined_genes[x][30] >= cut_off and combined_genes[x][29] >= cut_off and combined_genes[x][28] >= cut_off and combined_genes[x][27] >= cut_off and combined_genes[x][26] >= cut_off and combined_genes[x][25] >= cut_off and combined_genes[x][24] >= cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 17
        elif combined_genes[x][31] < cut_off and combined_genes[x][30] >= cut_off and combined_genes[x][29] >= cut_off and combined_genes[x][28] >= cut_off and combined_genes[x][27] >= cut_off and combined_genes[x][26] >= cut_off and combined_genes[x][25] >= cut_off and combined_genes[x][24] >= cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 16
        elif combined_genes[x][30] < cut_off and  combined_genes[x][29] >= cut_off and combined_genes[x][28] >= cut_off and combined_genes[x][27] >= cut_off and combined_genes[x][26] >= cut_off and combined_genes[x][25] >= cut_off and combined_genes[x][24] >= cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 15
        elif combined_genes[x][29] < cut_off and combined_genes[x][28] >= cut_off and combined_genes[x][27] >= cut_off and combined_genes[x][26] >= cut_off and combined_genes[x][25] >= cut_off and combined_genes[x][24] >= cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 14
        elif combined_genes[x][28] < cut_off and combined_genes[x][27] >= cut_off and combined_genes[x][26] >= cut_off and combined_genes[x][25] >= cut_off and combined_genes[x][24] >= cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 13
        elif combined_genes[x][27] < cut_off and combined_genes[x][26] >= cut_off and combined_genes[x][25] >= cut_off and combined_genes[x][24] >= cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 12
        elif combined_genes[x][26] < cut_off and combined_genes[x][25] >= cut_off and combined_genes[x][24] >= cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 11
        elif combined_genes[x][25] < cut_off and combined_genes[x][24] >= cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 10
        elif combined_genes[x][24] < cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 9
        elif combined_genes[x][23] < cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 8
        elif combined_genes[x][22] < cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 7
        elif combined_genes[x][21] < cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 6
        elif combined_genes[x][20] < cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 5
        elif combined_genes[x][19] < cut_off and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 4
        elif combined_genes[x][18] < cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 3
        elif combined_genes[x][17] < cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 2
        elif combined_genes[x][16] < cut_off and combined_genes[x][15] >= cut_off:
            length = 1
        elif combined_genes[x][15] < cut_off:
            length = 0
        else:
            "Broken ------------------> 19"
    if combined_genes[x][33] != 'x' and combined_genes[x][34] != 'x':
        i = 20
        if combined_genes[x][34] >= cut_off and combined_genes[x][33] >= cut_off and combined_genes[x][32] >= cut_off and combined_genes[x][31] >= cut_off and combined_genes[x][30] >= cut_off and combined_genes[x][29] >= cut_off and combined_genes[x][28] >= cut_off and combined_genes[x][27] >= cut_off and combined_genes[x][26] >= cut_off and combined_genes[x][25] >= cut_off and combined_genes[x][24] >= cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 20
        elif combined_genes[x][34] < cut_off and combined_genes[x][33] >= cut_off and combined_genes[x][32] >= cut_off and combined_genes[x][31] >= cut_off and combined_genes[x][30] >= cut_off and combined_genes[x][29] >= cut_off and combined_genes[x][28] >= cut_off and combined_genes[x][27] >= cut_off and combined_genes[x][26] >= cut_off and combined_genes[x][25] >= cut_off and combined_genes[x][24] >= cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 19
        elif combined_genes[x][33] < cut_off and combined_genes[x][32] >= cut_off and combined_genes[x][31] >= cut_off and combined_genes[x][30] >= cut_off and combined_genes[x][29] >= cut_off and combined_genes[x][28] >= cut_off and combined_genes[x][27] >= cut_off and combined_genes[x][26] >= cut_off and combined_genes[x][25] >= cut_off and combined_genes[x][24] >= cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 18
        elif combined_genes[x][32] < cut_off and combined_genes[x][31] >= cut_off and combined_genes[x][30] >= cut_off and combined_genes[x][29] >= cut_off and combined_genes[x][28] >= cut_off and combined_genes[x][27] >= cut_off and combined_genes[x][26] >= cut_off and combined_genes[x][25] >= cut_off and combined_genes[x][24] >= cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 17
        elif combined_genes[x][31] < cut_off and combined_genes[x][30] >= cut_off and combined_genes[x][29] >= cut_off and combined_genes[x][28] >= cut_off and combined_genes[x][27] >= cut_off and combined_genes[x][26] >= cut_off and combined_genes[x][25] >= cut_off and combined_genes[x][24] >= cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 16
        elif combined_genes[x][30] < cut_off and  combined_genes[x][29] >= cut_off and combined_genes[x][28] >= cut_off and combined_genes[x][27] >= cut_off and combined_genes[x][26] >= cut_off and combined_genes[x][25] >= cut_off and combined_genes[x][24] >= cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 15
        elif combined_genes[x][29] < cut_off and combined_genes[x][28] >= cut_off and combined_genes[x][27] >= cut_off and combined_genes[x][26] >= cut_off and combined_genes[x][25] >= cut_off and combined_genes[x][24] >= cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 14
        elif combined_genes[x][28] < cut_off and combined_genes[x][27] >= cut_off and combined_genes[x][26] >= cut_off and combined_genes[x][25] >= cut_off and combined_genes[x][24] >= cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 13
        elif combined_genes[x][27] < cut_off and combined_genes[x][26] >= cut_off and combined_genes[x][25] >= cut_off and combined_genes[x][24] >= cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 12
        elif combined_genes[x][26] < cut_off and combined_genes[x][25] >= cut_off and combined_genes[x][24] >= cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 11
        elif combined_genes[x][25] < cut_off and combined_genes[x][24] >= cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 10
        elif combined_genes[x][24] < cut_off and combined_genes[x][23] >= cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 9
        elif combined_genes[x][23] < cut_off and combined_genes[x][22] >= cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 8
        elif combined_genes[x][22] < cut_off and combined_genes[x][21] >= cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 7
        elif combined_genes[x][21] < cut_off and combined_genes[x][20] >= cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 6
        elif combined_genes[x][20] < cut_off and combined_genes[x][19] >= cut_off  and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 5
        elif combined_genes[x][19] < cut_off and combined_genes[x][18] >= cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 4
        elif combined_genes[x][18] < cut_off and combined_genes[x][17] >= cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 3
        elif combined_genes[x][17] < cut_off and combined_genes[x][16] >= cut_off and combined_genes[x][15] >= cut_off:
            length = 2
        elif combined_genes[x][16] < cut_off and combined_genes[x][15] >= cut_off:
            length = 1
        elif combined_genes[x][15] < cut_off:
            length = 0
        else:
            "Broken ------------------> 20"


    z = length #KEEP need to set z to length
    # starting to apply lag corrections and regression
    if combined_genes[x][13] == 0 : # this is left over from previous code, any gene longer than length 1nt should not trigger this if statment
        print("gene has no length (1 nt or less)")
        print(combined_genes[x][0])
        break
    if len(time_points_set) >= 2 and combined_genes[x][13] > 0 and 0 == (elong_rate*time_points[1]): #number based on 1 minute lag at elong_rate nt/sec
        if length >= 3 and length <= 20:
            add_4 = [] # create list
            add_4 = combined_genes[x][0:15] # add all info
            y = 15 + z # correct time point slicing
            add_4.append(combined_genes[x][15]) #moving the 0 point to where we think it should be based on TSS to AUG
            add_4.extend(combined_genes[x][16:y]) # losing first datapoint so change from -7 above to -6
            data_log_corrected = [] #dummy until do complete loop, then can index in with [x][-z:]
            # data_log_corrected.append(add_4) #dummy until do complete loop, then can index in with [x][-z:]
            data_log_corrected = add_4 #dummy until do complete loop, then can index in with [x][-z:]
            add_5 = [] # create list
            add_5 = combined_genes[x][0:15] # add gene info
            add_5.extend(time_points[0:z]) # changed from 0 above to 1, because losing first datapoint
            time_corrected = [] #dummy until do complete loop, then can index in with [x][-z:]
            # time_corrected.append(add_5) #dummy until do complete loop, then can index in with [x][-z:]
            time_corrected = add_5 #dummy until do complete loop, then can index in with [x][-z:]
            a = combined_genes[x][1]
            [b, c, d, e, f] = (stats.linregress(time_corrected[-z:], data_log_corrected[-z:])) #-z, to get all time points
            csv_results_table.append([a,b,c,d,np.log(2)/-b,length,0,combined_genes[x][(y-1)],combined_genes[x][(6)],combined_genes[x][(0)],combined_genes[x][(12)],i])
            if d > flag_r_value:
                csv_results_table[-1].append("R_value_FLAG")
            else:
                csv_results_table[-1].extend("x")
            if (np.log(2)/-b) > 2*time_points[z-1] or np.log(2)/-b < 0:
                csv_results_table[-1].append("Half_Life_Flag")
            else:
                csv_results_table[-1].extend("x")
            csv_results_table[-1].extend(combined_genes[x][15:15+len(time_points_set)])
            # print(a,',',b,',',c,',',d,',',np.log(2)/-b,',',length,',',-1,',',combined_genes[x][(y-1)],',',combined_genes[x][(6)],',',combined_genes[x][(0)],',',combined_genes[x][(12)],',',combined_genes[x][(-7)],',',combined_genes[x][(-6)],',',combined_genes[x][(-5)],',',combined_genes[x][(-4)],',',combined_genes[x][(-3)],',',combined_genes[x][(-2)],',',i)
            operon_count_total += 1
        elif length <= 2 and length >= 0:
            add_2 = []
            add_2 = combined_genes[x][0:15]
            y = 15 + z # correct time point slicing
            zero = [0]
            add_2.extend(zero)
            data_log_corrected = [] #dummy until do complete loop, then can index in with [x][-z:]
            # data_log_corrected.append(add_4) #dummy until do complete loop, then can index in with [x][-z:]
            data_log_corrected = add_2 #dummy until do complete loop, then can index in with [x][-z:]
            add_3 = []
            add_3 = combined_genes[x][0:15]
            add_3.extend(zero)
            time_corrected = [] #dummy until do complete loop, then can index in with [x][-z:]
            # time_corrected.append(add_5) #dummy until do complete loop, then can index in with [x][-z:]
            time_corrected = add_3 #dummy until do complete loop, then can index in with [x][-z:]
            operon_count_total += 1
            csv_results_table.append([add_2[1],'x','x','x','x',length,0,combined_genes[x][(y-1)],combined_genes[x][(6)],combined_genes[x][(0)],combined_genes[x][(12)],i,'x','x'])
            csv_results_table[-1].extend(combined_genes[x][15:15+len(time_points_set)])
            # print(add_2[1],',',add_2[-1],',','x',',','x',',','x',',',length,',',0,',',combined_genes[x][(y-1)],',',combined_genes[x][(6)],',',combined_genes[x][(0)],',',combined_genes[x][(12)],',',combined_genes[x][(-7)],',',combined_genes[x][(-6)],',',combined_genes[x][(-5)],',',combined_genes[x][(-4)],',',combined_genes[x][(-3)],',',combined_genes[x][(-2)],',',i)
        else:
            print(combined_genes[x][1], 2, z, length)
            print(combined_genes[x])
            break
        # if length >= 3 and length <= 6:
        #     add_4 = [] # create list
        #     add_4 = combined_genes[x][0:15] # add all info
        #     y = -7 + z # correct time point slicing
        #     add_4.extend(combined_genes[x][-7:y]) # add gene info
        #     data_log_corrected = [] #dummy until do complete loop, then can index in with [x][-z:]
        #     # data_log_corrected.append(add_4) #dummy until do complete loop, then can index in with [x][-z:]
        #     data_log_corrected = add_4 #dummy until do complete loop, then can index in with [x][-z:]
        #     add_5 = [] # create list
        #     add_5 = combined_genes[x][0:15] # add gene info
        #     add_5.extend(time_points[:z]) #running 0,1,3,9,...8,15
        #     time_corrected = [] #dummy until do complete loop, then can index in with [x][-z:]
        #     # time_corrected.append(add_5) #dummy until do complete loop, then can index in with [x][-z:]
        #     time_corrected = add_5 #dummy until do complete loop, then can index in with [x][-z:]
        #     a = combined_genes[x][1]
        #     [b, c, d, e, f] = (stats.linregress(time_corrected[-z:], data_log_corrected[-z:]))
        #     print a,',',b,',',c,',',d,',',np.log(2)/-b,',',length,',',0,',',combined_genes[x][(y-1)],',',combined_genes[x][(6)],',',combined_genes[x][(0)],',',combined_genes[x][(12)],',',combined_genes[x][(-7)],',',combined_genes[x][(-6)],',',combined_genes[x][(-5)],',',combined_genes[x][(-4)],',',combined_genes[x][(-3)],',',combined_genes[x][(-2)],',',i
        #     operon_count_total += 1
        # elif length <= 2 and length >= 0:
        #     add_2 = []
        #     add_2 = combined_genes[x][0:15]
        #     y = -7 + z # correct time point slicing
        #     zero = [0]
        #     add_2.extend(zero)
        #     data_log_corrected = [] #dummy until do complete loop, then can index in with [x][-z:]
        #     # data_log_corrected.append(add_4) #dummy until do complete loop, then can index in with [x][-z:]
        #     data_log_corrected = add_2 #dummy until do complete loop, then can index in with [x][-z:]
        #     add_3 = []
        #     add_3 = combined_genes[x][0:15]
        #     add_3.extend(zero)
        #     time_corrected = [] #dummy until do complete loop, then can index in with [x][-z:]
        #     # time_corrected.append(add_5) #dummy until do complete loop, then can index in with [x][-z:]
        #     time_corrected = add_3 #dummy until do complete loop, then can index in with [x][-z:]
        #     operon_count_total += 1
        #     print ""
        #     print
        #     print add_2[1],',',add_2[-1],',','x',',','x',',','x',',',length,',',0,',',combined_genes[x][(y-1)],',',combined_genes[x][(6)],',',combined_genes[x][(0)],',',combined_genes[x][(12)],',',combined_genes[x][(-7)],',',combined_genes[x][(-6)],',',combined_genes[x][(-5)],',',combined_genes[x][(-4)],',',combined_genes[x][(-3)],',',combined_genes[x][(-2)],',',i
        # else:
        #     print combined_genes[x][1], 1, combined_genes[x][-7]
        #     break
    if len(time_points_set) >= 2 and combined_genes[x][13] > 0 and combined_genes[x][13] < (elong_rate*time_points[1]) and (elong_rate*time_points[1]) != 0: #number based on 1 minute lag at elong_rate nt/sec
        if length >= 3 and length <= 20:
            add_4 = [] # create list
            add_4 = combined_genes[x][0:15] # add all info
            y = 15 + z # correct time point slicing
            add_4.append(combined_genes[x][15]) #moving the 0 point to where we think it should be based on TSS to AUG
            add_4.extend(combined_genes[x][16:y]) # losing first datapoint so change from -7 above to -6
            data_log_corrected = [] #dummy until do complete loop, then can index in with [x][-z:]
            # data_log_corrected.append(add_4) #dummy until do complete loop, then can index in with [x][-z:]
            data_log_corrected = add_4 #dummy until do complete loop, then can index in with [x][-z:]
            add_5 = [] # create list
            add_5 = combined_genes[x][0:15] # add gene info
            add_5.append(combined_genes[x][13]/(elong_rate)) #moving the 0 point to where we think it should be based on TSS to AUG
            add_5.extend(time_points[1:z]) # changed from 0 above to 1, because losing first datapoint
            time_corrected = [] #dummy until do complete loop, then can index in with [x][-z:]
            # time_corrected.append(add_5) #dummy until do complete loop, then can index in with [x][-z:]
            time_corrected = add_5 #dummy until do complete loop, then can index in with [x][-z:]
            a = combined_genes[x][1]
            [b, c, d, e, f] = (stats.linregress(time_corrected[-z:], data_log_corrected[-z:])) #-z, to get all time points
            csv_results_table.append([a,b,c,d,np.log(2)/-b,length,-1,combined_genes[x][(y-1)],combined_genes[x][(6)],combined_genes[x][(0)],combined_genes[x][(12)],i])
            if d > flag_r_value:
                csv_results_table[-1].append("R_value_FLAG")
            else:
                csv_results_table[-1].extend("x")
            if (np.log(2)/-b) > 2*time_points[z-1] or np.log(2)/-b < 0:
                csv_results_table[-1].append("Half_Life_Flag")
            else:
                csv_results_table[-1].extend("x")
            csv_results_table[-1].extend(combined_genes[x][15:15+len(time_points_set)])
            #print(a,',',b,',',c,',',d,',',np.log(2)/-b,',',length,',',-1,',',combined_genes[x][(y-1)],',',combined_genes[x][(6)],',',combined_genes[x][(0)],',',combined_genes[x][(12)],',',combined_genes[x][(-7)],',',combined_genes[x][(-6)],',',combined_genes[x][(-5)],',',combined_genes[x][(-4)],',',combined_genes[x][(-3)],',',combined_genes[x][(-2)],',',i)
            operon_count_total += 1
        elif length <= 2 and length >= 0:
            add_2 = []
            add_2 = combined_genes[x][0:15]
            y = 15 + z # correct time point slicing
            zero = [0]
            add_2.extend(zero)
            data_log_corrected = [] #dummy until do complete loop, then can index in with [x][-z:]
            # data_log_corrected.append(add_4) #dummy until do complete loop, then can index in with [x][-z:]
            data_log_corrected = add_2 #dummy until do complete loop, then can index in with [x][-z:]
            add_3 = []
            add_3 = combined_genes[x][0:15]
            add_3.extend(zero)
            time_corrected = [] #dummy until do complete loop, then can index in with [x][-z:]
            # time_corrected.append(add_5) #dummy until do complete loop, then can index in with [x][-z:]
            time_corrected = add_3 #dummy until do complete loop, then can index in with [x][-z:]
            operon_count_total += 1
            csv_results_table.append([add_2[1],'x','x','x','x',length,0,combined_genes[x][(y-1)],combined_genes[x][(6)],combined_genes[x][(0)],combined_genes[x][(12)],i,'x','x'])
            csv_results_table[-1].extend(combined_genes[x][15:15+len(time_points_set)])
            #print(add_2[1],',',add_2[-1],',','x',',','x',',','x',',',length,',',0,',',combined_genes[x][(y-1)],',',combined_genes[x][(6)],',',combined_genes[x][(0)],',',combined_genes[x][(12)],',',combined_genes[x][(-7)],',',combined_genes[x][(-6)],',',combined_genes[x][(-5)],',',combined_genes[x][(-4)],',',combined_genes[x][(-3)],',',combined_genes[x][(-2)],',',i)
        else:
            print(combined_genes[x][1], 2, z, length)
            print(combined_genes[x])
            break

    if len(time_points_set) >= 3 and combined_genes[x][13] >= (elong_rate*time_points[1]) and combined_genes[x][13] < (elong_rate*time_points[2]) and (elong_rate*time_points[1]) != 0:
        if length >= 4 and length <= 20:
            add_4 = [] # create list
            add_4 = combined_genes[x][0:15] # add all info
            y = 15 + z # correct time point slicing
            add_4.append((combined_genes[x][15]+combined_genes[x][16])/2) #avg. 0 and 1 point to move to where we think it should be based on TSS to AUG
            add_4.extend(combined_genes[x][17:y]) # losing two datapoint so change from -7 above to -5
            data_log_corrected = [] #dummy until do complete loop, then can index in with [x][-z:]
            # data_log_corrected.append(add_4) #dummy until do complete loop, then can index in with [x][-z:]
            data_log_corrected = add_4 #dummy until do complete loop, then can index in with [x][-z:]
            add_5 = [] # create list
            add_5 = combined_genes[x][0:15] # add gene info
            add_5.append(combined_genes[x][13]/(elong_rate)) #moving the avg. point to where we think it should be based on TSS to AUG
            add_5.extend(time_points[2:z]) # changed from 0 above to 2, because losing two datapoint
            time_corrected = [] #dummy until do complete loop, then can index in with [x][-z:]
            # time_corrected.append(add_5) #dummy until do complete loop, then can index in with [x][-z:]
            time_corrected = add_5 #dummy until do complete loop, then can index in with [x][-z:]
            a = combined_genes[x][1]
            [b, c, d, e, f] = (stats.linregress(time_corrected[-z+1:], data_log_corrected[-z+1:])) #-z+1, to index correctly from length
            csv_results_table.append([a,b,c,d,np.log(2)/-b,length,-2,combined_genes[x][(y-1)],combined_genes[x][(6)],combined_genes[x][(0)],combined_genes[x][(12)],i])
            if d > flag_r_value:
                csv_results_table[-1].append("R_value_FLAG")
            else:
                csv_results_table[-1].extend("x")
            if (np.log(2)/-b) > 2*time_points[z-1] or np.log(2)/-b < 0:
                csv_results_table[-1].append("Half_Life_Flag")
            else:
                csv_results_table[-1].extend("x")
            csv_results_table[-1].extend(combined_genes[x][15:15+len(time_points_set)])
            #print(a,',',b,',',c,',',d,',',np.log(2)/-b,',',length,',',-2,',',combined_genes[x][(y-1)],',',combined_genes[x][(6)],',',combined_genes[x][(0)],',',combined_genes[x][(12)],',',combined_genes[x][(-7)],',',combined_genes[x][(-6)],',',combined_genes[x][(-5)],',',combined_genes[x][(-4)],',',combined_genes[x][(-3)],',',combined_genes[x][(-2)],',',i)
        elif length <= 3 and length >= 0:
            add_2 = []
            add_2 = combined_genes[x][0:15]
            y = 15 + z # correct time point slicing
            zero = [0]
            add_2.extend(zero)
            data_log_corrected = [] #dummy until do complete loop, then can index in with [x][-z:]
            # data_log_corrected.append(add_4) #dummy until do complete loop, then can index in with [x][-z:]
            data_log_corrected = add_2 #dummy until do complete loop, then can index in with [x][-z:]
            add_3 = []
            add_3 = combined_genes[x][0:15]
            add_3.extend(zero)
            time_corrected = [] #dummy until do complete loop, then can index in with [x][-z:]
            # time_corrected.append(add_5) #dummy until do complete loop, then can index in with [x][-z:]
            time_corrected = add_3 #dummy until do complete loop, then can index in with [x][-z:]
            operon_count_total += 1
            csv_results_table.append([add_2[1],'x','x','x','x',length,0,combined_genes[x][(y-1)],combined_genes[x][(6)],combined_genes[x][(0)],combined_genes[x][(12)],i,'x','x'])
            csv_results_table[-1].extend(combined_genes[x][15:15+len(time_points_set)])
            #print(add_2[1],',',add_2[-1],',','x',',','x',',','x',',',length,',',0,',',combined_genes[x][(y-1)],',',combined_genes[x][(6)],',',combined_genes[x][(0)],',',combined_genes[x][(12)],',',combined_genes[x][(-7)],',',combined_genes[x][(-6)],',',combined_genes[x][(-5)],',',combined_genes[x][(-4)],',',combined_genes[x][(-3)],',',combined_genes[x][(-2)],',',i)
        else:
            print(combined_genes[x][1], 2)
            break
    if len(time_points_set) >= 4 and combined_genes[x][13] >= (elong_rate*time_points[2]) and combined_genes[x][13] < (elong_rate*time_points[3]) and (elong_rate*time_points[1]) != 0:
        if length >= 5 and length <= 20:
            add_4 = [] # create list
            add_4 = combined_genes[x][0:15] # add all info
            y = 15 + z # correct time point slicing
            add_4.append(sum(combined_genes[x][15:18])/3) #avg. 0, 1, 2 point to move to where we think it should be based on TSS to AUG
            add_4.extend(combined_genes[x][18:y]) # losing three datapoint so change from -7 above to -4
            data_log_corrected = [] #dummy until do complete loop, then can index in with [x][-z:]
            # data_log_corrected.append(add_4) #dummy until do complete loop, then can index in with [x][-z:]
            data_log_corrected = add_4 #dummy until do complete loop, then can index in with [x][-z:]
            add_5 = [] # create list
            add_5 = combined_genes[x][0:15] # add gene info
            add_5.append(combined_genes[x][13]/(elong_rate)) #moving the avg. point to where we think it should be based on TSS to AUG
            add_5.extend(time_points[3:z]) # changed from 0 above to 3, because losing three datapoint
            time_corrected = [] #dummy until do complete loop, then can index in with [x][-z:]
            # time_corrected.append(add_5) #dummy until do complete loop, then can index in with [x][-z:]
            time_corrected = add_5 #dummy until do complete loop, then can index in with [x][-z:]
            a = combined_genes[x][1]
            [b, c, d, e, f] = (stats.linregress(time_corrected[-z+2:], data_log_corrected[-z+2:])) #-z+2, to index correctly from length: missing 1st three time points
            csv_results_table.append([a,b,c,d,np.log(2)/-b,length,-3,combined_genes[x][(y-1)],combined_genes[x][(6)],combined_genes[x][(0)],combined_genes[x][(12)],i])
            if d > flag_r_value:
                csv_results_table[-1].append("R_value_FLAG")
            else:
                csv_results_table[-1].extend("x")
            if (np.log(2)/-b) > 2*time_points[z-1] or np.log(2)/-b < 0:
                csv_results_table[-1].append("Half_Life_Flag")
            else:
                csv_results_table[-1].extend("x")
            # csv_results_table[-1].extend(data_log_corrected[-z:])
            csv_results_table[-1].extend(combined_genes[x][15:15+len(time_points_set)])
            #print(a,',',b,',',c,',',d,',',np.log(2)/-b,',',length,',',-3,',',combined_genes[x][(y-1)],',',combined_genes[x][(6)],',',combined_genes[x][(0)],',',combined_genes[x][(12)],',',combined_genes[x][(-7)],',',combined_genes[x][(-6)],',',combined_genes[x][(-5)],',',combined_genes[x][(-4)],',',combined_genes[x][(-3)],',',combined_genes[x][(-2)],',',i)
        elif length <= 4 and length >= 0:
            add_2 = []
            add_2 = combined_genes[x][0:15]
            y = 15 + z # correct time point slicing
            zero = [0]
            add_2.extend(zero)
            data_log_corrected = [] #dummy until do complete loop, then can index in with [x][-z:]
            # data_log_corrected.append(add_4) #dummy until do complete loop, then can index in with [x][-z:]
            data_log_corrected = add_2 #dummy until do complete loop, then can index in with [x][-z:]
            add_3 = []
            add_3 = combined_genes[x][0:15]
            add_3.extend(zero)
            time_corrected = [] #dummy until do complete loop, then can index in with [x][-z:]
            # time_corrected.append(add_5) #dummy until do complete loop, then can index in with [x][-z:]
            time_corrected = add_3 #dummy until do complete loop, then can index in with [x][-z:]
            operon_count_total += 1
            csv_results_table.append([add_2[1],'x','x','x','x',length,0,combined_genes[x][(y-1)],combined_genes[x][(6)],combined_genes[x][(0)],combined_genes[x][(12)],i,'x','x'])
            csv_results_table[-1].extend(combined_genes[x][15:15+len(time_points_set)])
            #print(add_2[1],',',add_2[-1],',','x',',','x',',','x',',',length,',',0,',',combined_genes[x][(y-1)],',',combined_genes[x][(6)],',',combined_genes[x][(0)],',',combined_genes[x][(12)],',',combined_genes[x][(-7)],',',combined_genes[x][(-6)],',',combined_genes[x][(-5)],',',combined_genes[x][(-4)],',',combined_genes[x][(-3)],',',combined_genes[x][(-2)],',',i)
        else:
            print(combined_genes[x][1], 3)
            break

    if len(time_points_set) >= 5 and combined_genes[x][13] >= (elong_rate*time_points[3]) and combined_genes[x][13] < (elong_rate*time_points[4]) and (elong_rate*time_points[1]) != 0:
        if length >= 6 and length <= 20:
            add_4 = [] # create list
            add_4 = combined_genes[x][0:15] # add all info
            y = 15 + z # correct time point slicing
            add_4.append(sum(combined_genes[x][15:19])/4) #avg. 0, 1, 2, 4 point to move to where we think it should be based on TSS to AUG
            add_4.extend(combined_genes[x][19:y]) # losing three datapoint so change from -7 above to -4
            data_log_corrected = [] #dummy until do complete loop, then can index in with [x][-z:]
            # data_log_corrected.append(add_4) #dummy until do complete loop, then can index in with [x][-z:]
            data_log_corrected = add_4 #dummy until do complete loop, then can index in with [x][-z:]
            add_5 = [] # create list
            add_5 = combined_genes[x][0:15] # add gene info
            add_5.append(combined_genes[x][13]/(elong_rate)) #moving the avg. point to where we think it should be based on TSS to AUG
            add_5.extend(time_points[4:z]) # changed from 0 above to 4, because losing four datapoint
            time_corrected = [] #dummy until do complete loop, then can index in with [x][-z:]
            # time_corrected.append(add_5) #dummy until do complete loop, then can index in with [x][-z:]
            time_corrected = add_5 #dummy until do complete loop, then can index in with [x][-z:]
            a = combined_genes[x][1]
            [b, c, d, e, f] = (stats.linregress(time_corrected[-z+3:], data_log_corrected[-z+3:])) #-z+2, to drop the first two points
            csv_results_table.append([a,b,c,d,np.log(2)/-b,length,-4,combined_genes[x][(y-1)],combined_genes[x][(6)],combined_genes[x][(0)],combined_genes[x][(12)],i])
            if d > flag_r_value:
                csv_results_table[-1].append("R_value_FLAG")
            else:
                csv_results_table[-1].extend("x")
            if (np.log(2)/-b) > 2*time_points[z-1] or np.log(2)/-b < 0:
                csv_results_table[-1].append("Half_Life_Flag")
            else:
                csv_results_table[-1].extend("x")
            # csv_results_table[-1].extend(data_log_corrected[-z:])
            csv_results_table[-1].extend(combined_genes[x][15:15+len(time_points_set)])
            #print(a,',',b,',',c,',',d,',',np.log(2)/-b,',',length,',',-4,',',combined_genes[x][(y-1)],',',combined_genes[x][(6)],',',combined_genes[x][(0)],',',combined_genes[x][(12)],',',combined_genes[x][(-7)],',',combined_genes[x][(-6)],',',combined_genes[x][(-5)],',',combined_genes[x][(-4)],',',combined_genes[x][(-3)],',',combined_genes[x][(-2)],',',i)
        elif length <= 5 and length >= 0:
            add_2 = []
            add_2 = combined_genes[x][0:15]
            y = 15 + z  # correct time point slicing
            zero = [0]
            add_2.extend(zero)
            data_log_corrected = [] #dummy until do complete loop, then can index in with [x][-z:]
            # data_log_corrected.append(add_4) #dummy until do complete loop, then can index in with [x][-z:]
            data_log_corrected = add_2 #dummy until do complete loop, then can index in with [x][-z:]
            add_3 = []
            add_3 = combined_genes[x][0:15]
            add_3.extend(zero)
            time_corrected = [] #dummy until do complete loop, then can index in with [x][-z:]
            # time_corrected.append(add_5) #dummy until do complete loop, then can index in with [x][-z:]
            time_corrected = add_3 #dummy until do complete loop, then can index in with [x][-z:]
            operon_count_total += 1
            csv_results_table.append([add_2[1],'x','x','x','x',length,0,combined_genes[x][(y-1)],combined_genes[x][(6)],combined_genes[x][(0)],combined_genes[x][(12)],i,'x','x'])
            csv_results_table[-1].extend(combined_genes[x][15:15+len(time_points_set)])
            #print(add_2[1],',',add_2[-1],',','x',',','x',',','x',',',length,',',0,',',combined_genes[x][(y-1)],',',combined_genes[x][(6)],',',combined_genes[x][(0)],',',combined_genes[x][(12)],',',combined_genes[x][(-7)],',',combined_genes[x][(-6)],',',combined_genes[x][(-5)],',',combined_genes[x][(-4)],',',combined_genes[x][(-3)],',',combined_genes[x][(-2)],',',i)
        else:
            print(combined_genes[x][1], 3)
            break

    if len(time_points_set) >= 6 and combined_genes[x][13] >= (elong_rate*time_points[4]) and combined_genes[x][13] < (elong_rate*time_points[5]) and (elong_rate*time_points[1]) != 0:
        if length >= 7 and length <= 20:
            add_4 = [] # create list
            add_4 = combined_genes[x][0:15] # add all info
            y = 15 + z # correct time point slicing
            add_4.append(sum(combined_genes[x][15:20])/5) #avg. 0, 1, 2, 4 point to move to where we think it should be based on TSS to AUG
            add_4.extend(combined_genes[x][20:y]) # losing three datapoint so change from -7 above to -4
            data_log_corrected = [] #dummy until do complete loop, then can index in with [x][-z:]
            # data_log_corrected.append(add_4) #dummy until do complete loop, then can index in with [x][-z:]
            data_log_corrected = add_4 #dummy until do complete loop, then can index in with [x][-z:]
            add_5 = [] # create list
            add_5 = combined_genes[x][0:15] # add gene info
            add_5.append(combined_genes[x][13]/(elong_rate)) #moving the avg. point to where we think it should be based on TSS to AUG
            add_5.extend(time_points[5:z]) # changed from 0 above to 5, because losing four datapoint
            time_corrected = [] #dummy until do complete loop, then can index in with [x][-z:]
            # time_corrected.append(add_5) #dummy until do complete loop, then can index in with [x][-z:]
            time_corrected = add_5 #dummy until do complete loop, then can index in with [x][-z:]
            a = combined_genes[x][1]
            [b, c, d, e, f] = (stats.linregress(time_corrected[-z+4:], data_log_corrected[-z+4:])) #-z+2, to drop the first two points
            csv_results_table.append([a,b,c,d,np.log(2)/-b,length,-5,combined_genes[x][(y-1)],combined_genes[x][(6)],combined_genes[x][(0)],combined_genes[x][(12)],i])
            if d > flag_r_value:
                csv_results_table[-1].append("R_value_FLAG")
            else:
                csv_results_table[-1].extend("x")
            if (np.log(2)/-b) > 2*time_points[z-1] or np.log(2)/-b < 0:
                csv_results_table[-1].append("Half_Life_Flag")
            else:
                csv_results_table[-1].extend("x")
            # csv_results_table[-1].extend(data_log_corrected[-z:])
            csv_results_table[-1].extend(combined_genes[x][15:15+len(time_points_set)])
            #print(a,',',b,',',c,',',d,',',np.log(2)/-b,',',length,',',-5,',',combined_genes[x][(y-1)],',',combined_genes[x][(6)],',',combined_genes[x][(0)],',',combined_genes[x][(12)],',',combined_genes[x][(-7)],',',combined_genes[x][(-6)],',',combined_genes[x][(-5)],',',combined_genes[x][(-4)],',',combined_genes[x][(-3)],',',combined_genes[x][(-2)],',',i)

        elif length <= 6 and length >= 0:
            add_2 = []
            add_2 = combined_genes[x][0:15]
            y = 15 + z  # correct time point slicing
            zero = [0]
            add_2.extend(zero)
            data_log_corrected = [] #dummy until do complete loop, then can index in with [x][-z:]
            # data_log_corrected.append(add_4) #dummy until do complete loop, then can index in with [x][-z:]
            data_log_corrected = add_2 #dummy until do complete loop, then can index in with [x][-z:]
            add_3 = []
            add_3 = combined_genes[x][0:15]
            add_3.extend(zero)
            time_corrected = [] #dummy until do complete loop, then can index in with [x][-z:]
            # time_corrected.append(add_5) #dummy until do complete loop, then can index in with [x][-z:]
            time_corrected = add_3 #dummy until do complete loop, then can index in with [x][-z:]
            operon_count_total += 1
            csv_results_table.append([add_2[1],'x','x','x','x',length,0,combined_genes[x][(y-1)],combined_genes[x][(6)],combined_genes[x][(0)],combined_genes[x][(12)],i,'x','x'])
            csv_results_table[-1].extend(combined_genes[x][15:15+len(time_points_set)])
            #print(add_2[1],',',add_2[-1],',','x',',','x',',','x',',',length,',',0,',',combined_genes[x][(y-1)],',',combined_genes[x][(6)],',',combined_genes[x][(0)],',',combined_genes[x][(12)],',',combined_genes[x][(-7)],',',combined_genes[x][(-6)],',',combined_genes[x][(-5)],',',combined_genes[x][(-4)],',',combined_genes[x][(-3)],',',combined_genes[x][(-2)],',',i)
        else:
            print(combined_genes[x][1], 3)
            break

    if len(time_points_set) >= 7 and combined_genes[x][13] >= (elong_rate*time_points[5]) and combined_genes[x][13] < (elong_rate*time_points[6]) and (elong_rate*time_points[1]) != 0: # change 2 numbers
        if length >= 8 and length <= 20: # incease 1 number by +1
            add_4 = [] # create list
            add_4 = combined_genes[x][0:15] # add all info
            y = 15 + z # correct time point slicing
            add_4.append(sum(combined_genes[x][15:21])/6) #increase 2 number by +1
            add_4.extend(combined_genes[x][21:y]) #increase 1 number by +1
            data_log_corrected = []
            data_log_corrected = add_4
            add_5 = []
            add_5 = combined_genes[x][0:15]
            add_5.append(combined_genes[x][13]/(elong_rate))
            add_5.extend(time_points[6:z]) #increase 1 number by +1
            time_corrected = []
            time_corrected = add_5
            a = combined_genes[x][1]
            [b, c, d, e, f] = (stats.linregress(time_corrected[-z+5:], data_log_corrected[-z+5:])) #increase 2 number by +1, and reduce 1 below by -1
            csv_results_table.append([a,b,c,d,np.log(2)/-b,length,-6,combined_genes[x][(y-1)],combined_genes[x][(6)],combined_genes[x][(0)],combined_genes[x][(12)],i])
            if d > flag_r_value:
                csv_results_table[-1].append("R_value_FLAG")
            else:
                csv_results_table[-1].extend("x")
            if (np.log(2)/-b) > 2*time_points[z-1] or np.log(2)/-b < 0:
                csv_results_table[-1].append("Half_Life_Flag")
            else:
                csv_results_table[-1].extend("x")
            # csv_results_table[-1].extend(data_log_corrected[-z:])
            csv_results_table[-1].extend(combined_genes[x][15:15+len(time_points_set)])
            #print(a,',',b,',',c,',',d,',',np.log(2)/-b,',',length,',',-6,',',combined_genes[x][(y-1)],',',combined_genes[x][(6)],',',combined_genes[x][(0)],',',combined_genes[x][(12)],',',combined_genes[x][(-7)],',',combined_genes[x][(-6)],',',combined_genes[x][(-5)],',',combined_genes[x][(-4)],',',combined_genes[x][(-3)],',',combined_genes[x][(-2)],',',i)
        elif length <= 7 and length >= 0: #increase 1 number by +1
            add_2 = []
            add_2 = combined_genes[x][0:15]
            y = 15 + z
            zero = [0]
            add_2.extend(zero)
            data_log_corrected = []
            data_log_corrected = add_2
            add_3 = []
            add_3 = combined_genes[x][0:15]
            add_3.extend(zero)
            time_corrected = []
            time_corrected = add_3
            operon_count_total += 1
            csv_results_table.append([add_2[1],'x','x','x','x',length,0,combined_genes[x][(y-1)],combined_genes[x][(6)],combined_genes[x][(0)],combined_genes[x][(12)],i,'x','x'])
            csv_results_table[-1].extend(combined_genes[x][15:15+len(time_points_set)])
            #print(add_2[1],',',add_2[-1],',','x',',','x',',','x',',',length,',',0,',',combined_genes[x][(y-1)],',',combined_genes[x][(6)],',',combined_genes[x][(0)],',',combined_genes[x][(12)],',',combined_genes[x][(-7)],',',combined_genes[x][(-6)],',',combined_genes[x][(-5)],',',combined_genes[x][(-4)],',',combined_genes[x][(-3)],',',combined_genes[x][(-2)],',',i)
        else:
            print(combined_genes[x][1], 3)
            break

    if len(time_points_set) >= 8 and combined_genes[x][13] >= (elong_rate*time_points[6]) and combined_genes[x][13] < (elong_rate*time_points[7]) and (elong_rate*time_points[1]) != 0: # change 2 numbers
        if length >= 9 and length <= 20: # incease 1 number by +1
            add_4 = []
            add_4 = combined_genes[x][0:15]
            y = 15 + z
            add_4.append(sum(combined_genes[x][15:22])/7) #increase 2 number by +1
            add_4.extend(combined_genes[x][22:y]) #increase 1 number by +1
            data_log_corrected = []
            data_log_corrected = add_4
            add_5 = []
            add_5 = combined_genes[x][0:15]
            add_5.append(combined_genes[x][13]/(elong_rate))
            add_5.extend(time_points[7:z]) #increase 1 number by +1
            time_corrected = []
            time_corrected = add_5
            a = combined_genes[x][1]
            [b, c, d, e, f] = (stats.linregress(time_corrected[-z+6:], data_log_corrected[-z+6:])) #increase 2 number by +1, and reduce 1 below by -1
            csv_results_table.append([a,b,c,d,np.log(2)/-b,length,-7,combined_genes[x][(y-1)],combined_genes[x][(6)],combined_genes[x][(0)],combined_genes[x][(12)],i])
            if d > flag_r_value:
                csv_results_table[-1].append("R_value_FLAG")
            else:
                csv_results_table[-1].extend("x")
            if (np.log(2)/-b) > 2*time_points[z-1] or np.log(2)/-b < 0:
                csv_results_table[-1].append("Half_Life_Flag")
            else:
                csv_results_table[-1].extend("x")
            # csv_results_table[-1].extend(data_log_corrected[-z:])
            csv_results_table[-1].extend(combined_genes[x][15:15+len(time_points_set)])
            #print(a,',',b,',',c,',',d,',',np.log(2)/-b,',',length,',',-7,',',combined_genes[x][(y-1)],',',combined_genes[x][(6)],',',combined_genes[x][(0)],',',combined_genes[x][(12)],',',combined_genes[x][(-7)],',',combined_genes[x][(-6)],',',combined_genes[x][(-5)],',',combined_genes[x][(-4)],',',combined_genes[x][(-3)],',',combined_genes[x][(-2)],',',i)
        elif length <= 8 and length >= 0: #increase 1 number by +1
            add_2 = []
            add_2 = combined_genes[x][0:15]
            y = 15 + z
            zero = [0]
            add_2.extend(zero)
            data_log_corrected = []
            data_log_corrected = add_2
            add_3 = []
            add_3 = combined_genes[x][0:15]
            add_3.extend(zero)
            time_corrected = []
            time_corrected = add_3
            operon_count_total += 1
            csv_results_table.append([add_2[1],'x','x','x','x',length,0,combined_genes[x][(y-1)],combined_genes[x][(6)],combined_genes[x][(0)],combined_genes[x][(12)],i,'x','x'])
            csv_results_table[-1].extend(combined_genes[x][15:15+len(time_points_set)])
            #print(add_2[1],',',add_2[-1],',','x',',','x',',','x',',',length,',',0,',',combined_genes[x][(y-1)],',',combined_genes[x][(6)],',',combined_genes[x][(0)],',',combined_genes[x][(12)],',',combined_genes[x][(-7)],',',combined_genes[x][(-6)],',',combined_genes[x][(-5)],',',combined_genes[x][(-4)],',',combined_genes[x][(-3)],',',combined_genes[x][(-2)],',',i)
        else:
            print(combined_genes[x][1], 3)
            break

    if len(time_points_set) >= 9 and combined_genes[x][13] >= (elong_rate*time_points[7]) and combined_genes[x][13] < (elong_rate*time_points[8]) and (elong_rate*time_points[1]) != 0: # change 2 numbers
        if length >= 10 and length <= 20: # incease 1 number by +1
            add_4 = []
            add_4 = combined_genes[x][0:15]
            y = 15 + z
            add_4.append(sum(combined_genes[x][15:23])/8) #increase 2 number by +1
            add_4.extend(combined_genes[x][23:y]) #increase 1 number by +1
            data_log_corrected = []
            data_log_corrected = add_4
            add_5 = []
            add_5 = combined_genes[x][0:15]
            add_5.append(combined_genes[x][13]/(elong_rate))
            add_5.extend(time_points[8:z]) #increase 1 number by +1
            time_corrected = []
            time_corrected = add_5
            a = combined_genes[x][1]
            [b, c, d, e, f] = (stats.linregress(time_corrected[-z+7:], data_log_corrected[-z+7:])) #increase 2 number by +1, and reduce 1 below by -1
            csv_results_table.append([a,b,c,d,np.log(2)/-b,length,-8,combined_genes[x][(y-1)],combined_genes[x][(6)],combined_genes[x][(0)],combined_genes[x][(12)],i])
            if d > flag_r_value:
                csv_results_table[-1].append("R_value_FLAG")
            else:
                csv_results_table[-1].extend("x")
            if (np.log(2)/-b) > 2*time_points[z-1] or np.log(2)/-b < 0:
                csv_results_table[-1].append("Half_Life_Flag")
            else:
                csv_results_table[-1].extend("x")
            # csv_results_table[-1].extend(data_log_corrected[-z:])
            csv_results_table[-1].extend(combined_genes[x][15:15+len(time_points_set)])
            #print(a,',',b,',',c,',',d,',',np.log(2)/-b,',',length,',',-8,',',combined_genes[x][(y-1)],',',combined_genes[x][(6)],',',combined_genes[x][(0)],',',combined_genes[x][(12)],',',combined_genes[x][(-7)],',',combined_genes[x][(-6)],',',combined_genes[x][(-5)],',',combined_genes[x][(-4)],',',combined_genes[x][(-3)],',',combined_genes[x][(-2)],',',i)
        elif length <= 9 and length >= 0: #increase 1 number by +1
            add_2 = []
            add_2 = combined_genes[x][0:15]
            y = 15 + z
            zero = [0]
            add_2.extend(zero)
            data_log_corrected = []
            data_log_corrected = add_2
            add_3 = []
            add_3 = combined_genes[x][0:15]
            add_3.extend(zero)
            time_corrected = []
            time_corrected = add_3
            operon_count_total += 1
            csv_results_table.append([add_2[1],'x','x','x','x',length,0,combined_genes[x][(y-1)],combined_genes[x][(6)],combined_genes[x][(0)],combined_genes[x][(12)],i,'x','x'])
            csv_results_table[-1].extend(combined_genes[x][15:15+len(time_points_set)])
            #print(add_2[1],',',add_2[-1],',','x',',','x',',','x',',',length,',',0,',',combined_genes[x][(y-1)],',',combined_genes[x][(6)],',',combined_genes[x][(0)],',',combined_genes[x][(12)],',',combined_genes[x][(-7)],',',combined_genes[x][(-6)],',',combined_genes[x][(-5)],',',combined_genes[x][(-4)],',',combined_genes[x][(-3)],',',combined_genes[x][(-2)],',',i)
        else:
            print(combined_genes[x][1], 3)
            break

    if len(time_points_set) >= 10 and combined_genes[x][13] >= (elong_rate*time_points[8]) and combined_genes[x][13] < (elong_rate*time_points[9]) and (elong_rate*time_points[1]) != 0: #increase 2 number by +1
        if length >= 11 and length <= 20: # incease 1 number by +1
            add_4 = []
            add_4 = combined_genes[x][0:15]
            y = 15 + z
            add_4.append(sum(combined_genes[x][15:24])/9) #increase 2 number by +1
            add_4.extend(combined_genes[x][24:y]) #increase 1 number by +1
            data_log_corrected = []
            data_log_corrected = add_4
            add_5 = []
            add_5 = combined_genes[x][0:15]
            add_5.append(combined_genes[x][13]/(elong_rate))
            add_5.extend(time_points[9:z]) #increase 1 number by +1
            time_corrected = []
            time_corrected = add_5
            a = combined_genes[x][1]
            [b, c, d, e, f] = (stats.linregress(time_corrected[-z+8:], data_log_corrected[-z+8:])) #increase 2 number by +1, and reduce 1 below by -1
            csv_results_table.append([a,b,c,d,np.log(2)/-b,length,-9,combined_genes[x][(y-1)],combined_genes[x][(6)],combined_genes[x][(0)],combined_genes[x][(12)],i])
            if d > flag_r_value:
                csv_results_table[-1].append("R_value_FLAG")
            else:
                csv_results_table[-1].extend("x")
            if (np.log(2)/-b) > 2*time_points[z-1] or np.log(2)/-b < 0:
                csv_results_table[-1].append("Half_Life_Flag")
            else:
                csv_results_table[-1].extend("x")
            # csv_results_table[-1].extend(data_log_corrected[-z:])
            csv_results_table[-1].extend(combined_genes[x][15:15+len(time_points_set)])
            #print(a,',',b,',',c,',',d,',',np.log(2)/-b,',',length,',',-9,',',combined_genes[x][(y-1)],',',combined_genes[x][(6)],',',combined_genes[x][(0)],',',combined_genes[x][(12)],',',combined_genes[x][(-7)],',',combined_genes[x][(-6)],',',combined_genes[x][(-5)],',',combined_genes[x][(-4)],',',combined_genes[x][(-3)],',',combined_genes[x][(-2)],',',i)
        elif length <= 10 and length >= 0: #increase 1 number by +1
            add_2 = []
            add_2 = combined_genes[x][0:15]
            y = 15 + z
            zero = [0]
            add_2.extend(zero)
            data_log_corrected = []
            data_log_corrected = add_2
            add_3 = []
            add_3 = combined_genes[x][0:15]
            add_3.extend(zero)
            time_corrected = []
            time_corrected = add_3
            operon_count_total += 1
            csv_results_table.append([add_2[1],'x','x','x','x',length,0,combined_genes[x][(y-1)],combined_genes[x][(6)],combined_genes[x][(0)],combined_genes[x][(12)],i,'x','x'])
            csv_results_table[-1].extend(combined_genes[x][15:15+len(time_points_set)])
            #print(add_2[1],',',add_2[-1],',','x',',','x',',','x',',',length,',',0,',',combined_genes[x][(y-1)],',',combined_genes[x][(6)],',',combined_genes[x][(0)],',',combined_genes[x][(12)],',',combined_genes[x][(-7)],',',combined_genes[x][(-6)],',',combined_genes[x][(-5)],',',combined_genes[x][(-4)],',',combined_genes[x][(-3)],',',combined_genes[x][(-2)],',',i)
        else:
            print(combined_genes[x][1], 9)
            break

    if len(time_points_set) >= 11 and combined_genes[x][13] >= (elong_rate*time_points[9]) and combined_genes[x][13] < (elong_rate*time_points[10]) and (elong_rate*time_points[1]) != 0: #increase 2 number by +1
        if length >= 12 and length <= 20: # incease 1 number by +1
            add_4 = []
            add_4 = combined_genes[x][0:15]
            y = 15 + z
            add_4.append(sum(combined_genes[x][15:25])/10) #increase 2 number by +1
            add_4.extend(combined_genes[x][25:y]) #increase 1 number by +1
            data_log_corrected = []
            data_log_corrected = add_4
            add_5 = []
            add_5 = combined_genes[x][0:15]
            add_5.append(combined_genes[x][13]/(elong_rate))
            add_5.extend(time_points[10:z]) #increase 1 number by +1
            time_corrected = []
            time_corrected = add_5
            a = combined_genes[x][1]
            [b, c, d, e, f] = (stats.linregress(time_corrected[-z+9:], data_log_corrected[-z+9:])) #increase 2 number by +1, and reduce 1 below by -1
            csv_results_table.append([a,b,c,d,np.log(2)/-b,length,-10,combined_genes[x][(y-1)],combined_genes[x][(6)],combined_genes[x][(0)],combined_genes[x][(12)],i])
            if d > flag_r_value:
                csv_results_table[-1].append("R_value_FLAG")
            else:
                csv_results_table[-1].extend("x")
            if (np.log(2)/-b) > 2*time_points[z-1] or np.log(2)/-b < 0:
                csv_results_table[-1].append("Half_Life_Flag")
            else:
                csv_results_table[-1].extend("x")
            # csv_results_table[-1].extend(data_log_corrected[-z:])
            csv_results_table[-1].extend(combined_genes[x][15:15+len(time_points_set)])
            #print(a,',',b,',',c,',',d,',',np.log(2)/-b,',',length,',',-10,',',combined_genes[x][(y-1)],',',combined_genes[x][(6)],',',combined_genes[x][(0)],',',combined_genes[x][(12)],',',combined_genes[x][(-7)],',',combined_genes[x][(-6)],',',combined_genes[x][(-5)],',',combined_genes[x][(-4)],',',combined_genes[x][(-3)],',',combined_genes[x][(-2)],',',i)
        elif length <= 11 and length >= 0: #increase 1 number by +1
            add_2 = []
            add_2 = combined_genes[x][0:15]
            y = 15 + z
            zero = [0]
            add_2.extend(zero)
            data_log_corrected = []
            data_log_corrected = add_2
            add_3 = []
            add_3 = combined_genes[x][0:15]
            add_3.extend(zero)
            time_corrected = []
            time_corrected = add_3
            operon_count_total += 1
            csv_results_table.append([add_2[1],'x','x','x','x',length,0,combined_genes[x][(y-1)],combined_genes[x][(6)],combined_genes[x][(0)],combined_genes[x][(12)],i,'x','x'])
            csv_results_table[-1].extend(combined_genes[x][15:15+len(time_points_set)])
            #print(add_2[1],',',add_2[-1],',','x',',','x',',','x',',',length,',',0,',',combined_genes[x][(y-1)],',',combined_genes[x][(6)],',',combined_genes[x][(0)],',',combined_genes[x][(12)],',',combined_genes[x][(-7)],',',combined_genes[x][(-6)],',',combined_genes[x][(-5)],',',combined_genes[x][(-4)],',',combined_genes[x][(-3)],',',combined_genes[x][(-2)],',',i)
        else:
            print(combined_genes[x][1], 10)
            break

    if len(time_points_set) >= 12 and combined_genes[x][13] >= (elong_rate*time_points[10]) and combined_genes[x][13] < (elong_rate*time_points[11]) and (elong_rate*time_points[1]) != 0: #increase 2 number by +1
        if length >= 13 and length <= 20: # incease 1 number by +1
            add_4 = []
            add_4 = combined_genes[x][0:15]
            y = 15 + z
            add_4.append(sum(combined_genes[x][15:26])/11) #increase 2 number by +1
            add_4.extend(combined_genes[x][26:y]) #increase 1 number by +1
            data_log_corrected = []
            data_log_corrected = add_4
            add_5 = []
            add_5 = combined_genes[x][0:15]
            add_5.append(combined_genes[x][13]/(elong_rate))
            add_5.extend(time_points[11:z]) #increase 1 number by +1
            time_corrected = []
            time_corrected = add_5
            a = combined_genes[x][1]
            [b, c, d, e, f] = (stats.linregress(time_corrected[-z+10:], data_log_corrected[-z+10:])) #increase 2 number by +1, and reduce 1 below by -1
            csv_results_table.append([a,b,c,d,np.log(2)/-b,length,-11,combined_genes[x][(y-1)],combined_genes[x][(6)],combined_genes[x][(0)],combined_genes[x][(12)],i])
            if d > flag_r_value:
                csv_results_table[-1].append("R_value_FLAG")
            else:
                csv_results_table[-1].extend("x")
            if (np.log(2)/-b) > 2*time_points[z-1] or np.log(2)/-b < 0:
                csv_results_table[-1].append("Half_Life_Flag")
            else:
                csv_results_table[-1].extend("x")
            # csv_results_table[-1].extend(data_log_corrected[-z:])
            csv_results_table[-1].extend(combined_genes[x][15:15+len(time_points_set)])
            #print(a,',',b,',',c,',',d,',',np.log(2)/-b,',',length,',',-11,',',combined_genes[x][(y-1)],',',combined_genes[x][(6)],',',combined_genes[x][(0)],',',combined_genes[x][(12)],',',combined_genes[x][(-7)],',',combined_genes[x][(-6)],',',combined_genes[x][(-5)],',',combined_genes[x][(-4)],',',combined_genes[x][(-3)],',',combined_genes[x][(-2)],',',i)
        elif length <= 12 and length >= 0: #increase 1 number by +1
            add_2 = []
            add_2 = combined_genes[x][0:15]
            y = 15 + z
            zero = [0]
            add_2.extend(zero)
            data_log_corrected = []
            data_log_corrected = add_2
            add_3 = []
            add_3 = combined_genes[x][0:15]
            add_3.extend(zero)
            time_corrected = []
            time_corrected = add_3
            operon_count_total += 1
            csv_results_table.append([add_2[1],'x','x','x','x',length,0,combined_genes[x][(y-1)],combined_genes[x][(6)],combined_genes[x][(0)],combined_genes[x][(12)],i,'x','x'])
            csv_results_table[-1].extend(combined_genes[x][15:15+len(time_points_set)])
            #print(add_2[1],',',add_2[-1],',','x',',','x',',','x',',',length,',',0,',',combined_genes[x][(y-1)],',',combined_genes[x][(6)],',',combined_genes[x][(0)],',',combined_genes[x][(12)],',',combined_genes[x][(-7)],',',combined_genes[x][(-6)],',',combined_genes[x][(-5)],',',combined_genes[x][(-4)],',',combined_genes[x][(-3)],',',combined_genes[x][(-2)],',',i)
        else:
            print(combined_genes[x][1], 11)
            break
    if len(time_points_set) >= 13 and combined_genes[x][13] >= (elong_rate*time_points[11]) and combined_genes[x][13] < (elong_rate*time_points[12]) and (elong_rate*time_points[1]) != 0: #increase 2 number by +1
        if length >= 14 and length <= 20: # incease 1 number by +1
            add_4 = []
            add_4 = combined_genes[x][0:15]
            y = 15 + z
            add_4.append(sum(combined_genes[x][15:27])/12) #increase 2 number by +1
            add_4.extend(combined_genes[x][27:y]) #increase 1 number by +1
            data_log_corrected = []
            data_log_corrected = add_4
            add_5 = []
            add_5 = combined_genes[x][0:15]
            add_5.append(combined_genes[x][13]/(elong_rate))
            add_5.extend(time_points[12:z]) #increase 1 number by +1
            time_corrected = []
            time_corrected = add_5
            a = combined_genes[x][1]
            [b, c, d, e, f] = (stats.linregress(time_corrected[-z+11:], data_log_corrected[-z+11:])) #increase 2 number by +1, and reduce 1 below by -1
            csv_results_table.append([a,b,c,d,np.log(2)/-b,length,-12,combined_genes[x][(y-1)],combined_genes[x][(6)],combined_genes[x][(0)],combined_genes[x][(12)],i])
            if d > flag_r_value:
                csv_results_table[-1].append("R_value_FLAG")
            else:
                csv_results_table[-1].extend("x")
            if (np.log(2)/-b) > 2*time_points[z-1] or np.log(2)/-b < 0:
                csv_results_table[-1].append("Half_Life_Flag")
            else:
                csv_results_table[-1].extend("x")
            # csv_results_table[-1].extend(data_log_corrected[-z:])
            csv_results_table[-1].extend(combined_genes[x][15:15+len(time_points_set)])
            #print(a,',',b,',',c,',',d,',',np.log(2)/-b,',',length,',',-12,',',combined_genes[x][(y-1)],',',combined_genes[x][(6)],',',combined_genes[x][(0)],',',combined_genes[x][(12)],',',combined_genes[x][(-7)],',',combined_genes[x][(-6)],',',combined_genes[x][(-5)],',',combined_genes[x][(-4)],',',combined_genes[x][(-3)],',',combined_genes[x][(-2)],',',i)
        elif length <= 13 and length >= 0: #increase 1 number by +1
            add_2 = []
            add_2 = combined_genes[x][0:15]
            y = 15 + z
            zero = [0]
            add_2.extend(zero)
            data_log_corrected = []
            data_log_corrected = add_2
            add_3 = []
            add_3 = combined_genes[x][0:15]
            add_3.extend(zero)
            time_corrected = []
            time_corrected = add_3
            operon_count_total += 1
            csv_results_table.append([add_2[1],'x','x','x','x',length,0,combined_genes[x][(y-1)],combined_genes[x][(6)],combined_genes[x][(0)],combined_genes[x][(12)],i,'x','x'])
            csv_results_table[-1].extend(combined_genes[x][15:15+len(time_points_set)])
            #print(add_2[1],',',add_2[-1],',','x',',','x',',','x',',',length,',',0,',',combined_genes[x][(y-1)],',',combined_genes[x][(6)],',',combined_genes[x][(0)],',',combined_genes[x][(12)],',',combined_genes[x][(-7)],',',combined_genes[x][(-6)],',',combined_genes[x][(-5)],',',combined_genes[x][(-4)],',',combined_genes[x][(-3)],',',combined_genes[x][(-2)],',',i)
        else:
            print(combined_genes[x][1], 12)
            break

    if len(time_points_set) >= 14 and combined_genes[x][13] >= (elong_rate*time_points[12]) and combined_genes[x][13] < (elong_rate*time_points[13]) and (elong_rate*time_points[1]) != 0: #increase 2 number by +1
        if length >= 15 and length <= 20: # incease 1 number by +1
            add_4 = []
            add_4 = combined_genes[x][0:15]
            y = 15 + z
            add_4.append(sum(combined_genes[x][15:28])/13) #increase 2 number by +1
            add_4.extend(combined_genes[x][28:y]) #increase 1 number by +1
            data_log_corrected = []
            data_log_corrected = add_4
            add_5 = []
            add_5 = combined_genes[x][0:15]
            add_5.append(combined_genes[x][13]/(elong_rate))
            add_5.extend(time_points[13:z]) #increase 1 number by +1
            time_corrected = []
            time_corrected = add_5
            a = combined_genes[x][1]
            [b, c, d, e, f] = (stats.linregress(time_corrected[-z+12:], data_log_corrected[-z+12:])) #increase 2 number by +1, and reduce 1 below by -1
            csv_results_table.append([a,b,c,d,np.log(2)/-b,length,-13,combined_genes[x][(y-1)],combined_genes[x][(6)],combined_genes[x][(0)],combined_genes[x][(12)],i])
            if d > flag_r_value:
                csv_results_table[-1].append("R_value_FLAG")
            else:
                csv_results_table[-1].extend("x")
            if (np.log(2)/-b) > 2*time_points[z-1] or np.log(2)/-b < 0:
                csv_results_table[-1].append("Half_Life_Flag")
            else:
                csv_results_table[-1].extend("x")
            # csv_results_table[-1].extend(data_log_corrected[-z:])
            csv_results_table[-1].extend(combined_genes[x][15:15+len(time_points_set)])
            #print(a,',',b,',',c,',',d,',',np.log(2)/-b,',',length,',',-13,',',combined_genes[x][(y-1)],',',combined_genes[x][(6)],',',combined_genes[x][(0)],',',combined_genes[x][(12)],',',combined_genes[x][(-7)],',',combined_genes[x][(-6)],',',combined_genes[x][(-5)],',',combined_genes[x][(-4)],',',combined_genes[x][(-3)],',',combined_genes[x][(-2)],',',i)
        elif length <= 14 and length >= 0: #increase 1 number by +1
            add_2 = []
            add_2 = combined_genes[x][0:15]
            y = 15 + z
            zero = [0]
            add_2.extend(zero)
            data_log_corrected = []
            data_log_corrected = add_2
            add_3 = []
            add_3 = combined_genes[x][0:15]
            add_3.extend(zero)
            time_corrected = []
            time_corrected = add_3
            operon_count_total += 1
            csv_results_table.append([add_2[1],'x','x','x','x',length,0,combined_genes[x][(y-1)],combined_genes[x][(6)],combined_genes[x][(0)],combined_genes[x][(12)],i,'x','x'])
            csv_results_table[-1].extend(combined_genes[x][15:15+len(time_points_set)])
            #print(add_2[1],',',add_2[-1],',','x',',','x',',','x',',',length,',',0,',',combined_genes[x][(y-1)],',',combined_genes[x][(6)],',',combined_genes[x][(0)],',',combined_genes[x][(12)],',',combined_genes[x][(-7)],',',combined_genes[x][(-6)],',',combined_genes[x][(-5)],',',combined_genes[x][(-4)],',',combined_genes[x][(-3)],',',combined_genes[x][(-2)],',',i)
        else:
            print(combined_genes[x][1], 13)
            break

    if len(time_points_set) >= 15 and combined_genes[x][13] >= (elong_rate*time_points[13]) and combined_genes[x][13] < (elong_rate*time_points[14]) and (elong_rate*time_points[1]) != 0: #increase 2 number by +1
        if length >= 16 and length <= 20: # incease 1 number by +1
            add_4 = []
            add_4 = combined_genes[x][0:15]
            y = 15 + z
            add_4.append(sum(combined_genes[x][15:29])/14) #increase 2 number by +1
            add_4.extend(combined_genes[x][29:y]) #increase 1 number by +1
            data_log_corrected = []
            data_log_corrected = add_4
            add_5 = []
            add_5 = combined_genes[x][0:15]
            add_5.append(combined_genes[x][13]/(elong_rate))
            add_5.extend(time_points[14:z]) #increase 1 number by +1
            time_corrected = []
            time_corrected = add_5
            a = combined_genes[x][1]
            [b, c, d, e, f] = (stats.linregress(time_corrected[-z+13:], data_log_corrected[-z+13:])) #increase 2 number by +1, and reduce 1 below by -1
            csv_results_table.append([a,b,c,d,np.log(2)/-b,length,-14,combined_genes[x][(y-1)],combined_genes[x][(6)],combined_genes[x][(0)],combined_genes[x][(12)],i])
            if d > flag_r_value:
                csv_results_table[-1].append("R_value_FLAG")
            else:
                csv_results_table[-1].extend("x")
            if (np.log(2)/-b) > 2*time_points[z-1] or np.log(2)/-b < 0:
                csv_results_table[-1].append("Half_Life_Flag")
            else:
                csv_results_table[-1].extend("x")
            # csv_results_table[-1].extend(data_log_corrected[-z:])
            csv_results_table[-1].extend(combined_genes[x][15:15+len(time_points_set)])
            #print(a,',',b,',',c,',',d,',',np.log(2)/-b,',',length,',',-14,',',combined_genes[x][(y-1)],',',combined_genes[x][(6)],',',combined_genes[x][(0)],',',combined_genes[x][(12)],',',combined_genes[x][(-7)],',',combined_genes[x][(-6)],',',combined_genes[x][(-5)],',',combined_genes[x][(-4)],',',combined_genes[x][(-3)],',',combined_genes[x][(-2)],',',i)
        elif length <= 15 and length >= 0: #increase 1 number by +1
            add_2 = []
            add_2 = combined_genes[x][0:15]
            y = 15 + z
            zero = [0]
            add_2.extend(zero)
            data_log_corrected = []
            data_log_corrected = add_2
            add_3 = []
            add_3 = combined_genes[x][0:15]
            add_3.extend(zero)
            time_corrected = []
            time_corrected = add_3
            operon_count_total += 1
            csv_results_table.append([add_2[1],'x','x','x','x',length,0,combined_genes[x][(y-1)],combined_genes[x][(6)],combined_genes[x][(0)],combined_genes[x][(12)],i,'x','x'])
            csv_results_table[-1].extend(combined_genes[x][15:15+len(time_points_set)])
            #print(add_2[1],',',add_2[-1],',','x',',','x',',','x',',',length,',',0,',',combined_genes[x][(y-1)],',',combined_genes[x][(6)],',',combined_genes[x][(0)],',',combined_genes[x][(12)],',',combined_genes[x][(-7)],',',combined_genes[x][(-6)],',',combined_genes[x][(-5)],',',combined_genes[x][(-4)],',',combined_genes[x][(-3)],',',combined_genes[x][(-2)],',',i)
        else:
            print(combined_genes[x][1], 14)
            break

    if len(time_points_set) >= 16 and combined_genes[x][13] >= (elong_rate*time_points[14]) and combined_genes[x][13] < (elong_rate*time_points[15]) and (elong_rate*time_points[1]) != 0: #increase 2 number by +1
        if length >= 17 and length <= 20: # incease 1 number by +1
            add_4 = []
            add_4 = combined_genes[x][0:15]
            y = 15 + z
            add_4.append(sum(combined_genes[x][15:30])/15) #increase 2 number by +1
            add_4.extend(combined_genes[x][30:y]) #increase 1 number by +1
            data_log_corrected = []
            data_log_corrected = add_4
            add_5 = []
            add_5 = combined_genes[x][0:15]
            add_5.append(combined_genes[x][13]/(elong_rate))
            add_5.extend(time_points[15:z]) #increase 1 number by +1
            time_corrected = []
            time_corrected = add_5
            a = combined_genes[x][1]
            [b, c, d, e, f] = (stats.linregress(time_corrected[-z+14:], data_log_corrected[-z+14:])) #increase 2 number by +1, and reduce 1 below by -1
            csv_results_table.append([a,b,c,d,np.log(2)/-b,length,-15,combined_genes[x][(y-1)],combined_genes[x][(6)],combined_genes[x][(0)],combined_genes[x][(12)],i])
            if d > flag_r_value:
                csv_results_table[-1].append("R_value_FLAG")
            else:
                csv_results_table[-1].extend("x")
            if (np.log(2)/-b) > 2*time_points[z-1] or np.log(2)/-b < 0:
                csv_results_table[-1].append("Half_Life_Flag")
            else:
                csv_results_table[-1].extend("x")
            # csv_results_table[-1].extend(data_log_corrected[-z:])
            csv_results_table[-1].extend(combined_genes[x][15:15+len(time_points_set)])
            #print(a,',',b,',',c,',',d,',',np.log(2)/-b,',',length,',',-15,',',combined_genes[x][(y-1)],',',combined_genes[x][(6)],',',combined_genes[x][(0)],',',combined_genes[x][(12)],',',combined_genes[x][(-7)],',',combined_genes[x][(-6)],',',combined_genes[x][(-5)],',',combined_genes[x][(-4)],',',combined_genes[x][(-3)],',',combined_genes[x][(-2)],',',i)
        elif length <= 16 and length >= 0: #increase 1 number by +1
            add_2 = []
            add_2 = combined_genes[x][0:15]
            y = 15 + z
            zero = [0]
            add_2.extend(zero)
            data_log_corrected = []
            data_log_corrected = add_2
            add_3 = []
            add_3 = combined_genes[x][0:15]
            add_3.extend(zero)
            time_corrected = []
            time_corrected = add_3
            operon_count_total += 1
            csv_results_table.append([add_2[1],'x','x','x','x',length,0,combined_genes[x][(y-1)],combined_genes[x][(6)],combined_genes[x][(0)],combined_genes[x][(12)],i,'x','x'])
            csv_results_table[-1].extend(combined_genes[x][15:15+len(time_points_set)])
            #print(add_2[1],',',add_2[-1],',','x',',','x',',','x',',',length,',',0,',',combined_genes[x][(y-1)],',',combined_genes[x][(6)],',',combined_genes[x][(0)],',',combined_genes[x][(12)],',',combined_genes[x][(-7)],',',combined_genes[x][(-6)],',',combined_genes[x][(-5)],',',combined_genes[x][(-4)],',',combined_genes[x][(-3)],',',combined_genes[x][(-2)],',',i)
        else:
            print(combined_genes[x][1], 15)
            break

    if len(time_points_set) >= 17 and combined_genes[x][13] >= (elong_rate*time_points[15]) and combined_genes[x][13] < (elong_rate*time_points[16]) and (elong_rate*time_points[1]) != 0: #increase 2 number by +1
        if length >= 18 and length <= 20: # incease 1 number by +1
            add_4 = []
            add_4 = combined_genes[x][0:15]
            y = 15 + z
            add_4.append(sum(combined_genes[x][15:31])/16) #increase 2 number by +1
            add_4.extend(combined_genes[x][31:y]) #increase 1 number by +1
            data_log_corrected = []
            data_log_corrected = add_4
            add_5 = []
            add_5 = combined_genes[x][0:15]
            add_5.append(combined_genes[x][13]/(elong_rate))
            add_5.extend(time_points[16:z]) #increase 1 number by +1
            time_corrected = []
            time_corrected = add_5
            a = combined_genes[x][1]
            [b, c, d, e, f] = (stats.linregress(time_corrected[-z+15:], data_log_corrected[-z+15:])) #increase 2 number by +1, and reduce 1 below by -1
            csv_results_table.append([a,b,c,d,np.log(2)/-b,length,-16,combined_genes[x][(y-1)],combined_genes[x][(6)],combined_genes[x][(0)],combined_genes[x][(12)],i])
            if d > flag_r_value:
                csv_results_table[-1].append("R_value_FLAG")
            else:
                csv_results_table[-1].extend("x")
            if (np.log(2)/-b) > 2*time_points[z-1] or np.log(2)/-b < 0:
                csv_results_table[-1].append("Half_Life_Flag")
            else:
                csv_results_table[-1].extend("x")
            # csv_results_table[-1].extend(data_log_corrected[-z:])
            csv_results_table[-1].extend(combined_genes[x][15:15+len(time_points_set)])
            #print(a,',',b,',',c,',',d,',',np.log(2)/-b,',',length,',',-16,',',combined_genes[x][(y-1)],',',combined_genes[x][(6)],',',combined_genes[x][(0)],',',combined_genes[x][(12)],',',combined_genes[x][(-7)],',',combined_genes[x][(-6)],',',combined_genes[x][(-5)],',',combined_genes[x][(-4)],',',combined_genes[x][(-3)],',',combined_genes[x][(-2)],',',i)
        elif length <= 17 and length >= 0: #increase 1 number by +1
            add_2 = []
            add_2 = combined_genes[x][0:15]
            y = 15 + z
            zero = [0]
            add_2.extend(zero)
            data_log_corrected = []
            data_log_corrected = add_2
            add_3 = []
            add_3 = combined_genes[x][0:15]
            add_3.extend(zero)
            time_corrected = []
            time_corrected = add_3
            operon_count_total += 1
            csv_results_table.append([add_2[1],'x','x','x','x',length,0,combined_genes[x][(y-1)],combined_genes[x][(6)],combined_genes[x][(0)],combined_genes[x][(12)],i,'x','x'])
            csv_results_table[-1].extend(combined_genes[x][15:15+len(time_points_set)])
            #print(add_2[1],',',add_2[-1],',','x',',','x',',','x',',',length,',',0,',',combined_genes[x][(y-1)],',',combined_genes[x][(6)],',',combined_genes[x][(0)],',',combined_genes[x][(12)],',',combined_genes[x][(-7)],',',combined_genes[x][(-6)],',',combined_genes[x][(-5)],',',combined_genes[x][(-4)],',',combined_genes[x][(-3)],',',combined_genes[x][(-2)],',',i)
        else:
            print(combined_genes[x][1], 16)
            break

    if len(time_points_set) >= 18 and combined_genes[x][13] >= (elong_rate*time_points[16]) and combined_genes[x][13] < (elong_rate*time_points[17]) and (elong_rate*time_points[1]) != 0: #increase 2 number by +1
        if length >= 19 and length <= 20: # incease 1 number by +1
            add_4 = []
            add_4 = combined_genes[x][0:15]
            y = 15 + z
            add_4.append(sum(combined_genes[x][15:32])/17) #increase 2 number by +1
            add_4.extend(combined_genes[x][32:y]) #increase 1 number by +1
            data_log_corrected = []
            data_log_corrected = add_4
            add_5 = []
            add_5 = combined_genes[x][0:15]
            add_5.append(combined_genes[x][13]/(elong_rate))
            add_5.extend(time_points[17:z]) #increase 1 number by +1
            time_corrected = []
            time_corrected = add_5
            a = combined_genes[x][1]
            [b, c, d, e, f] = (stats.linregress(time_corrected[-z+16:], data_log_corrected[-z+16:])) #increase 2 number by +1, and reduce 1 below by -1
            csv_results_table.append([a,b,c,d,np.log(2)/-b,length,-17,combined_genes[x][(y-1)],combined_genes[x][(6)],combined_genes[x][(0)],combined_genes[x][(12)],i])
            if d > flag_r_value:
                csv_results_table[-1].append("R_value_FLAG")
            else:
                csv_results_table[-1].extend("x")
            if (np.log(2)/-b) > 2*time_points[z-1] or np.log(2)/-b < 0:
                csv_results_table[-1].append("Half_Life_Flag")
            else:
                csv_results_table[-1].extend("x")
            # csv_results_table[-1].extend(data_log_corrected[-z:])
            csv_results_table[-1].extend(combined_genes[x][15:15+len(time_points_set)])
            #print(a,',',b,',',c,',',d,',',np.log(2)/-b,',',length,',',-17,',',combined_genes[x][(y-1)],',',combined_genes[x][(6)],',',combined_genes[x][(0)],',',combined_genes[x][(12)],',',combined_genes[x][(-7)],',',combined_genes[x][(-6)],',',combined_genes[x][(-5)],',',combined_genes[x][(-4)],',',combined_genes[x][(-3)],',',combined_genes[x][(-2)],',',i)
        elif length <= 18 and length >= 0: #increase 1 number by +1
            add_2 = []
            add_2 = combined_genes[x][0:15]
            y = 15 + z
            zero = [0]
            add_2.extend(zero)
            data_log_corrected = []
            data_log_corrected = add_2
            add_3 = []
            add_3 = combined_genes[x][0:15]
            add_3.extend(zero)
            time_corrected = []
            time_corrected = add_3
            operon_count_total += 1
            csv_results_table.append([add_2[1],'x','x','x','x',length,0,combined_genes[x][(y-1)],combined_genes[x][(6)],combined_genes[x][(0)],combined_genes[x][(12)],i,'x','x'])
            csv_results_table[-1].extend(combined_genes[x][15:15+len(time_points_set)])
            #print(add_2[1],',',add_2[-1],',','x',',','x',',','x',',',length,',',0,',',combined_genes[x][(y-1)],',',combined_genes[x][(6)],',',combined_genes[x][(0)],',',combined_genes[x][(12)],',',combined_genes[x][(-7)],',',combined_genes[x][(-6)],',',combined_genes[x][(-5)],',',combined_genes[x][(-4)],',',combined_genes[x][(-3)],',',combined_genes[x][(-2)],',',i)
        else:
            print(combined_genes[x][1], 17)
            break

    if len(time_points_set) >= 19 and combined_genes[x][13] >= (elong_rate*time_points[17]) and combined_genes[x][13] < (elong_rate*time_points[18]) and (elong_rate*time_points[1]) != 0: #increase 2 number by +1
        if length <= 20: # incease 1 number by +1
            add_4 = []
            add_4 = combined_genes[x][0:15]
            y = 15 + z
            add_4.append(sum(combined_genes[x][15:33])/18) #increase 2 number by +1
            add_4.extend(combined_genes[x][33:y]) #increase 1 number by +1
            data_log_corrected = []
            data_log_corrected = add_4
            add_5 = []
            add_5 = combined_genes[x][0:15]
            add_5.append(combined_genes[x][13]/(elong_rate))
            add_5.extend(time_points[18:z]) #increase 1 number by +1
            time_corrected = []
            time_corrected = add_5
            a = combined_genes[x][1]
            [b, c, d, e, f] = (stats.linregress(time_corrected[-z+17:], data_log_corrected[-z+17:])) #increase 2 number by +1, and reduce 1 below by -1
            csv_results_table.append([a,b,c,d,np.log(2)/-b,length,-18,combined_genes[x][(y-1)],combined_genes[x][(6)],combined_genes[x][(0)],combined_genes[x][(12)],i])
            if d > flag_r_value:
                csv_results_table[-1].append("R_value_FLAG")
            else:
                csv_results_table[-1].extend("x")
            if (np.log(2)/-b) > 2*time_points[z-1] or np.log(2)/-b < 0:
                csv_results_table[-1].append("Half_Life_Flag")
            else:
                csv_results_table[-1].extend("x")
            # csv_results_table[-1].extend(data_log_corrected[-z:])
            csv_results_table[-1].extend(combined_genes[x][15:15+len(time_points_set)])
            #print(a,',',b,',',c,',',d,',',np.log(2)/-b,',',length,',',-18,',',combined_genes[x][(y-1)],',',combined_genes[x][(6)],',',combined_genes[x][(0)],',',combined_genes[x][(12)],',',combined_genes[x][(-7)],',',combined_genes[x][(-6)],',',combined_genes[x][(-5)],',',combined_genes[x][(-4)],',',combined_genes[x][(-3)],',',combined_genes[x][(-2)],',',i)
        elif length <= 18 and length >= 0: #increase 1 number by +1
            add_2 = []
            add_2 = combined_genes[x][0:15]
            y = 15 + z
            zero = [0]
            add_2.extend(zero)
            data_log_corrected = []
            data_log_corrected = add_2
            add_3 = []
            add_3 = combined_genes[x][0:15]
            add_3.extend(zero)
            time_corrected = []
            time_corrected = add_3
            operon_count_total += 1
            csv_results_table.append([add_2[1],'x','x','x','x',length,0,combined_genes[x][(y-1)],combined_genes[x][(6)],combined_genes[x][(0)],combined_genes[x][(12)],i,'x','x'])
            csv_results_table[-1].extend(combined_genes[x][15:15+len(time_points_set)])
            #print(add_2[1],',',add_2[-1],',','x',',','x',',','x',',',length,',',0,',',combined_genes[x][(y-1)],',',combined_genes[x][(6)],',',combined_genes[x][(0)],',',combined_genes[x][(12)],',',combined_genes[x][(-7)],',',combined_genes[x][(-6)],',',combined_genes[x][(-5)],',',combined_genes[x][(-4)],',',combined_genes[x][(-3)],',',combined_genes[x][(-2)],',',i)
        else:
            print(combined_genes[x][1], 18)
            break

    if len(time_points_set) >= 19 and combined_genes[x][13] == (elong_rate*time_points[18]) and (elong_rate*time_points[1]) != 0: #increase 2 number by +1
        if length <= 20: # incease 1 number by +1
            add_4 = []
            add_4 = combined_genes[x][0:15]
            y = 15 + z
            add_4.append(sum(combined_genes[x][15:33])/18) #increase 2 number by +1
            add_4.extend(combined_genes[x][33:y]) #increase 1 number by +1
            data_log_corrected = []
            data_log_corrected = add_4
            add_5 = []
            add_5 = combined_genes[x][0:15]
            add_5.append(combined_genes[x][13]/(elong_rate))
            add_5.extend(time_points[18:z]) #increase 1 number by +1
            time_corrected = []
            time_corrected = add_5
            a = combined_genes[x][1]
            [b, c, d, e, f] = (stats.linregress(time_corrected[-z+17:], data_log_corrected[-z+17:])) #increase 2 number by +1, and reduce 1 below by -1
            csv_results_table.append([a,b,c,d,np.log(2)/-b,length,-18,combined_genes[x][(y-1)],combined_genes[x][(6)],combined_genes[x][(0)],combined_genes[x][(12)],i])
            if d > flag_r_value:
                csv_results_table[-1].append("R_value_FLAG")
            else:
                csv_results_table[-1].extend("x")
            if (np.log(2)/-b) > 2*time_points[z-1] or np.log(2)/-b < 0:
                csv_results_table[-1].append("Half_Life_Flag")
            else:
                csv_results_table[-1].extend("x")
            # csv_results_table[-1].extend(data_log_corrected[-z:])
            csv_results_table[-1].extend(combined_genes[x][15:15+len(time_points_set)])
            #print(a,',',b,',',c,',',d,',',np.log(2)/-b,',',length,',',-18,',',combined_genes[x][(y-1)],',',combined_genes[x][(6)],',',combined_genes[x][(0)],',',combined_genes[x][(12)],',',combined_genes[x][(-7)],',',combined_genes[x][(-6)],',',combined_genes[x][(-5)],',',combined_genes[x][(-4)],',',combined_genes[x][(-3)],',',combined_genes[x][(-2)],',',i)
        elif length <= 18 and length >= 0: #increase 1 number by +1
            add_2 = []
            add_2 = combined_genes[x][0:15]
            y = 15 + z
            zero = [0]
            add_2.extend(zero)
            data_log_corrected = []
            data_log_corrected = add_2
            add_3 = []
            add_3 = combined_genes[x][0:15]
            add_3.extend(zero)
            time_corrected = []
            time_corrected = add_3
            operon_count_total += 1
            csv_results_table.append([add_2[1],'x','x','x','x',length,0,combined_genes[x][(y-1)],combined_genes[x][(6)],combined_genes[x][(0)],combined_genes[x][(12)],i,'x','x'])
            csv_results_table[-1].extend(combined_genes[x][15:15+len(time_points_set)])
            #print(add_2[1],',',add_2[-1],',','x',',','x',',','x',',',length,',',0,',',combined_genes[x][(y-1)],',',combined_genes[x][(6)],',',combined_genes[x][(0)],',',combined_genes[x][(12)],',',combined_genes[x][(-7)],',',combined_genes[x][(-6)],',',combined_genes[x][(-5)],',',combined_genes[x][(-4)],',',combined_genes[x][(-3)],',',combined_genes[x][(-2)],',',i)
        else:
            print(combined_genes[x][1], "18 equals")
            break

    if len(time_points_set) >= 19 and combined_genes[x][13] > (elong_rate*time_points[18]) and (elong_rate*time_points[1]) != 0: #increase 2 number by +1
        if length <= 20 and length >= 0:
            add_2 = []
            add_2 = combined_genes[x][0:15]
            y = 15 + z # correct time point slicing
            zero = [0]
            add_2.extend(zero)
            data_log_corrected = [] #dummy until do complete loop, then can index in with [x][-z:]
            # data_log_corrected.append(add_4) #dummy until do complete loop, then can index in with [x][-z:]
            data_log_corrected = add_2 #dummy until do complete loop, then can index in with [x][-z:]
            add_3 = []
            add_3 = combined_genes[x][0:15]
            add_3.extend(zero)
            time_corrected = [] #dummy until do complete loop, then can index in with [x][-z:]
            # time_corrected.append(add_5) #dummy until do complete loop, then can index in with [x][-z:]
            time_corrected = add_3 #dummy until do complete loop, then can index in with [x][-z:]
            operon_count_total += 1
            csv_results_table.append([add_2[1],'x','x','x','x',length,0,combined_genes[x][(y-1)],combined_genes[x][(6)],combined_genes[x][(0)],combined_genes[x][(12)],i,'x','x'])
            csv_results_table[-1].extend(combined_genes[x][15:15+len(time_points_set)])
            #print(add_2[1],',',add_2[-1],',','x',',','x',',','x',',',length,',',0,',',combined_genes[x][(y-1)],',',combined_genes[x][(6)],',',combined_genes[x][(0)],',',combined_genes[x][(12)],',',combined_genes[x][(-7)],',',combined_genes[x][(-6)],',',combined_genes[x][(-5)],',',combined_genes[x][(-4)],',',combined_genes[x][(-3)],',',combined_genes[x][(-2)],',',i)
        else:
            print(combined_genes[x][1], "end")
            break


with open(filename_output, "w") as csvfile:
    csvwriter = csv.writer(csvfile, delimiter="\t")

    for x in csv_results_table:
        csvwriter.writerow(x)

# x = "numpy_test.csv"
# csv_results_numpyarray = np.array([np.array(xi) for xi in csv_results_table])
# np.savetxt(x, csv_results_numpyarray,fmt='%5s', delimiter="\t")

# print elong_rate
