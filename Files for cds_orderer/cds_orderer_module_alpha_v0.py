from __future__ import print_function# inports python3's print function
if hasattr(__builtins__, 'raw_input'): #checks if python2
    input = raw_input # if python2 sets input function to python3's input function
import sys
# from operator import itemgetter
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

if __name__=="__main__":
    start_sites=load_predicted_start_sites(sys.argv[1]) #cds_V5_20190801_2_working_correct.csv

operon_plus = []
operon_minus = []
operon_header = []

for x in start_sites:
    if x[4] == "+":
        operon_plus.append(x)
    elif x[4] == "-":
        operon_minus.append(x)
    elif x[4] != "-" and x[4] != "+":
        operon_header.append(x)

operon_plus_sorted = []
operon_minus_sorted = []

operon_plus_sorted = sorted(operon_plus, key = lambda x: int(x[2]))
operon_minus_sorted = sorted(operon_minus, key = lambda x: int(x[2]), reverse=True)

start_sites_sorted = []
start_sites_sorted = operon_header + operon_plus_sorted + operon_minus_sorted




with open("output_ordered.csv", "w") as csvfile:
    csvwriter = csv.writer(csvfile, delimiter="\t")

    for x in start_sites_sorted:
        csvwriter.writerow(x)
