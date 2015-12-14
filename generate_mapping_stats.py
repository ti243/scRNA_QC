#!/usr/bin/python
import sys 
import numpy as np
import pandas as pd
import re
import matplotlib.pyplot as plt 
from os import listdir
from os.path import isfile, join
import seaborn as sb
import math
import scipy
import getopt
from collections import Counter, defaultdict, OrderedDict
from Bio import SeqIO
import glob
import fnmatch
import os
def generate_mappability_table(input_file, output_file, sample_name):

  
    multi_maps = defaultdict(str)
    sam_files = [input_file]
    mappability_freq_table = OrderedDict()

    #MAPPABILITY
    mappability_freq_table["total"] = 0 
    mappability_freq_table["mapped"] = 0
    mappability_freq_table["unmapped"] = 0
    mappability_freq_table["unique"] = 0
    mappability_freq_table["multi"] = 0
    
    #ERROR
    mappability_freq_table["perfect"] = 0
    mappability_freq_table["partly_perfect"] = 0
    mappability_freq_table["mapped_no_correct"] = 0
    
    for i in range(0,10):
        mappability_freq_table["S_"+str(i)] = 0
    mappability_freq_table["S_10+"] = 0
    mappability_freq_table["I"] = 0
    mappability_freq_table["D"] = 0
    mappability_freq_table["INDEL"] = 0
    reads = Counter()
    multi_reads = defaultdict(str)
    for sam in sam_files:
        print("RUN: " + sam)
        with open(sam) as f:

            for line in f:
                split=line.split("\t")
                if (not line.startswith("@PG") and not line.startswith("@HD") and not line.startswith("@SQ") and len(split) >= 10):
                    read_name=split[0]
                    chrom = split[2]
                    pos = split[3] 
                    errors=split[5]
                    read=split[9]
                
                    flagCode = int(split[1])

                    errors_a = list(errors)
                    number = ""
                    num = 0
                    error_table = defaultdict(int)
                   
                    name_and_flag = read_name

                    #CHECK IF READ MAPPED OR UNMAPPED 
                    #IT US UNMAPPED
                    if(flagCode & 0x0004 != 0):
                        mappability_freq_table["unmapped"] += 1
                        mappability_freq_table["total"] += 1
                        error_table["*"] += 1
                    #IT IS MAPPED
                    else:
                        if (flagCode & 0x0001 !=0):                             #This is paired end sequencing
                            if (flagCode & 0x0040 != 0):                        #1st read
                                name_and_flag += ";first"
                            if (flagCode & 0x0080 != 0):                        #2nd read
                                 name_and_flag  += ";second"

                        #CHECK IF READ UNIQUELY OR MULTI-MAPPED
                        if(read_name not in reads):
                            reads[name_and_flag] += 1
                            mappability_freq_table["unique"] += 1
                            mappability_freq_table["total"] += 1
                            mappability_freq_table["mapped"] += 1

                        else:
                            if(reads[name_and_flag] == 1):
                                mappability_freq_table["unique"] -= 1
                                mappability_freq_table["multi"] += 1
                            reads[name_and_flag] += 1


                        #PARSE MAPPING ERRORS         
                        for i in errors_a:
                            if (re.match("[0-9]",i)):
                                number +=(i)
                            elif(re.match("[A-Z]",i)):
                                num = int(number)
                                error_table[i] += num
                                number = ""
                        #print mappability_freq_table
                        #TABLE OF HOW MANY READS MAP PERFECT, PARTLY PERFECT, SUBSTITUINTS ETC
                        if("M" in  error_table and len(error_table)==1):
                            mappability_freq_table["perfect"] += 1
                        elif("M" in error_table and len(error_table) > 1):
                            mappability_freq_table["partly_perfect"] += 1
                        elif("M" not in error_table and "*" not in error_table):
                            mappability_freq_table["mapped_no_correct"] += 1
                        
                        if("S" in error_table):
                            if(int(error_table["S"]) < 10):
                                mappability_freq_table["S_"+str(error_table["S"])] += 1
                            else:
                                mappability_freq_table["S_10+"] += 1
                        elif("S" not in error_table):
                            mappability_freq_table["S_0"] += 1

                        if("I" in error_table):
                            mappability_freq_table["I"] += 1

                        if("D" in error_table):
                            mappability_freq_table["D"] += 1
                                        
                        if("I" in error_table or "D" in error_table):
                             mappability_freq_table["INDEL"] += 1

                    
    f = open(output_file,'w')

    o = sample_name
    for k,v in mappability_freq_table.items():
        print str(k) + ":" + str(v)
        o += "," + str(v)
    o += "\n" 
    f.write(o)
    f.close()

def main(argv):    
    input_file = ""
    output_dir = ""
    sample_name = ""
    try:    
        opts, args = getopt.getopt(argv, "i:o:n:")
    except getopt.GetoptError:    
        sys.exit(2)    
        usage()
    for opt, arg in opts:    
        if opt in ("-i"):    
            input_file = arg    
        elif opt in ("-o"): 
            output_dir = arg
        elif opt in ("-n"): 
            sample_name = arg
    if input_file != "" and  output_dir != "" and sample_name != "":  
        generate_mappability_table(input_file, output_dir, sample_name)    
    else: 
        print ("Argument missing")

main(sys.argv[1:])
