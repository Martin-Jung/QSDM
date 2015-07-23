#!/usr/bin/python

# Experiments in multiprocessing
#
# Test based on script to use a pool of workers to carry out some tasks
# written by R Oudkerk and published in the Python documentation for
# "multiprocessing".

# Peter D. Wilson
# 25 & 26 May 2010

import time
import random
import multiprocessing
import os
from subprocess import *
from multiprocessing import Process, Queue, current_process, freeze_support

#
# Function run by worker processes
#

def worker(input, output):
    for func, args in iter(input.get, 'STOP'):
        result = CallMaxent(args)
        output.put(result)

def CallMaxent(spp):
	spp2 = str(spp[0]).rstrip("\r\n") 
	#spp2 = spp.rstrip("\r\n")
	#spp2 = spp2.replace("_"," ")

	#log_file.write("*** Species = "+spp2+"\n")
	#print spp2

  # Location of samples file/s
	samplepath = '"' + "/mnt/T2/Species data/" + spp2 + "/" + spp2 + " global SWD.csv" + '"'
	#log_file.write("  samplepath = " + samplepath+"\n")

  # Location of environmental layers (current/native climate data files)
	envpath  = '"' + "/mnt/T2/Species data/" + spp2 + "/background1_SWD.csv" + '"'
	#log_file.write("  envpath = " + envpath+"\n")

 	# Location of projection layers for the current model
	projpath = '"' + "/mnt/T1/WorldClim/Current/5m_Global/mybioclim_WCTS_MXE" + '"'

	# Location of Maxent run output, sans enclosing quotes
	outpath  = "/mnt/T3/MaxEnt models/" + spp2 + "/CurrentReplicates"
						
	if not os.path.exists(outpath):
		os.mkdir(outpath)
		#print "Making outpath"
        
	#log_file.write("  outpath = " + outpath +"\n")
	
	# Add enclosing quotes for use as a parameter string
	outpath  = '"' + outpath + '"'
			
  # The Maxent launch command for this combination:
	maxent_cmd   = "java -mx1024m -jar /home/peterwil/MaxEnt3_3_2/maxent.jar -a -z replicates=10 -r -J" \
	+ " -e " + envpath + " -j " + projpath + " -s " + samplepath + " -o " + outpath
  #print maxent_cmd
	#log_file.write("  "+maxent_cmd+"\n")

	p = Popen(maxent_cmd, shell=True).wait()
	
	return spp2 + " processed"
	
def ProcessStuff(spp_list):
	print 'cpu_count() = %d\n' % multiprocessing.cpu_count()
	NUMBER_OF_PROCESSES = multiprocessing.cpu_count()
	TASKS = [(CallMaxent, (spp_list[i],)) for i in range(len(spp_list))]
	#TASKS2 = [(plus, (i, 8)) for i in range(10)]

    # Create queues
	task_queue = Queue()
	done_queue = Queue()

	# Submit tasks
	for task in TASKS:
		task_queue.put(task)

    # Start worker processes
	for i in range(NUMBER_OF_PROCESSES):
		Process(target=worker, args=(task_queue, done_queue)).start()

    # Get and print results
	print 'Unordered results:'
	for i in range(len(TASKS)):
		print '\t', done_queue.get()

    # Tell child processes to stop
	for i in range(NUMBER_OF_PROCESSES):
		task_queue.put('STOP')

if __name__ == '__main__':
	freeze_support()
	
	# Read the list of species to be processed in this run
	spp_file = open("/home/peterwil/Weed Project/EquisetumPatchList.txt",'r')
	spp_list = spp_file.readlines()
	num_spp = len(spp_list)
	print "List of "+str(len(spp_list))+" to be processed."
	
	ProcessStuff(spp_list)
