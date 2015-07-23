#!/usr/bin/python

# Experiments in multiprocessing
#
# Based on script to use a pool of workers to carry out some tasks
# written by R Oudkerk and published in the Python documentation for
# "multiprocessing".
#
# Runs projection of a set of 10 replicate MaxEnt models onto an
# ensemble of GCMs and future times
#
# Peter D. Wilson
# 11 & 16 July June 2010

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
        result = func(*args)
        output.put(result)

def CallMaxEnt(spp,epoch,model,replicate):
	#print "spp = "+spp
	spp2 = spp.replace(" ","_") 
	#print "spp2 = "+spp2
	#print epoch
	#print model
	#print replicate

	# Location of samples file/s
	lambdapath = '"' + "/mnt/T3/MaxEnt models/" + spp + "/CurrentReplicates/" + spp2 +"_"+replicate+".lambdas" + '"'
	#print lambdapath
	
 	# Location of projection layers for the current model
	projpath = '"' + "/mnt/T1/WC_5m_IPCC_AR4_Decade/WC_"+model+"_Decade/"+epoch+"/mybioclim_OZ_WCTS_MXE" + '"'
	#print projpath

	# Location of Maxent run output, sans enclosing quotes
	outpath  = "/mnt/T3/MaxEnt models/"+spp+"/"+model+"_"+epoch
	#print outpath
						
	#if not os.path.exists(outpath):
	#	os.mkdir(outpath)
	
	# Add enclosing quotes and file stub
	outpath  = '"' + outpath + "/" + spp2 + "_" + model+ "_" + epoch + "_" +replicate + '"'
	#print outpath
						
  # The Maxent launch command for this combination:
	maxent_cmd   = "java -mx1024m -cp /home/peterwil/MaxEnt3_3_2/maxent.jar density.Project " + lambdapath +" " + projpath + " " + outpath + " -a -r -z"

	p = Popen(maxent_cmd, shell=True).wait()
	
	return spp + " " + model + " " + epoch + " " + replicate + " processed"
	
def ProcessStuff(spp_list,epoch_list,model_list):
	print 'cpu_count() = %d\n' % multiprocessing.cpu_count()
	
	NUMBER_OF_PROCESSES = multiprocessing.cpu_count()
	task_queue = Queue()
	done_queue = Queue()
	
	for spp in spp_list:
		for model in model_list:
			for epoch in epoch_list:
				TASKS = [(CallMaxEnt,(spp.rstrip("\r\n"),epoch.rstrip("\r\n"),model.rstrip("\r\n"),str(i))) for i in range(10)]

				#print "Number of projections to be made = %d\n" % len(TASKS)
				#print TASKS
				print "   "+spp

				# Submit tasks
				for task in TASKS:
					#print task,"\n"
					task_queue.put(task)

				# Start worker processes
				for i in range(NUMBER_OF_PROCESSES):
					Process(target=worker, args=(task_queue, done_queue)).start()

				# Get and print results
				print 'Unordered results for '+spp.rstrip("\r\n")+':'
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
	
	model_file = open("/home/peterwil/MaxEnt Scripts/Model_list.txt",'r')
	model_list = model_file.readlines()
	num_models = len(model_list)
	#print "GCMs to be processed:"
	#print model_list

	epoch_file = open("/home/peterwil/MaxEnt Scripts/Epoch_list.txt",'r')
	epoch_list = epoch_file.readlines()
	num_epochs = len(epoch_list)
	#print "Epochs to be processed:"
	#print epoch_list
	
	ProcessStuff(spp_list,epoch_list,model_list)
