#!/usr/bin/env python

import os
import subprocess as sp
#import argparse

# ----------------------------------------------------------------

def check_queue_status(queue_name):
	"""
	Checks the queue status specified by queue_name
	
	queue_name: The name of the queue to check
	
	returns: True if running, False if Not Running
	"""

	cmd_args = ["/etc/init.d/%s" % queue_name, "status"]
	
	process = sp.Popen(cmd_args, stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE)
	response, err = process.communicate()

	# fetch the last element from the list
	response_str = response.strip().split(' ')[-1]
	
	# status initialization
	status = False

	if response_str.find("[Running]") != -1:
		status = True
	
	return status

# ----------------------------------------------------------------

def start_queue(queue_name, attempts = 6):
	"""
	Starts the queue specified by queue_name
	
	queue_name: The name of the queue to start
	attempts: The number of attempts to try and start the
	queue 

	returns: True on success, False on failure
	"""

	cmd_args = ["/etc/init.d/%s" % queue_name, "start"]
	
	queue_status = check_queue_status(queue_name)

	while not queue_status:
		process = sp.Popen(cmd_args, stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE)
		response, err = process.communicate()
	
		# fetch the last element from the list
        	response_str = response.split(' ')[-1]
		
		queue_status = check_queue_status(queue_name)
		
		# exit loop if status was successful
		if attempts == 0:
			break
			
		attempts -= 1

	# reached maximum attempts and queue status is still False	
	if queue_status is False and attempts==0:
		return False
	
	return True
	
# ----------------------------------------------------------------

def parse_arguments():
	"""
	Uses python's argparse to parse the command line arguments
	
	return: Argparse parser object
	"""

	# create a new argument parser object
    	parser = argparse.ArgumentParser(description='Update scores for new release')

    	# group required arguments together
    	req_args = parser.add_argument_group("required arguments")
    	req_args.add_argument('-q', help='A comma separated list of queues to watch',
                        type=list, required=True)
	req_args.add_argument('--attempts', help='Number of attempts to start the queue',
                        type=int, required=True, default=6)
	
	return parser
	
# ----------------------------------------------------------------

if __name__ == '__main__':
	
	# Make argparse available
	# create a new argument parser object
    	#parser = parse_arguments()
    	#args = parser.parse_args()

	queues = ["rfam_queue", "rfam_batch_queue"]

	for queue_name in queues:
		status = start_queue(queue_name, attempts = 6)

		if status is False:
			print "Post slack notification"
		
