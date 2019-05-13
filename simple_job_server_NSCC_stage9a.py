"""
This script could be deployed as a background task in a workstation that is used for mumax simulations.
It checks a server for input mx3 files, downloads them when there is an unused GPU, runs the simulations and uploads the results.
Optionally, one can specify certain text filters to be applied on the mx3, e.g. to access a local path.

This is modified from: https://jeffknupp.com/blog/2014/02/11/a-celerylike-python-task-queue-in-55-lines-of-code/
and https://medium.com/@shashwat_ds/a-tiny-multi-threaded-job-queue-in-30-lines-of-python-a344c3f3f7f0

best to insert the termination request txt file in the mx3running folder instead of mx3 to ensure fastest termination
"""

import zmq
import socket
import os, re, sys
import paramiko
import queue, threading
import subprocess
import time, datetime
from dataclasses import dataclass, field
import json
from argparse import ArgumentParser
import pandas as pd
import numpy as np

HOST_PORT = '127.0.0.1:5679'
TASK_SOCKET = zmq.Context().socket(zmq.REQ)
TASK_SOCKET.connect('tcp://' + HOST_PORT)
INPUT_PARAMS_FILENAME = 'input_parameters.json'


def update_obj_from_dict_recursively(some_obj, some_dict):
	"""
	Useful to convert nested json files into nested dataclasses

	Example
	-------
	@dataclass
	class InputOptions:
		path: str = None
		n1: float = 0
		some_thing: Opt2 = Opt2()

	@dataclass
	class Opt2:
		some_string:str = None
		some_num:float = 0.0

	...

	input_dict_list = UI_load_json_file()

	...

	in_opt = InputOptions()
	update_obj_from_dict_recursively (in_opt, input_dict)

	:param some_obj:
	:param some_dict:
	:return:
	"""
	for k in some_dict:
		tmp_v = some_dict[k]
		if isinstance(tmp_v, dict):
			tmp_obj = some_obj.__dict__[k]
			# check if tmp_obj is actually a dict
			if isinstance(tmp_obj, dict):
				# direct update
				some_obj.__dict__[k] = tmp_v
			else:
				# some other object
				update_obj_from_dict_recursively(tmp_obj, tmp_v)

		else:
			some_obj.__dict__.update({k:tmp_v})

@dataclass
class Parameters:
	cache_path: str = ''
	smi_path: str = 'C:\\Program Files\\NVIDIA Corporation\\NVSMI\\nvidia-smi'
	local_mx3path: str = ''
	ssh_hostname: str = 'astar.nscc.sg'
	mx3_filepath: str = ''
	mx3running_path: str = ''
	remote_datapath: str = ''
	remote_username: str = ''
	rsa_key_path: str = ''
	GPU_ids: [int] = None
	termination_datapath: str = ''
	mx3error_path: str = ''

	# used to scan an input mx3 file for keywords to replace
	replacement_dict: dict = field(default_factory=dict)

class Server(object):

	"""A remote task executor."""

	def __init__(self, host_port_in=HOST_PORT):

		self.params = Parameters()

		# Load parameters if available
		if os.path.isfile(INPUT_PARAMS_FILENAME):
			with open(INPUT_PARAMS_FILENAME) as data_file:
				params_dict = json.load(data_file)
				update_obj_from_dict_recursively(self.params, params_dict)

		#Initialize Server
		self.q = queue.Queue()
		self.host_port = host_port_in
		self._context = zmq.Context()
		self._socket = self._context.socket(zmq.REP)
		self.GPU_ids = self.params.GPU_ids # GPUs to use
		self.threads = []
		self.ssh_client = paramiko.SSHClient()
		self.ssh_client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
		self.rsakey= paramiko.RSAKey.from_private_key_file(self.params.rsa_key_path)
		self.remoteLog = self.params.remote_datapath +socket.gethostname() + '_log.txt'
		# self.runningfiles = []

		# start workers
		for i in self.GPU_ids:
			t = threading.Thread(target=self.worker, args=(i,))
			# daemon threads are killed automatically when main thread exits
			t.daemon = True
			t.start()
			self.threads.append(t)

		# initialise a dictionary mapping gpu number to filenames
		gpu_names = [num for num in self.GPU_ids]
		self.gpu_dict = {gpu: None for gpu in gpu_names}
		# print(self.gpu_dict)

		# initialise a dictionary mapping gpu number to its PID number
		self.gpu_pid_dict = {gpu: None for gpu in gpu_names}

		self.isterminate = False

	def run_server(self):
		"""Start listening for tasks."""

		self._socket.bind('tcp://' + self.host_port)
		self.PrintRemote('Simple Job Server (SJS) online. Waiting for jobs...\n')

		# filename = ''

		while True:
			#input_string = self._socket.recv_pyobj()
			#print('SJS received command %s.\n'%input_string)
			
			while self.No_work_and_GPU_Idle():
				print(self.AppendDateTime('GPU idle, attempting to check for files.\n'))
				# filename does NOT include local path
				filename = self.GetAFile()
			
				if filename == '': # or filename == KILL_STR:
					break

				elif filename.endswith('.mx3'):

					self.PrintRemote('File downloaded %s.\n'%filename)

					# filter file here
					local_path_and_filename = os.path.join(self.params.local_mx3path, filename)
					self.FilterFile(self.params.replacement_dict, local_path_and_filename, local_path_and_filename)

					# get the PIDs all the mumax3.exe files before a new one spawns
					original_pids = self.get_pid()

					# put in local running queue
					self.q.put(filename)

				# wait for 2 mins to allow the job to start utilisation of GPU
				# allow 1 minute for the mx3.exe instance to spawn 
				# for testing purposes reduce time to 30s
				time.sleep(60)

				# a mumax3.exe process has spawn
				# find the PID of the mumax3.exe that was just spawned
				new_pids = self.get_pid()
				spawned_pid = [x for x in new_pids if x not in original_pids]
				# print("Spawned PID is", spawned_pid)

				# get the gpu of the PID that was just spawned
				gpu_num = self.get_gpu(filename)
				self.gpu_pid_dict[gpu_num] = spawned_pid[0]
				# print("GPU is", gpu_num)

				# allow another minute before the server can check for any mx3 files to run
				time.sleep(60)

			# checks if there are txt or csv files in the mx3 running folder 
			# returns filename of txt/csv file if exists 
			end_filename = self.check_requests()

			if end_filename != "":
				print("End job(s) file has been detected in server")
				print('Attempting to read for end job(s) requests...\n')

				# read the files that are to be ended
				users, end_files, _ = self.read_file(end_filename)
				
				# check if the filename in txt file is the file that is currently running on the gpu
				for i in range(len(end_files)):
					user = users[i]
					filename = end_files[i]
					# loops over all the files in the gpus
					for key, runningfile in self.gpu_dict.items():
						if filename == runningfile:
							print("Received request to end job name:", runningfile)
							# turn on terminating switch
							self.isterminate = True 
							# endjob method here
							ret = self.stop_work(runningfile)

							if ret == 0:
								print('File:', runningfile, 'has successfully terminated\n')
								# remove the mx3 file from the mx3 running folder as well
								self.delete_request(runningfile)
								# remove this file from the dictionary
								self.gpu_dict[key] = None
								
							else:
								print("Request to end", runningfile, "is unsuccessful")
							
							# output txt file to mx3 output folder
							self.update_termination(user, filename, ret)

							# once all the file  from request list terminated, turn off the terminating switch
							# so that if a process is terminated not intentioanaly, output error message as usual
							self.isterminate = False

			# check the server every 2 mins
			# time.sleep(120)
			# for testing purposes reduce time to 10s
			time.sleep(120)

			# delete the txt/csv file after 2min, to ensure every computer has downloaded the file
			if end_filename != "":
				self.delete_request(end_filename)
				print('End request file:', end_filename, 'has been removed from mx3 running folder\n')

	def worker(self, GPU_id):
		self.PrintRemote('GPU %d started.\n' % GPU_id)

		while True:
			# get from job queue. will block until job is received
			filename = self.q.get()
			# insert the filename to the dictionary containing the gpu id
			self.gpu_dict[GPU_id] = filename
			# print(self.gpu_dict)

			self.PrintRemote('GPU %d received job %s.\n' % (GPU_id, filename))
			#  convert to tuple
			kwargs = {}

			self._do_work(filename, GPU_id, **kwargs)

			if not self.isterminate:
				self.PrintRemote('GPU %d completed job %s.\n' % (GPU_id, filename))

			# job completed
			self.q.task_done()

		# self.PrintRemote('GPU %d stopped.\n' % GPU_id)

	def _do_work(self, filename, GPU_id, **kwargs):

		"""
		1. Run the simulation
		2. Upload the results
		"""

		# sleep for a short while to ensure than the mx3 file has been written to disk
		time.sleep(0.5)
		local_path_and_filename = os.path.join(self.params.local_mx3path, filename)
		# runs the simulation
		retproc = subprocess.run(['mumax3','-cache', self.params.cache_path,'-gpu', '%d' % GPU_id, local_path_and_filename], capture_output=True)
		# check if the simulation ran with no error
		if retproc.returncode != 0:
			# if file is intentionally terminated, break out of the method so no error message printed out
			if self.isterminate == True:
				return
			
			# senario where the remote gpu is unable to run the file (syntax error etc)
			self.PrintRemote('***** Error in the mumax script %s, running on GPU %d. *****\nOutput: %s\nError msg: %s\n\n' % (filename, GPU_id, retproc.stdout, retproc.stderr))
			# move the file out of the mx3 running and into a folder called 'mx3 error'
			try:
				self.ssh_client.connect(self.params.ssh_hostname, username=self.params.remote_username, pkey = self.rsakey)
				self.ssh_client.exec_command('mv ' + self.params.mx3running_path+filename + ' ' + self.params.mx3error_path+filename)
				time.sleep(0.5)
				self.PrintRemote('Moved error file %s to mx3 error folder.\n'%filename)
				self.ssh_client.close()

			except:
				print("Unable to move error file from mx3running to mx3error ")

			return

		# upload the results
		self.UploadData(filename, False)
		try:
			self.ssh_client.connect(self.params.ssh_hostname, username=self.params.remote_username, pkey = self.rsakey)
			ftp_client = self.ssh_client.open_sftp()
			ftp_client.remove(self.params.mx3running_path + filename)
			ftp_client.close()
			self.ssh_client.close()
		except Exception as Argument:
			self.PrintRemote('***** Error in _do_work: Unable to connect to remove mx3 file in mx3running. *****\n Exception: %s'%Argument)
		
		self.PrintRemote('%s done.\n' % filename)

	def No_work_and_GPU_Idle(self):
		"""
		If there is no jobs in queue and any gpu usage is below 10%, counted as idle
		:return:
		"""
		if not self.q.empty():
			return False

		stdout = subprocess.run([self.params.smi_path,'--query-gpu=utilization.gpu','--format=csv'],capture_output=True).stdout
		gpu_usages = [int(s) for s in stdout.split() if s.isdigit()]
		for usage in gpu_usages:
			if usage < 10:
				return True
		return False
		
	def GetAFile(self):

		ftp_client = None
		try:
			newssh_client = self.GetNewSSHClient()
			self.ssh_client.connect(self.params.ssh_hostname, username=self.params.remote_username, pkey = self.rsakey)
			stdin,stdout,stderr = self.ssh_client.exec_command('ls ' + self.params.mx3_filepath)
			ftp_client = newssh_client.open_sftp()
			ftp_client.chdir(self.params.mx3_filepath)
		except:
			self.PrintRemote('***** Error in opening connection. Will return empty string as filename. *****\n')
			if not ftp_client is None:
				ftp_client.close()
			self.ssh_client.close()
			return ''
		
		# check if there are any files with priority
		for i in stdout.readlines():
			filename = i.rstrip()
			priority = re.findall('priority', filename)
			
			if priority != []:
				try:
					print('DISCOVERED FILE WITH PRIORITY')
					self.ssh_client.exec_command('mv ' + self.params.mx3_filepath+filename + ' ' + self.params.mx3running_path+filename)
					time.sleep(0.5)
					local_path_and_filename = os.path.join(self.params.local_mx3path,filename)
					ftp_client.get(self.params.mx3running_path+filename, local_path_and_filename)
					self.PrintRemote('Moved %s to local drive for running.\n'%filename)
					ftp_client.close()
					self.ssh_client.close()
					
					return filename
						
				except Exception as Argument:
					print('***** Error in trying to get a file. *****\n Exception: %s'%Argument)

		#  close the client first
		ftp_client.close()
		self.ssh_client.close()

		try:
			newssh_client = self.GetNewSSHClient()
			self.ssh_client.connect(self.params.ssh_hostname, username=self.params.remote_username, pkey = self.rsakey)
			stdin,stdout,stderr = self.ssh_client.exec_command('ls ' + self.params.mx3_filepath)
			ftp_client = newssh_client.open_sftp()
			ftp_client.chdir(self.params.mx3_filepath)
		except:
			self.PrintRemote('***** Error in opening connection. Will return empty string as filename. *****\n')
			if not ftp_client is None:
				ftp_client.close()
			self.ssh_client.close()
			return ''
			
		for i in stdout.readlines():
			filename = i.rstrip()
			self.PrintRemote('Found filename on server: %s\n'%filename)

			try:
				# fileObj = ftp_client.file(filename+'.lock',mode='x')   # create lock file for concurrency
				# TODO: Need to resolve concurrency issue. Check to see if move operation was successful
				self.ssh_client.exec_command('mv ' + self.params.mx3_filepath+filename + ' ' + self.params.mx3running_path+filename)
				time.sleep(0.5)
				local_path_and_filename = os.path.join(self.params.local_mx3path,filename)
				ftp_client.get(self.params.mx3running_path+filename, local_path_and_filename)

				# ftp_client.remove(mx3_filepath+filename)
				# fileObj.close()
				# ftp_client.remove(mx3_filepath+filename+'.lock')
				self.PrintRemote('Moved %s to local drive for running.\n'%filename)
				ftp_client.close()
				self.ssh_client.close()

				return filename

			except Exception as Argument:
				print('***** Error in trying to get a file. *****\n Exception: %s'%Argument)
		ftp_client.close()
		self.ssh_client.close()
		return ''

	# change to static method so that this can be used by others
	@staticmethod
	def FilterFile(replacement_dict: dict, old_path_and_filename, new_path_and_file_name):

		""" Open file and replace all keywords as indicated by replacement_dict. Note that the keys are treated as
		regular expressions."""

		if replacement_dict is not None:

			# set local path and filename
			# local_path_and_filename = os.path.join(params.local_mx3path, filename)

			with open(old_path_and_filename, 'r') as file:
				file_content = file.read()

			# do filtering here
			# NOTE: key is treated as a regular expression
			for key, val in replacement_dict.items():
				# give a lambda function that returns val, to avoid the extra processing of backslash escapes
				file_content = re.sub(key, lambda x: val, file_content)

			# write new contents to file
			# this discard all contents if the file already exist
			# use the linux newline standard
			with open(new_path_and_file_name, 'w', newline='\n') as file:
				file.write(file_content)
	
	def UploadData(self,filename,delete=False):
		try:
			newssh_client = self.GetNewSSHClient()
			ftp_client = newssh_client.open_sftp()

			# output folder name
			filename_out = os.path.splitext(filename)[0]+'.out'

			# TODO: check for existence of folder. Create a unique folder if needed
			remotedir = self.params.remote_datapath+filename_out
			ftp_client.mkdir(remotedir)
			# newssh_client.exec_command('chmod g+r ' + remotedir)
		except:
			self.PrintRemote('***** UploadData: Error in opening ssh connection. Files not uploaded. *****\n')
			return

		localdir = os.path.join(self.params.local_mx3path, filename_out)
		listoffiles = os.listdir(localdir)
		for file in listoffiles:
			ftp_client.put(os.path.join(localdir, file), remotedir+'/'+file)
			# newssh_client.exec_command('chmod g+r ' + remotedir+'/'+file)

			# should nv delete local copy programmatically

			# if delete:
			# 	os.remove(localdir+file)
			
		# if delete:
		# 	os.rmdir(localdir)
		# 	os.remove(self.params.local_mx3path+filename)
		
		# newssh_client.exec_command('chgrp -Rh 13000385' + remotedir)
		ftp_client.close()
		newssh_client.close()

		self.PrintRemote('Uploaded results of %s to remote location.\n' % filename)
		
	def PrintRemote(self,msg):
		# TODO: Is it possible to print to one central log file? Perhaps the program has to acquire a mutex before writing to log?
		# append date-time to msg
		msg = self.AppendDateTime(msg)

		print(msg)
		try:
			ssh_client = self.GetNewSSHClient()
			ftp_client = ssh_client.open_sftp()
			file = ftp_client.file(self.remoteLog,'a+')
			file.write(msg)
			file.close()
			ftp_client.close()
			ssh_client.close()
		except Exception as Argument:
			print('***** Error writing to remote log file. *****\n Exception: %s'%Argument)
			
	def GetNewSSHClient(self):
		ssh_client = paramiko.SSHClient()
		ssh_client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
		rsakey= paramiko.RSAKey.from_private_key_file(self.params.rsa_key_path)
		ssh_client.connect(self.params.ssh_hostname, username=self.params.remote_username, pkey = rsakey)
		return ssh_client

	@staticmethod
	def AppendDateTime(msg):
		return '[' + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+']:\t' + msg

	def check_requests(self):
		ftp_client = None

		# opens the remote directory
		try:
			self.ssh_client.connect(self.params.ssh_hostname, username=self.params.remote_username, pkey = self.rsakey)
			stdin,stdout,stderr = self.ssh_client.exec_command('ls ' + self.params.mx3running_path)
			ftp_client = self.ssh_client.open_sftp()
		except:
			self.PrintRemote('***** Error in opening connection. Will return empty string as filename. *****\n')
			if not ftp_client is None:
				ftp_client.close()
			self.ssh_client.close()
			return ''

		for i in stdout.readlines():
			filename = i.rstrip()

			# checks if the remote directory has a txt file
			if filename.endswith('.txt') or filename.endswith('.csv'):
				self.PrintRemote('Found filename on server: %s\n'%filename)
				try:
					# copies the txt file from the mx3 folder to the mx3 running path folder
					# DO NOT COPY ANY OTHER FILE TO LOCAL DIRECTORY 
					self.ssh_client.exec_command('mv ' + self.params.mx3running_path+filename)
					time.sleep(0.5)
					local_path_and_filename = os.path.join(self.params.local_mx3path,filename)
					# copies the txt file from the running path to the local directory
					ftp_client.get(self.params.mx3running_path+filename, local_path_and_filename)

					self.PrintRemote('Moved %s to local drive for running.\n'%filename)
					ftp_client.close()
					self.ssh_client.close()

					return filename

				except Exception as Argument:
					print('***** Error in trying to get a file. *****\n Exception: %s'%Argument)

		ftp_client.close()
		self.ssh_client.close()
		return ''

	def read_file(self, end_filename):
		local_path_and_filename = os.path.join(self.params.local_mx3path, end_filename)

		## if txt file, convert to csv file. if csv, stays csv
		csv_local_path_and_filename = local_path_and_filename

		try:
			request = pd.read_csv(csv_local_path_and_filename, delimiter = ',')
			users = request['users']
			filenames = request['filenames']
			status = request['status']

		## if there is formatting error, just return empty filename which does not exist
		except:
			print("Fomatting error in the txt file")
			users = ""
			filenames = ""
			status = ""
			
		return users, filenames, status
			
			
	def stop_work(self, filename):
		# get the gpu that is running this file
		gpu_num = self.get_gpu(filename)
		# get pid of the gpu that is running the end request file
		pid = self.gpu_pid_dict[gpu_num]

		# forcefully terminate mumax3.exe instance that is running that end request file
		access_func = "taskkill /F /PID " + pid
		return os.system(access_func)

	def update_termination(self, user, filename, status):
		local_path_and_filename = os.path.join(self.params.local_mx3path, filename)
		columns = ['user', 'filenames', 'status']

		# if successful termination, status column list as 't'
		if status == 0:
			data = [[user, filename, 't']]

		else:
			data = [[user, filename, 'f']]

		csvfile = pd.DataFrame(data = data, columns = columns)
		t_file = user + '_' + filename + '.csv'
		csvfile = csvfile.to_csv(os.path.join(self.params.local_mx3path, t_file), index = None, header=True)

		# upload the file to the remote directory 
		try:
			newssh_client = self.GetNewSSHClient()
			ftp_client = newssh_client.open_sftp()

		except:
			self.PrintRemote('***** UploadData: Error in opening ssh connection. Files not uploaded. *****\n')
			return

		ftp_client.put(os.path.join(self.params.local_mx3path, t_file), os.path.join(self.params.termination_datapath, t_file))
		ftp_client.close()
		self.ssh_client.close()

	# deleting any files from the mx3 running path (terminated mx3 file, request txt file etc.)
	def delete_request(self, filename):
		try:
			self.ssh_client.connect(self.params.ssh_hostname, username=self.params.remote_username, pkey = self.rsakey)
			stdin,stdout,stderr = self.ssh_client.exec_command('ls ' + self.params.mx3running_path)
			ftp_client = self.ssh_client.open_sftp()
		except:
			self.PrintRemote('***** Error in opening connection. Will return empty string as filename. *****\n')
			if not ftp_client is None:
				ftp_client.close()
			self.ssh_client.close()
			return ''

		ftp_client.remove(self.params.mx3running_path+filename)

		ftp_client.close()
		self.ssh_client.close()

	def get_gpu(self, filename):
		for gpu, running_filename in self.gpu_dict.items():
			if running_filename == filename:
				return gpu

	def get_pid(self, process_name = 'mumax3.exe'):
		return [item.split()[1] for item in os.popen('tasklist').read().splitlines()[4:] if process_name in item.split()]

def submit_local_job(input_str):
	"""Return the result of running the task *runnable* with the given
	arguments."""
	TASK_SOCKET.send_pyobj(input_str)
	results = TASK_SOCKET.recv_pyobj()
	return results

def main():
	# grab all the command line arguments
	# to start server, run this with no arguments
	
	# parser = ArgumentParser(description='Start or perform other actions on the job server.')
	# parser.add_argument('-a', '--add', dest='sim_file', help='add simulation to server', metavar='FILE')

	# args = parser.parse_args()

	# if not args:
	
	w = Server()
	w.run_server()

if __name__ == '__main__':

	main()
