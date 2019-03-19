"""
This script could be deployed as a background task in a workstation that is used for mumax simulations.
It checks a server for input mx3 files, downloads them when there is an unused GPU, runs the simulations and uploads the results.
Optionally, one can specify certain text filters to be applied on the mx3, e.g. to access a local path.

This is modified from: https://jeffknupp.com/blog/2014/02/11/a-celerylike-python-task-queue-in-55-lines-of-code/
and https://medium.com/@shashwat_ds/a-tiny-multi-threaded-job-queue-in-30-lines-of-python-a344c3f3f7f0
"""

import zmq
import socket
import os, re
import paramiko
import queue, threading
import subprocess
import time, datetime
from dataclasses import dataclass, field
import json
from argparse import ArgumentParser

HOST_PORT = '127.0.0.1:5679'
TASK_SOCKET = zmq.Context().socket(zmq.REQ)
TASK_SOCKET.connect('tcp://' + HOST_PORT)
KILL_STR = 'KILL.mx3'
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
	cache_path: str = 'D:\Xiaoye\Micromagnetics\Kernel_cache'
	smi_path: str = 'C:\\Program Files\\NVIDIA Corporation\\NVSMI\\nvidia-smi'
	local_mx3path: str = 'D:\\JF\\mx3\\'
	ssh_hostname: str = 'astar.nscc.sg'
	mx3_filepath: str = '/home/projects/13000385/mumax3/mx3/'
	mx3running_path: str = '/home/projects/13000385/mumax3/mx3running/'
	remote_datapath: str = '/home/projects/13000385/mumax3/mx3output/'
	remote_username: str = ''
	rsa_key_path: str = ''
	number_of_GPUs: int = 2

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
		self.GPU_ids = list(range(0, self.params.number_of_GPUs)) # GPUs to use
		self.threads = []
		self.ssh_client = paramiko.SSHClient()
		self.ssh_client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
		self.rsakey= paramiko.RSAKey.from_private_key_file(self.params.rsa_key_path)
		self.remoteLog = self.params.remote_datapath +socket.gethostname() + '_log.txt'

		# start workers
		for i in self.GPU_ids:
			t = threading.Thread(target=self.worker, args=(i,))
			t.daemon = True
			t.start()
			self.threads.append(t)

	def run_server(self):

		"""Start listening for tasks."""

		self._socket.bind('tcp://' + self.host_port)
		self.PrintRemote('Simple Job Server (SJS) online. Waiting for jobs...\n')

		filename = ''

		while True:
			#input_string = self._socket.recv_pyobj()
			#print('SJS received command %s.\n'%input_string)
			
			while self.No_work_and_GPU_Idle():
				print('GPU idle, attempting to check for files.\n')
				# filename does NOT include local path
				filename = self.GetAFile()
			
				if filename == '' or filename == KILL_STR:
					break

				elif filename.endswith('.mx3'):

					self.PrintRemote('File downloaded %s.\n'%filename)

					# filter file here
					self.FilterFile(filename)

					# put in local running queue
					self.q.put(filename)

				time.sleep(5)   # wait for some time to prevent instant checking
			
			if filename == KILL_STR:
				# empty queue
				while not self.q.empty():
					try:
						self.q.get(False)
					except queue.Empty:
						continue
					self.q.task_done()

				# kill workers
				for i in self.GPU_ids:
					self.q.put(filename)
				#self._socket.send_pyobj('')
				break
		
			# check the server every 2 mins
			time.sleep(120)

		# stopping
		self.PrintRemote('SJS stopping...\n')
		# block until all tasks are done
		for t in self.threads:
			t.join()

		self.PrintRemote('SJS stopped.\n')

	def worker(self, GPU_id):
		self.PrintRemote('GPU %d started.\n' % GPU_id)

		while True:
			# get from job queue. will block until job is received
			filename = self.q.get()
			self.PrintRemote('GPU %d received job %s.\n' % (GPU_id, filename))
			#  convert to tuple
			kwargs = {}

			if filename == KILL_STR:
				break

			else:
				self._do_work(filename, GPU_id, **kwargs)
				self.PrintRemote('GPU %d completed job %s.\n' % (GPU_id, filename))

			# job completed
			self.q.task_done()

		self.PrintRemote('GPU %d stopped.\n' % GPU_id)

	def _do_work(self, filename, GPU_id, **kwargs):

		"""
		1. Run the simulation
		2. Upload the results
		"""

		# sleep for a short while to ensure than the mx3 file has been written to disk
		time.sleep(0.5)
		local_path_and_filename = os.path.join(self.params.local_mx3path, filename)
		# runs the simulation
		retproc = subprocess.run(['mumax3','-cache', self.params.cache_path,'-gpu', '%d' % GPU_id, local_path_and_filename])
		# check if the simulation ran with no error
		if retproc.returncode != 0:
			self.PrintRemote('***** Error in the mumax script %s. *****\n' % filename)
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
			self.ssh_client.connect(self.params.ssh_hostname, username=self.params.remote_username, pkey = self.rsakey)
			stdin,stdout,stderr = self.ssh_client.exec_command('ls ' + self.params.mx3_filepath)
			ftp_client = self.ssh_client.open_sftp()
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

	def FilterFile(self, filename):

		""" Open file and replace all keywords as indicated by replacement_dict. Note that the keys are treated as
		regular expressions."""

		if self.params.replacement_dict is not None:

			# set local path and filename
			local_path_and_filename = os.path.join(self.params.local_mx3path, filename)

			with open(local_path_and_filename, 'r') as file:
				file_content = file.read()

			# do filtering here
			# NOTE: key is treated as a regular expression
			for key, val in self.params.replacement_dict.items():
				# give a lambda function that returns val, to avoid the extra processing of backslash escapes
				file_content = re.sub(key, lambda x: val, file_content)

			# delete file
			os.remove(local_path_and_filename)

			# write new contents to file
			with open(local_path_and_filename, 'w') as file:
				file.write(file_content)
	
	def UploadData(self,filename,delete=False):
		try:
			newssh_client = self.GetNewSSHClient()
			ftp_client = newssh_client.open_sftp()

			# output folder name
			filename_out = os.path.splitext(filename)[0]+'.out'

			remotedir = self.params.remote_datapath+filename_out
			ftp_client.mkdir(remotedir)
			# newssh_client.exec_command('chmod g+r ' + remotedir)
		except:
			self.PrintRemote('UploadData: Error in opening ssh connection. Files not uploaded.\n')
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
		# append date-time to msg
		msg = '[' + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+']:\t' + msg
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


def submit_local_job(input_str):
	"""Return the result of running the task *runnable* with the given
	arguments."""
	TASK_SOCKET.send_pyobj(input_str)
	results = TASK_SOCKET.recv_pyobj()
	return results

def kill():
	# kill the server by sending the KILL_STR
	submit_local_job(KILL_STR)

def check_running():
	# check if the server is running
	pass

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