# Modified from: https://jeffknupp.com/blog/2014/02/11/a-celerylike-python-task-queue-in-55-lines-of-code/
# and https://medium.com/@shashwat_ds/a-tiny-multi-threaded-job-queue-in-30-lines-of-python-a344c3f3f7f0

"""Broker-less distributed task queue."""

import zmq
import os
import paramiko
import queue, threading
import subprocess
import time
from argparse import ArgumentParser

HOST_PORT = '127.0.0.1:5679'
TASK_SOCKET = zmq.Context().socket(zmq.REQ)
TASK_SOCKET.connect('tcp://' + HOST_PORT)
KILL_STR = 'KILL.mx3'

cache_path = 'D:\Xiaoye\Micromagnetics\Kernel_cache'
smi_path = 'C:\\Program Files\\NVIDIA Corporation\\NVSMI\\nvidia-smi'
ssh_hostname = 'astar.nscc.sg'
mx3_filepath = '/home/projects/13000385/mumax3/mx3/'
mx3running_path = '/home/projects/13000385/mumax3/mx3running/'
local_mx3path = 'D:\\JF\\mx3\\'
remote_datapath = '/home/projects/13000385/mumax3/mx3output/'

class Server(object):
	"""A remote task executor."""

	def __init__(self, host_port_in=HOST_PORT):
		"""Initialize Server."""
		self.q = queue.Queue()
		self.host_port = host_port_in
		self._context = zmq.Context()
		self._socket = self._context.socket(zmq.REP)
		self.GPU_ids = [0,1] # GPUs to use
		self.threads = []
		self.ssh_client = paramiko.SSHClient()
		self.ssh_client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
		self.rsakey= paramiko.RSAKey.from_private_key_file('id_rsa')

		# start workers
		for i in self.GPU_ids:
			t = threading.Thread(target=self.worker, args=(i,))
			t.daemon = True
			t.start()
			self.threads.append(t)

	def run_server(self):
		"""Start listening for tasks."""
		self._socket.bind('tcp://' + self.host_port)
		print('Simple Job Server (SJS) online. Waiting for jobs...\n')

		while True:
			#input_string = self._socket.recv_pyobj()
			#print('SJS received command %s.\n'%input_string)
			
			while self.GPU_Idle():
				print('GPU idle, attempting to check for files.\n')
				filename = self.GetAFile()
			
				if filename == '' or filename == KILL_STR:
					break
				elif filename[len(filename)-3:len(filename)]=='mx3':
					print('File downloaded %s.\n'%filename)
					self.q.put(filename)
				time.sleep(5)   # wait for some time to prevent instant checking
			
			if filename == KILL_STR:
				# empty queue
				self.q.empty()

				# kill workers
				for i in self.GPU_ids:
					self.q.put(filename)
				#self._socket.send_pyobj('')
				break
		
			# check the server every 10 mins
			time.sleep(120)

		# stopping
		print('SJS stopping...\n')
		# block until all tasks are done
		for t in self.threads:
			t.join()

		print('SJS stopped.\n')

	def worker(self, GPU_id):
		print('GPU %d started.\n' % GPU_id)

		while True:
			# get from job queue
			input_string = self.q.get()
			print('GPU %d received job %s.\n' % (GPU_id, input_string))
			#  convert to tuple
			kwargs = {}

			if input_string == KILL_STR:
				break

			else:
				self._do_work(input_string, GPU_id, **kwargs)
				print('GPU %d completed job %s.\n' % (GPU_id, input_string))

		print('GPU %d stopped.\n' % GPU_id)

	def _do_work(self, mumax_file_str, GPU_id, **kwargs):
		"""Return the result of executing the given task."""
		# sleep for a short while to ensure than the mx3 file has been written to disk
		time.sleep(0.5)
		subprocess.run(['mumax3','-cache',cache_path,'-gpu', '%d' % GPU_id, mumax_file_str])
		self.UploadData(mumax_file_str,False)
	
		self.ssh_client.connect(ssh_hostname, username='kongjf', pkey = self.rsakey)
		ftp_client = self.ssh_client.open_sftp()
		ftp_client.remove(mx3running_path+mumax_file_str)
		ftp_client.close()
		
		print('%s done.\n' % mumax_file_str)
		
	# if any gpu usage is below 10%, counted as idle	
	def GPU_Idle(self):
		stdout = subprocess.run([smi_path,'--query-gpu=utilization.gpu','--format=csv'],capture_output=True).stdout
		gpu_usages = [int(s) for s in stdout.split() if s.isdigit()]
		for usage in gpu_usages:
			if usage < 10:
				return True
		return False
		
	def GetAFile(self):
		try:
			self.ssh_client.connect(ssh_hostname, username='kongjf', pkey = self.rsakey)
			stdin,stdout,stderr = self.ssh_client.exec_command('ls ' + mx3_filepath)
			ftp_client = self.ssh_client.open_sftp()
			ftp_client.chdir(mx3_filepath)
		except:
			print('Error in opening connection. Will return empty string as filename.\n')
			return ''
		
		for i in stdout.readlines():
			filename = i.rstrip()
			print('Filename on server: %s\n'%filename) 
			try:
				# fileObj = ftp_client.file(filename+'.lock',mode='x')   # create lock file for concurrency
				self.ssh_client.exec_command('mv ' + mx3_filepath+filename + ' ' + mx3running_path+filename)
				time.sleep(0.5)
				ftp_client.get(mx3running_path+filename,local_mx3path+filename)
				# ftp_client.remove(mx3_filepath+filename)
				# fileObj.close()
				# ftp_client.remove(mx3_filepath+filename+'.lock')
				print('Moved %s to local drive for running.\n'%filename)
				ftp_client.close()
				return filename
			except:
				print('Error in GetAFile\n')
		ftp_client.close()
		return ''
	
	def UploadData(self,filename,delete=False):
		self.ssh_client.connect(ssh_hostname, username='kongjf', pkey = self.rsakey)
		ftp_client = self.ssh_client.open_sftp()
		remotedir = remote_datapath+filename[0:-3]+'out'
		ftp_client.mkdir(remotedir)
		self.ssh_client.exec_command('chmod g+r ' + remotedir) 
		
		localdir = local_mx3path+filename[0:-3]+'out\\'
		listoffiles = os.listdir(localdir)
		for file in listoffiles:
			ftp_client.put(localdir+file,remotedir+'/'+file)
			self.ssh_client.exec_command('chmod g+r ' + remotedir+'/'+file) 
			if delete:
				os.remove(localdir+file)
			
		if delete:
			os.rmdir(localdir)
			os.remove(local_mx3path+filename)
			
		ftp_client.close()

		


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