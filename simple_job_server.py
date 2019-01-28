# Modified from: https://jeffknupp.com/blog/2014/02/11/a-celerylike-python-task-queue-in-55-lines-of-code/
# and https://medium.com/@shashwat_ds/a-tiny-multi-threaded-job-queue-in-30-lines-of-python-a344c3f3f7f0

"""Broker-less distributed task queue."""

import zmq
import queue, threading
import subprocess
import time
from argparse import ArgumentParser

HOST_PORT = '127.0.0.1:5679'
TASK_SOCKET = zmq.Context().socket(zmq.REQ)
TASK_SOCKET.connect('tcp://' + HOST_PORT)
KILL_STR = '#kill'

cache_path = 'D:\Xiaoye\Micromagnetics\Kernel_cache'

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
			input_string = self._socket.recv_pyobj()
			print('SJS received command %s.\n'%input_string)
			if input_string == KILL_STR:
				# empty queue
				self.q.empty()

				# kill workers
				for i in self.GPU_ids:
					self.q.put(input_string)
				self._socket.send_pyobj('')
				break

			else:
				#  put data received from socket into queue
				self.q.put(input_string)
				self._socket.send_pyobj('')

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