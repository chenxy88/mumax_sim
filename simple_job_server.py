# Modified from: https://jeffknupp.com/blog/2014/02/11/a-celerylike-python-task-queue-in-55-lines-of-code/
# and https://medium.com/@shashwat_ds/a-tiny-multi-threaded-job-queue-in-30-lines-of-python-a344c3f3f7f0

"""Broker-less distributed task queue."""

import zmq
import queue, threading
import subprocess
import time

HOST_PORT = '127.0.0.1:5678'
TASK_SOCKET = zmq.Context().socket(zmq.REQ)
TASK_SOCKET.connect('tcp://' + HOST_PORT)
KILL_STR = '#kill'

class Server(object):
	"""A remote task executor."""

	def __init__(self, host_port_in=HOST_PORT):
		"""Initialize Server."""
		self.q = queue.Queue()
		self.host_port = host_port_in
		self._context = zmq.Context()
		self._socket = self._context.socket(zmq.REP)
		self.num_GPUs = 1
		self.threads = []

		# start workers
		for i in range(self.num_GPUs):
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
			print('SJS received job %s.\n'%input_string)
			if input_string == KILL_STR:
				# empty queue
				self.q.empty()

				# kill workers
				for i in range(self.num_GPUs):
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

	def worker(self, num):
		print('Worker %d started.\n' %num)

		while True:
			# get from job queue
			input_string = self.q.get()
			print('GPU %d received job %s.\n' % (num,input_string))
			#  convert to tuple
			kwargs = {}

			if input_string == KILL_STR:
				break

			else:
				self._do_work(input_string, num, **kwargs)
				print('GPU %d completed job %s.\n' % (num,input_string))

		print('GPU %d stopped.\n' % num)

	def _do_work(self, mumax_file_str, gpu_number, **kwargs):
		"""Return the result of executing the given task."""
		# sleep for a short while to ensure than the mx3 file has been written to disk
		time.sleep(0.5)
		subprocess.run(['mumax3','-gpu', '%d'%gpu_number, mumax_file_str])


def submit_local_job(input_str):
	"""Return the result of running the task *runnable* with the given
	arguments."""
	TASK_SOCKET.send_pyobj(input_str)
	results = TASK_SOCKET.recv_pyobj()
	return results

def kill():
	# kill the server by sending the KILL_STR
	submit_local_job(KILL_STR)

if __name__ == '__main__':
	w = Server()
	w.run_server()