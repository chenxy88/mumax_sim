import subprocess, tkinter, time, queue, threading
from tkinter import filedialog

"""
Preamble
--------
This script allows the user to select a list of .mx3 files, and then run them sequentially. Compatible with multi-gpu workstations.


Usage
-----
1. Set the number of gpus and cache_path (optional)
2. Run script
3. Select mx3 files to run


Dependencies
------------
Tested on Python 3.7. No non-standard dependencies.

Good luck and code well!
"""

NUM_GPUS = 2
# optional
cache_path = r''

class Server(object):

	def __init__(self):
		"""Initialize Server."""
		self.q = queue.Queue()
		self.GPU_ids = range(NUM_GPUS)
		self.threads = []

	def start_workers(self):
		# start workers, one per GPU
		for i in self.GPU_ids:
			t = threading.Thread(target=self.worker, args=(i,))
			t.daemon = True
			t.start()
			self.threads.append(t)

	def run_server(self):
		# prompt user for simulations to run
		root = tkinter.Tk()
		input_files = filedialog.askopenfiles(title='Select mx3 files to run...')

		# close the Tkinter window
		root.destroy()

		if input_files == '':
			# user chose to cancel
			return

		for file in input_files:
			self.q.put(file.name)

		# start worker threads
		self.start_workers()

		for t in self.threads:
			t.join()

		print('All simulations completed!\n')

	def worker(self, GPU_id):
		print('GPU %d started.\n' % GPU_id)

		while not self.q.empty():
			# get from job queue
			input_string = self.q.get()
			print('GPU %d received job %s.\n' % (GPU_id, input_string))
			#  convert to tuple
			kwargs = {}

			self._do_work(input_string, GPU_id, **kwargs)
			print('GPU %d completed job %s.\n' % (GPU_id, input_string))

		print('GPU %d stopped.\n' % GPU_id)

	def _do_work(self, mumax_file_str, GPU_id, **kwargs):
		"""Return the result of executing the given task."""
		# sleep for a short while to ensure than the mx3 file has been written to disk
		time.sleep(0.5)
		if cache_path == '':
			cmd_list = ['mumax3', '-gpu', '%d' % GPU_id, mumax_file_str]
		else:
			cmd_list = ['mumax3','-cache',cache_path,'-gpu', '%d' % GPU_id, mumax_file_str]

		subprocess.run(cmd_list)


def main():
	w = Server()
	w.run_server()

if __name__ == '__main__':
	main()