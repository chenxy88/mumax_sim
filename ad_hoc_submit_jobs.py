"""
This script reads in a parameters file, a template sh file and a folder of mx3 files, then submit the mx3 files to NSCC server.
"""

import os, paramiko, time
import json, uuid
from copy import deepcopy
from dataclasses import dataclass, field
from tkinter import filedialog

from simple_job_server_NSCC import update_obj_from_dict_recursively, Server

@dataclass
class Sh_Replacement_Keys:
	name: str = ''
	walltime: str = ''
	project: str = ''
	mx3_file: str = ''


@dataclass
class Parameters:
	job_set_name: str = ''
	cache_path: str = ''
	local_path: str = ''
	ssh_hostname: str = 'astar.nscc.sg'
	remote_mx3_path: str = ''
	remote_username: str = ''
	rsa_key_path: str = ''
	walltime: str = '24:00:00'
	project: str = ''

	sh_replacement_keys: Sh_Replacement_Keys = Sh_Replacement_Keys()

	# used to scan an input mx3 file for keywords to replace
	replacement_dict: dict = field(default_factory=dict)


def submit_jobs_to_NSCC():

	params = Parameters()

	params_path_and_filename = filedialog.askopenfilename(title='Select parameters file')

	# Load parameters if available
	if params_path_and_filename != '' and os.path.isfile(params_path_and_filename):
		with open(params_path_and_filename) as data_file:
			params_dict = json.load(data_file)
			update_obj_from_dict_recursively(params, params_dict)

	else:
		# user cancelled, or file not found
		return

	# open local_mx3path and find all the mx3 files
	mx3_file_list = [file for file in os.listdir(params.local_path) if file.endswith('.mx3')]

	# check that mx3 files are found
	if len(mx3_file_list) == 0:
		print('No .mx3 files found in %s...\n' % params.local_path)
		return

	# open local_mx3path and find all the sh files
	sh_file_list = [file for file in os.listdir(params.local_path) if file.endswith('.sh')]

	# there should be only one sh template file
	assert (len(sh_file_list) == 1)
	sh_template_filename = sh_file_list[0]

	try:
		# setup ssh
		ssh_client = paramiko.SSHClient()
		ssh_client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
		rsakey = paramiko.RSAKey.from_private_key_file(params.rsa_key_path)

		ssh_client.connect(params.ssh_hostname, username=params.remote_username, pkey=rsakey)
		ftp_client = ssh_client.open_sftp()
		# update remote mx3 path
		params.remote_mx3_path = params.remote_mx3_path + '/' + params.job_set_name

		# check if remote folder exist
		try:
			ftp_client.chdir(params.remote_mx3_path)
			print('Remote directory %s already exist.'%params.remote_mx3_path)
		except:
			print('Remote directory %s does not exist. Making directory.'%params.remote_mx3_path)
			# create new remote folder
			ftp_client.mkdir(params.remote_mx3_path)

		# create a new folder in cache
		cache_folder_tmp = os.path.join(params.cache_path,str(uuid.uuid4()))
		os.makedirs(cache_folder_tmp)

	except Exception as Argument:
		print('Connection to remote location failed with exception: %s'%Argument)

		return

	# constructs replacement dict for sh file
	sh_replacement_dict_general = {
		params.sh_replacement_keys.name: params.job_set_name,
		params.sh_replacement_keys.project: params.project,
		params.sh_replacement_keys.walltime: params.walltime
	}

	# filters, copies to remote, then the job
	for filename in mx3_file_list:

		# filters mx3 file, and dump them in cache
		local_mx3_path_and_filename = os.path.join(cache_folder_tmp, filename)
		remote_mx3_path_and_filename = params.remote_mx3_path + '/' + filename
		Server.FilterFile(params.replacement_dict, os.path.join(params.local_path, filename), local_mx3_path_and_filename)

		# copies mx3 file to remote location
		ftp_client.put(local_mx3_path_and_filename, remote_mx3_path_and_filename)

		# constructs replacement dict for the specific sh file
		sh_replacement_dict = deepcopy(sh_replacement_dict_general)
		sh_replacement_dict[params.sh_replacement_keys.mx3_file] = remote_mx3_path_and_filename

		# filters sh file, and dump them in cache
		local_sh_path_and_filename = os.path.join(cache_folder_tmp, sh_template_filename)
		remote_sh_path_and_filename = params.remote_mx3_path + '/' + sh_template_filename
		Server.FilterFile(sh_replacement_dict, os.path.join(params.local_path, sh_template_filename), local_sh_path_and_filename)

		# copies sh file to remote location
		ftp_client.put(local_sh_path_and_filename, remote_sh_path_and_filename)

		# tries to submits the sh file using qsub
		try:
			# chan = ssh_client.get_transport().open_session()
			stdin, stdout, stderr = ssh_client.exec_command('qsub ' + remote_sh_path_and_filename)
			print('Submitted job: %s ' % filename)
			msg = stdout.channel.recv(4096).decode('ascii')
			exit_status = stdout.channel.recv_exit_status()
			print('Submitted job: %s with message: %s and exit status: %d'%(filename, msg, exit_status))

			if exit_status != 0:
				raise Exception("Error submitting job.")

		except Exception as Argument:
			print('Failed to submit job with exception: %s' % Argument)

	pass
	return


if __name__ == '__main__':
	submit_jobs_to_NSCC()
