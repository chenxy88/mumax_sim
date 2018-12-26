import os
import subprocess
import textwrap
from dataclasses import dataclass
import random
import string
import json
import math
from copy import deepcopy
import random as rand

from simple_job_server import submit_local_job

# This branch creates a list of simulations that calculates the M(H) loop of the sample
# For the first job, starting from random magnetisation, the system goes to the first magnetic field in the array and relaxes
# The magnetisation is saved for subsequent jobs. The magnetisation of the middle slide is saved to output.
# Subsequent jobs will load the magnetisation of the previous job, go to the next field value and relax.

# All other input parameters are expected to be single values instead of arrays

# constants
@dataclass
class Constants:
	full_mag_relaxed_str: str = 'full_mag_relaxed_'
	sliced_mag_relaxed_str: str = 'sliced_mag_relaxed_'
	sliced_mag_after_temp_str: str = 'sliced_mag_after_temp_'
	config_ovf_str: str = 'config.ovf'
	Bfield_mT_str: str = '_%.0fmT'
	forc_Hr_str: str = 'Hr_%.0fmT'


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
			update_obj_from_dict_recursively(tmp_obj, tmp_v)
		else:
			some_obj.__dict__.update({k:tmp_v})

def convert_obj_to_json_recursively(obj):

	# first check if obj has __dict__
	assert(hasattr(obj, '__dict__'))
	obj_dict = obj.__dict__

	json_str = '{'
	for ind,key in enumerate(obj_dict):
		if ind != 0:
			json_str+=','

		val = obj_dict[key]
		# if val is an object
		if hasattr(val, '__dict__'):
			json_str +=  '"' + key+'":' + convert_obj_to_json_recursively(val)
		else:
			json_str += '"' + key+'":' + json.dumps(val)

	json_str += '}'
	return json_str

# this function is the heart of one_sim
def outer_product_object_list(obj, start_ind = 0):
	# given an object which may contain some lists, return a list of objects that
	# holds individuals items of the list instead
	# if the object contains multiple lists, a list of objects equivalent to the outerproduct of
	# those lists are returned
	# this works with nest objects and but does not open/flatten nested list

	assert (hasattr(obj, '__dict__'))
	obj_dict = obj.__dict__
	dict_list = list(enumerate(obj_dict))

	# scan through each member of the object
	for ind in range(start_ind, len(dict_list)):

		key = dict_list[ind][1]
		val = obj_dict[key]

		# if found a nested object
		if hasattr(val, '__dict__'):
			# lets open up this object and explore

			sub_obj = outer_product_object_list(val)
			# if it turns out that this object contain lists
			if isinstance(sub_obj, list):
				# flatten any nested lists of sub objects
				sub_obj = flatten(sub_obj)
				result = [None] * len(sub_obj)

				for (ind_i, val_i) in enumerate(sub_obj):
					obj_tmp = deepcopy(obj)
					setattr(obj_tmp, key, val_i)
					result[ind_i] = outer_product_object_list(obj_tmp, ind+1) #plus one here is impt

				return result
			#
			# if the object does not contain any lists, don't bother with it
			else:
				pass

		# if found a list, stop and flatten list
		elif isinstance(val, list):

			result = [None]*len(val)

			for (ind_i, val_i) in enumerate(val):
				obj_tmp = deepcopy(obj)
				# replace list with a single value
				setattr(obj_tmp, key, val_i)
				result[ind_i] = outer_product_object_list(obj_tmp, ind+1) #plus one -> do not open up nested list

			return  result

	# if no objects or list are found, just return this object
	return obj

def flatten(some_list):
	"""Given a list, possibly nested to any level, return it flattened."""
	new_list = []
	for item in some_list:
		if isinstance(item, list):
			new_list.extend(flatten(item))
		else:
			new_list.append(item)
	return new_list



def load_json_file(input_file):
	"""
	load input json file
	:return:
	"""
	# root = tk.Tk()
	# input_file = filedialog.askopenfilenames(title='Select json file')

	data_loaded = None

	if input_file == '':
		return data_loaded

	converted_path = os.path.relpath(input_file)
	# if the user did selected something
	with open(converted_path) as data_file:
		data_loaded = json.load(data_file)

	return data_loaded

def find_FORC_prev_mag (path: str, ending: str):

	filename_list = os.listdir(path)


	# filter based on criteria
	for filename in filename_list:
		if filename.startswith(Constants.full_mag_relaxed_str) and  filename.endswith(ending):
			return filename

	return None

# simple helper classes
@dataclass
class Vector:
	x: float = 0
	y: float = 0
	z: float = 0

# material parameters
@dataclass
class MaterialParameters:
	landau_damping: float = 0.1
	mag_sat: float = 1 # magnetisation saturation in MA/m
	exchange: float = 10  # exchange in pJ/m
	dmi_bulk: float = 0 # bulk DMI mJ/m^2
	dmi_interface: float = 2 # interfacial DMI mJ/m^2
	anistropy_uni: float = 0.5  # 1st order uniaxial anisotropy constant, MJ/m^3. NOTE: This is NOT effective anistropy, which takes into account dipolar energy
	interlayer_exchange: float = 0 # interlayer exchange scaling

	def calc_effective_medium (self, prescaled_mat_params, scaling: float):
		# scale the effective medium, which includes the FM and spacer layers
		mu0 = 4e-7 * math.pi

		# no scaling
		self.landau_damping = prescaled_mat_params.landau_damping
		self.interlayer_exchange = prescaled_mat_params.interlayer_exchange

		# scaling according to Woo, S. et al. Nature Materials 15, 501 (2016).
		self.mag_sat = prescaled_mat_params.mag_sat * scaling
		self.exchange = prescaled_mat_params.exchange * scaling
		self.dmi_bulk = prescaled_mat_params.dmi_bulk * scaling
		self.dmi_interface = prescaled_mat_params.dmi_interface * scaling
		self.anistropy_uni = scaling*(prescaled_mat_params.anistropy_uni - mu0*1e6*prescaled_mat_params.mag_sat**2/2) + mu0*1e6*self.mag_sat**2/2


@dataclass
class GeometryParameter:
	z_fm_single_thickness: int = 1 # thickness in number of cells in z
	z_nm_single_thickness: int = 2 # thickness in number of cells in z
	z_layer_rep_num: int = 1 # number of repetitions
	z_single_rep_thickness: int = 3  # thickness of single repetition in number of cells in z

	z_cell_size: float = 1 # physical size of single cell in z, in nm

	phy_size: Vector = Vector(2048, 2048, 0) # in nm
	grid_cell_count: Vector = Vector(512, 512) # number of cells in simulation, should be power of 2
	pbc: Vector = Vector(0, 0, 0)

	effective_medium_scaling: float = 1 # scaling of material parameters = t_FM/t_repeat

	def calc_auto_parameters(self):
		self.z_single_rep_thickness = self.z_fm_single_thickness + self.z_nm_single_thickness
		self.grid_cell_count.z = self.z_single_rep_thickness * self.z_layer_rep_num - self.z_nm_single_thickness
		self.phy_size.z = self.grid_cell_count.z * self.z_cell_size


@dataclass
class SimulationMetadata:
	commit: str = '' # version control
	stage: int = 0
	loop: int = -1
	loop_start: int = 0
	sim_name: str = ''
	sim_name_full: str = ''
	walltime: str = '24:00:00'
	sim_id: str = ''
	output_dir: str = ''
	output_subdir: str = ''
	mumax_file: str = ''
	previous_jobid: str = ''
	config_ovf_path: str = ''
	config_ovf_name: str = Constants.config_ovf_str #'config.ovf'

	sh_file: str = ''
	# True: really doing simulations, False: debugging this script. This is set automatically by checking if not_production_run.txt exist in current dir
	production_run: bool = False
	# True: Running on local machine, False: Running on cluster e.g. NSCC
	local_run: bool = False
	# mumax_installed: bool = False
	project_code: str = 'Personal'

	def calc_auto_parameters(self, autoset_output_subdir:bool = True):

		if autoset_output_subdir:
			if self.output_dir == '':
				self.output_dir = os.path.join(os.getcwd(), self.sim_name)
			# convert to system specific paths
			# e.g. /scratch/users/astar/dsi/chenxy1/textures/mumax_sim_outputs
			self.output_dir = os.path.abspath(self.output_dir)
			self.output_subdir = os.path.join(self.output_dir, (self.sim_name + '_st%02d' % self.stage))

		self.sim_id = ''.join([random.choice(string.ascii_letters+string.digits) for ch in range(8)])
		self.sim_name_full = self.sim_name + '_st%02d_%02d' % (self.stage, self.loop)+'_'+self.sim_id

		# output subdirectory contains all the results from this series of M(H)
		# e.g. /scratch/users/astar/dsi/chenxy1/textures/mumax_sim_outputs/26m1_ar_st64
		self.mumax_file = os.path.join(self.output_subdir, self.sim_name_full + '.mx3')
		#self.sh_file = os.path.join(self.output_subdir, 'one_sim.sh')

		self.production_run = not os.path.isfile('./not_production_run.txt')  # this will be set automatically by checking if file not_production_run.txt exist in current dir
		# status, outputstr = subprocess.getstatusoutput('mumax3')
		# self.mumax_installed = status == 0

@dataclass
class TuningParameters:
	external_Bfield: float = 0
	# whether or not the first run of the M(H) loop starts with previous mag. Useful for minor loops.
	start_series_with_prev_mag: bool = False
	thermal_fluctuation: bool = False
	uniform_mag_initial: bool = True
	# whether or not previous config is copied out to common, and next run starts with prev config
	m_h_loop_run: bool = True
	# number of field points to calculate per run
	m_h_loop_points_per_run: int = 10
	temperature: float = 900
	# str is used in order to avoid outer_product of array
	# temperature_annealing_temp: str = ''
	# temperature_annealing_time: float = 1e-9
	temperature_run_time: float = 5e-10
	temperature_decay_time: float = 1e-9
	temperature_run_dt: float = 1e-15
	temperature_solver: int = 2
	mag_autosave_period: float = 0 # zero disable autosave
	table_autosave_period: float = 1e-11 # 100 to 1000 points for 1ns to 10ns run
	# forc related
	forc_run: bool = False
	# reverse field in T
	forc_Hr: [float] = 0
	# stopping field for minor loop, in T
	forc_H_stop: float = 0
	# field step in minor loop from Hr to Hstop, in T
	forc_H_step: float = 0.01
	# path of mag profile to load
	forc_starting_mag_path: str = ''

# all experimental parameters
@dataclass
class SimulationParameters:

	# SimulationMetadata
	sim_meta: SimulationMetadata = SimulationMetadata()
	# MaterialParameters
	mat: MaterialParameters = MaterialParameters()
	# MaterialParameters
	mat_scaled: MaterialParameters = MaterialParameters()
	# GeometryParameter
	geom: GeometryParameter = GeometryParameter()
	# TuningParameters
	tune: TuningParameters = TuningParameters()

	def generate_sims(self):
		# this is suppose to generate a parameter space defined by parameter arrays that are assumed to be orthogonal
		result = [outer_product_object_list(self)]
		return flatten(result)


def save_json_file(sim_param: SimulationParameters):

	if not os.path.exists(sim_param.sim_meta.output_subdir):
		os.makedirs(sim_param.sim_meta.output_subdir)

	json_file = open(os.path.join(sim_param.sim_meta.output_subdir, sim_param.sim_meta.sim_name_full +".json"), "w")
	json_str = convert_obj_to_json_recursively(sim_param)
	parsed = json.loads(json_str)
	json.dump(parsed, json_file, sort_keys=True, indent=2)

# Write PBS script for submission to queue
def writing_sh(sim_params: SimulationParameters, last_sim: bool = True):

	server_script = textwrap.dedent('''\
	#!/bin/bash
	#PBS -N %s
	#PBS -q gpu
	#PBS -l Walltime=%s
	#PBS -l select=1:ncpus=24:mem=24GB
	#PBS -P %s
	
	module load mumax
	
	''' % (sim_params.sim_meta.sim_name_full, sim_params.sim_meta.walltime, sim_params.sim_meta.project_code))

	# .out subsubdir autocreated by mumax
	out_subsubdir = os.path.join(sim_params.sim_meta.output_subdir, sim_params.sim_meta.sim_name_full + '.out')

	# move out the table and ovf (after relax) for easy harvesting
	server_script += textwrap.dedent('''\
	mumax3 %s
	mv -f  %s %s 
	''' % (sim_params.sim_meta.mumax_file,
		   os.path.join(out_subsubdir, '*.ovf'),
		   sim_params.sim_meta.output_subdir))

	# set special table copy destination for FORC, since there is an extra layer of folder

	if sim_params.tune.forc_run:
		table_dest_path_name = os.path.join(sim_params.sim_meta.output_dir,
											# (sim_params.sim_meta.sim_name + '_st%02d' % sim_params.sim_meta.stage),
											sim_params.sim_meta.sim_name_full + '_' + Constants.forc_Hr_str % (sim_params.tune.forc_Hr * 1e3)+ '.txt')
	else:
		table_dest_path_name = os.path.join(sim_params.sim_meta.output_subdir, sim_params.sim_meta.sim_name_full + '.txt')

	server_script += textwrap.dedent('''\
	cp -f  %s %s
	'''%(os.path.join(out_subsubdir, 'table.txt'), table_dest_path_name))

	# set up job chainning if this is an m_h loop
	if sim_params.tune.m_h_loop_run and not last_sim:
		# submit the next job
		next_sh_file = os.path.join(sim_params.sim_meta.output_subdir, 'one_sim_st%02d_%02d.sh' % (sim_params.sim_meta.stage, sim_params.sim_meta.loop+1))
		server_script = server_script + textwrap.dedent('''\
			qsub   %s
			''' % next_sh_file)

	# defining the location of the .mx3 script
	sim_params.sim_meta.sh_file = os.path.join(sim_params.sim_meta.output_subdir, 'one_sim_st%02d_%02d.sh' % (sim_params.sim_meta.stage, sim_params.sim_meta.loop))

	# write batch sh script if this is a cluster run
	if not sim_params.sim_meta.local_run:
		# opening and saving sh file
		executable_file = open(sim_params.sim_meta.sh_file, "w")
		executable_file.write(server_script)
		executable_file.close()

	else:
		# submit a local job to local server
		submit_local_job(sim_params.sim_meta.mumax_file)

# submit the job to PBS
def submit_sh(sim_params: SimulationParameters):
	# does not apply for local run
	if not sim_params.sim_meta.local_run and os.path.isfile(sim_params.sim_meta.sh_file):
		qsub_params = ['qsub',
					   '-o', sim_params.sim_meta.output_subdir,
					   '-e', sim_params.sim_meta.output_subdir,
					   sim_params.sim_meta.sh_file]

		# Run command and wait for it to complete
		subprocess.run(qsub_params) #, stdout=subprocess.PIPE).stdout.decode('utf-8')

# write mumax3 script
def writing_mumax_file(sim_params: SimulationParameters):

	middle_layer = (math.ceil(sim_params.geom.z_layer_rep_num / 2) - 1) * sim_params.geom.z_single_rep_thickness

	mumax_commands = textwrap.dedent('''\
	Mega :=1e6
	Pico :=1e-12
	Nano :=1e-9
	Mili :=1e-3
	
	// Micromagnetic variables
	Aex_var := %f*Pico  // Exchange in J/m^3
	Msat_var := %f*Mega  //Saturation magnetisation in A/m
	Ku1_var	:= %f*Mega  // Anistropy in J/m^3
	Dbulk_var  := %f*Mili  //Bulk DMI in J/m^2
	Dind_var  := %f*Mili  //Interfacial DMI in J/m^2
	
	// Setting micromagnetic parameters for Region 0: Non-magnetic
	Aex.SetRegion(0, 10*Pico)
	Msat.SetRegion(0, 0)
	Ku1.SetRegion(0, 0)
	Dbulk.SetRegion(0, 0)
	Dind.SetRegion(0, 0)
	
	// Setting micromagnetic parameters for Region 1: FM 1
	Aex.SetRegion(1, Aex_var)
	Msat.SetRegion(1, Msat_var)
	Ku1.SetRegion(1, Ku1_var)
	Dbulk.SetRegion(1, Dbulk_var)
	Dind.SetRegion(1, Dind_var)
	
	// Setting micromagnetic parameters for Region 2: FM 2
	Aex.SetRegion(2, Aex_var)
	Msat.SetRegion(2, Msat_var)
	Ku1.SetRegion(2, Ku1_var)
	Dbulk.SetRegion(2, Dbulk_var)
	Dind.SetRegion(2, Dind_var)
	
	// Micromagnetic parameters for all regions
	alpha  =%f		 // Damping
	AnisU = vector(0, 0, 1) //Uniaxial anisotropy direction 	
	
	// Physical size
	size_X	:=%f //sim_param.phy_size.x
	size_Y	:=%f
	size_Z	:=%f
	
	// Total number of simulations cells
	Nx	:=%.0f //sim_param.grid_size.x
	Ny	:=%.0f
	Nz	:=%.0f
	
	// PBC, if any
	PBC_x :=%.0f //sim_param.pbc.x
	PBC_y :=%.0f
	PBC_z :=%.0f

	z_single_rep_thickness := %.0f // thickness of single repetition in number of cells in z
	z_layer_rep_num := %.0f //this many repetitions
	num_of_regions := 2 // number of FM regions (exclude NM region 0)

	//SetPBC(PBC_x, PBC_y, PBC_z)
	SetGridsize(Nx, Ny, Nz)
	SetCellsize(size_X*Nano/Nx, size_Y*Nano/Ny, size_Z*Nano/Nz)

	//geometry
	for layer_number:=0; layer_number<Nz; layer_number+= z_single_rep_thickness {
		// set adjacent layers to be of different regions
		// so that we could set interlayer exchange coupling 
		// layer 1: FM1, layer 2: FM2
		defRegion(Mod(layer_number, num_of_regions)+1, layer(layer_number))
	}
	
	// interlayer exchange scaling
	ext_scaleExchange(1, 2, %f)	

	TableAdd(B_ext)
	TableAdd(E_Total)
	TableAdd(E_anis)
	TableAdd(E_demag)
	TableAdd(E_exch)
	TableAdd(E_Zeeman)
	TableAdd(ext_topologicalcharge)
	TableAdd(Temp)
	OutputFormat = OVF1_TEXT
	
	middle_layer := %d
	
	// save the middle slice of the config	
	AutoSave(CropLayer(m, middle_layer), %E) 
	tableautosave(%E)
	
	// define some variables here which may or may not be used later
	temperature_run_time := %E	
	mz := m.comp(2)
	// important to declare this as a float
	t_start := 0.0
	temperature_decay_time := %E
	temperature_val := %f
	
	''' % (sim_params.mat_scaled.exchange, sim_params.mat_scaled.mag_sat, sim_params.mat_scaled.anistropy_uni, sim_params.mat_scaled.dmi_bulk, sim_params.mat_scaled.dmi_interface,

		   sim_params.mat_scaled.landau_damping,

		   sim_params.geom.phy_size.x, sim_params.geom.phy_size.y, sim_params.geom.phy_size.z,

		   sim_params.geom.grid_cell_count.x, sim_params.geom.grid_cell_count.y, sim_params.geom.grid_cell_count.z,

		   sim_params.geom.pbc.x, sim_params.geom.pbc.y, sim_params.geom.pbc.z,

		   sim_params.geom.z_single_rep_thickness, sim_params.geom.z_layer_rep_num,

		   sim_params.mat_scaled.interlayer_exchange,

		   middle_layer,

		   sim_params.tune.mag_autosave_period, sim_params.tune.table_autosave_period,

		   sim_params.tune.temperature_run_time,

		   sim_params.tune.temperature_decay_time,

		   sim_params.tune.temperature))

	# set initial magnetisation
	if (sim_params.sim_meta.loop == 0 and not sim_params.tune.start_series_with_prev_mag) or not sim_params.tune.m_h_loop_run:

		if sim_params.tune.uniform_mag_initial:
			# start with uniform magnetisation for the first field
			mumax_commands += textwrap.dedent('''\
			// initialise with +z uniform mag since M(H) loop with start at saturation
			m.setRegion(0, Uniform(0, 0, 0))
			m.setRegion(1, Uniform(0, 0, 1))
			m.setRegion(2, Uniform(0, 0, 1))

			''')

		else:
			# start with random magnetisation for the first field
			mumax_commands += textwrap.dedent('''\
			// initialise with random magnetisation
			m.setRegion(0, Uniform(0, 0, 0))
			m.setRegion(1, RandomMagSeed(%d))
			m.setRegion(2, RandomMagSeed(%d))

			''' %(rand.randrange(0,2**32), rand.randrange(0,2**32)))

	else:
		# load previous magnetisation for subsequent fields
		# if config_ovf_path is not set, load from output_subdir
		if sim_params.sim_meta.config_ovf_path == '':
			sim_params.sim_meta.config_ovf_path = sim_params.sim_meta.output_subdir

		mumax_commands += textwrap.dedent('''\
		// load the previous magnetisation
		m.LoadFile(`%s`)
		
		''' % os.path.join(sim_params.sim_meta.config_ovf_path,
						   sim_params.sim_meta.config_ovf_name))

	# set magnetic field if external_Bfield is not a list
	if not isinstance(sim_params.tune.external_Bfield, list):
		mumax_commands = mumax_commands + textwrap.dedent('''\
		B_ext = vector(0, 0, %f) //in Teslas
		''' % sim_params.tune.external_Bfield)

	run_and_relax_commands = ''

	if sim_params.tune.m_h_loop_run:
		for ind, Bfield in enumerate(sim_params.tune.external_Bfield):

			# set Bfield for each field point
			run_and_relax_commands += textwrap.dedent('''\
			
			//---------------Set new field point for MH loop---------------
			B_ext = vector(0, 0, %f) //in Teslas
			'''%Bfield)

			# counter str is the loop number, followed by the field in mT
			counter_str =  Constants.Bfield_mT_str %(Bfield*1e3)

			# run thermal fluctuations
			if sim_params.tune.thermal_fluctuation:
				run_and_relax_commands += run_thermal_fluctuations_commands(sim_params, counter_str)

			# relax
			run_and_relax_commands += relax_commands(sim_params, counter_str)

		run_and_relax_commands += textwrap.dedent('''\

		// output final ovf to be loaded by the next run
		saveas(m,"%s")

		''' % Constants.config_ovf_str)

	else:
		run_and_relax_commands = relax_commands(sim_params)

	mumax_commands += run_and_relax_commands
	# opening and saving it
	mumax_file = open(sim_params.sim_meta.mumax_file, "w")
	mumax_file.write(mumax_commands)
	mumax_file.close()


def relax_commands(sim_params: SimulationParameters, counter:str = ''):
	commands = textwrap.dedent('''\
	
	// change back to normal settings
	SetSolver(3) // back to default solver for relax
	FixDt = 0 // turn off fixed time step
	Temp = 0 // turn off temperature
	relax()			// high-energy states best minimized by relax()
	
	// save only the middle layer
	saveas(CropLayer(m, middle_layer),"%s") 
	tablesave()
	
	''' % (Constants.sliced_mag_relaxed_str+sim_params.sim_meta.sim_name_full+counter))

	# do not save full mag for FORC run, will be way too large
	if not sim_params.tune.forc_run:
		commands += textwrap.dedent('''\
		saveas(m,"%s")	
		'''% (Constants.full_mag_relaxed_str+sim_params.sim_meta.sim_name_full+counter))
	return commands

def run_thermal_fluctuations_commands(sim_params: SimulationParameters, counter:str = ''):

	thermal_run_commands = textwrap.dedent('''\
	
	// apply a short burst of thermal fluctuations to allow the system to cross small energy barriers
	SetSolver(%d) // Solver for run with thermal fluctuations	
	ThermSeed(%d) // Set a random seed for thermal noise 
	
	FixDt = %E
	// record the start time
	t_start = t
	
	Temp = temperature_val*Exp(-(t-t_start)/temperature_decay_time)
	
	Run(temperature_run_time)		

	''' % (sim_params.tune.temperature_solver, rand.randrange(0, 2 ** 32),
		   sim_params.tune.temperature_run_dt))

	thermal_run_commands += textwrap.dedent('''\
	// save only the middle layer
	saveas(CropLayer(m, middle_layer),"%s") 

	''' % (Constants.sliced_mag_after_temp_str +sim_params.sim_meta.sim_name_full + counter))

	return thermal_run_commands

# generates list of sim_params based on input and then generate mx3 files
def process_sim_param (sim_params: SimulationParameters, forc_Hr_mag_name: str = ''):

	# introduce an extra layer of folder for FORC simulations
	if sim_params.tune.forc_run:

		sim_params.sim_meta.output_dir = sim_params.sim_meta.output_subdir
		sim_params.sim_meta.output_subdir = os.path.join(sim_params.sim_meta.output_subdir,
														 Constants.forc_Hr_str % (sim_params.tune.forc_Hr * 1e3) )
		sim_params.sim_meta.mumax_file = os.path.join(sim_params.sim_meta.output_subdir,
													  sim_params.sim_meta.sim_name_full + '.mx3')

	# if mh loop run, convert list of magnetic fields into a list of lists of magnetic fields
	# outer_product_object_list will not flatten nested list
	if sim_params.tune.m_h_loop_run:
		tmp_list_of_lists = []
		len_field_arr = len(sim_params.tune.external_Bfield)
		for i in range(int(math.ceil(len_field_arr / sim_params.tune.m_h_loop_points_per_run))):
			# segment of external_Bfield to put in nested list
			ind_start = i * sim_params.tune.m_h_loop_points_per_run
			ind_end = (i + 1) * sim_params.tune.m_h_loop_points_per_run
			# check for out of bounds
			if ind_end > len_field_arr:
				ind_end = len_field_arr
			# append to nested list
			tmp_list_of_lists.append(sim_params.tune.external_Bfield[ind_start:ind_end])

		sim_params.tune.external_Bfield = tmp_list_of_lists

	# generate list of simulations based on main input
	sim_params_list = sim_params.generate_sims()
	sim_params_list_len = len(sim_params_list)

	# this will check that the only array is that of external_Bfield
	if sim_params.tune.m_h_loop_run and len(sim_params_list) != len(sim_params.tune.external_Bfield):
		raise ValueError('There should not be other arrays than external_Bfield for M(H) loop simulations!')

	# if FORC, set the first sim_param in list to load from previous external mag
	if sim_params.tune.forc_run:
		sim_params_list[0].sim_meta.config_ovf_path = sim_params.tune.forc_starting_mag_path
		sim_params_list[0].sim_meta.config_ovf_name = forc_Hr_mag_name

	# create the subdirectory if it does not exist
	if not os.path.exists(sim_params.sim_meta.output_subdir):
		os.makedirs(sim_params.sim_meta.output_subdir)

	# save the main input json for record
	save_json_file(sim_params)

	# previous job id
	prev_jobid = ''
	# prev_sim_param = None

	# generate .mx3 and .sh files
	for loop, sim_params_i in enumerate(sim_params_list):
		# allows continuing from previous stage at a certain loop
		sim_params_i.sim_meta.loop = loop + sim_params.sim_meta.loop_start
		sim_params_i.sim_meta.previous_jobid = prev_jobid

		# generate new names and calc geometry
		sim_params_i.sim_meta.calc_auto_parameters(autoset_output_subdir=not sim_params.tune.forc_run)
		sim_params_i.geom.calc_auto_parameters()

		# do effective medium scaling
		sim_params_i.mat_scaled.calc_effective_medium(sim_params_i.mat, sim_params_i.geom.effective_medium_scaling)

		# save each individual simulation's json params
		save_json_file(sim_params_i)

		# these are written to the subdirectory
		writing_mumax_file(sim_params_i)
		writing_sh(sim_params_i, last_sim=loop == sim_params_list_len - 1)

		# save the prev_jobid
		prev_jobid = sim_params_i.sim_meta.previous_jobid
		# prev_sim_param = sim_params_i

	# submit .sh files to queue
	if sim_params.sim_meta.production_run:
		# check if this is a m_h loop run
		if not sim_params.tune.m_h_loop_run:
			# submit all jobs at once
			for loop, sim_params_i in enumerate(sim_params_list):
				submit_sh(sim_params_i)
		# m_h loop run, submit just the first job, and the jobs will chain
		else:
			submit_sh(sim_params_list[0])

#--------------- main ---------------#

def main():
	# --- CHOOSE INPUT OPTION --- #
	input_dict = load_json_file('../mumax_sim_inputs/input_parameters.json')

	sim_params = SimulationParameters()
	update_obj_from_dict_recursively(sim_params, input_dict)

	# generate simulation full name
	sim_params.sim_meta.calc_auto_parameters()

	# FORC logic here
	if sim_params.tune.forc_run:
		save_json_file(sim_params)
		# the only array should be forc_Hr
		# generate based on forc_Hr first
		sim_params_list_forc = sim_params.generate_sims()
		# this will check that the only array is that of forc_Hr
		if len(sim_params_list_forc) != len(sim_params.tune.forc_Hr):
			raise ValueError('There should not be arrays other than forc_Hr for FORC simulations!')

		# start main loop for FORC here
		# loop_start = 0
		for sim_params_Hr in sim_params_list_forc:
			sim_params_Hr.tune.external_Bfield = [x/1e3 for x in range(int (sim_params_Hr.tune.forc_Hr*1e3),
															int(sim_params_Hr.tune.forc_H_stop*1e3),
															int(sim_params_Hr.tune.forc_H_step*1e3))] + [sim_params_Hr.tune.forc_H_stop]
			# exclude the first point
			sim_params_Hr.tune.external_Bfield = sim_params_Hr.tune.external_Bfield[1:]
			# set certain compulsory settings
			sim_params_Hr.tune.start_series_with_prev_mag = True
			sim_params_Hr.tune.m_h_loop_run = True
			# search and set config_ovf_name
			forc_Hr_mag_name = find_FORC_prev_mag(sim_params_Hr.tune.forc_starting_mag_path,
													 ending=Constants.Bfield_mT_str %(sim_params_Hr.tune.forc_Hr*1e3)+'.ovf')
			if forc_Hr_mag_name is None:
				# prev mag not found
				raise ValueError('Previous magnetisation not found')

			# generate .mx3, write batch script and submit job
			process_sim_param(sim_params_Hr, forc_Hr_mag_name)

	else:
		# generate .mx3, write batch script and submit job
		process_sim_param(sim_params)


# run the main function
if __name__ == '__main__':
	main()
