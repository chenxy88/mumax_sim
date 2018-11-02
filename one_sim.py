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

# This branch creates a list of simulations that calculates the M(H) loop of the sample
# For the first job, starting from random magnetisation, the system goes to the first magnetic field in the array and relaxes
# The magnetisation is saved for subsequent jobs. The magnetisation of the middle slide is saved to output.
# Subsequent jobs will load the magnetisation of the previous job, go to the next field value and relax.

# All other input parameters are expected to be single values instead of arrays

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
		if (hasattr(val, '__dict__')):
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
	# this works with nest objects and but does NOT work with NESTED LISTS by design (let's not use nested lists)

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
					result[ind_i] = outer_product_object_list(obj_tmp, start_ind+1) #plus one here is impt

				return result

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
				result[ind_i] = outer_product_object_list(obj_tmp, start_ind+1) #plus one -> do not open up nested list

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
	config_ovf_name: str = 'config.ovf'

	sh_file: str = ''
	production_run: bool = False
	mumax_installed: bool = False
	project_code: string = 'Personal'

	def calc_auto_parameters(self):
		if self.output_dir == '':
			self.output_dir = os.path.join(os.getcwd(), self.sim_name)
		# convert to system specific paths
		# e.g. /scratch/users/astar/dsi/chenxy1/textures/mumax_sim_outputs
		self.output_dir = os.path.abspath(self.output_dir)

		self.sim_id = ''.join([random.choice(string.ascii_letters+string.digits) for ch in range(8)])
		self.sim_name_full = self.sim_name + '_st%02d_%02d' % (self.stage, self.loop)+'_'+self.sim_id

		# output subdirectory contains all the results from this series of M(H)
		# e.g. /scratch/users/astar/dsi/chenxy1/textures/mumax_sim_outputs/26m1_ar_st64
		self.output_subdir = os.path.join(self.output_dir, (self.sim_name + '_st%02d' % self.stage))
		self.mumax_file = os.path.join(self.output_subdir, self.sim_name_full + '.mx3')
		self.sh_file = os.path.join(self.output_subdir, 'one_sim.sh')

		self.production_run = not os.path.isfile('./not_production_run.txt')  # this will be set automatically by checking if file not_production_run.txt exist in current dir
		status, outputstr = subprocess.getstatusoutput('mumax3')
		self.mumax_installed = status == 0

@dataclass
class TuningParameters:
	external_Bfield: float = 0
	# whether or not the first run of the M(H) loop starts with previous mag. Useful for minor loops.
	start_series_with_prev_mag: bool = False
	thermal_fluctuation: bool = False
	uniform_mag_initial: bool = True
	# whether or not previous config is copied out to common, and next run starts with prev config
	m_h_loop_run: bool = True
	temperature: float = 300
	temperature_run_time: float = 5e-10
	temperature_run_dt: float = 1e-15
	temperature_stop_mz: float = 0
	mag_autosave_period: float = 0 # zero disable autosave
	table_autosave_period: float = 1e-11 # 100 to 1000 points for 1ns to 10ns run

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
def writing_sh(sim_param: SimulationParameters, prev_sim_param: SimulationParameters):
	server_script = textwrap.dedent('''\
	#!/bin/bash
	#PBS -N %s
	#PBS -q gpu
	#PBS -l Walltime=%s
	#PBS -l select=1:ncpus=24:mem=24GB
	#PBS -P %s
	module load mumax
	''' % (sim_param.sim_meta.sim_name_full, sim_param.sim_meta.walltime, sim_param.sim_meta.project_code))

	# move the previous config file out
	if sim_param.sim_meta.loop != 0 and sim_param.tune.m_h_loop_run:
		prev_output_config = os.path.join(prev_sim_param.sim_meta.output_subdir, prev_sim_param.sim_meta.sim_name_full + '.out',prev_sim_param.sim_meta.config_ovf_name)
		server_script = server_script + textwrap.dedent('''\
		
		cp -f %s %s
		''' % (prev_output_config, sim_param.sim_meta.output_subdir))

	server_script = server_script+ textwrap.dedent('''\
	mumax3 %s
	cp -f %s %s
	''' % (sim_param.sim_meta.mumax_file,
		   os.path.join(prev_sim_param.sim_meta.output_subdir,'table.txt'),
		   os.path.join(prev_sim_param.sim_meta.output_dir, prev_sim_param.sim_meta.sim_name_full +'_table.txt') ))
	# defining the location of the .mx3 script
	sh_file = sim_param.sim_meta.sh_file

	# opening and saving it
	executable_file = open(sh_file, "w")
	executable_file.write(server_script)
	executable_file.close()

	return 0

# submit the job to PBS
def submit_sh(sim_param: SimulationParameters):
	if os.path.isfile(sim_param.sim_meta.sh_file):
		if sim_param.sim_meta.loop == 0 or not sim_param.tune.m_h_loop_run:
			# first step
			qsub_params = ['qsub',
						   '-o',sim_param.sim_meta.output_subdir,
						   '-e',sim_param.sim_meta.output_subdir,
						   sim_param.sim_meta.sh_file]
			jobid_str = subprocess.run(qsub_params, stdout=subprocess.PIPE).stdout.decode('utf-8')
			# the previous_jobid is saved to sim_param, which is then extracted in the main loop
			sim_param.sim_meta.previous_jobid = jobid_str.split('.')[0]

		else:
			# subsequent steps will only start upon completion of previous
			qsub_params = ['qsub',
						   '-o', sim_param.sim_meta.output_subdir,
						   '-e', sim_param.sim_meta.output_subdir,
						   '-W','depend=afterok:%s'%sim_param.sim_meta.previous_jobid,
						   sim_param.sim_meta.sh_file]
			jobid_str = subprocess.run(qsub_params, stdout=subprocess.PIPE).stdout.decode('utf-8')
			sim_param.sim_meta.previous_jobid = jobid_str.split('.')[0]
	return 0

def run_n_convert_mumax(sim_param: SimulationParameters):
	if os.path.isfile(sim_param.sim_meta.mumax_file):
		os.system('mumax3 %s ' % sim_param.sim_meta.mumax_file)
		# move output file to current directory
		os.system('mv %s/*.ovf %s' % (os.path.join(sim_param.sim_meta.output_dir, sim_param.sim_meta.sim_name_full + '.out'), sim_param.sim_meta.output_dir))
		# os.system('mumax3-convert -vtk binary %s/*.ovf ' % (sim_param.sim_name_full + '.out'))

# write mumax3 script
def writing_mumax_file(sim_param: SimulationParameters):

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
	B_ext = vector(0, 0, %f) //in Teslas
	
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
	tableAdd(ext_topologicalcharge)
	OutputFormat = OVF1_TEXT
	
	''' % (sim_param.mat_scaled.exchange, sim_param.mat_scaled.mag_sat, sim_param.mat_scaled.anistropy_uni, sim_param.mat_scaled.dmi_bulk, sim_param.mat_scaled.dmi_interface,

		   sim_param.mat_scaled.landau_damping, sim_param.tune.external_Bfield,

		   sim_param.geom.phy_size.x, sim_param.geom.phy_size.y, sim_param.geom.phy_size.z,

		   sim_param.geom.grid_cell_count.x, sim_param.geom.grid_cell_count.y, sim_param.geom.grid_cell_count.z,

		   sim_param.geom.pbc.x, sim_param.geom.pbc.y, sim_param.geom.pbc.z,

		   sim_param.geom.z_single_rep_thickness, sim_param.geom.z_layer_rep_num,

		   sim_param.mat_scaled.interlayer_exchange))

	if (sim_param.sim_meta.loop == 0 and not sim_param.tune.start_series_with_prev_mag) or not sim_param.tune.m_h_loop_run:

		if sim_param.tune.uniform_mag_initial:
			# start with uniform magnetisation for the first field
			mumax_commands = mumax_commands + textwrap.dedent('''\
			// initialise with +z uniform mag since M(H) loop with start at saturation
			m.setRegion(0, Uniform(0, 0, 0))
			m.setRegion(1, Uniform(0, 0, 1))
			m.setRegion(2, Uniform(0, 0, 1))
			tablesave()

			''')

		else:
			# start with random magnetisation for the first field
			mumax_commands = mumax_commands+textwrap.dedent('''\
			// initialise with random magnetisation
			m.setRegion(0, Uniform(0, 0, 0))
			m.setRegion(1, RandomMagSeed(%d))
			m.setRegion(2, RandomMagSeed(%d))
			tablesave()			

			''' %(rand.randrange(0,2**32), rand.randrange(0,2**32)))

	else:
		# load previous magnetisation for subsequent fields
		mumax_commands = mumax_commands + textwrap.dedent('''\
		// load the previous magnetisation
		m.LoadFile("%s")
		tablesave()
		
		''' %(os.path.join(sim_param.sim_meta.output_subdir, sim_param.sim_meta.config_ovf_name)))

	middle_layer = (math.ceil(sim_param.geom.z_layer_rep_num / 2) - 1) * sim_param.geom.z_single_rep_thickness

	if sim_param.tune.thermal_fluctuation:
		mumax_commands = mumax_commands + textwrap.dedent('''\
		// apply a short burst of thermal fluctuations to allow the system to cross small energy barriers
		SetSolver(2) // Heun
		ThermSeed(%d) // Set a random seed for thermal noise 
		FixDt = %E
		Temp = %f
		temperature_run_time := %E		
		
		middle_layer := %d
		// save the middle slice of the config	
		AutoSave(CropLayer(m, middle_layer), %E) 
		tableautosave(%E)
		
		// decide whether to use autostop condition
		// set temperature_run_time to neg to use autostop condition
		if temperature_run_time > 0	{
			Run(temperature_run_time)
		} else {
			mz := m.comp(2)
			RunWhile(mz.average() > %f) 
		}		
		
		// save only the middle layer
		saveas(CropLayer(m, middle_layer),"%s") 
		
		// change back to normal settings
		SetSolver(5) // back to default solver
		FixDt = 0 // turn off fixed time step
		Temp = 0 // turn off temperature

		''' % (rand.randrange(0,2**32),
			   sim_param.tune.temperature_run_dt,
			   sim_param.tune.temperature,
			   sim_param.tune.temperature_run_time,
			   middle_layer,
			   sim_param.tune.mag_autosave_period,
			   sim_param.tune.table_autosave_period,
			   sim_param.tune.temperature_stop_mz,
			   'after_temp_'+sim_param.sim_meta.sim_name_full))

	# if production run, relax and save m
	if sim_param.sim_meta.production_run is True:
		mumax_commands = mumax_commands + textwrap.dedent('''\
		MinimizerStop = 1e-6
		relax()			// high-energy states best minimized by relax()
		saveas(m,"%s")
		// save only the middle layer
		saveas(CropLayer(m, %d),"%s") 
		tablesave()
		''' %(sim_param.sim_meta.config_ovf_name, middle_layer, sim_param.sim_meta.sim_name_full))

	# if not production run, just save inital mag
	else:
		mumax_commands = mumax_commands + textwrap.dedent('''\
		// save only the middle layer
		saveas(CropLayer(m, middle_layer),"%s") 
		'''%(sim_param.sim_meta.sim_name_full))

	# defining the location of the .mx3 script
	# executable = os.path.join(sim_param.sim_meta.output_dir, sim_param.sim_meta.sim_name_full + ".mx3")

	# opening and saving it
	mumax_file = open(sim_param.sim_meta.mumax_file, "w")
	mumax_file.write(mumax_commands)
	mumax_file.close()

	return 0

#--------------- main ---------------#

def main():
	# --- CHOOSE INPUT OPTION --- #
	input_dict = load_json_file('../mumax_sim_inputs/input_parameters.json')

	sim_params = SimulationParameters()
	update_obj_from_dict_recursively(sim_params, input_dict)

	# generate simulation full name
	sim_params.sim_meta.calc_auto_parameters()

	# generate list of simulations based on main input
	sim_params_list = sim_params.generate_sims()

	# this will check that the only array is that of external_Bfield
	if sim_params.tune.m_h_loop_run and  len(sim_params_list) != len(sim_params.tune.external_Bfield):
		raise ValueError('There should not be other arrays than external_Bfield for M(H) loop simulations!')

	# create the subdirectory if it does not exist
	if not os.path.exists(sim_params.sim_meta.output_subdir):
		os.makedirs(sim_params.sim_meta.output_subdir)

	# save the main input json for record
	save_json_file(sim_params)

	# previous job id
	prev_jobid = ''
	prev_sim_param = None

	for loop, sim_param_i in enumerate(sim_params_list):
		# allows continuing from previous stage at a certain loop
		sim_param_i.sim_meta.loop = loop + sim_params.sim_meta.loop_start
		sim_param_i.sim_meta.previous_jobid = prev_jobid

		# generate new names and calc geometry
		sim_param_i.sim_meta.calc_auto_parameters()
		sim_param_i.geom.calc_auto_parameters()

		# do effective medium scaling
		sim_param_i.mat_scaled.calc_effective_medium(sim_param_i.mat, sim_param_i.geom.effective_medium_scaling)

		# save each individual simulation's json params
		save_json_file(sim_param_i)

		# these are written to the subdirectory
		writing_mumax_file(sim_param_i)
		writing_sh(sim_param_i, prev_sim_param)

		if sim_param_i.sim_meta.production_run:
			# the real deal, write sh scripts
			submit_sh(sim_param_i)
		else:
			if sim_param_i.sim_meta.mumax_installed:
				run_n_convert_mumax(sim_param_i)

		# save the prev_jobid
		prev_jobid = sim_param_i.sim_meta.previous_jobid
		prev_sim_param = sim_param_i

	return

# run the main function
if __name__ == '__main__':
	main()
