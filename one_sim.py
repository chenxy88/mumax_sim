# everything is in IS
# say not to outdated not standard units of measurementent such as CGS
# eV is cool though
import os
import copy
import numpy as np
from textwrap import *
import random as rand
from dataclasses import dataclass


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

# simulation parameters

class params_2d:
	def __init__(self, x, y):
		self.x = x
		self.y = y


class params_3d(params_2d):
	def __init__(self, x, y, z):
		params_2d.__init__(self, x, y)
		self.z = z

class TextureType:
	neel = 'Neel'
	bloch = 'Bloch'

class DMI_strength:
	def __init__(self, bulk, interface):
		self.bulk = bulk
		self.interface = interface

# material parameters
class material_parameters:
	def __init__(self, landau_damping, mag_sat, exchange, dmi, anistropy_uni):
		self.landau_damping = landau_damping
		self.mag_sat = mag_sat  # magnetisation saturation
		self.exchange = exchange  # pJ/m
		self.dmi = dmi  # bulk/interfacial DMI mJ/m^2
		self.anistropy_uni = anistropy_uni  # 1st order uniaxial anisotropy constant, MJ/m^3


# all experimental parameters
@dataclass
class simulation_parameters:
	mat_1 = material_parameters(0.1, 0.648, 12.62, DMI_strength(0, 1.67), 0.488)  # 11a
	# landau_damping, mag_sat, exchange, (dmi_bulk, dmi_int), anistropy_uni

	skyrmion_size = 0  # in nm, diameter of the completely in-plane magnetisation boundary
	skyrmion_rand_pos_delta = params_2d(50, 75)
	skyrmion_scaling_factor = None
	texture_type = TextureType.neel
	skyrmion_chirality = 1
	num_sk_per_side = params_2d(4, 3)
	sk_spacing = params_2d(200, 250)  # in nm
	external_Bfield = 0.0  # in T
	pbc = params_3d(0, 0, 10)  # number of extra images, on both the positive and negative sides (e.g. 5 means 5 on +ve side and 5 on -ve side)

	# repeat structure in z
	z_fm_single_thickness = 1  # in nm, equals to magnetic layer thickness
	z_nm_single_thickness = 2  # in nm, equals to non-magnetic layer thickness
	z_layer_rep_num = 1  # this many repeats

	cell_size_z = 1  # in nm
	z_single_rep_thickness = z_fm_single_thickness + z_nm_single_thickness  # thickness of single repetition
	phy_size_z = z_single_rep_thickness * z_layer_rep_num
	grid_size_z = round(phy_size_z / cell_size_z)

	phy_size = params_3d(2048, 2048, phy_size_z)  # in nm
	grid_size = params_3d(1024, 1024, grid_size_z)  # number of cells in simulation, should be power of 2
	mesh_cell_size = None

	stage = 8
	loop = 0

	sim_name = 'chens_labyrinth'
	sim_name_full = ''

	walltime = '24:00:00'

	# AUTO SET FIELDS
	home_dir = os.path.join(os.getcwd(), sim_name)
	working_dir = ''
	results_dir_Data = ''
	results_dir_Plots = ''
	mumax_file = ''
	sh_file = ''
	production_run = not os.path.isfile(
		'./not_production_run.txt')  # this will be set automatically by checking if file not_production_run.txt exist in current dir

	update_sim_params()

# all experimental parameters
class simulation_parameters:
	def __init__(self):
		self.mat_1 = material_parameters(0.1, 0.648, 12.62, DMI_strength(0, 1.67), 0.488) #11a
		# landau_damping, mag_sat, exchange, (dmi_bulk, dmi_int), anistropy_uni

		self.skyrmion_size = 0  # in nm, diameter of the completely in-plane magnetisation boundary
		self.skyrmion_rand_pos_delta = params_2d(50,75)
		self.skyrmion_scaling_factor = None
		self.texture_type = TextureType.neel
		self.skyrmion_chirality = 1
		self.num_sk_per_side = params_2d(4, 3)
		self.sk_spacing = params_2d(200, 250)  # in nm
		self.external_Bfield = 0.0 # in T
		self.pbc = params_3d(0,0,10) #number of extra images, on both the positive and negative sides (e.g. 5 means 5 on +ve side and 5 on -ve side)

		# repeat structure in z
		self.z_fm_single_thickness = 1  # in nm, equals to magnetic layer thickness
		self.z_nm_single_thickness = 2  # in nm, equals to non-magnetic layer thickness
		self.z_layer_rep_num = 1  # this many repeats

		self.cell_size_z = 1  # in nm
		self.z_single_rep_thickness = self.z_fm_single_thickness + self.z_nm_single_thickness  # thickness of single repetition
		self.phy_size_z = self.z_single_rep_thickness * self.z_layer_rep_num
		self.grid_size_z = round(self.phy_size_z / self.cell_size_z)

		self.phy_size = params_3d(2048, 2048, self.phy_size_z)  # in nm
		self.grid_size = params_3d(1024, 1024, self.grid_size_z)  # number of cells in simulation, should be power of 2
		self.mesh_cell_size = None

		# stage 01: Try to achieve a stable labyrinthine state at zero field using PBC and a starting condition of randomly placed skyrmions
		# stage 02: Try to achieve a stable labyrinthine state at zero field using PBC and a starting condition of randomly placed skyrmions. DBulk instead of DInt.
		# stage 03: Setting DBulk only without setting DInt produced the runtime error "DMI: Msat, Aex, Dex must be uniform". Try setting DInt to zero. Also, streamlined the code to remove unnecessary material region definitions.
		# stage 04: The previous stage produced stripes at zero field, but since they evolved from skyrmions, the number of stripes domains are limited to the number of starting skyrmions (12). The simulation space is too small and confinement effect from the square geometry, too strong. Going to increase the simulation domain to 2048*2048 nm and start with random magnetisation. RESULTS: Dense forest of stripes and skymrions of a single chirality.
		# stage 05: Try to use experimental parameters of Co(0.8)/Pt(1)X20 of 11a, but with D = 0. RESULTS: Uniform magnetisation at 0.1 and 0.2 mT. At 0.05 mT, it is almost uniform, except for 2 domains, one of which is an antiskyrmion with a topological number of 13. Large patches of domains at 0 mT.
		# stage 06: Try a small interfacial DMI of 0.1-0.2. Bulk DMI should be zero in thin-films. Small fields of 0.05 to 0.1. 
		# stage 07: Try smaller fields of 0.01-0.025 T and larger interfacial DMI of 0.3 to 0.8. Also set random magnetisation with starting seed. RESULTS: Mostly Neel texture that has occassional nodes or kinks of Blochness. Topologically, they are still antiskyrmions. 
		# stage 08: Add in a little bulk DMI to try to stabilise Bloch textures. Keep DInt at 0.5 and external field at 0.015T.
		self.stage = 8
		self.loop = 0

		self.sim_name = 'chens_labyrinth'
		self.sim_name_full = ''

		self.walltime = '24:00:00'

		# AUTO SET FIELDS
		self.home_dir = os.path.join(os.getcwd(), self.sim_name)
		self.working_dir = ''
		self.results_dir_Data = ''
		self.results_dir_Plots = ''
		self.mumax_file = ''
		self.sh_file = ''
		self.production_run = not os.path.isfile('./not_production_run.txt') # this will be set automatically by checking if file not_production_run.txt exist in current dir

		self.update_sim_params()

	def update_sim_params(self):
		self.mesh_cell_size = params_3d(self.phy_size.x/self.grid_size.x, self.phy_size.y/self.grid_size.y, self.phy_size.z/self.grid_size.z)
		# 0.83255461115769769 is sqrt(ln(2))
		self.skyrmion_scaling_factor = self.skyrmion_size/(8*self.mesh_cell_size.x)/2/0.83255461115769769

		self.sim_name_full = self.sim_name + '_st%02d_%02d' %(self.stage, self.loop)
		self.mumax_file = os.path.join(self.home_dir, self.sim_name_full + ".mx3")
		self.sh_file = os.path.join(self.home_dir, 'ServerGPU.sh')

# parameters to loop through for a series of simulations
loop_params = simulation_parameters()
#----------TO EDIT FOR DIFFERENT SIMULATIONS----------#
loop_params.mat_1.exchange = [12.62]
loop_params.mat_1.dmi = [DMI_strength(0.25,0.5), DMI_strength(0.5,0.5), DMI_strength(0.75,0.5), DMI_strength(1,0.5), DMI_strength(1.5,0.5)]
loop_params.skyrmion_size = [50] # in nm
loop_params.texture_type = [TextureType.bloch]
loop_params.external_Bfield = [0.015]


# Write PBS script for submission to queue
def writing_sh(sim_param):
	server_script = dedent("""
	#!/bin/bash
	#PBS -N one_sim
	#PBS -q gpu
	#PBS -l Walltime=%s
	#PBS -l select=1:ncpus=24:mem=24GB
	#PBS -P Personal
	module load mumax

	mumax3 %s
	""" % (sim_param.walltime, sim_param.mumax_file)  # 00:00:30
						   )

	# defining the location of the .mx3 script
	executable = sim_param.sh_file

	# opening and saving it
	executable_file = open(executable, "w")
	executable_file.write(server_script)
	executable_file.close()

	return 0

# qsub the job
def submit_sh(sim_param):
	if os.path.isfile(sim_param.sh_file):
		os.system('qsub < %s ' % sim_param.sh_file)

	return 0

def run_n_convert_mumax(sim_param):
	if os.path.isfile(sim_param.mumax_file):
		os.system('mumax3 %s ' % sim_param.mumax_file)
		# move output file to current directory
		os.system('mv %s/*.ovf %s' % (os.path.join(sim_param.home_dir, sim_param.sim_name_full+ '.out'), sim_param.home_dir))
		# os.system('mumax3-convert -vtk binary %s/*.ovf ' % (sim_param.sim_name_full + '.out'))

# write mumax3 script
def writting_mx3(sim_param):
	if not os.path.exists(sim_param.working_dir):
		os.makedirs(sim_param.working_dir)
	print(sim_param.working_dir)

	if not os.path.exists(sim_param.working_dir):
		print(sim_param.working_dir + " doesn't exist")

	geometry_single_layer_cuboid = 'single_layer := cuboid(size_X*Nano, size_Y*Nano, z_fm_single_thickness*Nano)'

	geometry_single_layer_cylinder = 'single_layer := cylinder(size_X*Nano, z_fm_single_thickness*Nano)'

	geometry = dedent("""
		%s
		rep_layer := single_layer.repeat(0, 0, z_single_rep_thickness*Nano) // repeat once every this many times
		if Mod(z_layer_rep_num, 2) == 0{ //if even number of magnetic layers, translate a little bit
			rep_layer = rep_layer.transl(0,0,z_fm_single_thickness*Nano)
		}
		setgeom(rep_layer)
		//saveas(geom, "stack_geom")
	""" % (geometry_single_layer_cuboid))

	magnetic_layers = dedent("""
	nx_sk_per_side:= %.0f
	ny_sk_per_side:= %.0f
	sk_spacing_x := %f*Nano
	sk_spacing_y := %f*Nano
	sk_scale := %f
	sk_chirality := %d

	m = RandomMagSeed(%d)
	//m = uniform(0, 0, 1)
	xpos:= 0.0
	ypos:= 0.0
	""" % (sim_param.num_sk_per_side.x, sim_param.num_sk_per_side.y, 
		sim_param.sk_spacing.x, sim_param.sk_spacing.y, 
		sim_param.skyrmion_scaling_factor, sim_param.skyrmion_chirality,
		rand.randrange(0,2**32)))


	mumax_commands = dedent("""
	Mega :=1e6
	Pico :=1e-12
	Nano :=1e-9
	Mili :=1e-3

	alpha  =%f		 // Damping
	Aex = %f*Pico  // Exchange in J/m^3
	Msat = %f*Mega  //Saturation magnetisation in A/m

	Dbulk  = %f*Mili  //Bulk DMI in J/m^2
	Dind  = %f*Mili  //Interfacial DMI in J/m^2
	K1	:=%f*Mega  // Anistropy in J/m^3
	B_Max :=%f		 // BZ in T

	size_X	:=%f //sim_param.phy_size.x
	size_Y	:=%f
	size_Z	:=%f
	sk_size :=%f //sim_param.skyrmion_size, in nm, diameter of the completely in-plane magnetisation boundary

	Nx	:=%.0f //sim_param.grid_size.x
	Ny	:=%.0f
	Nz	:=%.0f
	
	PBC_x :=%.0f //sim_param.pbc.x
	PBC_y :=%.0f
	PBC_z :=%.0f

	z_fm_single_thickness := %f //in nm, equals to cell thickness
	z_single_rep_thickness := %f //in nm
	z_layer_rep_num := %.0f //this many repeats

	SetPBC(PBC_x, PBC_y, PBC_z)
	SetGridsize(Nx, Ny, Nz)
	SetCellsize(size_X*Nano/Nx, size_Y*Nano/Ny, size_Z*Nano/Nz)

	%s //geometry, use shape to define geometry: one big box, the same size as the simulation area

	
	Ku1	= K1
	AnisU = vector(0, 0, 1) //Uniaxial anisotropy direction 

	%s //magnetic_layers, M0_layers

	B_ext = vector(0, 0, B_Max) //in Teslas

	TableAdd(B_ext)
	TableAdd(E_Total)
	tableAdd(ext_topologicalcharge)
	OutputFormat = OVF1_TEXT
	saveas(m,"%s")
	""" % (sim_param.mat_1.landau_damping, sim_param.mat_1.exchange, sim_param.mat_1.mag_sat, sim_param.mat_1.dmi.bulk,  sim_param.mat_1.dmi.interface, sim_param.mat_1.anistropy_uni, sim_param.external_Bfield, 
		sim_param.phy_size.x, sim_param.phy_size.y,  sim_param.phy_size.z, 
		sim_param.skyrmion_size, 
		sim_param.grid_size.x, sim_param.grid_size.y,sim_param.grid_size.z, 
		sim_param.pbc.x, sim_param.pbc.y, sim_param.pbc.z, 
		sim_param.z_fm_single_thickness, sim_param.z_single_rep_thickness, sim_param.z_layer_rep_num,
		geometry, magnetic_layers, sim_param.sim_name_full + '_init'))

	if sim_param.production_run is True:
		mumax_commands = mumax_commands + dedent("""
			tablesave()
			MinimizerStop = 1e-6
			relax()			// high-energy states best minimized by relax()
			saveas(m,"%s")
			tablesave()
		""" % (sim_param.sim_name_full + '_relaxed.ovf'))

	# defining the location of the .mx3 script
	executable = os.path.join(sim_param.working_dir, sim_param.sim_name_full + ".mx3")

	# opening and saving it
	executable_file = open(executable, "w")
	executable_file.write(mumax_commands)
	executable_file.close()

	return 0


###################### complete sims ######################
###################### complete sims ######################
###################### complete sims ######################

def main():

	sim_param_i = simulation_parameters()
	sim_param_list = []

	if not os.path.exists(sim_param_i.home_dir):
		os.makedirs(sim_param_i.home_dir)

	parameter_list = []
	working_dir_list = []

	loop=0
	for Ex in loop_params.mat_1.exchange:  # pJ/m
		for D in loop_params.mat_1.dmi:
			for B_Max in loop_params.external_Bfield:
				for sk_size in loop_params.skyrmion_size:
					for sk_type in loop_params.texture_type:
						# set simulation params for single simulation
						sim_param_i.mat_1.exchange = Ex
						sim_param_i.mat_1.dmi = D
						sim_param_i.external_Bfield = B_Max
						sim_param_i.working_dir = sim_param_i.home_dir
						sim_param_i.skyrmion_size = sk_size
						sim_param_i.texture_type = sk_type
						sim_param_i.loop = loop

						# update sim name
						sim_param_i.update_sim_params()

						sim_param_list.append(copy.deepcopy(sim_param_i))

						loop = loop + 1

	for sim_param_i in sim_param_list:

		writting_mx3(sim_param_i)
		# TODO: Add control flow variable, and move production_run to control flow
		if sim_param_i.production_run is True:
			# the real deal, write sh scripts
			writing_sh(sim_param_i)
			submit_sh(sim_param_i)
		else:
			run_n_convert_mumax(sim_param_i)


# run the main function
if __name__ == '__main__':
	main()


