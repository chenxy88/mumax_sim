# everything is in IS
# say not to outdated not standard units of measurementent such as CGS
# eV is cool though
import os
import copy
import numpy as np
from textwrap import *
import random as rand


# simulation parameters

class params_2d:
	def __init__(self, x, y):
		self.x = x
		self.y = y


class params_3d(params_2d):
	def __init__(self, x, y, z):
		params_2d.__init__(self, x, y)
		self.z = z

class SkyrmionType:
	neel = 'NeelSkyrmion'
	bloch = 'BlochSkyrmion'

# material parameters
class material_parameters:
	def __init__(self, landau_damping, mag_sat, exchange, dmi, anistropy_uni):
		self.landau_damping = landau_damping
		self.mag_sat = mag_sat  # magnetisation saturation
		self.exchange = exchange  # pJ/m
		self.dmi = dmi  # mJ/m^2
		self.anistropy_uni = anistropy_uni  # 1st order uniaxial anisotropy constant, MJ/m^3


# all experimental parameters
class simulation_parameters:
	def __init__(self):
		self.mat_1 = material_parameters(0.1, 1.102, 11.36, 1.99, 0.812)
		# landau_damping, mag_sat, exchange, dmi, anistropy_uni

		self.skyrmion_size = 50  # in nm, diameter of the completely in-plane magnetisation boundary
		self.skyrmion_scaling_factor = None
		self.skyrmion_type = SkyrmionType.neel
		self.skyrmion_chirality = -1
		self.num_sk_per_side = params_2d(10,1)
		self.sk_spacing = params_2d(800, 800)  # in nm
		self.external_Bfield = 0.0

		# repeat structure in z
		self.z_fm_single_thickness = 1  # in nm, equals to magnetic layer thickness
		self.z_nm_single_thickness = 2  # in nm, equals to non-magnetic layer thickness
		self.z_layer_rep_num = 1  # this many repeats

		self.cell_size_z = 1  # in nm
		self.z_single_rep_thickness = self.z_fm_single_thickness + self.z_nm_single_thickness  # thickness of single repetition
		self.phy_size_z = self.z_single_rep_thickness * self.z_layer_rep_num
		self.grid_size_z = round(self.phy_size_z / self.cell_size_z)

		self.phy_size = params_3d(1024, 1024, self.phy_size_z)  # in nm
		self.grid_size = params_3d(2048, 2048, self.grid_size_z)  # number of cells in simulation, should be power of 2
		self.mesh_cell_size = None

		self.production_run = True  # whether this is just to test initialisation

		# stage 01: Just a bunch of random mag config, and relax.
		# stage 02: Random placement of a number of skrymions.
		self.stage = 2
		self.loop = 0

		self.sim_name = 'periodicity'
		self.sim_name_full = ''

		self.walltime = '24:00:00'

		self.home_dir = os.getcwd()
		self.working_dir = ''
		self.results_dir_Data = ''
		self.results_dir_Plots = ''
		self.exception_f = ''
		self.mumax_file = ''
		self.sh_file = ''

		self.update_sim_params()

	def update_sim_params(self):
		self.working_dir = self.home_dir
		self.mesh_cell_size = params_3d(self.phy_size.x/self.grid_size.x, self.phy_size.y/self.grid_size.y, self.phy_size.z/self.grid_size.z)
		# 0.83255461115769769 is sqrt(ln(2))
		self.skyrmion_scaling_factor = self.skyrmion_size/(8*self.mesh_cell_size.x)/2/0.83255461115769769

		self.sim_name_full = self.sim_name + '_st_%02d_loop_%02d' %(self.stage, self.loop)
		self.exception_f = os.path.join(self.home_dir, 'exception.txt')
		self.mumax_file = os.path.join(self.home_dir, self.sim_name_full + ".mx3")
		self.sh_file = os.path.join(self.home_dir, 'ServerGPU.sh')

# parameters to loop through for a series of simulations
loop_params = simulation_parameters()
#----------TO EDIT FOR DIFFERENT SIMULATIONS----------#
loop_params.mat_1.exchange = [11.36]
loop_params.mat_1.dmi = [1.99]
loop_params.skyrmion_size = [100] # in nm
loop_params.skyrmion_type = [SkyrmionType.neel]
loop_params.external_Bfield = [0.0]
# number of repeated simulations for random starting conditions
loop_params.repeated_sims = 2

# Write PBS script for submission to queue
def writting_sh(sim_param):
	server_script = dedent("""
	#!/bin/bash
	#PBS -N test0
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
		os.system('mv %s/*.ovf .' % (sim_param.sim_name_full + '.out'))
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

	M0_layers_alternate = dedent("""
	nx_sk_per_side:= %.0f
	ny_sk_per_side:= %.0f
	sk_spacing_x := %f
	sk_spacing_y := %f
	sk_scale := %f
	sk_chirality := %d

	DefRegion(1, cuboid(size_X*Nano, size_Y*Nano, size_Z*Nano))
	m = uniform(0, 0, 1)
	xpos:= 0.0
	ypos:= 0.0
	//m = RandomMagSeed(%d)
	
	""" % (sim_param.num_sk_per_side.x, sim_param.num_sk_per_side.y, sim_param.sk_spacing.x, sim_param.sk_spacing.y, sim_param.skyrmion_scaling_factor, sim_param.skyrmion_chirality, rand.randrange(0,2**32)))

	for i in range(0,sim_param.num_sk_per_side.x):
		M0_layers_alternate = M0_layers_alternate + dedent("""
			//Here comes the skyrmions
			xpos= %f*sk_spacing_x*Nano
			ypos= %f*sk_spacing_y*Nano
			m.setInShape(cylinder(sk_size*Nano*2, size_Z*Nano).transl(xpos, ypos, 0),%s(sk_chirality, -1).scale(sk_scale, sk_scale, 1).transl(xpos, ypos, 0))
		""" %(rand.random()-0.5, rand.random()-0.5, sim_param.skyrmion_type))

	mumax_commands = dedent("""

	Mega :=1e6
	Pico :=1e-12
	Nano :=1e-9
	Mili :=1e-3

	Damping  :=%f		 //
	Exchange :=%f*Pico  // in J/m^3
	Mag		:=%f*Mega  // in A/m
	D		  :=%f*Mili  // in J/m^2
	K1		 :=%f*Mega  // in J/m^3
	B_Max	 :=%f		 // BZ in T

	size_X	:=%f
	size_Y	:=%f
	size_Z	:=%f
	sk_size :=%f //in nm, diameter of the completely in-plane magnetisation boundary

	Nx	:=%.0f
	Ny	:=%.0f
	Nz	:=%.0f

	z_fm_single_thickness := %f //in nm, equals to cell thickness
	z_single_rep_thickness := %f //in nm
	z_layer_rep_num := %.0f //this many repeats

	SetGridsize(Nx, Ny, Nz)
	SetCellsize(size_X*Nano/Nx, size_Y*Nano/Ny, size_Z*Nano/Nz)

	%s //use shape to define geometry: one big box, the same size as the simulation area

	alpha = Damping
	Aex	= Exchange
	//Msat  = Mag
	Dind  = D	  //	Dbulk  = D
	//Ku1	= K1


	%s // M0_layers

	//anisU  = vector(0, 0, 1)

	Msat.SetRegion(0, 0)
	Msat.SetRegion(1, Mag)
	Ku1.SetRegion(0, 0)
	Ku1.SetRegion(1, K1)
	anisU.SetRegion(0, vector(0,0,0))
	anisU.SetRegion(1, vector(0,0,1))

	B_ext = vector(0, 0, B_Max) //in mT - doh not A/m

	TableAdd(B_ext)
	TableAdd(E_Total)
	tableAdd(ext_topologicalcharge)
	OutputFormat = OVF1_TEXT
	saveas(m.Comp(2),"%s")


	""" % (sim_param.mat_1.landau_damping, sim_param.mat_1.exchange, sim_param.mat_1.mag_sat, sim_param.mat_1.dmi, sim_param.mat_1.anistropy_uni, sim_param.external_Bfield, sim_param.phy_size.x, sim_param.phy_size.y,  sim_param.phy_size.z, sim_param.skyrmion_size, sim_param.grid_size.x, sim_param.grid_size.y,sim_param.grid_size.z, sim_param.z_fm_single_thickness, sim_param.z_single_rep_thickness, sim_param.z_layer_rep_num, geometry, M0_layers_alternate, sim_param.sim_name_full + '_before_relax'))

	if sim_param.production_run is True:
		mumax_commands = mumax_commands + dedent("""
			tablesave()
			MinimizerStop = 1e-6
			relax()			// high-energy states best minimized by relax()
			saveas(m.Comp(2),"%s")
			tablesave()
		""" % (sim_param.sim_name_full + '_after_relax.ovf'))

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


	exception_file = open(loop_params.exception_f, 'w')
	except_connection_list = []
	parameter_list = []
	working_dir_list = []

	for i in range(0,loop_params.repeated_sims):
		sim_param_i.loop = i
		# update sim name
		sim_param_i.update_sim_params()
		sim_param_list.append(copy.copy(sim_param_i))

	for sim_param_i in sim_param_list:

		writting_mx3(sim_param_i)
		# TODO: Add control flow variable, and move production_run to control flow
		if sim_param_i.production_run is True:
			# the real deal, write sh scripts
			writting_sh(sim_param_i)
			submit_sh(sim_param_i)
		else:
			run_n_convert_mumax(sim_param_i)


# run the main function
if __name__ == '__main__':
	main()


