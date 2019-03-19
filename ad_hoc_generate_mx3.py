"""
This script reads in a template text file and generates multiple mx3 text files.
"""

import textwrap
import os
import re
import random as rand

def FORC_cont_temp():

	local_run = True

	# Hr_list in mT
	Hr_list = list(range(-200, 0, 2))

	# sweep rate in T/s
	B_sweep_rate = 1e6
	B_step_size = 0.002
	# important to sweep pass zero, so that gradient at zero could be found accurately
	B_end = 0.02

	# ovf path
	ovf_path = r'D:\Skyrmions-data\FORC\Baby forc\st128'
	mx3_output_path = r'D:\Skyrmions-data\FORC\Baby forc\st129\mx3 files'
	ovf_file_prefix = 'xy_forc_st128_e309808aa0c543aba158bce2d7b33371_'

	sim_name = '33dFORC_st129'

	filename_list = os.listdir(ovf_path)

	if local_run:
		from simple_job_server import submit_local_job
	else:
		from one_sim import SimulationParameters, SimulationMetadata, TuningParameters, writing_sh, submit_sh
		loop_ind = 0


	for filename in filename_list:

		# Capture the Hr in mT
		search = re.search('full_mag_.+_(\+?-?\d+)mT\.ovf', filename)

		if search is None:
			continue

		Hr = int(search.group(1))

		if not Hr in Hr_list:
			continue

		# n is the ovf ind
		# n = int(round((MH_loop_start-Hr)/H_step))

		# ovf_file = os.path.join(ovf_path, filename)
		ovf_file = ovf_file_prefix+filename
		B_start = float(Hr)/1e3

		mx3_filename = sim_name + '_Hr_{0:04d}mT.mx3'.format(Hr)

		mumax_file_str = os.path.join(mx3_output_path, mx3_filename)

		if not local_run:
			sim_params = SimulationParameters(

				sim_meta=SimulationMetadata(
					sim_name_full=mx3_filename,
					walltime='24:00:00',
					project_code='13000385',
					output_dir=mx3_output_path,
					output_subdir=mx3_output_path,
					mumax_file=mumax_file_str,
					stage=71,
					loop=loop_ind,
					local_run=False),

				tune=TuningParameters(
					forc_run=False,
					m_h_loop_run=False
				))

			loop_ind+=1

		mumax_commands = textwrap.dedent('''\
		Mega :=1e6
		Pico :=1e-12
		Nano :=1e-9
		Mili :=1e-3
		
		// Micromagnetic variables
		Aex_var := 4.533333333*Pico  // Exchange in J/m^3
		Msat_var := 0.316519345*Mega  //Saturation magnetisation in A/m
		Ku1_var	:= 0.075923784*Mega  // Anistropy in J/m^3
		Dbulk_var  := 0.000000*Mili  //Bulk DMI in J/m^2
		Dind_var  := 0.62*Mili  //Interfacial DMI in J/m^2
				
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
		AnisU = vector(0, 0, 1) //Uniaxial anisotropy direction 	
		
		// Physical size
		size_X	:=1536.000000 //sim_param.phy_size.x
		size_Y	:=1536.000000
		size_Z	:=42.000000
		
		// Total number of simulations cells
		Nx	:=384 //sim_param.grid_size.x
		Ny	:=384
		Nz	:=14
		
		// PBC, if any
		PBC_x :=0 //sim_param.pbc.x
		PBC_y :=0
		PBC_z :=0
		
		z_single_rep_thickness := 1 // thickness of single repetition in number of cells in z
		z_layer_rep_num := 14 //this many repetitions
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
		ext_scaleExchange(1, 2, 0.000000)	
		
		TableAdd(B_ext)
		TableAdd(E_Total)
		TableAdd(E_anis)
		TableAdd(E_demag)
		TableAdd(E_exch)
		TableAdd(E_Zeeman)
		TableAdd(ext_topologicalcharge)
		TableAdd(Temp)
		TableAdd(LastErr)
		OutputFormat = OVF1_TEXT
		
		middle_layer := 6
		
		// define some variables here which may or may not be used later
		mz := m.comp(2)
		
		// constant temperature and damping throughout simulation
		Temp = 850.0
		alpha  =0.10000		 // Damping
		
		SetSolver(2) // Solver for run with thermal fluctuations	
		ThermSeed(%d) // Set a random seed for thermal noise 
		FixDt = 1.500000E-13
		
		// initialise with +z uniform mag since M(H) loop with start at saturation
		m.LoadFile(`%s`)
		
		// Sweeping B field
		B_start := %f
		B_end := %f
		// sweep rate in Tesla/s
		// -1e6 => 1mT/ns => 5mT/5ns
		B_sweep_rate := %E
		// take steps in field, say 2mT
		B_step_size := %f
		// field total num of steps
		B_total_steps := Ceil((B_end-B_start)/B_step_size)+1
		
		// define some variables here. Note these are all floats
		B_current := 0.0
		B_current_mT := 0.0
		t0 := 0.0
		
		// sim time per step in s
		sim_run_time_per_step := B_step_size / B_sweep_rate
		const_field_run_time := 1e-9
		
		// set autosaves	
		// in 0.1 mT step. typical: 0.01ns
		tableautosave_time := 1e-5/Abs(B_sweep_rate)
		tableautosave(tableautosave_time)
		
		// print values to check
		prt_str := sprintf("==> B_total_steps: %%f",B_total_steps)
		print(prt_str)
		prt_str = sprintf("==> sim_run_time_per_step: %%e", sim_run_time_per_step)
		print(prt_str)
		
		// Big loop
		for ind := 0; ind < B_total_steps; ind++ {
			
			// initial stablisation run at constant field
			B_current = B_start + B_step_size*ind
			B_ext = vector(0,0,B_current)
			
			// print values to check
			prt_str = sprintf("==> B_current: %%f. Running const_field_run_time",B_current)
			print(prt_str)	
			
			Run(const_field_run_time)
			
			// output mag	
			B_current_mT = B_current*1e3
			saveas(CropLayer(m, middle_layer), sprintf("%%04d_sliced_mag_%s_%%.0fmT", ind, B_current_mT ) ) 
			//saveas(m, sprintf("%%04d_full_mag_%s_%%.0fmT", ind, B_current_mT) )
			
			// sweep to next field
			t0 = t
			B_ext = vector(0, 0, B_current + B_sweep_rate*(t-t0))
			
			// print values to check
			prt_str = sprintf("==> B_current: %%f. Sweeping to next field value",B_current)
			print(prt_str)	
			
			Run(sim_run_time_per_step)
		}
		
		// save final ovf file
		saveas(m, sprintf("final_full_mag_%s_%%.0fmT", B_current_mT))
	
		''' %(rand.randrange(0, 2 ** 32), ovf_file, B_start, B_end, B_sweep_rate, B_step_size, sim_name, sim_name, sim_name))

		mumax_file = open(mumax_file_str, "w")
		mumax_file.write(mumax_commands)
		mumax_file.close()

		# submit job to local server
		# if local_run:
		# 	submit_local_job(mumax_file_str)
		# else:
		# 	writing_sh(sim_params)
			#submit_sh(sim_params)

def relax_many():

	from one_sim import SimulationParameters, SimulationMetadata, TuningParameters,  writing_sh, submit_sh

	# ovf path
	ovf_path = r'/scratch/users/astar/imre/chenxy14/Micromagnetics/FORC/mumax_sim_outputs/11dFORC_st65/11dFORC_st65.out'
	mx3_path = r'/scratch/users/astar/imre/chenxy14/Micromagnetics/FORC/mumax_sim_outputs/11dFORC_st71'
	# ovf_path = r'D:\Skyrmions-data\FORC\Baby forc\st61-65\st65\full ovf'
	# mx3_path = r'C:\Users\Xiaoye\Tmp\st71 test'
	sim_name = '11dFORC_st71'

	# find all ovf files in ovf_path
	ovf_list = [file for file in os.listdir(ovf_path) if file.endswith('.ovf')]

	for ovf_file in ovf_list:
		# abs path and filename
		ovf_full_name = os.path.join(ovf_path, ovf_file)
		# index of ovf file
		n_ind = int(re.search('m(\d+)\.ovf', ovf_file).group(1))
		# calc field from n index
		Bfield = 0.250 - 0.005*float(n_ind)

		Bfield_str = 'H{0:04d}mT'.format(int(round(Bfield*1e3)))
		mx3_filename = sim_name + '_' + Bfield_str + '.mx3'

		mumax_file_full = os.path.join(mx3_path, mx3_filename)

		sim_params = SimulationParameters(
			sim_meta=SimulationMetadata(
				sim_name_full=sim_name + '_' + Bfield_str,
				walltime = '24:00:00',
				project_code='13000385',
				output_dir=mx3_path,
				output_subdir=mx3_path,
				mumax_file=mumax_file_full,
				stage=71,
				loop=n_ind,
				local_run=False
			),
			tune=TuningParameters(
				forc_run=False,
				m_h_loop_run=False
			)
		)

		mumax_commands = textwrap.dedent('''\
			Mega :=1e6
			Pico :=1e-12
			Nano :=1e-9
			Mili :=1e-3

			// Micromagnetic variables
			Aex_var := 3.999960*Pico  // Exchange in J/m^3
			Msat_var := 0.367330*Mega  //Saturation magnetisation in A/m
			Ku1_var	:= 0.101102*Mega  // Anistropy in J/m^3
			Dbulk_var  := 0.000000*Mili  //Bulk DMI in J/m^2
			Dind_var  := 0.733326*Mili  //Interfacial DMI in J/m^2

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
			AnisU = vector(0, 0, 1) //Uniaxial anisotropy direction 	

			// Physical size
			size_X	:=1536.000000 //sim_param.phy_size.x
			size_Y	:=1536.000000
			size_Z	:=60.000000

			// Total number of simulations cells
			Nx	:=384 //sim_param.grid_size.x
			Ny	:=384
			Nz	:=20

			// PBC, if any
			PBC_x :=0 //sim_param.pbc.x
			PBC_y :=0
			PBC_z :=0

			z_single_rep_thickness := 1 // thickness of single repetition in number of cells in z
			z_layer_rep_num := 20 //this many repetitions
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
			ext_scaleExchange(1, 2, 0.000000)	

			TableAdd(B_ext)
			TableAdd(E_Total)
			TableAdd(E_anis)
			TableAdd(E_demag)
			TableAdd(E_exch)
			TableAdd(E_Zeeman)
			TableAdd(ext_topologicalcharge)
			TableAdd(Temp)
			OutputFormat = OVF1_TEXT

			middle_layer := 9

			// define some variables here which may or may not be used later
			mz := m.comp(2)
			
			B_ext = vector(0, 0, %f) //in Teslas

			// initialise with +z uniform mag since M(H) loop with start at saturation
			m.LoadFile(`%s`)
			
			// save mag before and after relax
			tablesave()

			// relax
			relax()
			
			tablesave()
			
			// save output
			// save only the middle layer
			saveas(CropLayer(m, middle_layer),"%s") 
			// output final ovf to be loaded by the next run
			saveas(m,"%s")

				''' % (Bfield, ovf_full_name,
					   'sliced_mag_relaxed_'+Bfield_str,'full_mag_relaxed_'+Bfield_str))

		mumax_file = open(mumax_file_full, "w")
		mumax_file.write(mumax_commands)
		mumax_file.close()

		# submit job to NSCC server
		writing_sh(sim_params)
		# submit_sh(sim_params)


if __name__ == '__main__':
	# FORC for continuous temp
	FORC_cont_temp()
	# relax individual temp results
	# relax_many()