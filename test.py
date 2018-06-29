import json
import tkinter as tk
from tkinter import filedialog
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

def UI_load_json_file():
	"""
	load input json file
	:return:
	"""
	# root = tk.Tk()
	files = filedialog.askopenfilenames(title='Select json file')

	data_loaded = []
	if files == '':
		return data_loaded

	for file in files:
		# if the user did selected something
		with open(file) as data_file:
			data_loaded.append(json.load(data_file))

	return data_loaded


@dataclass
class InputOptions:
	path: str = None
	n1: float = 0
	some_thing: Opt2 = Opt2()

@dataclass
class Opt2:
	some_string:str = None
	some_num:float = 0.0


def analysis_main():

	# --- CHOOSE INPUT OPTION --- #
	input_dict_list = UI_load_json_file()

	if not input_dict_list:
		# user chose to cancel
		return

	# --- BIG OUTSIDE LOOP THRU LIST OF INPUT JSON FILES --- #
	for input_dict in input_dict_list:
		in_opt = InputOptions()
		update_obj_from_dict_recursively (in_opt, input_dict)
		# in_opt = InputOptions(**input_dict)
		# in_opt.__dict__.update({'n1':-9.99})
		a = 3

# run the main function
if __name__ == '__main__':
	analysis_main()
	# test()
