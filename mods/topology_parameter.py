#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import re
import os



# =============== variables =============== #
RE_QUOTE = re.compile(r"[\'\"]")
RE_COMMENT = re.compile(r"^[\s\t]*;")
RE_DIRECTIVE = re.compile(r"\[ (.+) \]")
RE_POSRES = re.compile(r"([\s\t]+\d+){2}([\s\t]+\d+(?:\.\d+)?){3}")



# =============== function =============== #
def get_filepath(search_path, library, dirname=None):
	"""
	Function to find files and return their contents

	Args:
		search_path (str): search dir
		library (str): library path
		dirname (str, optional): source dir path (file with include line) (Default: None)

	Returns:
		str: parameter file path
	"""
	search_file = os.path.basename(search_path)

	library_paths = [x for x in library]
	if dirname == "":
		dirname = "."
	if dirname is not None:
		# when include source dir path, place the path at the top
		library_paths.insert(0, dirname)

	matched_files = []
	for library_path in library_paths:
		# loop for library path (search library)
		matched_files.extend([os.path.join(root, file) for root, dirs, files in os.walk(library_path, followlinks = True) for file in files if search_path in os.path.join(root, file)])
		if len(matched_files) != 0:
			break

	if len(matched_files) == 0:
		# when no file
		sys.stderr.write("ERROR: include topology file not found ({0}).\n".format(search_path))
		sys.exit(1)
	if 2 <= len(matched_files):
		# when multiple files are found
		sys.stderr.write("{0} is found in:\n".format(search_path))
		for idx, value in enumerate(matched_files):
			sys.stderr.write("{0:>3} {1}\n".format(idx, value))
		sys.stderr.write("Select file: ")
		sys.stderr.flush()

		user = int(sys.stdin.readline().strip())
		if user < len(matched_files):
			search_file = matched_files[user]
		else:
			sys.stderr.write("ERROR: Undefined value.\n")
			sys.exit(1)
	else:
		search_file = matched_files[0]

	return search_file


def clean_lines(lines):
	"""
	Function to deal with empty line

	Args:
		lines (list): line list

	Returns:
		list
	"""
	results = []

	flag_first = True
	flag_empty = False
	flag_comment = False
	idx_first = 0
	idx_end = 0
	cnt_result_idx = -1

	for line_val in lines:
		if len(line_val.strip()) == 0:
			# when empty line
			if flag_first:
				# when first empty line, skip
				idx_first = cnt_result_idx + 1
			elif flag_empty == False:
				# permit one empty line in sequential empty lines
				flag_empty = True
				results.append(line_val)
				cnt_result_idx += 1

		else:
			if RE_COMMENT.search(line_val):
				# when comment line
				if flag_comment == False:
					idx_end = cnt_result_idx + 1
					flag_comment = True
			else:
				flag_comment = False

			flag_first = False
			flag_empty = False
			results.append(line_val)
			cnt_result_idx += 1

	if flag_comment == False:
		idx_end = len(results)

	if len(results) == 0:
		return results
	else:
		if results[-1] == "":
			idx_end = -1

		return results[idx_first:idx_end]



# =============== class =============== #
class TopologyParameter:
	""" Topology Data class """
	def __init__(self, input_file=None, library=[], nmol=1, posres=[1000, 1000, 1000], prefix=""):
		self._n_mol = 0			# number of molecules
		self._library = []		# library with parameters
		self._posres = []		# position restraint energy
		self._prefix = ""		# prefix for output

		self._defaults = []			# defaults directive
		self._parameters = []		# force field parameter
		self._molecule_names = []	# molecule names
		self._molecules_info = []	# molecule information (molecule index)
		self._system_name = []		# system information
		self._system_mol = []		# molecule information in system

		# init
		self.set_library(library)
		self.set_n_mol(nmol)
		self.set_prefix(prefix)
		self.set_posres(posres)

		lines = self._load_data(input_file)
		self._parse_topology(lines)


	@property
	def defaults(self):
		return self._defaults

	@property
	def parameters(self):
		return self._parameters

	@property
	def molecule_names(self):
		return self._molecule_names

	@property
	def system_name(self):
		return self._system_name

	@property
	def system_mol(self):
		return self._system_mol


	def _load_data(self, input_file):
		"""
		Method to read file

		Args:
			input_file (str): source file

		Returns:
			None or contents of file
		"""
		if input_file is not None:
			return self._load_data_sub(input_file)
		else:
			return None


	def _load_data_sub(self, input_file):
		"""
		Method to read include line

		Args:
			input_file (str): source file

		Returns:
			list: contents
		"""
		lines = []
		with open(input_file, "r") as obj_input:
			for line_val in obj_input:
				if "#include" in line_val:
					# at include line

					# search for include file
					dirname = os.path.dirname(input_file)
					include_file = line_val.replace("#include", "")
					include_file = RE_QUOTE.sub("", include_file)
					include_file = get_filepath(include_file.strip(), self._library, dirname)
					lines.extend(self._load_data_sub(include_file))
					lines.append("\n")

				else:
					# at other lines
					lines.append(line_val)
		return lines


	def _parse_topology(self, lines):
		"""
		Method to parse parameters

		Args:
			lines (list): contents
		"""
		flag_directive = 0
		molecule_name = ""
		temporary_info = []
		flag_posres = False

		for line_val in lines:
			if RE_DIRECTIVE.search(line_val):
				# identify directive
				directive = RE_DIRECTIVE.search(line_val).group(1)

				if directive == "molecules":
					# at molecule information in system
					flag_directive = 5

				elif directive in "system":
					# at system information
					flag_directive = 4

				elif directive == "moleculetype":
					# each molecule information
					flag_directive = 3
					temporary_info = ["[ moleculetype ]\n"]

				elif directive == "defaults":
					# at defaults directive
					flag_directive = 1
					self._defaults.append(line_val)

				elif 1 <= flag_directive <= 2:
					# at force field
					flag_directive = 2
					self._parameters.append(line_val)

				elif flag_directive == 3:
					# other directives in system information
					if directive == "position_restraints":
						# 原子拘束の場合
						flag_posres = True
					else:
						flag_posres = False
					self._molecules_info[-1].append(line_val)

			elif 0 <= flag_directive <= 1:
				self._defaults.append(line_val)

			elif flag_directive == 2:
				self._parameters.append(line_val)

			elif flag_directive == 3:
				if len(temporary_info) != 0:
					# when moleculetype to newly register
					if  not RE_COMMENT.search(line_val):
						# not comment
						molecule_name = line_val.strip().split()[0]
						self._molecule_names.append(molecule_name)
						self._molecules_info.append(temporary_info)
						self._molecules_info[-1].append(line_val)
						temporary_info = []
					else:
						# comment (label) after directive
						temporary_info.append(line_val)

				else:
					# when moleculetype which has already registered
					if flag_posres and RE_POSRES.search(line_val):
						# at posres
						datas = line_val.strip().split()
						datas[2:5] = self._posres
						self._molecules_info[-1].append("{0[0]:>6} {0[1]:>5} {0[2]:>7} {0[3]:>7} {0[4]:>7}\n".format(datas))
					else:
						self._molecules_info[-1].append(line_val)

			elif flag_directive == 4:
				self._system_name.append(line_val)

			elif flag_directive == 5:
				if len(line_val.strip()) == 0:
					continue

				if RE_COMMENT.search(line_val):
					self._system_mol.append(line_val)
				else:
					mol_info = line_val.strip().split()
					mol_info[1] = int(mol_info[1]) * self._n_mol
					self._system_mol.append("{0[0]:<15} {0[1]:>7}\n".format(mol_info))

		# delete head empty line
		self._defaults = clean_lines(self._defaults)
		self._parameters = clean_lines(self._parameters)
		self._molecule_names = clean_lines(self._molecule_names)
		for idx, info in enumerate(self._molecules_info):
			self._molecules_info[idx] = clean_lines(info)
		self._system_name = clean_lines(self._system_name)
		self._system_mol = clean_lines(self._system_mol)


	def set_library(self, library):
		"""
		Method to set library path list for searching parameters

		Args:
			library (list): library path list for searching parameters

		Returns:
			self
		"""
		self._library = library
		return self


	def set_n_mol(self, n_mol):
		"""
		Method to set number of molecules

		Args:
			nmol (int): number of molecules

		Returns:
			self
		"""
		self._n_mol = n_mol
		return self


	def set_posres(self, posres):
		"""
		Method to set force energy for position restraint

		Args:
			posres (list): force energy list

		Returns:
			self
		"""
		self._posres = posres
		return self


	def set_prefix(self, prefix):
		"""
		Method to set output prefix

		Args:
			prefix (str): output prefix

		Returns:
			self
		"""
		self._prefix = prefix
		return self


	def _merge_parameters(self, parameters1, parameters2):
		"""
		Method to merge two parameters

		Args:
			parameters1 (TopologyParameter object): TopologyParameter object
			parameters2 (TopologyParameter object): TopologyParameter object

		Returns:
			list: merged parameters
		"""
		parameter_values = []
		for line in parameters1:
			if RE_DIRECTIVE.search(line):
				# at directive
				parameter_values.append([line, []])
			else:
				# at parameter line
				parameter_values[-1][1].append(line)

		directives = [v[0] for v in parameter_values]

		directive_idx = 0
		for line in parameters2:
			if RE_DIRECTIVE.search(line):
				# at directive
				if line in directives:
					# determine directive_idx
					directive_idx = directives.index(line)
				else:
					# when directive is not registered
					directive_idx += 1
					directives.insert(directive_idx, line)
					parameter_values.insert(directive_idx, [line, []])
			else:
				# at parameter line
				if line not in parameter_values[directive_idx][1]:
					# when parameter is not exist, register parameter (remove duplicate)
					parameter_values[directive_idx][1].append(line)

		results = []
		for parameters in parameter_values:
			results.append(parameters[0])
			results += parameters[1]

		return results


	def merge_topology(self, obj_molecule):
		"""
		Method to uptake other topology file (additional molecule)

		Args:
			obj_molecule (Molecule object): Molecule object

		Returns:
			self
		"""
		self._parameters = self._merge_parameters(self._parameters, obj_molecule.parameters)
		self._molecule_names.extend(obj_molecule.molecule_names)
		self._molecules_info.extend(obj_molecule.get_molecules_infos())
		for line in obj_molecule.system_mol:
			if not RE_COMMENT.search(line):
				self._system_mol.append(line)
		return self


	def get_molecules_infos(self, idx = None):
		"""
		Method to return molecule informations (each molecule)

		Args:
			idx (int, optional): molecule index (Default: None)

		Returns:
			list
		"""
		if idx is None:
			return self._molecules_info
		elif type(idx) == int:
			return self._molecules_info[idx]
		else:
			sys.stderr.write("ERROR: method error: get_molecules_infos() was given something without integer.\n")
			sys.exit(1)


	def write(self, output):
		"""
		Method to write out file

		Args:
			output (str): output file

		Returns:
			self
		"""
		# 出力
		with open(output, "w") as obj_output_top:
			for line in self._defaults:
				obj_output_top.write(line)

			parameter_name = "{0}_{1}.itp".format(self._prefix, "parameters")
			obj_output_top.write("; include parameter file\n")
			obj_output_top.write("#include \"{0}\"\n".format(parameter_name))

			with open(parameter_name, "w") as obj_output_itp:
				for line in self.parameters:
					obj_output_itp.write(line)
			sys.stderr.write("{0} was created.\n".format(parameter_name))
			obj_output_top.write("\n")

			obj_output_top.write("; include molecules information\n")
			for idx, molecule_name in enumerate(self.molecule_names):
				itp_name = "{0}_{1:02d}_{2}.itp".format(self._prefix, idx, molecule_name)
				obj_output_top.write("#include \"{0}\"\n".format(itp_name))
				with open(itp_name, "w") as obj_output_itp:
					for line in self.get_molecules_infos(idx):
						obj_output_itp.write(line)
				sys.stderr.write("{0} was created.\n".format(itp_name))
			obj_output_top.write("\n")

			obj_output_top.write("[ system ]\n")
			for line in self._system_name:
				obj_output_top.write(line)
			obj_output_top.write("\n")

			obj_output_top.write("[ molecules ]\n")
			for line in self._system_mol:
				obj_output_top.write(line)
			obj_output_top.write("\n")

		return self
