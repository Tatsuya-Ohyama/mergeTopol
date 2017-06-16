#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
mergeTopol.py

merge forcefields for Gromacs

"""

import sys, os, re, signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

import argparse
from py_module_basic import basic
import tempfile
import shutil


# =============== variables =============== #
re_include = re.compile(r"#include")
re_wsp = re.compile(r"[\s\t]+")
re_quote = re.compile(r"[\'\"]")
re_comment = re.compile(r"^[\s\t]*;")
re_ffcomment = re.compile(r"^\*")
re_empty = re.compile(r"^[\s\t]*\n$")
re_include = re.compile(r"#include ")
re_directive = re.compile(r"\[ (.+) \]")
re_posres = re.compile(r"([\s\t]+\d+){2}([\s\t]+\d+(?:\.\d+)?){3}")

# =============== class =============== #
class TopologyParameter:
	""" トポロジーデータクラス """
	def __init__(self, input_file, nmol = 1):
		self.__nmol = nmol			# 分子数

		self.__defaults = []		# defaults directive
		self.__parameters = []		# 力場パラメータを格納
		self.__molecules_name = []	# 各分子の名前を格納
		self.__molecules_info = []	# 各分子の情報を格納 (molecules 情報のインデックスに対応)
		self.__system_name = []		# 系情報を格納
		self.__system_mol = []		# 系内の分子情報を格納

		# ファイル読み込み
		lines = self.__load_data(input_file)
		with open("temp", "w") as obj_output:
			for line in lines:
				obj_output.write(line)
		self.__parse_topology(lines)


	def __load_data(self, input_file):
		""" ファイルを読み込むメソッド """
		# ファイル読み込み
		basic.check_exist(input_file, 2)
		return self.__load_data_sub(input_file)


	def __load_data_sub(self, input_file):
		""" include 行を読み込むメソッド (再帰的に利用可能) """
		basic.check_exist(input_file, 2)

		lines = []
		with open(input_file, "r") as obj_input:
			for line in obj_input:
				if re_include.search(line):
					# include 行がある場合

					# include ファイルの探索
					dirname = os.path.dirname(input_file)
					include_file = re_include.sub("", line)
					include_file = re_quote.sub("", include_file)
					include_file = get_filepath(include_file.rstrip("\n"), dirname)
					lines.extend(self.__load_data_sub(include_file))
					lines.append("\n")

				else:
					# その他の行
					lines.append(line)
		return lines


	def __parse_topology(self, lines):
		""" パラメータの構文解析をするメソッド """
		cnt_line = 0
		flag_directive = 0
		molecule_name = ""
		temporary_info = []
		flag_posres = False

		for line in lines:
			cnt_line += 1
			if re_directive.search(line):
				# directive の特定
				directive = re_directive.search(line).group(1)

				if directive == "molecules":
					# 系内の分子情報への切り替わりの場合
					flag_directive = 5

				elif directive in "system":
					# 系情報への切り替わりの場合
					flag_directive = 4

				elif directive == "moleculetype":
					# 各分子情報への切り替わり、あるいは既に各分子情報の場合
					flag_directive = 3
					temporary_info = ["[ moleculetype ]\n"]

				elif directive == "defaults":
					# defaults の場合
					flag_directive = 1
					self.__defaults.append(line)

				elif 1 <= flag_directive <= 2:
					# 力場パラメータの場合
					flag_directive = 2
					self.__parameters.append(line)

				elif flag_directive == 3:
					# 分子情報のその他の directive の場合
					if directive == "position_restraints":
						# 原子拘束の場合
						flag_posres = True
					else:
						flag_posres = False
					self.__molecules_info[-1].append(line)

			elif 0 <= flag_directive <= 1: # and (not re_comment.search(line) or not re_ffcomment.search(line)):
				self.__defaults.append(line)

			elif flag_directive == 2:
				self.__parameters.append(line)

			elif flag_directive == 3:
				if len(temporary_info) != 0:
					# 新しく登録する moleculetype の場合
					if  not re_comment.search(line):
						# コメントでない場合
						molecule_name = re_wsp.split(line.strip())[0]
						self.__molecules_name.append(molecule_name)
						self.__molecules_info.append(temporary_info)
						self.__molecules_info[-1].append(line)
						temporary_info = []
					else:
						# directive 直後のコメント (ラベル) の場合
						temporary_info.append(line)

				else:
					# 既に登録済みの moleculetype の場合
					if flag_posres and re_posres.search(line):
						# posres の場合
						datas = re_wsp.split(line.strip())
						datas[2] = args.posres[0]
						datas[3] = args.posres[1]
						datas[4] = args.posres[2]
						self.__molecules_info[-1].append("{0[0]:>6} {0[1]:>5} {0[2]:>7} {0[3]:>7} {0[4]:>7}\n".format(datas))
					else:
						self.__molecules_info[-1].append(line)

			elif flag_directive == 4:
				self.__system_name.append(line)

			elif flag_directive == 5:
				if re_comment.search(line):
					self.__system_mol.append(line)
				else:
					mol_info = re_wsp.split(line.strip())
					mol_info[1] = int(mol_info[1]) * self.__nmol
					self.__system_mol.append("{0[0]:<15} {0[1]:>7}\n".format(mol_info))

		# 先頭の空行を削除
		self.__defaults = clean_lines(self.__defaults)
		self.__parameters = clean_lines(self.__parameters)
		self.__molecules_name = clean_lines(self.__molecules_name)
		for idx, info in enumerate(self.__molecules_info):
			self.__molecules_info[idx] = clean_lines(info)
		self.__system_name = clean_lines(self.__system_name)
		self.__system_mol = clean_lines(self.__system_mol)


	def __merge_parameters(self, parameters1, parameters2):
		""" 2 つのパラメータをマージするメソッド """
		directives1 = []
		directives2 = []
		parameter_values1 = {}
		parameter_values2 = {}
		directive = ""

		# パラメータのラベル付
		for line in parameters1:
			if re_directive.search(line):
				# directive の場合
				directive = line
				directives1.append(line)
				parameter_values1[directive] = []
			else:
				parameter_values1[directive].append(line)

		for line in parameters2:
			if re_directive.search(line):
				# directive の場合
				directive = line
				directives2.append(line)
				parameter_values2[directive] = []
			else:
				parameter_values2[directive].append(line)

		# directive のマージ
		for idx, directive in enumerate(directives2):
			if directive not in directives1:
				# 追加パラメータの場合
				pre2 = directives2[idx - 1]
				next2 = ""
				if idx + 1 != len(directives2):
					next2 = directives2[idx + 1]

				idx_pre2 = None
				idx_next2 = None
				if pre2 in directives1:
					idx_pre2 = directives1.index(pre2)
				if next2 in directives1:
					idx_next2 = directives1.index(next2)

				if idx_pre2 is None and idx_next2 is None:
					# directives1 に含まれない場合
					directives1.append(directive)
				elif idx_pre2 is None:
					# 直前の directive が不明の場合、直後の directive の前に挿入する
					directives1.insert(idx_next2, directive)
				elif idx_next2 is None:
					# 直後の directive が不明の場合、直前の directive の後に挿入する
					directives1.insert(idx_pre2 + 1, directive)

		# value のマージ
		results = []
		for directive in directives1:
			results.append(directive)
			if directive in parameter_values1.keys():
				if parameter_values1[directive][-1] in ["", "\n"]:
					# 末尾に空行がある場合
					results.extend(parameter_values1[directive][: -1])
				else:
					results.extend(parameter_values1[directive])
			if directive in parameter_values2.keys():
				results.extend(parameter_values2[directive])
			results.append("\n")

		return results


	def merge_topology(self, obj_molecule):
		""" 他のトポロジーファイルを取り込むメソッド (追加分子など) """
		self.__parameters = self.__merge_parameters(self.__parameters, obj_molecule.get_parameters())
		self.__molecules_name.extend(obj_molecule.get_molecules_name())
		self.__molecules_info.extend(obj_molecule.get_molecules_info())
		for line in obj_molecule.get_system_mol():
			if not re_comment.search(line):
				self.__system_mol.append(line)

	def get_defaults(self):
		""" defaults directive を返すメソッド """
		return self.__defaults

	def get_parameters(self):
		""" 力場パラメータ情報を返すメソッド """
		return self.__parameters

	def get_molecules_name(self):
		""" 分子名リストを返すメソッド """
		return self.__molecules_name

	def get_molecules_info(self, idx = None):
		""" 分子情報 (各分子) を返すメソッド """
		if idx is None:
			return self.__molecules_info
		elif type(idx) == int:
			return self.__molecules_info[idx]
		else:
			sys.stderr.write("ERROR: method error: get_molecules_info() was given something without integer.\n")
			sys.exit(1)

	def get_system_name(self):
		""" 系の情報を返すメソッド """
		return self.__system_name

	def get_system_mol(self):
		""" 系内の分子情報を返すメソッド """
		return self.__system_mol

	def write(self, output):
		""" ファイルに書き出すメソッド """
		if args.flag_overwrite == False:
			basic.check_overwrite(output)

		# 出力
		with open(output, "w") as obj_output_top:
			for line in self.__defaults:
				obj_output_top.write(line)

			parameter_name = "{0}_{1}.itp".format(args.prefix, "parameters")
			obj_output_top.write("; include parameter file\n")
			obj_output_top.write("#include \"{0}\"\n".format(parameter_name))

			with open(parameter_name, "w") as obj_output_itp:
				for line in self.get_parameters():
					obj_output_itp.write(line)
			sys.stderr.write("{0} was created.\n".format(parameter_name))
			obj_output_top.write("\n")

			obj_output_top.write("; include molecules information\n")
			for idx, molecule_name in enumerate(self.get_molecules_name()):
				itp_name = "{0}_moleculeinfo_{1:02d}_{2}.itp".format(args.prefix, idx, molecule_name)
				obj_output_top.write("#include \"{0}\"\n".format(itp_name))
				with open(itp_name, "w") as obj_output_itp:
					for line in self.get_molecules_info(idx):
						obj_output_itp.write(line)
				sys.stderr.write("{0} was created.\n".format(itp_name))
			obj_output_top.write("\n")

			obj_output_top.write("[ system ]\n")
			for line in self.__system_name:
				obj_output_top.write(line)
			obj_output_top.write("\n")

			obj_output_top.write("[ molecules ]\n")
			for line in self.__system_mol:
				obj_output_top.write(line)
			obj_output_top.write("\n")


# =============== functions =============== #
def get_filepath(search_path, dirname = None):
	""" ライブラリ内に存在するファイルを検索し、ファイルの内容を返す関数 """
	search_file = os.path.basename(search_path)

	library_paths = args.library
	if dirname is not None:
		library_paths.insert(0, dirname)

	matched_files = []
	for library_path in library_paths:
		# ライブラリを検索
		for (root, dirs, files) in os.walk(library_path):
			# ライブラリ内のファイルを検索
			matched_files.extend([os.path.join(root, x) for x in files if search_path in os.path.join(root, x)])
		if len(matched_files) != 0:
			break

	if len(matched_files) == 0:
		# ファイルが存在しない場合
		sys.stderr.write("ERROR: include topology file not found ({0}).\n".format(search_file))
		sys.exit(1)
	if 2 <= len(matched_files):
		# 複数のファイルが見つかった場合
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
	""" 空行を処理する関数 """
	results = []

	flag_first = True
	flag_empty = False
	flag_comment = False
	idx_first = 0
	idx_end = 0
	cnt_result_idx = -1

	for line in lines:
		if line == "\n" or line == "":
			# 空行
			if flag_first:
				# 最初の空行はスキップ
				idx_first = cnt_result_idx + 1
			elif flag_empty == False:
				# 1 つ目の空行を認める
				flag_empty = True
				results.append(line)
				cnt_result_idx += 1

		else:
			if re_comment.search(line):
				# コメント行
				if flag_comment == False:
					idx_end = cnt_result_idx + 1
					flag_comment = True
			else:
				flag_comment = False

			flag_first = False
			flag_empty = False
			results.append(line)
			cnt_result_idx += 1

	if flag_comment == False:
		idx_end = len(results)
	if results[-1] == "":
		idx_end = -1

	return results[idx_first : idx_end]


# =============== main =============== #
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description = "mergeTopol.py - merge forcefield for Gromacs", formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument("-p", dest = "top", metavar = "TOP", required = True, help = ".top file for input")
	parser.add_argument("-o", dest = "output", metavar = "TOP", required = True, help = ".top file for output")
	parser.add_argument("-l", dest = "library", metavar = "GMX_LIB_PATH", nargs = "+", default = "[\".\"]", help = "force field path (Default: current directory only)")
	parser.add_argument("-b", dest = "prefix", metavar = "prefix", default = "mergeTopol", help = "prefix for outputing file (Default: output filename)")
	parser.add_argument("-a", dest = "add_molecule", metavar = "TOP", nargs = "*", default = [], help = "add other topology")
	parser.add_argument("-n", dest = "mol", metavar = "NUM", nargs = "*", type = int, default = [], help = "the number of added other topology (with -a)")
	parser.add_argument("-r", dest = "posres", metavar = "VECTOR", nargs = 3, type = int, default = [1000, 1000, 1000],help = "Force constants for position restraints (kJ/mol nm^2) (Default: [1000, 1000, 1000])")
	parser.add_argument("-O", dest = "flag_overwrite", action = "store_true", default = False, help = "overwrite forcibly")
	args = parser.parse_args()

	# ファイルやディレクトリの確認
	basic.check_exist(args.top, 2)
	for elem in args.library:
		basic.check_exist(elem, 3)

	# 追加分子指定の確認
	if len(args.add_molecule) != len(args.mol):
		sys.stderr.write("ERROR: the number of files and mol number for adding molecules are different.\n")
		sys.exit(1)

	# top ファイル
	topology = TopologyParameter(args.top)

	# 追加分子の追加
	for idx in range(len(args.add_molecule)):
		# 追加分子
		new_molecule = TopologyParameter(args.add_molecule[idx], args.mol[idx])
		topology.merge_topology(new_molecule)

	topology.write(args.output)
