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
re_empty = re.compile(r"^[\s\t]*\n$")
re_atomtypes = re.compile(r"\[ atomtypes \]")
re_moleculetype = re.compile(r"\[ moleculetype \]")
re_molecules = re.compile(r"\[ molecules \]")
re_system = re.compile(r"\[ system \]")
re_posres = re.compile(r"\[ position_restraints \]")
re_posres_main = re.compile(r"([\s\t]+\d+){2}([\s\t]+\d+(?:\.\d+)?){3}")


# =============== functions =============== #
# トポロジーファイルの読み込み
def load_file(input_file, output_file, library):
	global atomtypes
	global molecules
	global new_topol_paths

	with tempfile.NamedTemporaryFile(mode = "w", prefix = ".mergeTopol_", delete = False) as obj_output:
		tempfile_name = obj_output.name
		flag_read = 0
		with open(input_file, "r") as obj_input:
			for line in obj_input:
				if re_include.search(line):
					# include 行がある場合
					path = re_include.sub("", line).strip()
					path = re_quote.sub("", path)

					flag_found = 0
					for lib_path in library:
						# include パスの検索
						lib_path = os.path.join(lib_path, path)
						if os.path.isfile(lib_path):
							# パスが存在する場合

							# 該当ファイルをカレントディレクトリにコピー
							new_path = basename_output + "_" + os.path.basename(path)
							if args.flag_overwrite == False:
								basic.check_overwrite(new_path)
							shutil.copyfile(lib_path, new_path)
							sys.stderr.write(" * %s was created\n" % new_path)

							tmp_library = library
							if os.path.abspath(os.path.dirname(lib_path)) not in tmp_library:
								# 参照しているファイルと同じレベルのパスが登録されていない場合、登録する
								tmp_library.insert(0, os.path.abspath(os.path.dirname(lib_path)))

							# 再帰的に include 内容の読み込み
							load_file(lib_path, new_path, tmp_library)

							flag_found = 1
							obj_output.write("#include \"%s\"\n" % new_path)
							break

					if flag_found == 0:
						sys.stderr.write("ERROR: library not found (%s)\n" % path)
						sys.exit(1)
					continue

				elif re_atomtypes.search(line):
					# [ atomtypes ]
					flag_read = 1

				elif re_system.search(line):
					# [ system ] (トポロジーを追加する場合)
					for new_topol in new_topol_paths:
						obj_output.write("; Include topology for new molecules added by mergeTopol.py\n")
						obj_output.write("#include \"%s\"\n" % new_topol)
					obj_output.write("\n")
					new_topol_paths = []

				elif re_molecules.search(line):
					# [ molecules ]
					flag_read = 3

				elif re_posres.search(line):
					# [ position_restrains ]
					flag_read = 5

				elif flag_read == 1:
					# [ atomtypes ] の処理 (atomtypes を追加する場合)
					if re_empty.search(line):
						for atomtype in atomtypes:
							obj_output.write(atomtype)
						atomtypes = []

				elif flag_read == 3:
					# [ molecules ] の処理 (分子を追加する場合)
					if re_empty.search(line):
						for molecule in molecules:
							obj_output.write(molecule)
						molecules = []

				elif flag_read == 5:
					# [ position_restraints ] の処理 (拘束条件を変更する場合)
					if not re_comment.search(line) and re_posres_main.search(line) and args.posres != None:
						datas = re_wsp.split(line.strip())
						line = "%5d %5d %5d %5d %5d\n" % (int(datas[0]), int(datas[1]), args.posres[0], args.posres[1], args.posres[2])

				obj_output.write(line)


		if flag_read == 1 and len(atomtypes) != 0:
			# [ atomtypes ] の処理が未消化の場合
			for atomtype in atomtypes:
				obj_output.write(atomtype)
			atomtypes = []

		elif flag_read == 3 and len(molecules) != 0:
			# [ molecules ] の処理が未消化の場合
			for molecule in molecules:
				obj_output.write(molecule)
			molecules = []

	shutil.move(tempfile_name, output_file)


# 追加用 top ファイルの読み込み
def load_new_topol(input_file):
	global atomtypes
	global molecules
	global new_topol_paths

	flag_read = 0
	with tempfile.NamedTemporaryFile(mode = "w", prefix = ".mergeTopol_", delete = False) as obj_output:
		tempfile_name = obj_output.name
		with open(input_file, "r") as obj_input:
			for line in obj_input:
				if re_atomtypes.search(line):
					# [ atomtypes ]
					flag_read = 1
				elif re_moleculetype.search(line):
					# [ moleculetype ]
					flag_read = 5
					obj_output.write(line)
				elif re_system.search(line):
					# [ system ]
					flag_read = 0
				elif re_molecules.search(line):
					# [ molecules ]
					flag_read = 3
				elif flag_read == 1:
					# [ atomtypes ] の処理
					if re_empty.search(line):
						# 空行
						flag_read = 0
					elif not re_comment.search(line):
						# 通常行
						atomtypes.append(line)
				elif flag_read == 5:
					obj_output.write(line)
				elif flag_read == 3:
					if not re_comment.search(line):
						# 通常行
						molecules.append(line)

	output_file = basename_output + "_" + os.path.basename(input_file)
	output_file = re.sub(r"\.top", ".itp", output_file)
	if args.flag_overwrite == False:
		basic.check_overwrite(output_file)
	shutil.move(tempfile_name, output_file)
	new_topol_paths.append(output_file)


# =============== main =============== #
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description = "mergeTopol.py - merge forcefield for Gromacs", formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument("-p", dest = "top", metavar = "TOP", required = True, help = ".top file for input")
	parser.add_argument("-o", dest = "output", metavar = "TOP", required = True, help = ".top file for output")
	parser.add_argument("-l", dest = "library", metavar = "GMX_LIB_PATH", nargs = "+", default = "[\".\"]", help = "force field path (Default: current directory only)")
	parser.add_argument("-b", dest = "basename", metavar = "prefix", help = "prefix for outputing file (Default: output filename)")
	parser.add_argument("-a", dest = "add", metavar = "TOP", nargs = "*", default = [], help = "add other topology")
	parser.add_argument("-n", dest = "mol", metavar = "NUM", nargs = "*", default = [], help = "the number of added other topology (with -a)")
	parser.add_argument("-r", dest = "posres", metavar = "VECTOR", nargs = 3, type = int, help = "Force constants for position restraints (kJ/mol nm^2)")
	parser.add_argument("-O", dest = "flag_overwrite", action = "store_true", default = False, help = "overwrite forcibly")
	args = parser.parse_args()

	# ファイルやディレクトリの確認
	basic.check_exist(args.top, 2)
	for elem in args.library:
		basic.check_exist(elem, 3)

	# 出力ファイルの接頭辞の決定
	basename_output = ""
	if args.basename != None:
		# basename が指定されている場合
		basename_output = args.basename
	else:
		# basename が指定されていない場合
		basename_output = args.output
		pos = 0
		if "/" in basename_output:
			try:
				pos = basename_output.rindex("/")
			except ValueError:
				# 同じ階層のファイルの場合
				pos = 0
			basename_output = basename_output[pos + 1 : ]

		if "." in basename_output:
			# 拡張子前のドットがある場合
			try:
				pos = basename_output.rindex(".")
			except ValueError:
				pos = len(basename_output)
			basename_output = basename_output[0 : pos]

	atomtypes = []
	molecules = []
	new_topol_paths = []
	if len(args.add) != 0:
		# 分子を追加する場合
		count = 0
		for new_topol in args.add:
			load_new_topol(new_topol)
			if len(args.mol) != 0:
				data = molecules[count].strip()
				datas = re_wsp.split(data)
				datas[1] = args.mol[count]
				molecules[count] = "%-15s %5s\n" % (datas[0], datas[1])
			count += 1

	if args.flag_overwrite == False:
		basic.check_overwrite(args.output)
	load_file(args.top, args.output, args.library)
