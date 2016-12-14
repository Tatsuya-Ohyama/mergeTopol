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
re_quote = re.compile(r"[\'\"]")


# =============== functions =============== #
# トポロジーファイルの読み込み
def load_file(input_file, output_file, library):
	with tempfile.NamedTemporaryFile(mode = "w", prefix = ".mergeTopol_", delete = False) as obj_output:
		tempfile_name = obj_output.name
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

							tmp_library = library
							if os.path.abspath(os.path.dirname(lib_path)) not in tmp_library:
								# 参照しているファイルと同じレベルのパスが登録されていない場合、登録する
								tmp_library.insert(0, os.path.abspath(os.path.dirname(lib_path)))

							# 再帰的に include 内容の読み込み
							print("!!!!!", new_path)
							load_file(lib_path, new_path, tmp_library)

							flag_found = 1
							obj_output.write("#include \"%s\"\n" % new_path)
							break


					if flag_found == 0:
						sys.stderr.write("ERROR: library not found (%s)\n" % path)
						sys.exit(1)
				else:
					# その他の行
					obj_output.write(line)

	shutil.move(tempfile_name, output_file)


# =============== main =============== #
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description = "mergeTopol.py - merge forcefield for Gromacs", formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument("-p", dest = "top", metavar = "TOP", required = True, help = ".top file for input")
	parser.add_argument("-o", dest = "output", metavar = "TOP", required = True, help = ".top file for output")
	parser.add_argument("-l", dest = "library", metavar = "GMX_LIB_PATH", nargs = "+", default = "[\".\"]", help = "force field path (Default: current directory only)")
	parser.add_argument("-b", dest = "basename", metavar = "prefix", help = "prefix for outputing file (Default: output filename)")
	parser.add_argument("-O", dest = "flag_overwrite", action = "store_true", default = False, help = "overwrite forcibly")
	args = parser.parse_args()

	basic.check_exist(args.top, 2)

	for elem in args.library:
		basic.check_exist(elem, 3)

	basename_output = ""
	if args.basename != None:
		# basename が指定されている場合
		basename_output = args.basename
	else:
		# basename が指定されていない場合
		basename_output = args.output
		try:
			pos = basename_output.rindex("/")
		except ValueError:
			# 同じ階層のファイルの場合
			pos = 0
		basename_output = basename_output[pos + 1 : ]

		try:
			pos = basename_output.rindex("\.")
		except ValueError:
			# 拡張子前のドットがある場合
			basename_output = basename_output[0 : pos]

	output = basename_output + ".top"
	if args.flag_overwrite == False:
		basic.check_overwrite(output)
	load_file(args.top, output, args.library)
