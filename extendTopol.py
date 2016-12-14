#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
extendTopol.py

extend include statement for Gromacs

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
# load_file
def load_file(path, library):
	lines = []

	with open(path, "r") as obj_input:
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
						lines.append("\n; <--------------- extend %s --------------->\n" % path)
						tmp_library = library
						if os.path.abspath(os.path.dirname(lib_path)) not in tmp_library:
							# 参照しているファイルと同じレベルのパスが登録されていない場合、登録する
							tmp_library.insert(0, os.path.abspath(os.path.dirname(lib_path)))

						# 再帰的に include 内容の読み込み
						include_file = load_file(lib_path, tmp_library)

						# 追加する
						lines.extend(include_file)

						flag_found = 1
						break
				if flag_found == 0:
					sys.stderr.write("ERROR: library not found (%s)\n" % path)
					sys.exit(1)
			else:
				# その他の行
				lines.append(line)

	return lines

# =============== main =============== #
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description = "extendTopol - expand include statement for Gromacs", formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument("-p", dest = "top", metavar = "topol", required = True, help = ".top file for input")
	parser.add_argument("-l", dest = "library", metavar = "GMX_LIB_PATH", nargs = "+", default = "[\".\"]", help = "force field path (Default: current directory only)")
	parser.add_argument("-o", dest = "output", metavar = "TOP", required = True, help = "output")
	parser.add_argument("-O", dest = "flag_overwrite", action = "store_true", default = False, help = "overwrite forcibly (Default: False)")
	args = parser.parse_args()

	basic.check_exist(args.top, 2)

	for elem in args.library:
		basic.check_exist(elem, 3)

	lines = load_file(args.top, args.library)

	if args.flag_overwrite == False:
		basic.check_overwrite(args.output)

	with open(args.output, "w") as obj_output:
		for line in lines:
			obj_output.write(line)
