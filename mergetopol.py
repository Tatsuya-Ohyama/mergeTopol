#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
mergetopol

merge forcefields for gromacs

"""

import sys, os, re, signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

import argparse
from py_module_basic import basic
import textwrap

# =============== functions =============== #
# トポロジーファイルの読み込み
def load_topol(input_file, flag_system = False):
	flag_read = 0
	re_empty = re.compile(r"^[\s\t]*\n$")
	re_comment = re.compile(r"^[\s\t]*;")
	re_default = re.compile(r"\[ defaults \]")
	re_atomtypes = re.compile(r"\[ atomtypes \]")
	re_moleculetypes = re.compile(r"\[ moleculetype \]")
	re_molecules = re.compile(r"\[ molecules \]")
	re_system = re.compile(r"\[ system \]")

	global defaults
	global atomtypes
	global molecules
	global others

	with open(input_file, "r") as obj_input:
		for line in obj_input:
			if re_comment.search(line):
				continue
			elif re_atomtypes.search(line):
				# [ atomtypes ] directive
				flag_read = 1
			elif re_default.search(line):
				flag_read = 4
			elif re_system.search(line):
				# [ system ] directive
				if flag_system == True:
					flag_read = 2
				else:
					flag_read = 0
			elif re_moleculetypes.search(line):
				# [ moleculetype ] directive
				flag_read = 5
				others.append(line)
			elif re_molecules.search(line):
				# [ molecules ] directive
				flag_read = 3

			elif flag_read == 1:
				# atomtypes の登録
				if re_empty.search(line):
					flag_read = 5
				else:
					atomtypes.append(line)

			elif flag_read == 4:
				defaults.append(line)
				flag_read = 0

			elif flag_read == 2:
				# system の登録
				global system
				system = line
				flag_read = 5

			elif flag_read == 3:
				# molecules の登録
				if re_empty.search(line):
					flag_read = 5
				else:
					molecules.append(line)

			elif flag_read == 5:
				# その他の登録
				others.append(line)




# =============== main =============== #
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description = "mergetopol.py - merge forcefields for gromacs", formatter_class=argparse.RawTextHelpFormatter)

	subparser = parser.add_subparsers(help = "sub-command")
	subparser.required = True

	# トポロジーの追記
	add_topol = subparser.add_parser("add", help = "add topology to existing topology")
	add_topol.set_defaults(func = "add_topol")
	add_topol.add_argument("-p", dest = "top", metavar = "TOP", required = True, help = ".top file for input")
	add_topol.add_argument("-f", dest = "forcefield_path", metavar = "forcefield_path", nargs = "+", required = True, help = "forcefield path (file (.itp or .top) or directory)")
	add_topol.add_argument("-o", dest = "output", metavar = "TOP", required = True, help = ".top file for output")
	add_topol.add_argument("-O", dest = "flag_overwrite", action = "store_true", default = False, help = "overwrite forcibly (Default: False)")

	# トポロジーの作成
	create_topol = subparser.add_parser("create", help = "create topology from .gro and forcefields")
	create_topol.set_defaults(func = "create_topol")
	create_topol.add_argument("-c", dest = "gro", metavar = "GRO", required = True, help = ".gro file for input")
	create_topol.add_argument("-f", dest = "forcefield_path", metavar = "forcefield_path", nargs = "+", required = True, help = "forcefield path (file (.itp or .top) or directory)")
	create_topol.add_argument("-o", dest = "output", metavar = "TOP", required = True, help = ".top file for output")

	args = parser.parse_args()

	defaults = []
	atomtypes = ["[ atomtypes ]\n", ";name   bond_type     mass     charge   ptype   sigma         epsilon       Amb\n"]
	molecules = ["[ molecules ]\n", "; Compound        nmols\n"]
	system = ""
	others = []

	if args.func == "add_topol":
		# トポロジーの追記
		basic.check_exist(args.top, 2)
		load_topol(args.top, True)

		for file in args.forcefield_path:
			if basic.check_exist(file, 2):
				# ファイルの場合
				load_topol(file)

			elif basic.check_exist(file, 3):
				# ディレクトリの場合
				pass



	elif args.func == "create_topol":
		# トポロジーの作成
		basic.check_exist(args.gro, 2)


	if args.flag_overwrite == False:
		basic.check_overwrite(args.output)


	tmp = []
	for i in defaults:
		if i not in tmp:
			tmp.append(i)
	defaults = tmp

	if len(defaults) != 1:
		sys.stderr.write("WARNING: different default properties in these topology. Which use?\n")
		count = 1
		for i in default:
			sys.stderr.write(" %d %s\n" % (count, i))
			count += 1
		user = sys.stdin.readline()
		default = defaults[int(user) - 1]
	else:
		default = defaults[0]
	print(textwrap.dedent("""\
	[ defaults ]
	; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ
	%s
	""" % default))
	for elem in atomtypes:
		print(elem, end = "")
	print("")
	for elem in others:
		print(elem, end = "")
	print("[ system ]")
	print(system)
	for elem in molecules:
		print(elem, end = "")
