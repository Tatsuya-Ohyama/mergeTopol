#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
mergeTopol.py

merge forcefields for Gromacs

"""

import sys, signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

import argparse
from basic_func import check_exist, check_overwrite

from classes.topology_parameter import TopologyParameter


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
	check_exist(args.top, 2)
	for elem in args.library:
		check_exist(elem, 3)

	# 追加分子指定の確認
	if len(args.add_molecule) != len(args.mol):
		sys.stderr.write("ERROR: the number of files and mol number for adding molecules are different.\n")
		sys.exit(1)
	for add_molecule in args.add_molecule:
		check_exist(add_molecule, 2)

	# top ファイル
	topology = TopologyParameter(args.top, library = args.library, posres = args.posres, prefix = args.prefix)

	# 追加分子の追加
	for idx in range(len(args.add_molecule)):
		# 追加分子
		new_molecule = TopologyParameter(args.add_molecule[idx], library = args.library, nmol = args.mol[idx], posres = args.posres, prefix = args.prefix)
		topology.merge_topology(new_molecule)

	if args.flag_overwrite == False:
		check_overwrite(args.output)
	topology.write(args.output)
