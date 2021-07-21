#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
mergeTopol.py

merge forcefields for Gromacs

"""

import sys, signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

import argparse

from mods.func_prompt_io import check_exist, check_overwrite
from mods.topology_parameter import TopologyParameter



# =============== main =============== #
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="mergeTopol.py - merge forcefield for Gromacs", formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument("-p", dest="TOP_FILE", metavar="INPUT.top", required=True, help="gromacs topology file (Input)")
	parser.add_argument("-o", dest="OUTPUT_FILE", metavar="OUTPUT.top", required=True, help="gromacs topology file (Output)")
	parser.add_argument("-l", dest="LIBRARY_DIR", metavar="GMX_LIB_PATH", nargs="+", default="[\".\"]", help="force field path (Default: current directory only)")
	parser.add_argument("-b", dest="OUTPUT_PREFIX", metavar="OUTPUT_PREFIX", default="mergeTopol", help="prefix for outputing file (Default: output filename)")
	parser.add_argument("-a", dest="INSERT_MOL_FILE", metavar="INSERT.top", nargs="*", default=[], help="gromacs topology file for inserted molecules")
	parser.add_argument("-n", dest="N_MOL", metavar="NUM", nargs="*", type=int, default=[], help="the number of added other topology (with -a)")
	parser.add_argument("-r", dest="POSRES", metavar="POSRES_FORCE", nargs=3, type=int, default=[1000, 1000, 1000],help="Force constants for position restraints (kJ/mol nm^2) (Default: [1000, 1000, 1000])")
	parser.add_argument("-O", dest="FLAG_OVERWRITE", action="store_true", default=False, help="overwrite forcibly")
	args = parser.parse_args()


	# check files and directories
	check_exist(args.TOP_FILE, 2)
	for elem in args.LIBRARY_DIR:
		check_exist(elem, 3)

	# check additional molecule
	if len(args.INSERT_MOL_FILE) != len(args.N_MOL):
		sys.stderr.write("ERROR: the number of files and mol number for adding molecules are different.\n")
		sys.exit(1)
	for add_molecule in args.INSERT_MOL_FILE:
		check_exist(add_molecule, 2)

	# read topology file
	topology = TopologyParameter(args.TOP_FILE, library=args.LIBRARY_DIR, posres=args.POSRES, prefix=args.OUTPUT_PREFIX)

	# add additional molecules
	for idx in range(len(args.INSERT_MOL_FILE)):
		# at additional molecule
		new_molecule = TopologyParameter(args.INSERT_MOL_FILE[idx], library=args.LIBRARY_DIR, nmol=args.N_MOL[idx], posres=args.POSRES, prefix=args.OUTPUT_PREFIX)
		topology.merge_topology(new_molecule)

	if args.FLAG_OVERWRITE == False:
		check_overwrite(args.OUTPUT_FILE)
	topology.write(args.OUTPUT_FILE)
