#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys, re, os


# =============== variables =============== #
re_include = re.compile(r"#include")
re_wsp = re.compile(r"[\s\t]+")
re_quote = re.compile(r"[\'\"]")
re_comment = re.compile(r"^[\s\t]*;")
re_ffcomment = re.compile(r"^\*")
re_empty = re.compile(r"^[\s\t]*\n$")
re_directive = re.compile(r"\[ (.+) \]")
re_posres = re.compile(r"([\s\t]+\d+){2}([\s\t]+\d+(?:\.\d+)?){3}")



# =============== function =============== #
def get_filepath(search_path, library, dirname = None):
	"""
	ライブラリ内に存在するファイルを検索し、ファイルの内容を返す関数
	@param search_path: 検索するパス
	@param library: ライブラリパス
	@param dirname: 読み込み元 (include 行を持つファイル) のディレクトリパス (Default: None)
	@return パラメータファイルのパス
	"""
	search_file = os.path.basename(search_path)

	library_paths = [x for x in library]
	if dirname == "":
		dirname = "."
	if dirname is not None:
		# 読み込み元のディレクトリパスがある場合、そのパスを先頭に配置
		library_paths.insert(0, dirname)

	matched_files = []
	for library_path in library_paths:
		# ライブラリを検索
		matched_files.extend([os.path.join(root, file) for root, dirs, files in os.walk(library_path) for file in files if search_path in os.path.join(root, file)])
		if len(matched_files) != 0:
			break

	if len(matched_files) == 0:
		# ファイルが存在しない場合
		sys.stderr.write("ERROR: include topology file not found ({0}).\n".format(search_path))
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
	"""
	空行を処理する関数
	@param lines: テキストリスト
	@return 空行を除去したテキストリスト
	"""
	results = []

	flag_first = True
	flag_empty = False
	flag_comment = False
	idx_first = 0
	idx_end = 0
	cnt_result_idx = -1

	for line in lines:
		if len(line.strip()) == 0:
			# 空行の場合
			if flag_first:
				# 最初の空行の場合、スキップ
				idx_first = cnt_result_idx + 1
			elif flag_empty == False:
				# 1 つ目の空行を認める
				flag_empty = True
				results.append(line)
				cnt_result_idx += 1

		else:
			if re_comment.search(line):
				# コメント行の場合
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

	if len(results) == 0:
		return results
	else:
		if results[-1] == "":
			idx_end = -1

		return results[idx_first:idx_end]



# =============== class =============== #
class TopologyParameter:
	""" トポロジーデータクラス """
	def __init__(self, input_file = None, library = [], nmol = 1, posres = [1000, 1000, 1000], prefix = ""):
		self._nmol = 0			# 分子数
		self._library = []		# パラメータが格納されているライブラリ
		self._posres = []		# 拘束エネルギー
		self._prefix = ""		# 出力ファイル用の接頭辞

		self._defaults = []			# defaults directive
		self._parameters = []		# 力場パラメータを格納
		self._molecules_name = []	# 各分子の名前を格納
		self._molecules_info = []	# 各分子の情報を格納 (molecules 情報のインデックスに対応)
		self._system_name = []		# 系情報を格納
		self._system_mol = []		# 系内の分子情報を格納

		# initiation
		self.set_library(library)
		self.set_nmol(nmol)
		self.set_prefix(prefix)
		self.set_posres(posres)

		lines = self._load_data(input_file)
		self._parse_topology(lines)


	def _load_data(self, input_file):
		"""
		ファイルを読み込むメソッド
		@param input_file: 読み込むファイルパス
		@return: 読み込んだ内容
		"""
		if input_file is not None:
			return self._load_data_sub(input_file)
		else:
			return None


	def _load_data_sub(self, input_file):
		"""
		include 行を読み込むメソッド (再帰的に利用可能)
		@param input_file: 読み込むファイルパス
		"""
		lines = []
		with open(input_file, "r") as obj_input:
			for line in obj_input:
				if re_include.search(line):
					# include 行がある場合

					# include ファイルの探索
					dirname = os.path.dirname(input_file)
					include_file = re_include.sub("", line)
					include_file = re_quote.sub("", include_file)
					include_file = get_filepath(include_file.strip(), self._library, dirname)
					lines.extend(self._load_data_sub(include_file))
					lines.append("\n")

				else:
					# その他の行
					lines.append(line)
		return lines


	def _parse_topology(self, lines):
		""" パラメータの構文解析をするメソッド """
		flag_directive = 0
		molecule_name = ""
		temporary_info = []
		flag_posres = False

		for line_val in lines:
			if re_directive.search(line_val):
				# directive の特定
				directive = re_directive.search(line_val).group(1)

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
					self._defaults.append(line_val)

				elif 1 <= flag_directive <= 2:
					# 力場パラメータの場合
					flag_directive = 2
					self._parameters.append(line_val)

				elif flag_directive == 3:
					# 分子情報のその他の directive の場合
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
					# 新しく登録する moleculetype の場合
					if  not re_comment.search(line_val):
						# コメントでない場合
						molecule_name = re_wsp.split(line_val.strip())[0]
						self._molecules_name.append(molecule_name)
						self._molecules_info.append(temporary_info)
						self._molecules_info[-1].append(line_val)
						temporary_info = []
					else:
						# directive 直後のコメント (ラベル) の場合
						temporary_info.append(line_val)

				else:
					# 既に登録済みの moleculetype の場合
					if flag_posres and re_posres.search(line_val):
						# posres の場合
						datas = re_wsp.split(line_val.strip())
						datas[1:4] = self._posres
						self._molecules_info[-1].append("{0[0]:>6} {0[1]:>5} {0[2]:>7} {0[3]:>7} {0[4]:>7}\n".format(datas))
					else:
						self._molecules_info[-1].append(line_val)

			elif flag_directive == 4:
				self._system_name.append(line_val)

			elif flag_directive == 5:
				if re_comment.search(line_val):
					self._system_mol.append(line_val)
				else:
					mol_info = re_wsp.split(line_val.strip())
					mol_info[1] = int(mol_info[1]) * self._nmol
					self._system_mol.append("{0[0]:<15} {0[1]:>7}\n".format(mol_info))

		# 先頭の空行を削除
		self._defaults = clean_lines(self._defaults)
		self._parameters = clean_lines(self._parameters)
		self._molecules_name = clean_lines(self._molecules_name)
		for idx, info in enumerate(self._molecules_info):
			self._molecules_info[idx] = clean_lines(info)
		self._system_name = clean_lines(self._system_name)
		self._system_mol = clean_lines(self._system_mol)


	def set_library(self, library):
		"""
		パラメータ検索用のライブラリパスリストを設定するメソッド
		@param library: パラメータ検索用ライブラリパスリスト
		@return self
		"""
		self._library = library
		return self


	def set_nmol(self, nmol):
		"""
		分子数を設定するメソッド
		@param mol: 分子数
		@return self
		"""
		self._nmol = nmol
		return self


	def set_posres(self, posres):
		"""
		位置拘束用の拘束エネルギーを設定するメソッド
		@param posres: 拘束エネルギーリスト [x, y, z]
		@return self
		"""
		self._posres = posres
		return self


	def set_prefix(self, prefix):
		"""
		出力用の接頭辞を設定するメソッド
		@param prefix: 接頭辞
		@return self
		"""
		self._prefix = prefix
		return self


	def _merge_parameters(self, parameters1, parameters2):
		"""
		2 つのパラメータをマージするメソッド
		@param parameter1: obj_TopologyParameter
		@param parameter2: obj_TopologyParameter
		return: マージしたパラメータ
		"""
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
		self._parameters = self._merge_parameters(self._parameters, obj_molecule.get_parameters())
		self._molecules_name.extend(obj_molecule.get_molecules_name())
		self._molecules_info.extend(obj_molecule.get_molecules_info())
		for line in obj_molecule.get_system_mol():
			if not re_comment.search(line):
				self._system_mol.append(line)

	def get_defaults(self):
		""" defaults directive を返すメソッド """
		return self._defaults

	def get_parameters(self):
		""" 力場パラメータ情報を返すメソッド """
		return self._parameters

	def get_molecules_name(self):
		""" 分子名リストを返すメソッド """
		return self._molecules_name

	def get_molecules_info(self, idx = None):
		""" 分子情報 (各分子) を返すメソッド """
		if idx is None:
			return self._molecules_info
		elif type(idx) == int:
			return self._molecules_info[idx]
		else:
			sys.stderr.write("ERROR: method error: get_molecules_info() was given something without integer.\n")
			sys.exit(1)

	def get_system_name(self):
		""" 系の情報を返すメソッド """
		return self._system_name

	def get_system_mol(self):
		""" 系内の分子情報を返すメソッド """
		return self._system_mol

	def write(self, output):
		""" ファイルに書き出すメソッド """
		# 出力
		with open(output, "w") as obj_output_top:
			for line in self._defaults:
				obj_output_top.write(line)

			parameter_name = "{0}_{1}.itp".format(self._prefix, "parameters")
			obj_output_top.write("; include parameter file\n")
			obj_output_top.write("#include \"{0}\"\n".format(parameter_name))

			with open(parameter_name, "w") as obj_output_itp:
				for line in self.get_parameters():
					obj_output_itp.write(line)
			sys.stderr.write("{0} was created.\n".format(parameter_name))
			obj_output_top.write("\n")

			obj_output_top.write("; include molecules information\n")
			for idx, molecule_name in enumerate(self.get_molecules_name()):
				itp_name = "{0}_{1:02d}_{2}.itp".format(self._prefix, idx, molecule_name)
				obj_output_top.write("#include \"{0}\"\n".format(itp_name))
				with open(itp_name, "w") as obj_output_itp:
					for line in self.get_molecules_info(idx):
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


# =============== main =============== #
# if __name__ == '__main__':
# 	main()
