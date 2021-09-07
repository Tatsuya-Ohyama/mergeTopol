# MergeTopol

## 概要
Gromacs の力場をマージするプログラム


## 使用方法
```sh
$ mergeTopol.py [-h] -p INPUT.top -o OUTPUT.top [-l GMX_LIB_PATH [GMX_LIB_PATH ...]] [-b OUTPUT_PREFIX] [-a INSERT.top [INSERT.top ...]] [-n NUM [NUM ...]] [-r POSRES_FORCE POSRES_FORCE POSRES_FORCE] [-O]
```

* `-h`, `--help`
	: ヘルプメッセージを表示して終了する。
* `-p INPUT.top`
	: Gromacs トポロジーファイル (Input)
* `-o OUTPUT.top`
	: Gromacs トポロジーファイル (Output)
* `-l GMX_LIB_PATH [GMX_LIB_PATH ...]`
	: 力場パス (Default: `.`) (Ubuntu の apt でインストールした場合: `/usr/share/gromacs/top`)
* `-b OUTPUT_PREFIX`
	: 出力ファイルの接頭辞 (Default: output filename)
* `-a INSERT.top [INSERT.top ...]`
	: 追加するリガンドの力場ファイル (`gmx insert-molecules` で追加した分子の力場ファイル)
* `-n NUM [NUM ...]`
	: 追加するリガンドの数 (`-a` オプションが必要)
* `-r POSRES_FORCE POSRES_FORCE POSRES_FORCE`
	: 追加したリガンドの位置拘束でかける力 (kJ/mol nm^2^) (Default: [1000, 1000, 1000])
* `-O`
	: プロンプトを出さずに上書きする。


## 動作要件
* Python3


## License
The MIT License (MIT)

Copyright (c) year name


## Authors
* Tatsuya Ohyama


## ChangeLog
### Ver. 1.21 (2021-09-07)
* github で公開した。
