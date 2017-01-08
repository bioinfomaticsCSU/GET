License
=========

Copyright (C) 2014 Jianxin Wang(jxwang@mail.csu.edu.cn), Junwei Luo(luojunwei@csu.edu.cn)

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.

Jianxin Wang(jxwang@mail.csu.edu.cn), Junwei Luo(luojunwei@csu.edu.cn)
School of Information Science and Engineering
Central South University
ChangSha
CHINA, 410083


GET
=================
1) Introduction

	GET is a gap-filling evaluation tool.
	The input includes the original scaffolds, the gap-filled scaffolds and the reference genome. 

2) Before installing and running
	
	Users should install MUMmer and add it to your PATH. 
	Users can download MUMmer from http://mummer.sourceforge.net/

3) Installing

	GET is written C++ and therefore will require a machine with GNU C++ pre-installed.
	Create a main directory (eg:GET). Copy all source code to this directory.
	Type "make all".

4) Running

	Run command line:  
	perl evaluate_gap_filling.pl original_scaffolds.fa gap-filled_scaffolds.fa reference_genome.fa output_directory 

5) Output:

	The final evaluation result is named "result.txt" in the output_directory.

