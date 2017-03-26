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
	perl evaluate_gap_filling.pl -s original_scaffolds.fa -f gap-filled_scaffolds.fa -r reference_genome.fa -o output_directory [options]

	-s <scaffold file>
		Mandatory parameter. It points out the original scaffold file.
	-f <gap-filled scaffold file>
		Mandatory parameter. It points out the gap-filled scaffold file.
	-r <reference genome file>
		Mandatory parameter. It points out the reference genome file.
	-o <output directory>
		Mandatory parameter. It points out the output directory.
	-c <mininum contig length>
		Optional parameter. When the length of a sub-contig is shorter than this mininum value, it will be ingored as gap region. Default value 100.
	-m <minimum distance>
		Optional parameter. If the left-most (or right-most) aligned segment of a sub-contig aligns over this minimum value from the starting (or ending) position of the sub-contig, the sub-contig will be regarded as non-aligned when it is considered as a right (or left) flanking sub-contig of a gap. Default value 200.
	-d <mininum distance>
		Optional parameter. It is used to define Relocation gap. Default value 3000.
	-a <mininum times>
	Optional parameter. It is also used to define Relocation gap. Default value 3.

5) Output:

	The final evaluation result is named "result.txt" in the output directory.
	The first column is the gap type, it shows the evaluation mertics from second column to nine column based on each gap type.
	The last row shows the evaluation metrics based on all gap types.
	An example of result.txt:

Gap Type  | Gap Count  | Reference Gap Length  | Gap-filled Length  | Match Number  | Mis-match Number  | Precision  | Recall  | F1-score
--------- | --------  | --------  | --------  | --------  | --------  | --------  | --------  | --------
Missing Gap  | 0  | 0  | 0  | 0  | 0  | -  | -  | -
Translocation Gap  | 0  | 0  | 0  | 0  | 0  | -  | -  | -
Inversion Gap  | 11  | 154  | 139  | 135  | 4  | 0.971  | 0.877  | 0.922
Relocation Gap  | 17  | 2784  | 2345  | 2311  | 34  | 0.986  | 0.830  | 0.901
Normal Gap  | 97  | 31777  | 13522  | 13341  | 181  | 0.987  | 0.420  | 0.589
All Gap  | 125  | 34715  | 16006  | 15787  | 219  | 0.986  | 0.455  | 0.623


	The gap locations on the original scaffold, gap-filled scaffold and the reference genome are showed in the file "infor.txt".  
	The gap references and gap-filled regions of all gaps are shown in the files "gap_in_reference.txt" and "gap_in_filling.txt" respectively.