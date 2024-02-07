# MixIntCount



## Building MixIntCount

* Make sure that g++, boost, glpk, Armadillo are installed on your machine.

Then execute: 

	sh build.sh



Some quick tests:

	./mixIntCount benchmarks/cube_10.in

	./mixIntCount -pv benchmarks/random/4_1_1000.in

	./mixIntCount -e benchmarks/random/4_1_1000.in



For details of options, please refer to the help menu of mixIntCount by:

	./mixIntCount -h





### Quick guide for building on Ubuntu



Execute:

	sudo apt-get install g++

	sudo apt-get install libglpk-dev

	sudo apt-get install libboost-dev

	sudo apt-get install libarmadillo-dev

	sh build.sh

	

If you encounter any problem with Armadillo library, you can try to compile and install it (http://arma.sourceforge.net/).






## Benchmarks



* random:

	Basic benchmarks of random polytopes P(m, n, i, l) for Table 1, where m = n ranges from 4 to 15, i in {1, 2, 3}, and l ranges from 10 to 10000.

	Benchmarks are in name "n_i_l.in".



* random_n5, random_n10, random_n15:

	Additional benchmarks of random polytopes P(m, n, i, l) for Figure 5, where m = n = 5, 10, 15, i is fixed to 1, and l ranges from 10 to 10000.

	Benchmarks are in name "n_i_l.in".



* random_v2l:

	Additional benchmarks of random polytopes P(m, n, i, l) for Figure 6, where m = n = i ranges from 3 to 8, and l ranges from 10 to 10000.

	Benchmarks are in name "n_i_l.in", "n_i_l.in2" and "n_i_l.smt2", they are inputs for tool MixIntCount, barvinok and Vol2Lat respectively.





## Data



* mic_ran.xlsx

	Full results of Table 1 is in 'comparison' sheet in mic_ran.xlsx.

	

* mic_ran2.xlsx

	Raw data of Figure 2, three sheets are results of n = 5, 10, 15 respectively.

	Figure 2 contains all rows of results in these sheets.

	

* v2l_mic.xlsx

	Raw data of Figure 3, sheet 'v2l' and 'mic' are results of Vol2Lat and MixIntCount respectively.

	Figure 3 contains all rows of results in these sheets.






