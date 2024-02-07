
Folder "mixintcount" contains the main program of MixIntCount and benchmarks on random polytopes. For more details, please refer to "mixintcount/README.md".



Folder "VOL2LAT_MIC" is the combination of MixIntCount with VOL2LAT in order to compare MixIntCount with Vol2Lat on benchmarks generated from program analysis. For more details about compilation, please refer to "VOL2LAT_MIC/README.md".





# MixIntCount



## Building MixIntCount

* Make sure that g++, boost, glpk, Armadillo are installed on your machine.

* Then execute: 

	cd mixintcount

	sh build.sh



Some quick tests:

	cd mixintcount

	./mixIntCount benchmarks/cube_10.in

	./mixIntCount -pv benchmarks/random/4_1_1000.in

	./mixIntCount -e benchmarks/random/4_1_1000.in



For details of options, please refer to the help menu of mixIntCount by:

	./mixIntCount -h





### Quick guide for building on Ubuntu



Execute:

	cd mixintcount

	sudo apt-get install g++

	sudo apt-get install libglpk-dev

	sudo apt-get install libboost-dev

	sudo apt-get install libarmadillo-dev

	sh build.sh

	

If you encounter any problem with Armadillo library, you can try to compile and install it (http://arma.sourceforge.net/).






## Benchmarks



* mixintcount/benchmarks/random:

	Basic benchmarks of random polytopes P(m, n, i, l) for Table 1, where m = n ranges from 4 to 15, i in {1, 2, 3}, and l ranges from 10 to 10000.

	Benchmarks are in name "n_i_l.in".



* mixintcount/benchmarks/random_n5, random_n10, random_n15:

	Additional benchmarks of random polytopes P(m, n, i, l) for Figure 5, where m = n = 5, 10, 15, i is fixed to 1, and l ranges from 10 to 10000.

	Benchmarks are in name "n_i_l.in".



* mixintcount/benchmarks/random_v2l:

	Additional benchmarks of random polytopes P(m, n, i, l) for Figure 6, where m = n = i ranges from 3 to 8, and l ranges from 10 to 10000.

	Benchmarks are in name "n_i_l.in", "n_i_l.in2" and "n_i_l.smt2", they are inputs for tool MixIntCount, barvinok and Vol2Lat respectively.


* VOL2LAT_MIC/benchmarks.zip

	Application benchmarks generated from program analysis.




## Data



* mixintcount/mic_ran.xlsx

	Full results of Table 1 is in 'comparison' sheet in mic_ran.xlsx.

	

* mixintcount/mic_ran2.xlsx

	Raw data of Figure 2, three sheets are results of n = 5, 10, 15 respectively.

	Figure 2 contains all rows of results in these sheets.

	

* mixintcount/v2l_mic.xlsx

	Raw data of Figure 3, sheet 'v2l' and 'mic' are results of Vol2Lat and MixIntCount respectively.

	Figure 3 contains all rows of results in these sheets.


* VOL2LAT_MIC/results/

	
	Raw data of Table 2, each file correspond to the results on a family of benchmarks.





