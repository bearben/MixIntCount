#!/bin/sh
ADDR=$(cd "$(dirname "$0")"; pwd)
echo "Working directory: $ADDR"

if [ ! -d ${ADDR}/bin/ ]; then
	mkdir ${ADDR}/bin
fi

#vinci
echo "Cheking and compiling Vinci ..."
if [ ! -d ${ADDR}/vinci-1.0.5/ ]; then
	unzip vinci-1.0.5.zip
	cd ${ADDR}/vinci-1.0.5
	make
	cd ../
	mv ${ADDR}/vinci-1.0.5/vinci ${ADDR}/bin/vinci
fi

#volce main program
echo "Compiling mixIntCount ..."
cd ${ADDR}
make depend
make
