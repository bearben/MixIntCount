#!/bin/bash

contains()
  case "$1" in
    (*"$2") true;;
    (*) false;;
  esac

function read_dir() {
for file in `ls $1`

do

    if [ -d $1"/"$file ] 

    then

    read_dir $1"/"$file

    else

		if contains $file ".in" ;then
   			echo $1"/"$file
    		./timeout -m 16000000 -t 3600 ./mixIntCount $2 $1"/"$file
		fi

    fi

done

}

read_dir $1 $2
