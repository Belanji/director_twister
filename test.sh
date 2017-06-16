#!/bin/bash 

rm director_twister
gcc director.c -lm -lgsl -lgslcblas -o director_twister

homeDirr=$(pwd)
export wup=( 10 50 100 200 1000 )
export wdown=( 10 50 100 200 1000 )

for wd in ${wdown[@]}
do
    if [ ! -d w_down=${wd} ]
    then

	mkdir w_down=${wd}
	cd w_down=${wd}

    else

	rm  -r w_down=${wd}
	mkdir w_down=${wd}
	cd w_down=${wd}

    fi


    
    for wu in ${wup[@]}
    do
	
	    
	mkdir w_top=${wu}
	cd w_top=${wu}
	    
	
	../../director_twister $wd $wu 
	cd ..
	
    done

    cd ..

done

#gnuplot -e "wa='${ww[*]}'"  plot.gp
#evince watop_const.pdf &
