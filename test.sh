#!/bin/bash 

rm director_twister
icc director.c -lm -lgsl -lgslcblas -o director_twister

homeDirr=$(pwd)
export wup=( 10.0 50.0 100.0 200.0 1000.0 )
export wdown=( 10.0 50.0 100.0 200.0 1000.0 )

#export wup=( 0.0 0.1 1.0 10.0 100.0 )  
#export wdown=( 0.0 0.1 1.0 10.0 100.0 )


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

gnuplot   plot.gp

