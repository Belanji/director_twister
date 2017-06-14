
set terminal pdfcairo dashed enhanced size 3.5 , 2.5 font "Helvetica,13"
set key bottom left font "Helvetica,12"  samplen 2


w_top='10 50 100 200 1000'
w_down=" 10 50 100 200 1000 "

set pointsize 0.3
set ylabel "sin({/Symbol F}(Z=1/2,t))"
set xlabel "t^*"

do for [ wu in w_top ] { 

    set output "w_top=".wu.".pdf"
    
    plot for [wd in w_down] "./w_down=".wd."/w_top=".wu."/middle_sin.dat" u 1:2 every 2 w l lw 2.0 t "w_{-}=".wd

    unset output
    
}

