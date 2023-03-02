swig -lua -c++ src/gnuplot.i
gcc -fpermissive -O2 -fPIC -shared -o plot.so src/gnuplot_wrap.cxx src/gnuplot_i.c -lstdc++ -lm -lluajit
