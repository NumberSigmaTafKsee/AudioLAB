swig -Iinclude -lua -c++ plot.i
gcc -Iinclude -O2 -fPIC -shared -o plot.so plot_wrap.cxx include/gnuplot_i.c -lstdc++ -lm -lluajit
