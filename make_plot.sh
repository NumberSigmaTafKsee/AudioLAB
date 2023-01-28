swig -lua -c++ -Iinclude plot.i
gcc -Wfatal-errors -Iinclude -fpermissive -O2 -fPIC -shared -o plot.so plot_wrap.cxx -lstdc++ -lm -lluajit
