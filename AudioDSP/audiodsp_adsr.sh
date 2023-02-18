swig -lua -c++ -I../include audiodsp_adsr.i
gcc -fmax-errors=1 -fopenmp -pthread -std=c++17 -I. -I../include -I.-O2 -fPIC -march=native -mavx2 -shared -fopenmp -pthread -o audiodsp_adsr.so \
audiodsp_adsr_wrap.cxx -lstdc++ -lm -lluajit
