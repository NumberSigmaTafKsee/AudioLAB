swig -lua -c++ -I/usr/local/include std_arrayfire.i
gcc -fmax-errors=1 -I/usr/local/include -L/usr/local/lib -O2 -fPIC -march=native -mavx2 -shared -o af.so std_arrayfire_wrap.cxx -lstdc++ -lm -lluajit -lafcpu -laf
