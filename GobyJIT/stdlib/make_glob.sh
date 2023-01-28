swig -lua -c++ glob.i
gcc -fmax-errors=1 -O2 -fPIC -shared -o Glob.so glob_wrap.cxx -lstdc++ -lm -lluajit
