swig -lua -c++ -Iinclude audiodsp_cppfilters.i
gcc -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o cppfilters.so audiodsp_cppfilters_wrap.cxx -lstdc++ -lm -lluajit
