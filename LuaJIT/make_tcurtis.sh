swig -lua -c++ -Iinclude src/TCurtisVCF.i
gcc -Iinclude -fmax-errors=1 -O2 -fPIC -march=native -mavx2 -shared -o TCurtisVCF.so src/TCurtisVCF_wrap.cxx -lstdc++ -lm -lluajit
