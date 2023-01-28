swig -lua -c++ -Iinclude audiodsp_cdca.i
gcc -fmax-errors=1 -std=c++17 -I. -Iinclude -I.-O2 -fPIC -march=native -mavx2 -shared -o cdca.so audiodsp_cdca_wrap.cxx  -lstdc++ -lm -lluajit
