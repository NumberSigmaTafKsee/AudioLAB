swig -lua -c++ -Iinclude audiodsp.i
gcc -fmax-errors=1 -std=c++17 -I. -Iinclude -I/usr/local/include/lilv-0 -c audiodsp_wrap.cxx
