swig -lua -c++ -Iinclude audiodsp_chamberlin_svf_filter.i
gcc -fmax-errors=1 -std=c++17 -I. -Iinclude -I. -O2 -fPIC -march=native -mavx2 -shared -o audiodsp_chamberlin_svf_filter.so \
audiodsp_chamberlin_svf_filter_wrap.cxx -lstdc++ -lm -lluajit
