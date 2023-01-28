    swig -lua -c++ -Iinclude VirtualAnalogDiodeLadderFilter.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VirtualAnalogDiodeLadderFilter.so VirtualAnalogDiodeLadderFilter_wrap.cxx -lstdc++ -lm -lluajit    
    