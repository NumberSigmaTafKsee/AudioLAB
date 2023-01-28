    swig -lua -c++ -Iinclude VirtualAnalogStateVariableFilter.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VirtualAnalogStateVariableFilter.so VirtualAnalogStateVariableFilter_wrap.cxx -lstdc++ -lm -lluajit    
    