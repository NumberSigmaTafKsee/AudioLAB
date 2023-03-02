    swig -lua -c++ -Iinclude VAKorg35HPFFilter.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VAKorg35HPFFilter.so VAKorg35HPFFilter_wrap.cxx -lstdc++ -lm -lluajit    
    