    swig -lua -c++ -Iinclude EMAFilter.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o EMAFilter.so EMAFilter_wrap.cxx -lstdc++ -lm -lluajit    
    