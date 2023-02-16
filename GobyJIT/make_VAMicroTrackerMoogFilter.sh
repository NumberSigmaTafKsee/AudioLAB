    swig -lua -c++ -Iinclude VAMicroTrackerMoogFilter.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VAMicroTrackerMoogFilter.so VAMicroTrackerMoogFilter_wrap.cxx -lstdc++ -lm -lluajit    
    