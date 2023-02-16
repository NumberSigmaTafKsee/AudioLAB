    swig -lua -c++ -Iinclude LV2Plugin.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o LV2Plugin.so LV2Plugin_wrap.cxx -lstdc++ -lm -lluajit    
    