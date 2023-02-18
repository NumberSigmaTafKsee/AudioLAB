    swig -lua -c++ -I../include -I../include LV2Plugin.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o LV2Plugin.so LV2Plugin_wrap.cxx -lstdc++ -lm -lluajit 
    