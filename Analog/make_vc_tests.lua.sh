    swig -lua -c++ -I../include -I../include/Analog vc_tests.lua.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o vc_tests.lua.so vc_tests.lua_wrap.cxx -lstdc++ -lm -lluajit 
    