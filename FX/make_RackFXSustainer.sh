    swig -lua -c++ -I../include -I../include RackFXSustainer.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o RackFXSustainer.so RackFXSustainer_wrap.cxx -lstdc++ -lm -lluajit 
    