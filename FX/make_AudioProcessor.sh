    swig -lua -c++ -I../include -I../include AudioProcessor.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o AudioProcessor.so AudioProcessor_wrap.cxx -lstdc++ -lm -lluajit 
    