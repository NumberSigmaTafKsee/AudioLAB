    swig -python -c++ -I../include -I../include PyVAKorg35HPFilter.cpp.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -I/usr/local/include/python3.9 -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o _PyVAKorg35HPFilter.cpp.so PyVAKorg35HPFilter.cpp_wrap.cxx -lstdc++ -lm -L/usr/local/lib -lpython3.9 -lutil -lrt -ldl
    