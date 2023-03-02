    swig -python -c++ -I../include -I../include PyVADiodeLadderFilter.cpp.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -I/usr/local/include/python3.9 -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o _PyVADiodeLadderFilter.cpp.so PyVADiodeLadderFilter.cpp_wrap.cxx -lstdc++ -lm -L/usr/local/lib -lpython3.9 -lutil -lrt -ldl
    