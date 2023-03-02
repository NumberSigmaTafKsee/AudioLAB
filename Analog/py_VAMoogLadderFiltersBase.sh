    swig -python -c++ -I../include -I../include PyVAMoogLadderFiltersBase.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -I/usr/local/include/python3.9 -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o _PyVAMoogLadderFiltersBase.so PyVAMoogLadderFiltersBase_wrap.cxx -lstdc++ -lm -L/usr/local/lib -lpython3.9 -lutil -lrt -ldl
    