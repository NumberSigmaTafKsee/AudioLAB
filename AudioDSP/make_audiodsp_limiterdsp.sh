    swig -lua -c++ -Iinclude audiodsp_limiterdsp.swg
    gcc -Wfatal-errors -std=c++17 -Iinclude -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o audiodsp_limiterdsp.so audiodsp_limiterdsp_wrap.cxx -lstdc++ -lm -lluajit 
    