    swig -lua -c++ -Iinclude audiodsp_svflineartrapoptimized2.swg
    gcc -Wfatal-errors -std=c++17 -Iinclude -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o audiodsp_svflineartrapoptimized2.so audiodsp_svflineartrapoptimized2_wrap.cxx -lstdc++ -lm -lluajit 
    