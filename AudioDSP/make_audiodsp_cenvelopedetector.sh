    swig -lua -c++ -Iinclude audiodsp_cenvelopedetector.swg
    gcc -Wfatal-errors -std=c++17 -Iinclude -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o audiodsp_cenvelopedetector.so audiodsp_cenvelopedetector_wrap.cxx -lstdc++ -lm -lluajit 
    