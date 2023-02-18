    swig -lua -c++ -I../include audiodsp_faustfx.swg
    gcc -Wfatal-errors -std=c++17 -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o audiodsp_faustfx.so audiodsp_faustfx_wrap.cxx \
    -lstdc++ -lm -lluajit -lfaustwithllvm -lpthread -lrt -ldl -lLLVM-10 -lz -lcurses 
    
