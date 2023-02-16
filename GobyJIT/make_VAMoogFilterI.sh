    swig -lua -c++ -Iinclude VAMoogFilterI.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VAMoogFilterI.so VAMoogFilterI_wrap.cxx -lstdc++ -lm -lluajit    
    