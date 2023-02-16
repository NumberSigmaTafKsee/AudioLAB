    swig -lua -c++ -Iinclude AmplifiersUdo1.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o AmplifiersUdo1.so AmplifiersUdo1_wrap.cxx -lstdc++ -lm -lluajit    
    