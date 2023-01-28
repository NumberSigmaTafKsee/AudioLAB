    swig -lua -c++ -Iinclude AdaptiveFilterProcessors.i
    gcc -std=c++17 -I. -I../ -I../include -O2 -fPIC -march=native -mavx2 -shared -o AdaptiveFilterProcessors.so AdaptiveFilterProcessors_wrap.cxx -lstdc++ -lm -lluajit    
    
