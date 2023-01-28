    swig -lua -c++ -Iinclude RackFXRecChord.cpp.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o RackFXRecChord.cpp.so RackFXRecChord.cpp_wrap.cxx -lstdc++ -lm -lluajit    
    