swig -lua -c++ filesystem.i 
gcc -Iinclude/Std -fmax-errors=1 -std=c++2a -shared -fPIC -O2 -o filesystem.so filesystem_wrap.cxx -lstdc++fs -lstdc++ -lluajit
