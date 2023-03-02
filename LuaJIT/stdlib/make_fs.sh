swig -lua -c++ -I../include/Std filesys.i 
gcc -I../include/Std -fmax-errors=1 -std=c++17 -shared -fPIC -O2 -ofs.so filesys_wrap.cxx -lstdc++ -lstdc++fs -lluajit
