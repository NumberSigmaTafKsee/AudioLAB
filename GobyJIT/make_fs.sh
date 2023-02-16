swig -lua -c++ -I../include/Std src/filesys.i 
gcc -I../include/Std -fmax-errors=1 -std=c++17 -shared -fPIC -O2 -ofs.so src/filesys_wrap.cxx -lstdc++ -lstdc++fs -lluajit
