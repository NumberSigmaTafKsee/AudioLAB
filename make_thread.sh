swig -lua -c++ src/thread.i
gcc -fmax-errors=1  -O2 -fPIC -shared -o thread.so src/thread_wrap.cxx -lstdc++ -lluajit -lpthread
