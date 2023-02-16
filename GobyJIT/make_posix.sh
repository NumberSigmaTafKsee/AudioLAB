swig -Iinclude -lua -c++ posix.i
gcc -Iinclude -O2 -fPIC -shared -o posix6.so posix_wrap.cxx lib/libposix.a -lstdc++ -lm -lluajit -lrt -lpthread
