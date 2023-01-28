swig -Iinclude -lua -c++ MQ.i
gcc -Iinclude -O2 -fPIC -shared -o MQ.so MQ_wrap.cxx -lstdc++ -lluajit -lrt -lpthread
