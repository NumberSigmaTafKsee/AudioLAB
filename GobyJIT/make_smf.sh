swig -lua -c++ -I/usr/include src/smf.i
gcc -I/usr/include/glib-2.0 -I/usr/lib/x86_64-linux-gnu/glib-2.0/include  -O2 -fPIC -march=native -mavx2 -shared -o smf.so src/smf_wrap.cxx -lstdc++ -lm -lluajit -lsmf
