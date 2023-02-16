swig -Icppy3/cppy3 -lua -c++ jam.i
gcc -fmax-errors=1  -I/usr/local/include/python3.9  -Icppy3/cppy3 -O2 -fPIC -shared -o jam.so jam_wrap.cxx libcppy3.a -lm -lluajit -lpython3.9 -lrt -lutil -lpthread -lstdc++
