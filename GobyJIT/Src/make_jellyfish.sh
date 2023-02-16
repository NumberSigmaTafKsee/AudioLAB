swig -Icppy3/cppy3 -lua -c++ jellyfish.i
gcc -fmax-errors=1  -I/usr/local/include/python3.9  -O2 -mavx2 -march=native -mfma -fPIC -shared -o jellyfish.so jellyfish_wrap.cxx libcppy3.a -lm -lluajit -lpython3.9 -lrt -lutil -lpthread -lstdc++
