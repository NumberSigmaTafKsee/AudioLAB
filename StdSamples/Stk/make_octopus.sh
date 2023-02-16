#swig -I/usr/local/include/octave-7.2.0 -lua -c++ octopus.i
#PKG_CONFIG_PATH=~/pkgconfig gcc `pkg-config --cflags octave octinterp plot` -O2 -fPIC -mavx2 -march=native -shared -o octopus.so octopus_wrap.cxx -lstdc++ -lm  -lluajit `pkg-config --libs plot octave octinterp`

gcc -std=gnu++11 -pthread -fopenmp -fmax-errors=1  -I/usr/local/include/octave-7.2.0 \
     -O2 -fPIC -march=native -mavx2 \
     -o octopus  Octopus.cpp \
     -lstdc++ -lm -lluajit   -L/usr/local/lib/octave/7.2.0 -loctave -loctgui -loctinterp
