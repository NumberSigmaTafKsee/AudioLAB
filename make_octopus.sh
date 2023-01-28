#swig -I/usr/local/include/octave-7.2.0 -lua -c++ octopus.i
#PKG_CONFIG_PATH=~/pkgconfig gcc `pkg-config --cflags octave octinterp plot` -O2 -fPIC -mavx2 -march=native -shared -o octopus.so octopus_wrap.cxx -lstdc++ -lm  -lluajit `pkg-config --libs plot octave octinterp`

gcc -Wfatal-errors -std=gnu++11 -pthread -fopenmp -fmax-errors=1 -Iinclude  -I/usr/local/include/octave-7.2.0 \
     -DMKL_ILP64  -m64  -I"${MKLROOT}/include" -O2 -fPIC -march=native -mavx2 \
     -o octopus  Octopus.cpp \
     -lstdc++ -lm -lluajit -lfftw3 -lfftw3f -lpthread -lsndfile -L/usr/local/lib/octave/7.2.0 -loctave -loctgui -loctinterp \
     -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_rt -lmkl_intel_ilp64 -lmkl_core -lmkl_intel_thread -liomp5 -lpthread -lm -ldl

swig -lua -c++ -I/usr/local/include/octave-7.3.0 octopus.i
gcc -I/usr/local/include/octave-7.3.0 -Wfatal-errors -O2 -fPIC -mavx2 -march=native -mfma -shared -o octopus.so octopus_wrap.cxx -lstdc++ -lm -lluajit -L/usr/local/lib/octave/7.3.0 -loctave -loctinterp
