swig -Ilibtempo/libtempo -lua -c++ libtempo.i
gcc -fopenmp -Isigpack-1.2.7/sigpack -Ilibtempo/libtempo  -DMKL_ILP64  -m64  -I"${MKLROOT}/include" -O2 -fPIC -mavx2 -march=native -shared -o libtempo.so libtempo_wrap.cxx lib/liblibtempo.a -lstdc++ -lm -lluajit -larmadillo -lfftw3  -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl
