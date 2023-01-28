swig -Isigpack -python -c++ sigpack.i
gcc -DMKL_ILP64  -m64 -L/usr/local/lib -I/usr/local/include/python3.9 -Isigpack  -I"${MKLROOT}/include" -fmax-errors=1 -O2 -fPIC -march=native -mavx2 -shared -o _sigpack.so sigpack_wrap.cxx -lstdc++ -lm -lpython3.9 -lfftw3 -larmadillo  -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl
