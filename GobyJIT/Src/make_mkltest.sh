gcc -fopenmp -pthread -fmax-errors=1 -I. -Iinclude/MKL -DMKL_ILP64  -m64  -I"${MKLROOT}/include" -O2 -fPIC -march=native -mavx2 -o test_mkl test_mkl.cpp -lstdc++ -lluajit -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lcblas -lmkl_intel_thread -lmkl_rt -lmkl_core -lmkl_intel_ilp64 -liomp5 -lpthread -lm -ldl

