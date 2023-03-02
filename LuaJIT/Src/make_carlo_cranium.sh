gcc -fopenmp -pthread -fmax-errors=1 -I. -Iinclude/MKL -DMKL_ILP64  -m64  -I"${MKLROOT}/include" -O2 -fPIC -march=native -mavx2 -o carlo_cranium carlo_cranium.cpp -lstdc++ -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lcblas -lmkl_rt -lmkl_intel_ilp64 -lmkl_core -lmkl_intel_thread -liomp5 -lpthread -lm -ldl

