gcc -fcf-protection=none -fmax-errors=1 -std=c++17 -Iinclude  -O2 -march=native -mavx2 -fopenmp -o samples test_samples.cpp -lstdc++ -lm -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lcblas -lmkl_rt -lmkl_intel_ilp64 -lmkl_core -lmkl_intel_thread -liomp5 -lpthread -lm -ldl