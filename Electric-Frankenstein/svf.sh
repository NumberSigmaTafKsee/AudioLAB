gcc -DUSE_MKL -fmax-errors=1 -DMKL_ILP64-m64 -std=c++17 -mavx2 -mfma  -I"${MKLROOT}/include" -I../include -I../StdSamples -O2  \
-march=native -mavx2 -fopenmp -pthread -o svf electric_svf.cpp -lstdc++ -L/usr/local/cuda/lib64 \
-L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_rt -lmkl_intel_ilp64 -lmkl_core -lmkl_intel_thread -liomp5 -lpthread -lm -ldl -lfftw3 -lfftw3f -lsndfile -lboost_math_c99 -lluajit
