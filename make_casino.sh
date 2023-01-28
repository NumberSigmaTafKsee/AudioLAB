swig -lua -c++ -IDSP/Carlo -Iinclude casino.i
gcc -fopenmp -pthread -std=c++17 -fmax-errors=1 -I. -Iinclude -IDSP/Carlo -DMKL_ILP64-m64 -mavx2 -mfma  \
    -DMKL_Complex8="std::complex<float>" -DMKL_Complex16="std::complex<double>" \
    -ICarlo -I"${MKLROOT}/include" -O2 -fPIC -march=native -shared -mavx2 -mfma \
    -o casino.so casino_wrap.cxx \
    -lstdc++ -lluajit -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_rt -lmkl_intel_ilp64 -lmkl_core -lmkl_intel_thread -lsvml -lippvm -lippcore -lipps -liomp5 -lpthread -lm -ldl -lfftw3 -lfftw3f -lsndfile 

