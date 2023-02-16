gcc -fopenmp -pthread -std=c++17 -fmax-errors=1 -IDSP -Iinclude -Iinclude/MKL -DMKL_ILP64-m64 -mavx2 -mfma  -I"${MKLROOT}/include" -O2 -fPIC -march=native \
    -o casino DSP/Carlo/carlo_casino.cpp lib/libaudiofft.a \
    -lstdc++ -lluajit -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_rt -lmkl_intel_ilp64 -lmkl_core -lmkl_intel_thread -lsvml -lippvm -lippcore -lipps -liomp5 -lpthread -lm -ldl -lfftw3 -lfftw3f -lsndfile 

