swig -octave -c++ -Iinclude -IKits/AudioDSP oct_stdsamples.i
mkoctfile -Wfatal-errors -fopenmp -pthread -fmax-errors=1 -I. \
-Iinclude -IKits/AudioDSP -std=c++17 -DMKL_ILP64  -m64  -I"${MKLROOT}/include"  \
-O2 -fPIC -march=native -mavx2 -shared -o stdsamples oct_stdsamples_wrap.cxx -lstdc++ -lluajit  \
-L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_rt -lmkl_intel_ilp64 \
-lmkl_core -lmkl_intel_thread -lippcore -lipps -liomp5 -lpthread -lm -ldl -lfftw3 -lfftw3f -lsndfile -lsamplerate

