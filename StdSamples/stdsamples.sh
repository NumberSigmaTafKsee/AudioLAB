swig -lua -c++ -I../include stdsamples.i
gcc -Wfatal-errors -fopenmp -pthread -fmax-errors=1 -I. \
-I../include -std=c++17 -DMKL_ILP64  -m64  -I"${MKLROOT}/include"  \
-O2 -fPIC -march=native -mavx2 -shared -o stdsamples.so stdsamples_wrap.cxx \
-lstdc++ -lluajit  \
-L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_rt -lmkl_intel_ilp64 \
-lmkl_core -lmkl_intel_thread -lippcore -lipps -liomp5 -lpthread -lm -ldl -lfftw3 -lfftw3f -lsndfile

