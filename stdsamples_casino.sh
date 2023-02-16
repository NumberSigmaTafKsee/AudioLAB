swig -lua -c++ -IDSP/Carlo -Iinclude audiodsp_casino.i
gcc -Wfatal-errors -Iinclude -IDSP/Carlo -O2 -fPIC -march=native -mavx2 -mfma -fopenmp -pthread -shared \
-DMKL_ILP64-m64 -DMKL_Complex8="std::complex<float>" -DMKL_Complex16="std::complex<double>" -IDSP/Carlo -I"${MKLROOT}/include" \
-o audiodsp_casino.so audiodsp_casino_wrap.cxx -lstdc++ -lm -lluajit \
-L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_rt -lmkl_intel_ilp64 -lmkl_core -lmkl_intel_thread -lsvml -lippvm -lippcore -lipps -liomp5 -lpthread -lm -ldl -lfftw3 -lfftw3f -lsndfile 
