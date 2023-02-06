swig -w -lua -c++ -Iinclude -Iinclude/Eigen eigen.i
#gcc -fmax-errors=1 -Iinclude -O2 -fPIC -march=native -shared -oeigen.so eigen_wrap.cxx -lstdc++ -lluajit -L/usr/local/cuda/lib64 -lnvblas -lopenblas -llapacke

gcc -fmax-errors=1 -Iinclude -Iinclude/Eigen -DMKL_ILP64-m64 -mavx2 -mfma  -I"${MKLROOT}/include" \
     -O2 -march=native -mavx2 -shared -fPIC -o eigen.so eigen_wrap.cxx -lstdc++ -lm -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_rt -lmkl_intel_ilp64 -lmkl_core -lmkl_intel_thread -lsvml -lippvm -lippcore -lipps -liomp5 -lpthread -lm -ldl -lfftw3 -lfftw3f -lsndfile 
