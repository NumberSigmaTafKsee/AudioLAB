swig -lua -c++ -IMkl -IStd -Icppmkl mkl.i
gcc -fopenmp -pthread -fmax-errors=1 -I. -IStd -IMkl -I/usr/local/include/python3.9 -DMKL_ILP64  -m64  -I"${MKLROOT}/include" -O2 -fPIC -march=native -mavx2 -shared -o mkl.so mkl_wrap.cxx -lstdc++ -lluajit -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_rt -lmkl_intel_ilp64 -lmkl_core -lmkl_intel_thread -liomp5 -lpthread -lm -ldl

