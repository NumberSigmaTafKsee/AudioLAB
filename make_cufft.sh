swig -lua -c++ cufft.i
nvcc --ptxas-options=-v --compiler-options '-fPIC -fmax-errors=1' -arch=sm_61 -shared -o cufft.so cufft_wrap.cxx -lstdc++ -lm -lluajit -lcudart -lcufft -lcufftw
