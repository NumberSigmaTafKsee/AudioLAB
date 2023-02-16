clear
echo 'compiling'
nvcc  -use_fast_math -arch=sm_61 -gencode=arch=compute_61,code=sm_61 --compiler-options '-fPIC -fmax-errors=1' -c viper_vector.cu viper_matrix.cu viper_complex_vector.cu viper_complex_matrix.cu

