octave_test;
v = octave_test.double_vector(10);
r = zeros(1,256);
v = octave_test.vectorize_double(r)
r = octave_test.vectorize(v)