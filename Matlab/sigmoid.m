function x = sigmoid(a)
    x = 1.0 ./ (1.0 + exp(-a));
endfunction