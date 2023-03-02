function x = sigmoid_grad(m)
    x = sigmoid(m) .* (1 .- sigmoid(m));
endfunction