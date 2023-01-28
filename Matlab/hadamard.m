function x = hadamard(op,a,b)
    if(op == "plus") x = hadamard_plus(a,b);
    elseif(op == "minus") x = hadamard_minus(a,b);
    elseif(op == "times") x = hadamard_times(a,b);
    elseif(op == "divide") x = hadamard_div(a,b);
    elseif(op == "add") x =a+b;
    elseif(op == "sub") x =a-b;
    elseif(op == "mul") x = a*b;
    elseif(op == "div") x =a/b;
    else x = a .* 0
    endif    
endfunction

function x = hadamard_plus(a,b)
    x = a.+b;
endfunction
function x = hadamard_minus(a,b)
    x = a.-b;
endfunction
function x = hadamard_times(a,b)
    x = a.*b;
endfunction
function x = hadamard_div(a,b)
    x = a./b;
endfunction