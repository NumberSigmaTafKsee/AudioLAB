function x = cwiseops(op,a,b)        
    if(strcmp(op,"plus"))      
        x = a .+ b;
    elseif(strcmp(op,"minus")) 
        x = a .- b;
    elseif(strcmp(op,"times"))   
        x = a .* b;
    elseif(strcmp(op,"divide"))  
        x = a ./ b;
    elseif(strcmp(op,"add"))
        x = a + b;
    elseif(strcmp(op,"sub"))
        x = a - b;
    elseif(strcmp(op,"mul"))
        x = a * b;
    elseif(strcmp(op,"div"))
        x = a / b;
    else 
        display("error");
    endif    
endfunction