function m = cwisemax(a,max)
    for i= 1:size(a,1);
        for j=1:size(a,2)
            if a(i,j) < max
                a(i,j) = max;
            end
        end
    end
endfunction