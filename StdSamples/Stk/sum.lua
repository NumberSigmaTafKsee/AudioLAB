require('vf')
a = vf.Matrix(3,3)
for i=0,2 do
for j=0,2 do
a[i][j] = i*3 + j 
end 
end 
a:upload_device()
a:sigmoid()
a:download_host()
a:print()
a:sigmoid_deriv()
a:print()
