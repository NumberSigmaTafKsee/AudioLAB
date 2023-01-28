require('vf')
a = vf.Matrix(3,3)
a[0][0]=1
a[0][1]=2
a[0][2]=3
a[1][0]=4
a[1][1]=5
a[1][2]=6
a[2][0]=7
a[2][1]=8
a[2][2]=9
a:upload_device()
a:print()
for i=0,2 do 
	for j=0,2 do
		a[i][j]=i*3+j 
	end
end 
a:upload_device()
a:print()
