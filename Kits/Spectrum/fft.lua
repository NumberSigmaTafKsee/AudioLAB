fft=require('luakissfft')
r = fft.table_to_vector({1,2,3,4,5})
for i=0,r:size()-1 do print(r[i]) end
x = fft.vector_to_table(r)
table.foreach(x,print)
