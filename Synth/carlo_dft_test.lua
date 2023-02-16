require('casino')

-- there isn't a real dft as the way ipp does it is retarded
-- to do real make a complex and set imaginary to 0

a = casino.complex_vector(16)
b = casino.complex_vector(16)
c = casino.complex_vector(16)

for i=1,16 do 
    a[i] = casino.complex(i,0)
end

dft = casino.CDFT32(16)
dft:Forward(a:data(),b:data())

for i = 1,16 do
    io.write('('..b[i]:real()..','..b[i]:imag()..')')
end
print()

dft:Inverse(b:data(),c:data())
for i = 1,16 do
    io.write('('..c[i]:real()..','..c[i]:imag()..')')
end
print()


for i=1,16 do 
    a[i] = casino.complex(i,0)
end


casino.dft(dft,a:data(),b:data())

for i = 1,16 do
    io.write('('..b[i]:real()..','..b[i]:imag()..')')
end
print()

casino.idft(dft,b:data(),c:data())
for i = 1,16 do
    io.write('('..c[i]:real()..','..c[i]:imag()..')')
end
print()

a = casino.vector(16)
c = casino.vector(16)

for i=1,16 do 
    a[i] = i*i
end
casino.dft(dft,a:data(),b:data())

for i = 1,16 do
    io.write('('..b[i]:real()..','..b[i]:imag()..')')
end
print()

casino.idft(dft,b:data(),c:data())
for i = 1,16 do
    io.write(c[i]..',')
end
print()
