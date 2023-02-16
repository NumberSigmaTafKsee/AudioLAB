require('stdsamples')
a = stdsamples.sample_vector(16)
b = stdsamples.sample_vector(16)
a:fill(10)
b:fill(20)
c = a + b
c:print()
c = a - b
c:print()
c = a * b
c:print()
c = a / b
c:print()
c = a ^ 2
c:print()

