require('se')

v = se.CreateRowVectorFloat({1,2,3})
v:print()

-- m = se.FloatRowVector({1,2,3})

print("Array")
a = se.CreateArrayFloat({1,2,3,4,5,6,7,8,9})
a:print()

c = se.CreateColVectorFloat({1,2,3})
c:print() 


m = se.CreateMatrixFloat({{1,2,3},{4,5,6},{7,8,9}})
m:print()


--doesnt' work yet
--v = se.FloatRowVector({1,2,3})
