require('tcc')
t = tcc.TinyCC()
t:AddLibraryPath("tcc")
t:SetOutputType(tcc.TCC_OUTPUT_MEMORY)
s = 'int main() { return 10; }';
t:CompileString(s)
t:Relocate(tcc.TCC_RELOCATE_AUTO)
x = t:Run()
print(x)
x=t:Exec("main")
print(x)

t:New()
t:SetOutputType(tcc.TCC_OUTPUT_MEMORY)
s = '#include<tcclib.h>\nint test(void * x) { return 25; }\nint main() { return test(NULL); }';
r = t:CompileString(s)
print(r)
t:Relocate(tcc.TCC_RELOCATE_AUTO)
t:Run()
x = t:Exec("test")
print(x)
x = tcc.tcc_get_symbol(t.state,"test")
print(x)

