#include "cppy3.hpp"

/*
require('jellyfish')
vm = jellyfish.PythonVM()
v = jellyfish.eval("1+2")
x = v:toDouble()
print(x)
jellyfish.execScript("test.py")
*/
namespace cppy3
{
    inline Var execScript(const std::string& path) {
    std::ifstream t(path);
    
    std::string script((std::istreambuf_iterator<char>(t)),
                    std::istreambuf_iterator<char>());
    return exec(script.c_str());
    }
}

int main(int argc, char * argv[])
{
    auto vm = cppy3::PythonVM();
    if(argc > 1)
        cppy3::execScript(argv[1]);    
}