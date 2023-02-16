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

int main(int argc, char * argv)
{
    auto vm = cppy3::PythonVM();
    auto v = cppy3::eval("1+2");
    double x = v.toDouble();
    std::cout << x << std::endl;
    cppy3::List l(PyList_New(10));
    l.append(v.data());
    cppy3::execScript("test.py");
    cppy3::arguments args;
    args.push_back(l);
    cppy3::call("foo",args);
}