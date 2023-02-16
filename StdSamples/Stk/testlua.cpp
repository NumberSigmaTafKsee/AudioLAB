#include "include/LuaJIT.hpp"

using namespace Lua;

int main()
{
    LuaJIT lua("test.lua");
    LuaValues values;
    LuaValue a;
    LuaValue b;
    a.setDouble(10.0);
    b.setDouble(20.0);
    values.push_back(a);
    values.push_back(b);
    LuaFunction test("foo");
    values = test(lua,values,1);    
    std::cout << values[0].d << std::endl;
    lua.PushString("Hello World");
    std::cout << lua.asString() << std::endl;

    
    LuaVariable v(lua,"x");
    v = 10.0;
    v.print();
    std::cout << v.asDouble() << std::endl;
    v = true;
    v.print();
    std::cout << v.asBool() << std::endl;
    v = "Hello";
    v.print();
    
}