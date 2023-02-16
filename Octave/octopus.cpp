#include "audiodsp_lua.hpp"
#include "octopus.hpp"
#include "Threads.hpp"

using namespace Octopus;

volatile bool bLoaded=false;

void* octave_thread(void * p)
{
    Octopus::OctaveInterpreter interpreter;
    bLoaded = true;
    while(1)
    {
        sleep(10);
    }
}

extern Octopus::OctopusValueList* convert_octave_value_list(lua_State * L);

int lua_GetInterpreter(lua_State * L)
{
    //lua_pushlightuserdata(L,&interpreter);
    return 0;
}
int lua_Plot(lua_State * L)
{           
    Octopus::OctopusValueList* v = convert_octave_value_list(L);
    lua_pop(L,1);    
    try {    
        octave::feval("plot",v->vlist,0);     
        ValueList x;        
        octave::feval("pause",x,0);
    } catch(std::runtime_error & e) {
        std::cout << e.what() << std::endl;
    }
    return 0;
}
int lua_FreqzPlot(lua_State * L)
{           
    Octopus::OctopusValueList* v = convert_octave_value_list(L);
    lua_pop(L,1);    
    try {    
        ValueList r,l;
        r = octave::feval("freqz",v->vlist,2);
        l(0) = r(1);
        l(1) = r(0);
        octave::feval("freqz_plot",l,0);        
        octave::feval("pause",v->vlist,0);
    } catch(std::runtime_error & e) {
        std::cout << e.what() << std::endl;
    }
    return 0;
}
extern "C" int luaopen_octopus(lua_State * L);

int main(int argc, char * argv[])
{
    Thread oct(octave_thread,nullptr);
    while(!bLoaded) sleep(1);
    try {
        Lua::LuaJIT interp;     
        luaopen_octopus(interp.L);
        interp.CreateCFunction("getInterpreter",lua_GetInterpreter);   
        interp.CreateCFunction("Plot",lua_Plot);   
        interp.CreateCFunction("FreqzPlot",lua_FreqzPlot);   
        if( interp.DoFile(argv[1]) != 0) {
            std::cout << lua_tostring(interp.L,-1) << std::endl;
        }
    }
    catch(std::runtime_error & what)
    {
        std::cout << what.what() << std::endl;
    }
    catch(...)
    {
    }
}