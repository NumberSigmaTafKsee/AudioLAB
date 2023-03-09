#pragma once

#include <lua.hpp>
#include <iostream>
#include <string>
#include <vector>

using std::string;
using std::vector;
using std::cin;

namespace Lua
{
    struct LuaJIT
    {
        int status;
        lua_State *L;

        LuaJIT() {
            L = luaL_newstate(); // open Lua
            if (!L) {
                    std::cout << "Error creating lua\n";
                    exit(-1);
            }

            luaL_openlibs(L); // load Lua libraries  		
        }
        LuaJIT(const char * script) {
            L = luaL_newstate(); // open Lua
            if (!L) {
                    std::cout << "Error creating lua\n";
                    exit(-1);
            }

            luaL_openlibs(L); // load Lua libraries  
            status = luaL_loadfile(L, script);  // load Lua script
            int ret = lua_pcall(L, 0, 0, 0); // tell Lua to run the script
            if (ret != 0) {
                std::cout << "Lua Error: " << lua_tostring(L, -1) << std::endl;
                exit(-1);
            }
        }
        ~LuaJIT()
        {
            if(L) lua_close(L);
        }

        int DoFile(const string& filename)
        {
            return luaL_dofile(L, filename.c_str());
        }
        int LoadFile(const string& filename)
        {
            int status = luaL_loadfile(L, filename.c_str());
            if (status == 0)
            {
                lua_setglobal(L, filename.c_str());
            }
            return status;
        }

        int operator()(const string& cmd) {
            return DoCmd(cmd);
        }

        int DoCmd(const string& cmd) {
            return luaL_dostring(L,cmd.c_str());
        }

        int Call(const string& func) {
            GetGlobal(func);
            PCall();
            return 0;
        }
        void pcall(int n, int r)
        {
            if(lua_pcall(L, n,r, 0) != 0)
                std::cout << lua_tostring(L,-1) << std::endl;
        }
        int PCall()
        {
            return lua_pcall(L, 0, LUA_MULTRET, 0);
        }

        void NewTable(const string& name)
        {
            lua_newtable(L);
            lua_setglobal(L, name.c_str());
        }

        void GetGlobal(const string& name)
        {
            lua_getglobal(L, name.c_str());
        }
        void SetGlobal(const string& name)
        {
            lua_setglobal(L, name.c_str());
        }

        void PushNumber(double number)
        {
            lua_pushnumber(L, number);
        }
        void PushBool(bool b)
        {
            lua_pushboolean(L, b);
        }
        void PushString(const string& str)
        {
            lua_pushstring(L, str.c_str());
        }
		void PushLightUserData(void * ptr) {
			lua_pushlightuserdata(L,ptr);
		}
        void SetTableValue(double index, double value)
        {
            PushNumber(index);
            PushNumber(value);
            lua_rawset(L, -3);
        }
        void SetTableValue(double index, bool value)
        {
            PushNumber(index);
            PushBool(value);
            lua_rawset(L, -3);
        }
        void SetTableValue(double index, const string& value)
        {
            PushNumber(index);
            PushString(value);
            lua_rawset(L, -3);
        }    
        double getTableDouble(double index)
        {
            PushNumber(index);
            lua_rawget(L,-3);
            double r = lua_tonumber(L,-1);
            Pop(1);
            return r;
        }
        string getTableString(double index)
        {
            PushNumber(index);
            lua_rawget(L,-3);
            std::string r = lua_tostring(L,-1);
            Pop(1);
            return r;
        }
        bool getTableBool(double index)
        {
            PushNumber(index);
            lua_rawget(L,-3);
            bool r = lua_toboolean(L,-1);
            Pop(1);
            return r;
        }

        void Pop(int number)
        {
            lua_pop(L, number);
        }

        void CreateCFunction(const string& name, lua_CFunction function)
        {
            lua_pushcfunction(L, function);
            lua_setglobal(L, name.c_str());
        }

        double asNumber()
        {
            return lua_tonumber(L, -1);
        }
        bool asBool()
        {
            return lua_toboolean(L, -1);
        }
        std::string asString()
        {
            std::string r = (lua_tostring(L, -1));            
            return r;
        }

        bool isNumber() { 
            return lua_isnumber(L,-1);
        }
        bool isBool() { 
            return lua_isboolean(L,-1);
        }
        bool isString() {             
            return lua_isstring(L,-1);
        }

        void rawget(int stack=-3) {
            lua_rawget(L,stack);
        }
        void rawgeti(int i, int stack=-2) {
            lua_rawgeti(L,stack,i);
        }
        void rawset(int stack=-3) {
            lua_rawset(L,stack);
        }
        void rawseti(int i,int stack=-2) {
            lua_rawseti(L,stack,i);
        }
    };




    struct LuaValue
    {
        
        double d;
        std::string s;
        bool b;
    
        enum {
            DOUBLE,
            STRING,
            BOOL,
        };
        int type = DOUBLE;

        LuaValue() {
            d = 0;
            type = DOUBLE;
        }    
        LuaValue(const LuaValue & v) {
            *this = v;
        }
        ~LuaValue() {
        }

        void setDouble(double v) {
            d = v;
            type = DOUBLE;
        }
        void setString(const string& v) {
            s = v;
            type = STRING;
        }
        void setBool(bool v) {
            b = v;
            type = BOOL;
        }
        LuaValue& operator = (const LuaValue & v) {
            if(type == DOUBLE) d = v.d;
            else if(type == BOOL) b = v.b;
            else s = v.s;
            type = v.type;
            return *this;
        }
        double asDouble() {
            if(type==DOUBLE) return d;
            if(type==BOOL) return b;
            return std::stod(s);
        }
        std::string asString() {
            if(type == DOUBLE) {
                char temp[256];
                sprintf(temp,"%f",d);
                return std::string(temp);
            }
            if(type == BOOL) return (b == true)? "True":"False";
            return s;
        }
        void GetValue(LuaJIT & lua) {
            if(lua.isNumber()) setDouble(lua.asNumber());
            else if(lua.isBool()) setBool(lua.asBool());
            else if(lua.isString()) setString(lua.asString());
        }
        bool asBool() {
            if(type == DOUBLE) return d;
            if(type == STRING) return false;
            return b;
        }
        void push(LuaJIT & lua) {
            if(type == DOUBLE) lua.PushNumber(d);
            else if(type == BOOL) lua.PushBool(b);
            else lua.PushString(s.c_str());
        }
        
        void print() {
            switch(type)
            {
                case DOUBLE: std::cout << "<DOUBLE>"; break;
                case BOOL: std::cout << "<BOOL>"; break;
                case STRING: std::cout << "<STRING>"; break;
            }
            if(type == DOUBLE) std::cout << this->d;
            else if(type == BOOL) std::cout << this->b;
            else std::cout << this->s;
            std::cout << std::endl;
        }
    };

    using LuaValues = std::vector<LuaValue>;

    struct LuaFunction
    {
        std::string func;
        LuaFunction(const std::string& name) : func(name)
        {

        }
        LuaValues operator()(LuaJIT& lua, int num_outputs=0)
        {
            LuaValues values;
            LuaValue value;
            lua.GetGlobal(func.c_str());
            lua.pcall(0,num_outputs);
            if(num_outputs > 0)
            {
                for(size_t i = 0; i < num_outputs; i++)
                {
                    if(lua.isNumber())
                    {
                        value.d = lua.asNumber();                    
                        value.type = LuaValue::DOUBLE;
                        values.push_back(value);
                        lua.Pop(1);
                    }
                    else if(lua.isBool())
                    {
                        value.b = lua.asBool();
                        value.type = LuaValue::BOOL;
                        values.push_back(value);
                        lua.Pop(1);
                    }
                    else if(lua.isString())
                    {
                        value.s = lua.asString();
                        value.type = LuaValue::STRING;
                        values.push_back(value);
                        lua.Pop(1);
                    }
                    else {
                        throw std::runtime_error("Unknown Lua Type on stack\n");
                    }
                }
            }
            return values;
        }
        LuaValues operator()(LuaJIT& lua,LuaValues& input, int num_outputs=0)
        {
            LuaValues values;
            LuaValue  value;
            lua.GetGlobal(func.c_str());
            for(size_t i = 0; i < input.size(); i++) input[i].push(lua);
            lua.pcall(input.size(),num_outputs);
            if(num_outputs > 0)
            {
                for(size_t i = 0; i < num_outputs; i++)
                {
                    if(lua.isNumber())
                    {
                        value.d = lua.asNumber();                    
                        value.type = LuaValue::DOUBLE;
                        values.push_back(value);
                        lua.Pop(1);
                    }
                    else if(lua.isBool())
                    {
                        value.b = lua.asBool();
                        value.type = LuaValue::BOOL;
                        values.push_back(value);
                        lua.Pop(1);
                    }
                    else if(lua.isString())
                    {
                        value.s = lua.asString();
                        value.type = LuaValue::STRING;
                        values.push_back(value);
                        lua.Pop(1);
                    }
                    else {
                        throw std::runtime_error("Unknown Lua Type on stack\n");
                    }
                }
            }
            return values;
        }
    };

    struct LuaVariable
    {
        LuaValue value;
        std::string name;
        LuaJIT& lua;

        void setValue() {
            if(value.type == LuaValue::DOUBLE) lua.PushNumber(value.d);
            else if(value.type == LuaValue::BOOL) lua.PushBool(value.b);
            else lua.PushString(value.s);
            lua.SetGlobal(name);            
        }
        void updateValue() {
            lua.GetGlobal(name);
            if(lua.isNumber()) value.setDouble(lua.asNumber());
            else if(lua.isBool()) value.setBool(lua.asBool());
            else value.setString(lua.asString());            
            lua.Pop(1);
        }


        LuaVariable(LuaJIT & l, const std::string& t) : lua(l),name(t) {}

        LuaVariable& operator = (double v) {
            value.setDouble(v);
            setValue();
            return *this;
        }
        LuaVariable& operator = (bool v) {
            value.setBool(v);
            setValue();
            return *this;
        }
        LuaVariable& operator = (const string& v) {            
            value.setString(v);            
            setValue();
            return *this;
        }        
        LuaVariable& operator = (const char * v) {            
            value.setString(std::string(v));            
            setValue();
            return *this;
        }        
        LuaVariable& operator = (const LuaVariable& v) {
            value = v.value;
            name = v.name;
            lua = v.lua;
            setValue();
            return *this;
        }
        LuaVariable& operator = (const LuaValue& v) {
            value = v;
            setValue();
            return *this;
        }

        double asDouble()  { 
            updateValue();
            return value.d;
        }
        std::string asString() {
            updateValue();
            return value.s;
        }
        bool asBool()  {
            updateValue();
            return value.b;
        }
        void print() {
            value.print();
        }
    };
    
    struct LuaTableValue
    {            
        LuaValue value;
        int index;
        std::string table;
        LuaJIT& lua;

        LuaTableValue(LuaJIT& luajit, const std::string &name, int i) : lua(luajit), table(name), index(i) {}
        
        
        void setValue() {
            lua.GetGlobal(table);
            if(value.type == LuaValue::DOUBLE) lua.PushNumber(value.d);
            else if(value.type == LuaValue::BOOL) lua.PushBool(value.b);
            else lua.PushString(value.s);            
            lua.rawseti(index,-2); 
            lua.Pop(1);               
        }
        void updateValue() {
            lua.GetGlobal(table);
            lua.rawgeti(index,-1);        
            value.GetValue(lua);
            lua.Pop(1);                      
        }
        
        LuaTableValue& operator = (double v) {
            value.setDouble(v);
            setValue();
            return *this;
        }
        LuaTableValue& operator = (bool v) {
            value.setBool(v);
            setValue();
            return *this;
        }
        LuaTableValue& operator = (const string& v) {            
            value.setString(v);            
            setValue();
            return *this;
        }        
        LuaTableValue& operator = (const char * v) {            
            value.setString(std::string(v));            
            setValue();
            return *this;
        }        
        LuaTableValue& operator = (const LuaVariable& v) {
            value = v.value;            
            lua = v.lua;
            setValue();
            return *this;
        }
        LuaTableValue& operator = (const LuaValue& v) {
            value = v;
            setValue();
            return *this;
        }
    };

    struct LuaTableKey
    {            
        LuaValue value;
        std::string index;
        std::string table;
        LuaJIT& lua;

        LuaTableKey(LuaJIT& luajit, const std::string &name, const std::string& i) : lua(luajit), table(name), index(i) {}
        
        
        void setValue() {
            lua.GetGlobal(table);
            lua.PushString(index);
            if(value.type == LuaValue::DOUBLE) lua.PushNumber(value.d);
            else if(value.type == LuaValue::BOOL) lua.PushBool(value.b);
            else lua.PushString(value.s);            
            lua.rawset(-3); 
            lua.Pop(1);               
        }
        void updateValue() {
            lua.GetGlobal(table);
            lua.PushString(index);
            lua.rawgeti(-2);        
            value.GetValue(lua);
            lua.Pop(1);                      
        }
        
        LuaTableKey& operator = (double v) {
            value.setDouble(v);
            setValue();
            return *this;
        }
        LuaTableKey& operator = (bool v) {
            value.setBool(v);
            setValue();
            return *this;
        }
        LuaTableKey& operator = (const string& v) {            
            value.setString(v);            
            setValue();
            return *this;
        }        
        LuaTableKey& operator = (const char * v) {            
            value.setString(std::string(v));            
            setValue();
            return *this;
        }        
        LuaTableKey& operator = (const LuaVariable& v) {
            value = v.value;            
            lua = v.lua;
            setValue();
            return *this;
        }
        LuaTableKey& operator = (const LuaValue& v) {
            value = v;
            setValue();
            return *this;
        }
    };

    struct LuaTable
    {
        LuaJIT& lua;
        std::string name;

        LuaTable(LuaJIT & l, const std::string& t) : lua(l),name(t) {}

        LuaTableValue operator[](size_t i) {
            LuaTableValue v(lua,name,i);            
            v.updateValue();
            return v;
        }        
        LuaTableKey operator[](const std::string& i) {
            LuaTableKey v(lua,name,i);            
            v.updateValue();
            return v;
        }        
    };
}            
