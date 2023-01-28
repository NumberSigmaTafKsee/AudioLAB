%{
#include <string>
#include <algorithm>
#include <iostream>
%}

namespace std
{
    %naturalvar string;

    %typemap(in,checkfn="lua_isstring") string
    %{$1.assign(lua_tostring(L,$input),lua_rawlen(L,$input));%}

    %typemap(out) string
    %{ lua_pushlstring(L,$1.data(),$1.size()); SWIG_arg++;%}

    %typemap(in,checkfn="lua_isstring") const string& ($*1_ltype temp)
    %{temp.assign(lua_tostring(L,$input),lua_rawlen(L,$input)); $1=&temp;%}

    %typemap(out) const string&
    %{ lua_pushlstring(L,$1->data(),$1->size()); SWIG_arg++;%}

    // for throwing of any kind of string, string ref's and string pointers
    // we convert all to lua strings
    %typemap(throws) string, string&, const string&
    %{ lua_pushlstring(L,$1.data(),$1.size()); SWIG_fail;%}

    %typemap(throws) string*, const string*
    %{ lua_pushlstring(L,$1->data(),$1->size()); SWIG_fail;%}

    %typecheck(SWIG_TYPECHECK_STRING) string, const string& {
    $1 = lua_isstring(L,$input);
    }

    %typemap(in) string &INPUT=const string &;
    %typemap(in, numinputs=0) string &OUTPUT ($*1_ltype temp)
    %{ $1 = &temp; %}
    %typemap(argout) string &OUTPUT
    %{ lua_pushlstring(L,$1->data(),$1->size()); SWIG_arg++;%}
    %typemap(in) string &INOUT =const string &;
    %typemap(argout) string &INOUT = string &OUTPUT;

    class string
    {
    public:
        string();
        string(const char * s);
        string(const string& s);

        %extend {
            char __getitem__(size_t i) { return (*$self)[i]; }
            void __setitem__(size_t i, char c) { (*$self)[i] = c; }
            const char* __str__() { return $self->c_str(); }
            /*
            void mutate() { size_t n = randint(0,$self->size()); (*$self)[n] = randchar(); }
            void uniform_mutate(float p = 0.001) {
                for(size_t i = 0; i < $self->size(); i++)
                    if(flip(p)) (*$self)[i] = randchar();
            }
            void crossover(const string & mom, const string & dad) {
                assert(mom.size() == dad.size());
                resize(mom.size());
                for(size_t i = 0; i < mom.size(); i++)
                    ($self)[i] = flip()? mom[i]:dad[i];
            }
            */
            int32_t to_int32(int base =10) { return std::stoi($self->c_str(),nullptr,base); }
            int64_t to_int64(int base =10) { return std::stoll($self->c_str(),nullptr,base); }
            float   to_float()  { return std::stof(*$self); }
            double  to_double() { return std::stod(*$self); }

            void reverse() { std::reverse($self->begin(),$self->end()); }
            void sort() { std::sort($self->begin(),$self->end()); }
            void shuffle() { std::random_shuffle($self->begin(),$self->end()); }

            void getline() {
                std::getline(std::cin, *$self);
            }
            void fill(size_t i, char c) {
                $self->resize(i);
                for(size_t n = 0; n < i; n++) (*$self)[i] = c;
            }
        }

        const char* data() const;
        size_t size() const;
        bool empty();
        void clear();

        //void insert(size_t i, char c);
        //void insert(size_t i, size_t n, char c);
        
        //void insert(size_t i, const string& c);
                
        void erase(size_t i, size_t n =1);
        void push_back(size_t c);
        void pop_back();
        //void append(char c, size_t count=1);
        string replace(size_t pos, size_t count, const string &s);
        void resize(size_t i);                        
        string substr(size_t i, size_t n);

        string& operator = (const string & s);                
        //bool operator == (const string & b);

        size_t find(const string &s, size_t pos=0);
        size_t rfind(const string &s, size_t pos=0);
        int compare(const string & b);
    };
}