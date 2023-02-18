%{
#include <stack>
%}

namespace std
{
    template<typename T>
    class stack
    {
    public:

        stack();
        stack(const stack&);

        stack& operator=(const stack& s);

        T& top();
        bool empty();
        size_t size();
        void push(const T& value);
        void pop();        
    };
}