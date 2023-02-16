%{
#include <queue>

%}

namespace std {
    template<class T>
    struct queue
    {
        Queue();
        Queue(const Queue & que);

        bool empty();
        size_t size();
        T& front();
        T& back();
        void push(const T val);
        void pop();
    };
}    