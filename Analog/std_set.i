%{
#include <set>
%}

namespace std 
{

    template<typename T>
    class set
    {
    public:
        set();
        set(const set& other);

        set<T>& operator = (const set& s);

        bool empty() const;
        size_t size() const;
        void clear();

        %extend {
            void insert(const T & val) {
                $self->insert(val);
            }
            size_t erase(const T & val) {
                return $self->erase(val);
            }
            bool find(const T & val) {
                return $self->find(val) != $self->end();
            }
        }
    }
}