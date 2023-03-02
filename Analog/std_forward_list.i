%{
#include <forward_list>
#include <vector>
%}

namespace std {
    template<class T>
    struct forward_list
    {

        forward_list(size_t size = 1);
        forward_list(size_t n, const T & value);

        %extend {        
            T& get(size_t index)
            {
                assert(index < $self->size());
                typename std::list<T>::iterator it = $self->begin();
                for(size_t i = 0; i < index; i++) it++;
                return *it;
            }
            void set(size_t index, const T val)
            {
                assert(index < $self->size());
                typename std::list<T>::iterator it = $self->begin();
                for(size_t i = 0; i < index; i++) it++;
                *it = val;
            }
            T    __getitem__(size_t key) { return $self->get(key); }
            void __setitem__(const size_t key, const T val) { $self->set(key,val); }

            void insert_vector(size_t index, const std::vector<T> & a)
            {
                assert(index < $self->size());    
                typename std::list<T>::iterator it = $self->begin();        
                for(size_t i = 1; i < index; i++) it++;
                $self->insert(it,a.v.begin(),a.v.end());
            }   
            void insert_at(size_t index, const T val)
            {        
                assert(index < $self->size());    
                typename std::list<T>::iterator it = $self->begin();        
                for(size_t i = 1; i < index; i++) it++;
                $self->insert(it,val);
            }
            T remove_at(size_t index)
            {
                assert(index < $self->size());
                typename std::list<T>::iterator it = $self->begin();
                for(size_t i = 1; i < index; i++) it++;
                T r = *it;
                $self->erase(it);
                return r;
            }                     
        }

        void push_front(const T val);        
        void pop_front();
        bool    empty();
        size_t  max_size();
        void unique();
        T front();
        void resize(size_t s);
        void clear();
        size_t size();
        void swap(forward_list & x);
        void reverse();
    };
}    
