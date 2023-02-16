%{
#include <list>
%}

%inline %{
  namespace std {
    template<typename T>
        struct list_iterator 
        {
            typename std::list<T>::iterator _iter;
            std::list<T> _list;;

            map_iterator(const std::list<T> & l, const typename std::list<T>::iterator & i) {
                _list = l;
                _iter = i;
            }

            void next() { if(_iter != _map.end()) _iter++; }
            void prev() { if(_iter != _map.begin()) _iter--; }
            Data& get() { return *_iter;}

            //void insert(const Data& value) { _map.insert(_iter,value); }
            //bool find(const Key& key) { return _map.find(key) != _map.end(); }   
            
        };
    }
%}
namespace std {

    typename<T>
    class list 
    {
    publc:

        list();
        list(const list& l);

         %extend {
                T __getitem(size_t i) {
                    iterator it;
                    size_t p=0;
                    for(it = begin(); it != end(); it++)
                    {
                        if(p == i) break;                        
                        p++;
                    }
                    return *it;
                }
                void __setitem(size_t i, const T value)
                {
                    iterator it;
                    size_t p = 0;
                    for(it = begin(); it != end(); it++)
                    {
                        if(p == i) break;                        
                        p++;
                    }
                    insert(it,value);
                }

                std::list_iterator<Key,Data> get_begin() {
                    return std::list_iterator<T>(*$self,$self->begin());
                }            
                std::list_iterator<Key,Data> get_end() {
                    return std::list_iterator<T>(*$self,$self->end());
                }           

                list<T>& __add__(const list& b) { std::copy(b.begin(),b.end(),$self->end(); return *$self; }
                list<T>& __sub__(const list& b) { 
                    for(iterator i = b.begin(); i != b.end(); i++) $self->remove(*i);
                    return *$self; 
                }
                bool __eq__(const list& b) { return *$self == b; } 
                bool __lt__(const list& b) { return *$self == b; }
                bool __le__(const list& b) { return *$self <= b; }
            }

            list<T>& operator = (const list& m);
                    
            T& at(Key & key);
            bool empty();
            size_t size() const;
            void resize(size_t size);

            T& front() const;
            T& back() const;

            void push_front(const T& value);
            void push_back(const T& value);

            void pop_front();
            void pop_back();

            void insert(const_iterator position, const T& value);
            void insert(const_iterator position, size_t n, const T& value);

            void clear();            
            
            void erase(const_iterator &pos);
            
            void splice(const_iterator &pos, list& x);
            void merge(list& x);
            void sort();
            void reverse() noexcept;
            
            void remove(const T& val);            
    };

}


template<typename T>
void copy(std::list<T>& dst, std::list<T> & src) {
    std::copy(src.begin(),src.end(),dst.begin());
}
template<typename T>
void copy(std::list_iterator<T>& d1, std::list_iterator<T> & s1, std::list_iterator<T> & s2) {
    std::copy(s1._iter,s2._iter,d1._iter);
}