
%{
#include <map>
%}

%inline %{
  namespace std {
    template<typename Key,  typename Data>
        struct map_iterator 
        {
            typename std::map<Key,Data>::iterator _iter;
            std::map<Key,Data> _map;

            map_iterator(const std::map<Key,Data> & m, const typename std::map<Key,Data>::iterator & i) {
                _map = m;
                _iter = i;
            }

            void next() { if(_iter != _map.end()) _iter++; }
            void prev() { if(_iter != _map.begin()) _iter--; }
            Data& get() { return *_iter;}

            void insert(const Data& value) { _map.insert(_iter,value); }
            bool find(const Key& key) { return _map.find(key) != _map.end(); }   
            
        };
    }
%}
namespace std {
        
        
        template<typename Key, typename Data>
        class map
        {
        public:
            map();
            map(const map& m);

            Data& operator[](const Key& key);

            %extend {
                Data __getitem(const Key & key) {
                    return (*$self)[key];
                }
                void __setitem(const Key & key, const Data value)
                {
                    (*$self)[key] = value;
                }

                std::map_iterator<Key,Data> get_begin() {
                    return std::map_iterator<Key,Data>(*$self,$self->begin());
                }            
                std::map_iterator<Key,Data> get_end() {
                    return std::map_iterator<Key,Data>(*$self,$self->end());
                }            
            }


            map<Key,Data>& operator = (const map& m);
                        
            Data& at(Key & key);
            bool empty();
            size_t size() const;
            void clear();            
            void erase(const Key & key);
            size_t count(const Key & key) const;
        };        
}

template<typename Key, typename Data>
void copy(std::map<Key,Data>& dst, std::map<Key,Data> & src) {
    std::copy(src.begin(),src.end(),dst.begin());
}
template<typename Key, typename Data>
void copy(std::map_iterator<Key,Data>& d1, std::map_iterator<Key,Data> & s1, std::map_iterator<Key,Data> & s2) {
    std::copy(s1._iter,s2._iter,d1._iter);
}