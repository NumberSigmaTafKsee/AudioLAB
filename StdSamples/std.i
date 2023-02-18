%module std
%{
#include <vector>
#include <list>
#include <map>
%}

%include "stdint.i"
%include "std_math.i"
%include "std_vector.i"
%include "std_map.i"
%include "std_string.i"
%include "lua_fnptr.i"

%inline %{
    bool operator < (SWIGLUA_REF a, SWIGLUA_REF b) {
        return a.ref < b.ref;
    }
    bool operator == (SWIGLUA_REF a, SWIGLUA_REF b) {
        return a.L == b.L && a.ref == b.ref;
    }
%}

%template(lua_map) std::map<std::string,SWIGLUA_REF>;
%template(lua_vector) std::vector<SWIGLUA_REF>;

%inline %{
    template<class T>
    struct CircularBuffer : public std::vector<T>
    {
        size_t r,w;
        CircularBuffer(size_t n) {
            std::vector<T>::resize(n);
            r = 0;
            w = 0;
        }
        void insert(T data) {
            (*this)[w++] = data;
            w = w % std::vector<T>::size();
        }
        T read() {
            T x = (*this[r++])
            r = r % std::vector<T>::size();
            return x;
        }
    };

    template<class Data>
    struct List : public std::list<Data>
    {
        List() = default;

        size_t count() { return this->size(); }

        Data operator[](size_t i) {
            return getAt(i);
        }
        Data __getitem__(size_t i) { return getAt(i); }
        void __setitem__(size_t i, Data d) { setAt(i,d); }

        Data getAt(size_t cnt) {
            auto it = this->begin();
            if(cnt >= count()) cnt = count()-1;
            for(size_t i = 1; i <  cnt; i++)
            {
                it++;
            }
            return (*it);
        }
        void setAt(size_t idx, Data data) {
            auto it = this->begin();
            if(idx >= count()) return;
            for(size_t i = 1; i <  idx; i++)
            {
                it++;
            }
            (*it) = data;
        }

        void pushBack(Data d) {
            this->push_back(d);
        }
        Data popBack() {
            Data x = this->back();
            this->pop_back();
            return x;
        }
        void pushFront(Data d) {
            this->push_front(d);
        }
        Data popFront() {
            Data x = this->front();
            this->pop_front();
            return x;        
        }
        void insertAt(size_t idx, Data d) 
        {
            auto it = this->begin();
            if(idx >= count()) return;
            for(size_t i = 1; i <  idx; i++)
            {
                it++;
            }
            this->insert(it,d);
        }
        void removeAt(size_t idx) {
            auto it = this->begin();
            if(idx >= count()) return;
            for(size_t i = 1; i <  idx; i++)
            {
                it++;
            }
            this->erase(it);
        }
        List<Data> split(size_t idx) {
            List<Data> splitlist;
            auto it = this->begin();
            if(idx >= count()) return splitlist;
            for(size_t i = 1; i <  idx; i++)
            {
                it++;
            }
            std::copy(this->begin(),it,splitlist.begin());
            this->erase(this->begin(),it);
            return splitlist;
        }
        void mergeBack(const List<Data> & l) {
            std::copy(l.begin(),l.end(),this->end());
        }
        void mergeFront(const List<Data> & l) {
            std::copy(l.begin(),l.end(),this->begin());
        }
    };
    template<class Key,class Data>
    struct Dict : public std::map<Key,Data>
    {
        Data  null;

        Dict() = default;

        Data operator[](const Key & k) {
            if(this->find(k) == this->end()) return null;
            return (*this)[k];
        }
        bool findKey(const Key & k) {
            return (this->find(k) != this->end());
            
        }
        void insert(const Key & k, const Data & d) {
            (*this)[k] = d;
        }
        void remove(const Key & k) {
            auto i = this->find(k);
            if(i != this->end()) {
                this->erase(i);
            }
        }
        bool isNull(const Data & d) {
            return d == null;
        }
    };

%}

%template(lua_list) List<SWIGLUA_REF>;
%template(lua_dict) Dict<std::string,SWIGLUA_REF>;