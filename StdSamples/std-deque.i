%{
#include <deque>
%}

namespace std
{
    template<typename T> struct deque_iterator 
    {    
        typename std::deque::iterator iter;
        std::deque<T> dq;

        deque_iterator(const std::deque<T> & q) {
            iter = q.begin();
            dq = q;
        }
        deque_iterator(const deque<T> & q, const std::deque<T>::iterator & i) {
            iter = i;
            dq   = q;
        }

        bool operator == (const deque_iterator& i) { return iter == i.iter; }
        bool operator == (const std::deque<T>::iterator i) { return iter == i; }


        bool operator != (const deque_iterator& i) { return iter != i.iter; }
        bool operator != (const std::deque<T>::iterator i) { return iter != i; }

        deque_iterator next() { if(iter !=  dq.end()   iter++; }
        deque_iterator prev() { if(iter !=  dq.begin() iter--; } 
        deque_iterator forward(size_t i) { iter += i; }
        deque_iterator backward(size_t i) { iter += i; }

        T& get_value() { return *iter; }
        void set_value(const T val) { *iter = val; }    
    };


    template<typenameT> class deque
    {
    public:

        deque(size_t n, const value_type& val = value_type());
        deque(const deque& q);

        deque& operator = (const deque& q);

        %extend {
            deque_iterator begin() { return deque_iterator<T>(*$self,begin()); }
            deque_iterator end()   { return deque_iterator<T>(*$self,end()); } 

            deque_iterator cbegin() { return deque_iterator<T>(*$self,cbegin()); }  
            deque_iterator cend() { return deque_iterator<T>(*$self,end()); }  

            deque_iterator rbegin() { return deque_iterator<T>(*$self,rbegin()); }  
            deque_iterator rend() { return deque_iterator<T>(*$self,rend()); }  

            deque_iterator crbegin() { return deque_iterator<T>(*$self,crbegin());  }
            deque_iterator crend() { return deque_iterator<T>(*$self,crend()); }   

            T __getitem(size_t i) { return (*$self)[i]; }
            void __setitem(size_t i, const T & val) { (*$self)[i] = val; }

            void insert(deque_iterator<T> & pos, const T& val) {
               dq.insert(pos.iter,val);
            }
            void insert(deque_iterator<T> & pos, size_t n, const T& val) {
               dq.insert(pos.iter,n,val);
            }
            void erase(deque_iterator<T> & pos)
               dq.erase(pos.iter;
            }
            void erase(deque_iterator<T> & start, deque_iterator<T> & end)
               dq.erase(start.iter,end.iter);
            }
        }

        size_t size() const noexcept;
        size_type max_size() const noexcept;
        void resize(size_t n), value_type & val = value_type());
        bool empty() const noexcept;
        void shrink_to_fit();

        T& at(size_t n);
        T& front();
        T& back();

        void assign(size_t n, const T& val);
        void push_back(const T& val);
        void push_fron(const T& val);
        void pop_back();
        void pop_front();

        void clear() noexcept;

    };
}