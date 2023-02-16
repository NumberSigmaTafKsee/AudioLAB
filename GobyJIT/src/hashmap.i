%module hashmap
%{
#include "rigtorp/HashMap.h"
%}

%include "std_string.i"
%include "lua_fnptr.i"

namespace std
{
    template<class T1, class T2>
    struct pair
    {
        T1 first;
        T2 second;
        pair();
        pair(const T1&, const T2&);

        pair& operator = (const pair& p);
        void swap(pair&);
    };
}

%inline %{

    struct Ref
    {
        SWIGLUA_REF ref;

        Ref() = default;
        Ref(SWIGLUA_REF r) : ref(r) {}

        Ref& operator = (const Ref& r) { 
            ref = r.ref;
            return *this;
        }        
    };
    // it doesnt make much sense but to stop it from complaining
    bool operator < (const Ref & a, const Ref & b) {
        return a.ref.ref < b.ref.ref;
    }
    bool operator == (const Ref & a, const Ref & b) {
        return a.ref.L == b.ref.L && a.ref.ref == b.ref.ref;
    }
%}

namespace rigtorp {
    template <typename Key, typename T, typename Hash = std::hash<Key>,
            typename KeyEqual = std::equal_to<void>,
            typename Allocator = std::allocator<std::pair<Key, T>>>
    class HashMap {
    public:
        using key_type = Key;
        using mapped_type = T;
        using value_type = std::pair<Key, T>;
        using size_type = size_t;
        using hasher = Hash;
        using key_equal = KeyEqual;
        using allocator_type = Allocator;
        using reference = value_type &;
        using const_reference = const value_type &;
        using buckets = std::vector<value_type, allocator_type>;

        template <typename ContT, typename IterVal> struct hm_iterator;
        using iterator = hm_iterator<HashMap, value_type>;
        using const_iterator = hm_iterator<const HashMap, const value_type>;  

    public:
    HashMap(size_type bucket_count, key_type empty_key,
            const allocator_type &alloc = allocator_type());

    HashMap(const HashMap &other, size_type bucket_count);

    allocator_type get_allocator() const noexcept;
    

    iterator begin() noexcept;
    const_iterator begin() const noexcept;

    const_iterator cbegin() const noexcept;
    iterator end() noexcept;

    const_iterator end() const noexcept;
    const_iterator cend() const noexcept;

    bool empty() const noexcept;
    size_type size() const noexcept;
    size_type max_size() const noexcept;
    void clear() noexcept;

    std::pair<iterator, bool> insert(const value_type &value);
    std::pair<iterator, bool> insert(value_type &&value);

    
    //void erase(iterator it);
    size_type erase(const key_type &key);

    //template <typename K> size_type erase(const K &x) { return erase_impl(x); }

    void swap(HashMap &other) noexcept;

    mapped_type &at(const key_type &key);

    //template <typename K> mapped_type &at(const K &x) { return at_impl(x); }

    const mapped_type &at(const key_type &key) const;
    //template <typename K> const mapped_type &at(const K &x) const;

    %extend {
        mapped_type __getitem__(const key_type & key) { return (*$self)[key];}
        void __setitem__(const key_type & key, const mapped_type v) { (*$self)[key] = v; }
    }

    size_type count(const key_type &key) const;
    template <typename K> size_type count(const K &x) const;
    iterator find(const key_type &key);
    template <typename K> iterator find(const K &x);
    const_iterator find(const key_type &key) const;
    template <typename K> const_iterator find(const K &x) const;
    size_type bucket_count() const noexcept;
    size_type max_bucket_count() const noexcept;

    void rehash(size_type count);
    void reserve(size_type count);
    
    hasher hash_function() const;
    key_equal key_eq() const;

    };

   template <typename ContT, typename IterVal> struct hm_iterator 
   {
        using difference_type = std::ptrdiff_t;
        using value_type = IterVal;
        using pointer = value_type *;
        using reference = value_type &;
        using iterator_category = std::forward_iterator_tag;

        bool operator==(const hm_iterator &other) const;
        bool operator!=(const hm_iterator &other) const;
        hm_iterator &operator++();

        reference operator*() const;
        pointer operator->() const;
    };

} // namespace rigtorp

%template(lua_hashmap) rigtorp::HashMap<std::string,SWIGLUA_REF>;