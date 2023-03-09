template<class Data>
struct List : public std::list<Data>
{
    List() = default;

    size_t count() { return this->size(); }

    Data& operator[](size_t i) {
        return getAt(i);
    }
    Data& getAt(size_t cnt) {
        auto it = this->begin();
        if(cnt >= count()) cnt = count()-1;
        for(size_t i = 0; i <  cnt; i++)
        {
            it++;
        }
        return (*it);
    }
    void setAt(size_t idx, const Data & data) {
        auto it = this->begin();
        if(cnt >= count()) cnt = count()-1;
        for(size_t i = 0; i <  cnt; i++)
        {
            it++;
        }
        (*it) = data;
    }

    void pushBack(const Data & d) {
        this->push_back(d);
    }
    Data popBack() {
        Data x = this->back();
        this->pop_back();
        return x;
    }
    void pushFront(const Data & d) {
        this->push_front(d);
    }
    Data popFront() {
        Data x = this->front();
        this->pop_front();
        return x;        
    }
    void insertAt(size_t idx, const Data & d) 
    {
        auto it = this->begin();
        if(cnt >= count()) cnt = count()-1;
        for(size_t i = 0; i <  cnt; i++)
        {
            it++;
        }
        this->insert(it,d);
    }
    void removeAt(size_t idx) {
        auto it = this->begin();
        if(cnt >= count()) cnt = count()-1;
        for(size_t i = 0; i <  cnt; i++)
        {
            it++;
        }
        this->erase(it);
    }
    List<Data> split(size_t idx) {
        List<Data> splitlist;
        auto it = this->begin();
        if(cnt >= count()) cnt = count()-1;
        for(size_t i = 0; i <  cnt; i++)
        {
            it++;
        }
        std::copy(this->begin(),it,splitlist.begin());
        this->erase(this->begin(),it);
        return splitlist
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

    Data& operator[](const Key & k) {
        if(!find(key)) return null;
        return (*this)[key];
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

struct SampleQueue
{
    std::queue<double> q;

    SampleQueue() = default;
    SampleQueue(const SampleQueue & que) { q = que.q; }
    ~SampleQueue()= default;

    SampleQueue& operator = (const SampleQueue & que) {
        q = que.q;
        return *this;
    }
    bool isEmpty() const { return q.empty(); }
    size_t getSize() const { return q.size(); }

    void push(double v)
    {
        q.push(v);
    }
    double pop(double v) {
        double r = q.front();
        q.pop();
        return r;
    }
};
