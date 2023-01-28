#pragma once

#include <iostream>
#include <cassert>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <condition_variable>
#include <memory>
#include <mutex>
#include <unordered_map>
#include <deque>
#include <unordered_map>
#include <future>
#include <thread>

using namespace Eigen;

// A matrix of ints with a dynamic size, Use it when the size is not known
typedef Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> matXXi;

// enum for operation types
enum class operationType
{
    NA,
    addition,
    negative,
    multiply,
    dotproduct,
    sigmoid,
    log,
    sum
};

// Node types
enum nodeType
{
    NA,
    minimizer,
    operation,
    variable,
    placeholder
};

/* a wrapper class to avoid data races when accessing data 
from raw pointers by different threads. It is move/copy constructable.
Pointer access is atomic. */
template <typename T>
class Locking_ptr
{
public:
    Locking_ptr(T *ptr);
    Locking_ptr(T *ptr, std::mutex *mtx);
    Locking_ptr(Locking_ptr<T> const &other);
    Locking_ptr(Locking_ptr<T> &&other);

    ~Locking_ptr();

    Locking_ptr<T> &operator=(Locking_ptr<T> &&other);
    Locking_ptr<T> &operator=(Locking_ptr<T> const &other);
    T *operator->() const;
    T &operator*() const;
    bool operator==(Locking_ptr<T> const &rhs);

    void reset(T *ptr);
    T *get() const;

protected:
    std::atomic<std::mutex *> Mtx_;

private:
    T *ptr_;
};


/* a wrapper class to avoid data races when accessing data 
from shared smart pointers by different threads. 
It is move/copy constructable. Pointer access is atomic. */
template <typename T>
class Locking_shared_ptr
{
public:
    Locking_shared_ptr(T *ptr);
    Locking_shared_ptr(T *ptr, std::mutex *mtx);
    Locking_shared_ptr(std::shared_ptr<T> ptr);

    Locking_shared_ptr(Locking_shared_ptr<T> const &other);
    Locking_shared_ptr(Locking_shared_ptr<T> &&other);

    ~Locking_shared_ptr();

    Locking_shared_ptr<T> &operator=(Locking_shared_ptr<T> &&other);
    Locking_shared_ptr<T> &operator=(Locking_shared_ptr<T> const &other);
    T *operator->() const;
    T &operator*() const;
    bool operator==(Locking_shared_ptr<T> const &rhs);

    void reset(T *ptr);
    T *get() const;

protected:
    std::atomic<std::mutex *> Mtx_;

private:
    std::shared_ptr<T> ptr_;
};

/* a wrapper class to avoid data races when accessing data 
from unique smart pointers by different threads. 
It is only movable. Pointer access is atomic. */
template <typename T>
class Locking_unique_ptr
{
public:
    Locking_unique_ptr(T *ptr);
    Locking_unique_ptr(T *ptr, std::mutex *mtx);
    Locking_unique_ptr(std::unique_ptr<T> &&ptr);
    Locking_unique_ptr(Locking_unique_ptr<T> &&other);

    ~Locking_unique_ptr();

    Locking_unique_ptr<T> &operator=(Locking_unique_ptr<T> &&other);
    T *operator->() const;
    T &operator*() const;
    bool operator==(Locking_unique_ptr<T> const &rhs);

    void reset(T *ptr);
    T *get() const;

protected:
    std::atomic<std::mutex *> Mtx_;

private:
    Locking_unique_ptr<T> &operator=(Locking_unique_ptr<T> const &other) = delete;
    Locking_unique_ptr(Locking_unique_ptr<T> const &other) = delete;
    std::unique_ptr<T> ptr_;
};

// Base node class
class BaseNode
{
public:
    virtual ~BaseNode(){};
    void addInputs(BaseNode *n);
    void eraseInput(BaseNode *n);
    void addConsumers(BaseNode *n);
    void eraseConsumer(BaseNode *n);
    void setName(std::string n);

    // get output value of this node
    template <typename T>
    Locking_shared_ptr<T> getValue();

    // get total gradient from node's consumer
    template <typename T>
    T getGradient();

    // set gradient from consumer
    template <typename T>
    void setGrad(T t);

    // make this abstract base class
    virtual void clearGrads() = 0;
    virtual void compute() = 0;
    virtual void gradient() = 0;

    nodeType getNodeType();
    operationType getOperationType();
    std::vector<Locking_ptr<BaseNode>> &getConsumers();
    std::vector<Locking_ptr<BaseNode>> &getInputs();
    std::string getName();

    // keep the size of consumers as an atomic data
    std::atomic_int consumerSize_{0};
    std::mutex Mtx_;     // for ptrs
    std::mutex nodeMtx_; // for data

protected:
    std::string _name = " ";
    nodeType _nType;       // node type
    operationType _opType; // type if node is operation

private:
    std::vector<Locking_ptr<BaseNode>> _consumers = {}; // parent nodes
    std::vector<Locking_ptr<BaseNode>> _inputs = {};    // child nodes
};

/* Class for nodes of the computational graph; 
each node is one of the three:
 - an operation
 - a variable
 - a placeholder
*/
template <typename T>
class Node : public BaseNode
{
public:
    ~Node(){};
    Locking_shared_ptr<T> getValue();
    T getGradient();

    void setValue(T &&t);
    void setGrad(T t);
    void clearGrads();

private:
    std::condition_variable cond_;
    // ouput might be shared
    Locking_shared_ptr<T> _output = Locking_shared_ptr<T>(nullptr, &(this->Mtx_));
    std::vector<Locking_shared_ptr<T>> _grad;

    std::atomic<bool> _dataAvailable{false};
    std::atomic<bool> _gradientAvailable{false};
};

// A class for variables of type T
template <typename T>
class Variable : public Node<T>
{
public:
    Variable(T &&a);
    Variable(Variable<T> &&other);
    Variable<T> &operator=(Variable<T> &&other);

    void compute();
    void gradient();
    void updateValue(float lr);

private:
    Variable(Variable<T> const &other) = delete;
    Variable<T> &operator=(Variable<T> const &other) = delete;
};

// A class for placeholders for values of type T
template <typename T>
class Placeholder : public Node<T>
{
public:
    Placeholder(std::string n);

    void compute();
    void gradient();
};

// forward declaration
template <typename T>
class Minimizer;

class GradientDescentOptimizer
{
public:
    GradientDescentOptimizer(float lr);

    // compute gradients
    void computeGradients(BaseNode *loss);

    // get node list doing level order traversal
    std::vector<Locking_ptr<BaseNode>> getNodeQueue(BaseNode *loss);

    template <typename T>
    Minimizer<T> minimize(BaseNode *loss);

    float learningRate_;
};
class Session
{
public:
    // Runs calculation of the node and returns the output value for the node;
    // Takes input data for placeholders with an unordered map using placeholder's name
    template <typename T>
    void Run(BaseNode *n, std::unordered_map<std::string, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> *> feed);
    std::vector<Locking_ptr<BaseNode>> getNodesList();
    // Perform post-order traversal
    void updateNodesList(BaseNode *n);

private:
    std::vector<Locking_ptr<BaseNode>> _nodesList = {};
};

// Class for operations with data value of type T
template <typename T>
class Operation : public Node<T>
{
public:
    Operation();
    virtual void compute() = 0;
    virtual void gradient() = 0;
};

template <typename T>
class UnaryOperation : public Operation<T>
{
public:
    UnaryOperation(BaseNode *rhs);
};

template <typename T>
class BinaryOperation : public Operation<T>
{
public:
    BinaryOperation(BaseNode *lhs, BaseNode *rhs);
};

// ---- Operations ----

// addition operation with T return type value, T1 and T2 input type value
template <typename T, typename T1, typename T2>
class Add : public BinaryOperation<T>
{
public:
    Add(BaseNode &a, BaseNode &b);
    Add(BaseNode *a, BaseNode *b);

    void compute();
    void gradient();
};

// negative operation
template <typename T>
class Negative : public UnaryOperation<T>
{
public:
    Negative(BaseNode &a);
    Negative(BaseNode *a);

    void compute();
    void gradient();
};

// Elementwise multiplication
template <typename T, typename T1, typename T2>
class Multiply : public BinaryOperation<T>
{
public:
    Multiply(BaseNode &a, BaseNode &b);
    Multiply(BaseNode *a, BaseNode *b);

    void compute();
    void gradient();
};

// Matrix multiply
template <typename T, typename T1, typename T2>
class MatMultiply : public BinaryOperation<T>
{
public:
    MatMultiply(BaseNode &a, BaseNode &b);
    MatMultiply(BaseNode *a, BaseNode *b);

    void compute();
    void gradient();
};

// Vector dot product
template <typename T, typename T1, typename T2>
class Dot : public BinaryOperation<T>
{
public:
    Dot(BaseNode &a, BaseNode &b);
    Dot(BaseNode *a, BaseNode *b);

    void compute();
    void gradient();
};

// Element-wise sigmoid function
template <typename T>
class Sigmoid : public UnaryOperation<T>
{
public:
    Sigmoid(BaseNode &a);
    Sigmoid(BaseNode *a);

    void compute();
    void gradient();
};

// Element-wise Log operation
template <typename T>
class Log : public UnaryOperation<T>
{
public:
    Log(BaseNode &a);
    Log(BaseNode *a);

    void compute();
    void gradient();
};

// Reduce sum operation
template <typename T>
class Sum : public UnaryOperation<T>
{
public:
    Sum(BaseNode &a, int axis);
    Sum(BaseNode *a, int axis);

    void compute();
    void gradient();

private:
    //  axis 0 is columnwise, axis 1 is rowwise
    int _axis = 0;
};

/* Minimization Operation. Minimization node doesn't have any inputs or consumers. 
It's only move constructable. */
template <typename T>
class Minimizer : public Operation<T>
{
public:
    Minimizer(GradientDescentOptimizer *grd, BaseNode *loss);
    Minimizer(Minimizer<T> &&other);
    Minimizer<T> &operator=(Minimizer<T> &&other);

    void compute();
    void gradient();

private:
    // copy/assgn deleted
    Minimizer(Minimizer<T> &other) = delete;
    Minimizer<T> &operator=(Minimizer<T> &other) = delete;
    // cash the optimization class to access methods
    Locking_ptr<GradientDescentOptimizer> grdOpt_ = Locking_ptr<GradientDescentOptimizer>(nullptr, &(this->Mtx_));
    Locking_ptr<BaseNode> loss_ = Locking_ptr<BaseNode>(nullptr, &(this->Mtx_));
    std::mutex gMtx_;
    float learningRate_;
};

// a container class that ownes the nodes
class Graph
{
public:
    // add nodes of unary operation to graph
    template <template <typename> class U, typename T>
    BaseNode *addNodeOne(std::unique_ptr<U<T>> n);

    // add nodes of binary operation to graph
    template <template <typename, typename, typename> class U, typename T, typename T1, typename T2>
    BaseNode *addNodeTwo(std::unique_ptr<U<T, T1, T2>> n);

    // get nodes from graph
    std::vector<std::unique_ptr<BaseNode>> &getNodes();

private:
    std::vector<std::unique_ptr<BaseNode>> _baseNodes;
};

class NN
{
public:
    NN();

    // variable node for learnable parameters of the model
    template <typename T>
    BaseNode *variable(Matrix<T, Dynamic, Dynamic> &&t);

    // placeholder node for contant prameters of the node that needs to be fed later with data
    template <typename T>
    BaseNode *placeholder(std::string n);

    // add operation node
    template <typename T>
    BaseNode *add(BaseNode *a, BaseNode *b);

    // negative operation node
    template <typename T>
    BaseNode *negative(BaseNode *a);

    // element-wise multiplication operation
    template <typename T>
    BaseNode *multiply(BaseNode *a, BaseNode *b);

    // matrix multiplication node
    template <typename T>
    BaseNode *matmultiply(BaseNode *a, BaseNode *b);

    // dot product operation node
    template <typename T>
    BaseNode *dot(BaseNode *a, BaseNode *b);

    // sigmoid operation node
    template <typename T>
    BaseNode *sigmoid(BaseNode *a);

    // element-wise log operation node
    template <typename T>
    BaseNode *log(BaseNode *a);

    // reduce sum operation node
    template <typename T>
    BaseNode *sum(BaseNode *a, int axis);

    // run session for operation node 
    template <typename T>
    void run(BaseNode *n, std::unordered_map<std::string, Matrix<T, Dynamic, Dynamic> *> feed);

    // helper function to check gradient calculations
    template <typename T>
    void checkAllGradient(BaseNode *loss, std::unordered_map<std::string, Matrix<T, Dynamic, Dynamic> *> feed);

private:
    // factory method to create unary operations
    template <typename T, template <typename> class U>
    BaseNode *UnaryOperation(BaseNode *a);

    // factory method to create binary operations
    template <typename T, template <typename, typename, typename> class U>
    BaseNode *BinaryOperation(BaseNode *a, BaseNode *b);

    // helper function to check gradient of a node numerically
    template <typename T>
    void checkGradient(BaseNode *n, BaseNode *loss, std::unordered_map<std::string, Matrix<T, Dynamic, Dynamic> *> feed);

    // helper function to swap nodes in and out of graph
    void swapNodes(BaseNode *a, BaseNode *b);

    Graph _graph;
    Session _session;
};

NN::NN()
{
    _graph = Graph();
    _session = Session();
};

template <typename T>
void NN::run(BaseNode *n, std::unordered_map<std::string, Matrix<T, Dynamic, Dynamic> *> feed)
{
    _session.Run<T>(n, feed);
}

template <typename T, template <typename> class U>
BaseNode *NN::UnaryOperation(BaseNode *a)
{
    // Create a shared_ptr of new U class
    auto v = std::unique_ptr<U<Matrix<T, Dynamic, Dynamic>>>(new U<Matrix<T, Dynamic, Dynamic>>(a));
    // Copy the shared_ptr into graph to make sure its life is extended as necessary
    return _graph.addNodeOne<U, Matrix<T, Dynamic, Dynamic>>(std::move(v));
}

template <typename T, template <typename, typename, typename> class U>
BaseNode *NN::BinaryOperation(BaseNode *a, BaseNode *b)
{
    // Create a shared_ptr of new U class
    auto c = std::unique_ptr<U<Matrix<T, Dynamic, Dynamic>, Matrix<T, Dynamic, Dynamic>, Matrix<T, Dynamic, Dynamic>>>(new U<Matrix<T, Dynamic, Dynamic>, Matrix<T, Dynamic, Dynamic>, Matrix<T, Dynamic, Dynamic>>(a, b));
    // Copy the shared_ptr into graph
    return _graph.addNodeTwo<U, Matrix<T, Dynamic, Dynamic>, Matrix<T, Dynamic, Dynamic>, Matrix<T, Dynamic, Dynamic>>(std::move(c));
}

template <typename T>
BaseNode *NN::variable(Matrix<T, Dynamic, Dynamic> &&t)
{
    //std::cout << "nn variable constructor" << std::endl;
    auto v = std::unique_ptr<Variable<Matrix<T, Dynamic, Dynamic>>>(new Variable<Matrix<T, Dynamic, Dynamic>>(std::move(t)));
    return _graph.addNodeOne<Variable, Matrix<T, Dynamic, Dynamic>>(std::move(v));

}

template <typename T>
BaseNode *NN::placeholder(std::string n)
{
    //std::cout << "nn placeholder constructor" << std::endl;
    auto plc = std::unique_ptr<Placeholder<Matrix<T, Dynamic, Dynamic>>>(new Placeholder<Matrix<T, Dynamic, Dynamic>>(n));
    return _graph.addNodeOne<Placeholder, Matrix<T, Dynamic, Dynamic>>(std::move(plc));

}

template <typename T>
BaseNode *NN::add(BaseNode *a, BaseNode *b)
{

    return BinaryOperation<T, Add>(a, b);
}

template <typename T>
BaseNode *NN::negative(BaseNode *a)
{
    return UnaryOperation<T, Negative>(a);
}

template <typename T>
BaseNode *NN::multiply(BaseNode *a, BaseNode *b)
{
    return BinaryOperation<T, Multiply>(a, b);
}

template <typename T>
BaseNode *NN::matmultiply(BaseNode *a, BaseNode *b)
{
    return BinaryOperation<T, MatMultiply>(a, b);
}

template <typename T>
BaseNode *NN::dot(BaseNode *a, BaseNode *b)
{
    return UnaryOperation<T, Dot>(a);
}

template <typename T>
BaseNode *NN::sigmoid(BaseNode *a)
{
    return UnaryOperation<T, Sigmoid>(a);
}

template <typename T>
BaseNode *NN::log(BaseNode *a)
{
    return UnaryOperation<T, Log>(a);
}

template <typename T>
BaseNode *NN::sum(BaseNode *a, int axis)
{
    return UnaryOperation<T, Sum>(a);
}

template <typename T>
void NN::checkAllGradient(BaseNode *loss, std::unordered_map<std::string, Matrix<T, Dynamic, Dynamic> *> feed)
{
    _session.updateNodesList(loss);
    auto nodes = _session.getNodesList();
    for (auto &n : nodes)
    {
        checkGradient(n.get(), loss, feed);
    }
    std::cout << "All gradients are correct!\n";
}

// Check node gradient numerically
template <typename T>
void NN::checkGradient(BaseNode *n, BaseNode *loss, std::unordered_map<std::string, Matrix<T, Dynamic, Dynamic> *> feed)
{
    typedef Matrix<T, Dynamic, Dynamic> matxxT;
    // define epsilon and threshold
    long double E = 1.0e-14;
    long double err = 1.0e-3;
    matxxT grad;
    // create +E and -E values
    matxxT value1 = n->getValue<Matrix<T, Dynamic, Dynamic>>()->array() + E;
    matxxT value2 = n->getValue<Matrix<T, Dynamic, Dynamic>>()->array() - E;
    // check if n is loss node
    if (n == loss)
    {
        // get numerical gradient directly
        grad = (value1 - value2) / (2 * E);
    }
    else
    {
        // otherwise create variable nodes and run the session on nodes
        Variable<matxxT> newNodeP(std::move(value1));
        Variable<matxxT> newNodeN(std::move(value2));
        // swap the node with the new variable node
        swapNodes(&newNodeP, n);
        // compute value of loss
        _session.Run(loss, feed);
        matxxT outP = *loss->getValue<matxxT>();
        //std::cout << "Loss+:" << outP << std::endl;
        // swap the node with other new node
        swapNodes(&newNodeN, &newNodeP);
        _session.Run(loss, feed);
        matxxT outN = *loss->getValue<matxxT>();
        //std::cout << "Loss-:" << outN << std::endl;
        // swap the node back in and compute the graph
        swapNodes(n, &newNodeN);
        // find numerical gradient and check the node gradient
        grad = (outP - outN) / (2 * E);
    }
    // As E is added to each element of the matrix we consider the difference of the sum of matrices
    auto er = grad.sum() - n->getGradient<matxxT>().sum();
    // check if the error is within the threshold
    assert(er < err);
    std::cout << "Numerical gradient: " << grad.sum() << std::endl;
    std::cout << "Node gradient: " << n->getGradient<matxxT>().sum() << std::endl;
}

// swap nodes in and out of computational graph
void NN::swapNodes(BaseNode *in, BaseNode *out)
{
    std::vector<Locking_ptr<BaseNode>> consumers = out->getConsumers();
    // only end node (i.e. loss) has no consumers
    for (auto cns : consumers)
    {
        // remove other node from its consumers' input
        cns->eraseInput(out);
        // add this node to consumer's input
        cns->addInputs(in);
        // remove consumers of other node
        out->eraseConsumer(cns.get());
        // add consumers to new node
        in->addConsumers(cns.get());
    }
}

GradientDescentOptimizer::GradientDescentOptimizer(float lr) : learningRate_(lr) {}

void GradientDescentOptimizer::computeGradients(BaseNode *loss)
{
    Locking_ptr<BaseNode> lss{loss};
    // get node queue in level order traversal 
    std::vector<Locking_ptr<BaseNode>> nodes = getNodeQueue(loss);
    // store ftrs to wait for them later
    std::vector<std::future<void>> ftrs;
    // clear gradients from previous epoch
    std::for_each(nodes.begin(), nodes.end(), [](Locking_ptr<BaseNode> n){ n->clearGrads(); });    
    // calculate gradients
    for (auto &node : nodes)
    {
        // calculate node gradients 
        ftrs.emplace_back(std::async(std::launch::async, [node] { node->gradient(); }));
    }
    // wait for results
    for_each(ftrs.begin(), ftrs.end(), [](std::future<void> &ftr) { ftr.wait(); });
    //std::cout << "gradients computed! " << std::endl;
}

template <typename T>
Minimizer<T> GradientDescentOptimizer::minimize(BaseNode *loss)
{
    // Instantiate a minimzer object and return it
    return Minimizer<T>(this, loss);
}

std::vector<Locking_ptr<BaseNode>> GradientDescentOptimizer::getNodeQueue(BaseNode *loss)
{
    // Do level-order traversal 
    // create a deque 
    std::deque<Locking_ptr<BaseNode>> nodeQueue;
    // create a vector of nodes to return the nodes
    std::vector<Locking_ptr<BaseNode>> nodesList;
    // create a map for exitence ckeck in constant time
    std::unordered_map<BaseNode *, bool> visitedNodes;
    nodeQueue.push_front(Locking_ptr<BaseNode>(loss));
    while (!nodeQueue.empty())
    {
        // get the front element
        Locking_ptr<BaseNode> node = nodeQueue.front();
        // cash in node list
        nodesList.push_back(node);
        // set node to visited
        visitedNodes[node.get()] = true;
        // remove the visited node from queue
        nodeQueue.pop_front();
        // get the inputs
        auto nodes = node->getInputs();
        // go through all inputs of the node
        for (auto &n : nodes)
        {
            // check if the node is visited before
            if (visitedNodes[n.get()] != true)
            {
                // if node not visited add to queue
                nodeQueue.push_back(n);
            }
        }
    }
    // return the node list
    return nodesList;
}

template <template <typename> class U, typename T>
BaseNode* Graph::addNodeOne(std::unique_ptr<U<T>> n)
{
    // move the node to the list
    _baseNodes.push_back(std::move(n));
    // return the moved node
    return _baseNodes.back().get();
}

template <template <typename, typename, typename> class U, typename T, typename T1, typename T2>
BaseNode* Graph::addNodeTwo(std::unique_ptr<U<T, T1, T2>> n)
{
    // move the node to the list
    _baseNodes.push_back(std::move(n));
    // return the moved node
    return _baseNodes.back().get();
}

std::vector<std::unique_ptr<BaseNode>> &Graph::getNodes()
{
    return _baseNodes;
}

template <typename T>
Locking_ptr<T>::Locking_ptr(T *ptr) : ptr_(ptr), Mtx_(&(ptr->Mtx_)) {}

template <typename T>
Locking_ptr<T>::Locking_ptr(T *ptr, std::mutex *mtx) : ptr_(ptr), Mtx_(mtx) {}

template <typename T>
Locking_ptr<T>::~Locking_ptr() {}

template <typename T>
Locking_ptr<T>::Locking_ptr(Locking_ptr<T> const &other)
{
    //std::cout << "Locking_ptr copy constructor..." << std::endl;
    Mtx_.store(other.Mtx_.load());
    std::unique_lock<std::mutex> lc1(*(Mtx_.load()));
    ptr_ = other.ptr_;
}

template <typename T>
Locking_ptr<T>::Locking_ptr(Locking_ptr<T> &&other)
{
    //std::cout << "Locking_ptr move constructor..." << std::endl;
    Mtx_.store(other.Mtx_.load());
    std::unique_lock<std::mutex> lc1(*Mtx_);
    ptr_ = std::move(other.ptr_);
    other.Mtx_.store(nullptr);
}

template <typename T>
Locking_ptr<T> &Locking_ptr<T>::operator=(Locking_ptr<T> &&other)
{
    Mtx_.store(other.Mtx_.load());
    std::unique_lock<std::mutex> lc1(*Mtx_);
    //std::cout << "Locking_ptr move assignment constructor..." << std::endl;
    if (this != &other)
    {
        ptr_ = std::move(other.ptr_);
        other.Mtx_.store(nullptr);
    }
    return *this;
}

template <typename T>
Locking_ptr<T> &Locking_ptr<T>::operator=(Locking_ptr<T> const &other)
{
    Mtx_.store(other.Mtx_.load());
    std::unique_lock<std::mutex> lc1(*Mtx_);
    //std::cout << "Locking_ptr copy assignment constructor..." << std::endl;
    if (this != &other)
    {
        ptr_ = other.get();
    }
    return *this;
}

template <typename T>
T *Locking_ptr<T>::operator->() const
{
    std::unique_lock<std::mutex> lc(*Mtx_);
    return ptr_;
}

template <typename T>
T &Locking_ptr<T>::operator*() const
{
    std::unique_lock<std::mutex> lc(*Mtx_);
    return *ptr_;
}

template <typename T>
bool Locking_ptr<T>::operator==(Locking_ptr<T> const &rhs)
{
    return this->get() == rhs.get();
}

template <typename T>
void Locking_ptr<T>::reset(T *ptr)
{
    std::unique_lock<std::mutex> lc(*Mtx_);
    ptr_ = ptr;
}

template <typename T>
T *Locking_ptr<T>::get() const
{
    std::unique_lock<std::mutex> lc(*Mtx_);
    return ptr_;
}
//-------- shared locking ptr ---------

template <typename T>
Locking_shared_ptr<T>::Locking_shared_ptr(T *ptr) : ptr_(ptr), Mtx_(&(ptr->Mtx_)){}

template <typename T>
Locking_shared_ptr<T>::Locking_shared_ptr(T *ptr, std::mutex *mtx) : ptr_(ptr), Mtx_(mtx) {}

template <typename T>
Locking_shared_ptr<T>::Locking_shared_ptr(std::shared_ptr<T> ptr)
{
    Mtx_.store((&ptr->Mtx_).load());
    ptr_ = ptr;
}

template <typename T>
Locking_shared_ptr<T>::~Locking_shared_ptr() { }

template <typename T>
Locking_shared_ptr<T>::Locking_shared_ptr(Locking_shared_ptr<T> const &other)
{
    //std::cout << "Locking_shared_ptr copy constructor..." << std::endl;
    Mtx_.store(other.Mtx_.load());
    std::unique_lock<std::mutex> lc1(*Mtx_);
    ptr_ = other.ptr_;
}

template <typename T>
Locking_shared_ptr<T>::Locking_shared_ptr(Locking_shared_ptr<T> &&other)
{
    //std::cout << "Locking_shared_ptr move constructor..." << std::endl;
    Mtx_.store(other.Mtx_.load());
    std::unique_lock<std::mutex> lc1(*Mtx_);
    ptr_ = std::move(other.ptr_);
    other.Mtx_.store(nullptr);
}

template <typename T>
Locking_shared_ptr<T> &Locking_shared_ptr<T>::operator=(Locking_shared_ptr<T> &&other)
{
    Mtx_.store(other.Mtx_.load());
    std::unique_lock<std::mutex> lc1(*Mtx_);
    //std::cout << "Locking_shared_ptr move assignment constructor..." << std::endl;
    if (this != &other)
    {
        ptr_ = std::move(other.ptr_);
        other.Mtx_.store(nullptr);
    }
    return *this;
}

template <typename T>
Locking_shared_ptr<T> &Locking_shared_ptr<T>::operator=(Locking_shared_ptr<T> const &other)
{
    Mtx_.store(other.Mtx_.load());
    std::unique_lock<std::mutex> lc1(*Mtx_);
    //std::cout << "Locking_shared_ptr copy assignment constructor..." << std::endl;
    if (this != &other)
    {
        ptr_ = other.ptr_;
    }
    return *this;
}

template <typename T>
T *Locking_shared_ptr<T>::operator->() const
{
    std::unique_lock<std::mutex> lc(*Mtx_);
    return ptr_.get();
}

template <typename T>
T &Locking_shared_ptr<T>::operator*() const
{
    std::unique_lock<std::mutex> lc(*Mtx_);
    return *ptr_;
}

template <typename T>
bool Locking_shared_ptr<T>::operator==(Locking_shared_ptr<T> const &rhs)
{
    return this->get() == rhs.get();
}

template <typename T>
void Locking_shared_ptr<T>::reset(T *ptr)
{
    std::unique_lock<std::mutex> lc(*Mtx_);
    ptr_.reset(ptr);
}

template <typename T>
T *Locking_shared_ptr<T>::get() const
{
    std::unique_lock<std::mutex> lc(*Mtx_);
    return ptr_.get();
}

//-------- locking unique ptr ---------

template <typename T>
Locking_unique_ptr<T>::Locking_unique_ptr(T *ptr) : ptr_(ptr), Mtx_(&(ptr->Mtx_)){}

template <typename T>
Locking_unique_ptr<T>::Locking_unique_ptr(T *ptr, std::mutex *mtx) : ptr_(ptr), Mtx_(mtx) {}

template <typename T>
Locking_unique_ptr<T>::Locking_unique_ptr(std::unique_ptr<T> &&ptr)
{
    Mtx_.store((&ptr->Mtx_).load());
    std::unique_lock<std::mutex> lc1(*Mtx_);
    ptr_ = std::move(ptr);
}

template <typename T>
Locking_unique_ptr<T>::~Locking_unique_ptr() { }

template <typename T>
Locking_unique_ptr<T>::Locking_unique_ptr(Locking_unique_ptr<T> &&other)
{
    //std::cout << "Locking_shared_ptr move constructor..." << std::endl;
    Mtx_.store(other.Mtx_.load());
    std::unique_lock<std::mutex> lc1(*Mtx_);
    ptr_ = std::move(other.ptr_);
    other.Mtx_.store(nullptr);
}

template <typename T>
Locking_unique_ptr<T> &Locking_unique_ptr<T>::operator=(Locking_unique_ptr<T> &&other)
{
    Mtx_.store(other.Mtx_.load());
    std::unique_lock<std::mutex> lc1(*Mtx_);
    //std::cout << "Locking_shared_ptr move assignment constructor..." << std::endl;
    if (this != &other)
    {
        ptr_ = std::move(other.ptr_);
        other.Mtx_.store(nullptr);
    }
    return *this;
}

template <typename T>
T *Locking_unique_ptr<T>::operator->() const
{
    std::unique_lock<std::mutex> lc(*Mtx_);
    return ptr_.get();
}

template <typename T>
T &Locking_unique_ptr<T>::operator*() const
{
    std::unique_lock<std::mutex> lc(*Mtx_);
    return *ptr_;
}

template <typename T>
bool Locking_unique_ptr<T>::operator==(Locking_unique_ptr<T> const &rhs)
{
    return this->get() == rhs.get();
}

template <typename T>
void Locking_unique_ptr<T>::reset(T *ptr)
{
    std::unique_lock<std::mutex> lc(*Mtx_);
    ptr_.reset(ptr);
}

template <typename T>
T *Locking_unique_ptr<T>::get() const
{
    std::unique_lock<std::mutex> lc(*Mtx_);
    return ptr_.get();
}


// --- BaseNode ---

void BaseNode::setName(std::string n)
{
    std::lock_guard<std::mutex> lck(nodeMtx_);
    _name = n;
}

void BaseNode::addInputs(BaseNode *n)
{
    std::lock_guard<std::mutex> lck1(nodeMtx_);
    // Look if a node is previously replaced with nullptr
    for (int i = 0; i < _inputs.size(); i++)
    {
        if (_inputs[i].get() == nullptr)
        {
            // Replace BaseNode* with nullptr
            _inputs[i].reset(n);
            return;
        }
    }
    _inputs.push_back(Locking_ptr<BaseNode>(n));
}

void BaseNode::eraseInput(BaseNode *n)
{
    std::lock_guard<std::mutex> lck(nodeMtx_);
    // remove the input node but keep the place
    for (int i = 0; i < _inputs.size(); i++)
    {
        if (_inputs[i] == n)
        {
            // use nullptr as a placeholder
            _inputs[i].reset(nullptr);
        }
    }
}

void BaseNode::addConsumers(BaseNode *n)
{
    std::lock_guard<std::mutex> lck(nodeMtx_);
    // remove consumer but keep the place
    for (int i = 0; i < _consumers.size(); i++)
    {
        // check if there is deleted consumer
        if (_consumers[i].get() == nullptr)
        {
            // replace the deleted consumer in place
            _consumers[i].reset(n);
            return;
        }
    }
    // add consumer and increment the size
    _consumers.push_back(Locking_ptr<BaseNode>(n));
    // Increment consumer size
    consumerSize_++;
}

void BaseNode::eraseConsumer(BaseNode *n)
{
    std::lock_guard<std::mutex> lck(nodeMtx_);
    // remove consumer but keep the place
    for (int i = 0; i < _consumers.size(); i++)
    {
        // find the node
        if (_consumers[i] == n)
        {
            // use nullptr as a placeholder
            _consumers[i].reset(nullptr);
        }
    }
}

template <typename T>
Locking_shared_ptr<T> BaseNode::getValue()
{
    std::unique_lock<std::mutex> lck(nodeMtx_);
    auto node = static_cast<Node<T> *>(this);
    lck.unlock();
    return node->getValue();
}

template <typename T>
T BaseNode::getGradient()
{
    std::unique_lock<std::mutex> lck(nodeMtx_);
    auto node = static_cast<Node<T> *>(this);
    lck.unlock();
    return node->getGradient();
}

template <typename T>
void BaseNode::setGrad(T t)
{
    std::unique_lock<std::mutex> lck(nodeMtx_);
    auto node = static_cast<Node<T> *>(this);
    lck.unlock();
    node->setGrad(t);
}

std::string BaseNode::getName()
{
    std::lock_guard<std::mutex> lck(nodeMtx_);
    return _name;
}

std::vector<Locking_ptr<BaseNode>> &BaseNode::getInputs()
{
    std::lock_guard<std::mutex> lck(nodeMtx_);
    return _inputs;
}

std::vector<Locking_ptr<BaseNode>> &BaseNode::getConsumers()
{
    std::lock_guard<std::mutex> lck(nodeMtx_);
    return _consumers;
}

nodeType BaseNode::getNodeType()
{
    std::lock_guard<std::mutex> lck(nodeMtx_);
    return _nType;
}

operationType BaseNode::getOperationType()
{
    std::lock_guard<std::mutex> lck(nodeMtx_);
    return _opType;
}

// --- Node  ---

template <typename T>
Locking_shared_ptr<T> Node<T>::getValue()
{
    std::lock_guard<std::mutex> lck(nodeMtx_);
    //std::cout << "Variable get value..." << std::endl;
    assert (_dataAvailable.load());
    //std::cout << "Output get: " << *_output << ", size: " << (*_output).rows() << "," << (*_output).cols() << std::endl;
    return _output;

}

template <typename T>
T Node<T>::getGradient()
{
    std::unique_lock<std::mutex> lck1(nodeMtx_);
    //std::cout << "Get gradient ...\n";
    //std::cout << "Thread ID: " << std::this_thread::get_id() << std::endl;
    // Initialize node's gradient
    T grad;
    // check if node has a consumer; consumerSize_ is atomic
    int cnsSize = this->consumerSize_.load();
    if (cnsSize > 0)
    {
        // check if gradient data is available
        if ((this->_gradientAvailable.load()))
        {
            //std::cout << "No wait required; Thread ID: " << std::this_thread::get_id() << std::endl;
            grad.setZero(_grad[0]->rows(), _grad[0]->cols());
            //  get total derivative
            for (auto &&g : _grad)
            {
                // sum gradients to get total
                grad += *g;
            }
            //std::cout << "Total gradient get: " << grad << ", size: " << grad.rows() << "," << grad.cols() << std::endl;
            return grad;
        }
        else
        {
            //  wait until gradient is available
            cond_.wait(lck1, [this]() { return this->_gradientAvailable.load(); });
            //std::cout << "Notified Thread ID: " << std::this_thread::get_id() << std::endl;
            grad.setZero(_grad[0]->rows(), _grad[0]->cols());
            //  get total derivative
            for (auto &&g : _grad)
            {
                // sum gradients to get total
                grad += *g;
            }
            //std::cout << "Total gradient get: " << grad << ", size: " << grad.rows() << "," << grad.cols() << std::endl;
            return grad;
        }
    }
    else
    {
        //return 1s  of size output if there is no consumer
        //std::cout << "No consumer" << std::endl;
        grad.setOnes(_output->rows(), _output->cols());
        return grad;
    }
}

template <typename T>
void Node<T>::setValue(T &&t)
{
    std::lock_guard<std::mutex> lck(nodeMtx_);
    _dataAvailable.store(true);
    _output.reset(new T(t));
    //std::cout << "Output set: " << *_output << ", size: " << (*_output).rows() << "," << (*_output).cols() << std::endl;
}

template <typename T>
void Node<T>::setGrad(T t)
{
    std::lock_guard<std::mutex> lck1(nodeMtx_);
    //std::cout << "Gradient set: " << t << ", size: " << t.rows() << "," << t.cols() << std::endl;
    // gradient and value must have same dimensions
    assert(t.cols() == _output->cols() or t.rows() == _output->rows());
    // add gradient
    _grad.push_back(Locking_shared_ptr<T>((new T(t)), &(this->Mtx_)));
    // get the number of consumer of the node; consumerSize_ is atomic
    int cnsSize = this->consumerSize_.load();
    // check if gradient of all consumers are set; use >= as a node might not have a consumer
    if (_grad.size() >= cnsSize)
    {
        // set flag to true
        _gradientAvailable.store(true);
        // notify all threads waiting for this data
        cond_.notify_all();
    }
}

template <typename T>
void Node<T>::clearGrads()
{
    // lock for all _grad read and write
    std::lock_guard<std::mutex> lck1(nodeMtx_);
    _gradientAvailable.store(false);
    // reset gradients
    _grad.clear();
}

// --- Variable ---

template <typename T>
Variable<T>::Variable(T &&a)
{
    std::cout << "Variable contructor ..." << std::endl;
    this->_nType = nodeType::variable;
    this->_opType = operationType::NA;
    // set value locks the node, no need to create a lock
    this->setValue(std::move(a));
}

template <typename T>
Variable<T>::Variable(Variable<T> &&other)
{
    std::cout << "Variable move contructor ..." << std::endl;
    // move
    this->_nType = nodeType::variable;
    this->_opType = operationType::NA;
    // set value and get value locks the node, no need to create a lock
    T val = (&other)->getValue();
    this->setValue(std::move(*val));
}

template <typename T>
Variable<T> &Variable<T>::operator=(Variable<T> &&other)
{
    std::cout << "Variable move assignment contructor ..." << std::endl;
    if (this != &other)
    {
        // move
        this->_nType = nodeType::variable;
        this->_opType = operationType::NA;
        // set value and get value locks the node, no need to create a lock
        T val = (&other)->getValue();
        this->setValue(*val);
    }
    return *this;
}

template <typename T>
void Variable<T>::compute() { return; }

template <typename T>
void Variable<T>::gradient()
{
    //std::cout << "Variable gradient ..." << std::endl;
}

template <typename T>
void Variable<T>::updateValue(float lr)
{
    //std::cout << "Variable update value ..." << std::endl;
    // variable has only one input gradient
    T grad = this->getGradient();
    Locking_shared_ptr<T> output = this->getValue();
    // update variable values based on learning rate and gradient
    this->setValue(output->array() - (grad.array() * lr));
}

// --- Placeholder ---

template <typename T>
Placeholder<T>::Placeholder(std::string n)
{
    this->_nType = nodeType::placeholder;
    this->setName(n);
}

template <typename T>
void Placeholder<T>::compute(){ return; }

template <typename T>
void Placeholder<T>::gradient(){ return; }


// --- Operation ---

template <typename T>
Operation<T>::Operation()
{
    this->_nType = nodeType::operation;
}

// --- UnaryOperation ---

template <typename T>
UnaryOperation<T>::UnaryOperation(BaseNode *rhs)
{
    // use Locking_ptr<BaseNode> to cast to BaseNode
    Locking_ptr<BaseNode> rhsptr(rhs);
    Locking_ptr<BaseNode> ptrthis(this);
    // add inpupts
    ptrthis->addInputs(rhs);
    // add consumers
    rhsptr->addConsumers(this);
}

// --- BinaryOperation ---

template <typename T>
BinaryOperation<T>::BinaryOperation(BaseNode *lhs, BaseNode *rhs)
{
    // use Locking_ptr<BaseNode> to cast to BaseNode
    Locking_ptr<BaseNode> rhsptr(rhs);
    Locking_ptr<BaseNode> lhsptr(lhs);
    Locking_ptr<BaseNode> ptrthis(this);
    // add inputs
    ptrthis->addInputs(lhs);
    ptrthis->addInputs(rhs);
    // add consumers
    lhsptr->addConsumers(this);
    rhsptr->addConsumers(this);
}

// --- add operation ---

template <typename T, typename T1, typename T2>
Add<T, T1, T2>::Add(BaseNode &a, BaseNode &b) : BinaryOperation<T>(&a, &b)
{
    this->_opType = operationType::addition;
}

template <typename T, typename T1, typename T2>
Add<T, T1, T2>::Add(BaseNode *a, BaseNode *b) : BinaryOperation<T>(a, b)
{
    this->_opType = operationType::addition;
}

template <typename T, typename T1, typename T2>
void Add<T, T1, T2>::compute()
{
    //std::cout << "Compute add operation ..." << std::endl;
    std::vector<Locking_ptr<BaseNode>> inputs = this->getInputs();
    // get value of input nodes
    Locking_shared_ptr<T1> A = inputs[0]->getValue<T1>();
    Locking_shared_ptr<T2> B = inputs[1]->getValue<T2>();
    // broadcast column or row vectors if both inputs have same number of columns
    if (A->rows() != B->rows() & A->cols() == B->cols())
    {
        if (B->rows() == 1)
        {
            // add B to each row of A
            this->setValue(A->rowwise() + B->row(0));
        }
        else if (A->rows() == 1)
        {
            // Add A to each row of B
            this->setValue(B->rowwise() + A->row(0));
        }
    }
    // broadcast column or row vectors if both inputs have same number of rows
    else if (A->cols() != B->cols() & A->rows() == B->rows())
    {
        if (B->cols() == 1)
        {
            // add B to each column of A
            this->setValue(A->colwise() + B->col(0));
        }
        else if (A->cols() == 1)
        {
            // add A to each colun of B
            this->setValue(B->colwise() + A->col(0));
        }
    }
    // if A is scalar
    else if (A->cols() == 1 & A->rows() == 1)
    {
        this->setValue((*A)(0) + B->array());
    }
    // if B is scalar
    else if (B->cols() == 1 & B->rows() == 1)
    {
        this->setValue((*B)(0) + A->array());
    }
    else
    {
        // they are same size so element-wise addition without broadcasting
        this->setValue(A->array() + B->array());
    }
}

template <typename T, typename T1, typename T2>
void Add<T, T1, T2>::gradient()
{
    //std::cout << "Compute Add operation geradient ..." << std::endl;
    T grad = this->getGradient();
    std::vector<Locking_ptr<BaseNode>> inputs = this->getInputs();
    Locking_shared_ptr<T1> A = inputs[0]->getValue<T1>();
    Locking_shared_ptr<T2> B = inputs[1]->getValue<T2>();
    /* Check for broadcasting. If Gradient is larger than A, then A was broadcasted.
     Broadcasted variable is as though it has that many consumers. So the gradient is 
     the total gradient  (the sum of gradients in the broadcasted direction) */
    if (grad.cols() > A->cols() or grad.rows() > A->rows())
    {
        T g;
        g.setOnes(A->rows(), A->cols());
        if (A->rows() == 1 & A->cols() == 1)
        {
            // if A is scalar
            auto a = grad.sum();
            T gr(1, 1);
            gr << a;
            inputs[0]->setGrad<T>(gr);
        }
        else if (A->rows() == 1)
        {
            // broadcasted in columns direction
            T gr = g * grad;
            inputs[0]->setGrad<T>(gr);
        }
        else if (A->cols() == 1)
        {
            // broadcasted in rows direction
            T gr = grad * g;
            inputs[0]->setGrad<T>(gr);
        }
    }
    else
    {
        // No broadcasting; Input gradient is the same as output gradient
        inputs[0]->setGrad<T>(grad);
    }
    // If Gradient is larger than B, then B was broadcasted
    if (grad.cols() > B->cols() or grad.rows() > B->rows())
    {
        T g;
        g.setOnes(B->rows(), B->cols());
        if (B->rows() == 1 & B->cols() == 1)
        {
            // if B is scalar
            auto a = grad.sum();
            T gr(1, 1);
            gr << a;
            inputs[1]->setGrad<T>(gr);
        }
        else if (B->rows() == 1)
        {
            // broadcasted in columns direction
            T gr = g * grad;
            inputs[1]->setGrad<T>(gr);
        }
        else if (B->cols() == 1)
        {
            // broadcasted in rows direction
            T gr = grad * g;
            inputs[1]->setGrad<T>(gr);
        }
    }
    else
    {
        // No broadcasting; Input gradient is the same as output gradient
        inputs[1]->setGrad<T>(grad);
    }
}

// --- negative operation---

template <typename T>
Negative<T>::Negative(BaseNode &a) : UnaryOperation<T>(&a)
{
    this->_opType = operationType::negative;
}

template <typename T>
Negative<T>::Negative(BaseNode *a) : UnaryOperation<T>(a)
{
    this->_opType = operationType::negative;
}

template <typename T>
void Negative<T>::compute()
{
    //std::cout << "Compute negative operation ..." << std::endl;
    Locking_ptr<BaseNode> ptrthis(this);
    this->setValue(-(*(ptrthis->getInputs()[0]->getValue<T>())));
}

template <typename T>
void Negative<T>::gradient()
{
    // get inputs of this node; it only has one input
    std::vector<Locking_ptr<BaseNode>> inputs = this->getInputs();
    //std::cout << "Compute negative operation geradient ..." << std::endl;
    inputs[0]->setGrad<T>(-(this->getGradient()));
}

// --- Multiply Operation ---

template <typename T, typename T1, typename T2>
Multiply<T, T1, T2>::Multiply(BaseNode &a, BaseNode &b) : BinaryOperation<T>(&a, &b)
{
    this->_opType = operationType::multiply;
}

template <typename T, typename T1, typename T2>
Multiply<T, T1, T2>::Multiply(BaseNode *a, BaseNode *b) : BinaryOperation<T>(a, b)
{
    this->_opType = operationType::multiply;
}

template <typename T, typename T1, typename T2>
void Multiply<T, T1, T2>::compute()
{
    //std::cout << "Compute multiplication operation..." << std::endl;
    std::vector<Locking_ptr<BaseNode>> inputs = this->getInputs();
    // multiplication of scalar and matrix
    Locking_shared_ptr<T1> A = inputs[0]->getValue<T1>();
    Locking_shared_ptr<T2> B = inputs[1]->getValue<T2>();
    // perform matrix multiplication
    this->setValue(A->array() * B->array());
}

template <typename T, typename T1, typename T2>
void Multiply<T, T1, T2>::gradient()
{
    //std::cout << "Compute multiplication operation gradient..." << std::endl;
    // get output gradient from consumer
    T G = this->getGradient();
    // get inputs of this node
    std::vector<Locking_ptr<BaseNode>> inputs = this->getInputs();
    Locking_shared_ptr<T1> A = inputs[0]->getValue<T1>();
    Locking_shared_ptr<T2> B = inputs[1]->getValue<T2>();
    // calculate and set gradient for first input "A"
    inputs[0]->setGrad<T>(G.array() * B->array());
    // calculate and set gradient for first input "B"
    inputs[1]->setGrad<T>(G.array() * A->array());
}

// --- MatMultiply Operation ---

template <typename T, typename T1, typename T2>
MatMultiply<T, T1, T2>::MatMultiply(BaseNode &a, BaseNode &b) : BinaryOperation<T>(&a, &b)
{
    this->_opType = operationType::multiply;
}

template <typename T, typename T1, typename T2>
MatMultiply<T, T1, T2>::MatMultiply(BaseNode *a, BaseNode *b) : BinaryOperation<T>(a, b)
{
    this->_opType = operationType::multiply;
}

template <typename T, typename T1, typename T2>
void MatMultiply<T, T1, T2>::compute()
{
    //std::cout << "Compute matrix multiplication operation..." << std::endl;
    std::vector<Locking_ptr<BaseNode>> inputs = this->getInputs();
    // get input node values
    Locking_shared_ptr<T1> A = inputs[0]->getValue<T1>();
    Locking_shared_ptr<T2> B = inputs[1]->getValue<T2>();
    // perform matrix multiplication
    this->setValue((*A) * (*B));
}

template <typename T, typename T1, typename T2>
void MatMultiply<T, T1, T2>::gradient()
{
    //std::cout << "Compute matrix multiplication operation gradient..." << std::endl;
    // get output gradient from consumer
    T G = this->getGradient();
    // get inputs of this node
    std::vector<Locking_ptr<BaseNode>> inputs = this->getInputs();
    Locking_shared_ptr<T1> A = inputs[0]->getValue<T1>();
    Locking_shared_ptr<T2> B = inputs[1]->getValue<T2>();
    // calculate and set gradient for first input "A"
    T C = G * B->transpose();
    inputs[0]->setGrad<T>(C);
    // calculate and set gradient for second input "B"
    T D = A->transpose() * G;
    inputs[1]->setGrad<T>(D);
}

// --- DotProduct ---

template <typename T, typename T1, typename T2>
Dot<T, T1, T2>::Dot(BaseNode &a, BaseNode &b) : BinaryOperation<T>(&a, &b)
{
    this->_opType = operationType::dotproduct;
}

template <typename T, typename T1, typename T2>
Dot<T, T1, T2>::Dot(BaseNode *a, BaseNode *b) : BinaryOperation<T>(a, b)
{
    this->_opType = operationType::dotproduct;
}

template <typename T, typename T1, typename T2>
void Dot<T, T1, T2>::compute()
{
    //std::cout << "Compute dot product operation ..." << std::endl;
    Locking_ptr<BaseNode> ptrthis(this);
    this->setValue(ptrthis->getInputs()[0]->getValue<T1>().dot(ptrthis->getInputs()[1]->getValue<T2>()));
}

template <typename T, typename T1, typename T2>
void Dot<T, T1, T2>::gradient()
{
    //std::cout << "Compute dot product operation gradient..." << std::endl;
    // get output gradient from consumer
    T G = this->getGradient();
    // get inputs of this node
    std::vector<Locking_ptr<BaseNode>> inputs = this->getInputs();
    Locking_shared_ptr<T1> A = inputs[0]->getValue<T1>();
    Locking_shared_ptr<T2> B = inputs[1]->getValue<T2>();
    // calculate and set gradient for first input "A"
    T C = G * B->transpose();
    inputs[0]->setGrad<T>(C);
    // calculate and set gradient for first input "B"
    T D = A->transpose() * G;
    inputs[1]->setGrad<T>(D);
}

// --- Sigmoid ---

template <typename T>
Sigmoid<T>::Sigmoid(BaseNode &a) : UnaryOperation<T>(&a)
{
    this->_opType = operationType::sigmoid;
}

template <typename T>
Sigmoid<T>::Sigmoid(BaseNode *a) : UnaryOperation<T>(a)
{
    this->_opType = operationType::sigmoid;
}

template <typename T>
void Sigmoid<T>::compute()
{
    Locking_ptr<BaseNode> ptrthis(this);
    //std::cout << "Compute sigmoid operation ..." << std::endl;
    this->setValue(1 / (1 + exp(-(ptrthis->getInputs()[0]->getValue<T>()->array()))));
}

template <typename T>
void Sigmoid<T>::gradient()
{
    //std::cout << "Compute sigmoid gradient..." << std::endl;
    // Cast this to BaseNode and wrap around with a locking class
    Locking_ptr<BaseNode> ptrthis(this);
    // get inputs of this node
    std::vector<Locking_ptr<BaseNode>> inputs = this->getInputs();
    // get output gradient from consumer
    T G = this->getGradient();
    // get sigmoid value
    Locking_shared_ptr<T> sig = ptrthis->getValue<T>();
    // compute gradient
    T grad = G.array() * sig->array() * (1 - sig->array());
    // set gradient for input
    inputs[0]->setGrad<T>(grad);
}

// --- log ---

template <typename T>
Log<T>::Log(BaseNode &a) : UnaryOperation<T>(&a)
{
    this->_opType = operationType::log;
}

template <typename T>
Log<T>::Log(BaseNode *a) : UnaryOperation<T>(a)
{
    this->_opType = operationType::log;
}

template <typename T>
void Log<T>::compute()
{
    // Cast to this to BaseNode and wrape around with a lock
    Locking_ptr<BaseNode> ptrthis(this);
    //std::cout << "Compute log operation ..." << std::endl;
    this->setValue(log(ptrthis->getInputs()[0]->getValue<T>()->array()));
}

template <typename T>
void Log<T>::gradient()
{
    //std::cout << "Compute log gradient..." << std::endl;
    // get output gradient from consumer
    T G = this->getGradient();
    // copy inputs of this node to local variable to avoid data
    std::vector<Locking_ptr<BaseNode>> inputs = this->getInputs();
    // get log input value
    Locking_shared_ptr<T> log = inputs[0]->getValue<T>();
    // compute gradient; elementwise division
    T grad = G.array() / log->array();
    // set gradient for input
    inputs[0]->setGrad<T>(grad);
}

// --- Sum ---

template <typename T>
Sum<T>::Sum(BaseNode &a, int axis) : UnaryOperation<T>(&a), _axis(axis)
{
    this->_opType = operationType::sum;
}

template <typename T>
Sum<T>::Sum(BaseNode *a, int axis) : UnaryOperation<T>(a), _axis(axis)
{
    this->_opType = operationType::sum;
}

template <typename T>
void Sum<T>::compute()
{
    Locking_ptr<BaseNode> ptrthis(this);
    //std::cout << "Compute Sum operation ..." << std::endl;
    if (_axis == 0)
    {
        // if axis = 0 then sum colwise
        this->setValue(ptrthis->getInputs()[0]->getValue<T>().colwise().sum());
    }
    else if (_axis == 1)
    {
        // if axis = 1 then sum rowwise
        this->setValue(ptrthis->getInputs()[0]->getValue<T>().rowwise().sum());
    }
}

template <typename T>
void Sum<T>::gradient()
{
    //std::cout << "Compute sum operation gradient..." << std::endl;
    // get output gradient from consumer
    T G = this->getGradient();
    T g;
    // get inputs of this node
    std::vector<Locking_ptr<BaseNode>> inputs = this->getInputs();
    Locking_shared_ptr<T> A = inputs[0]->getValue<T>();
    if (G.rows() == 1)
    {
        g = G.replicate(A->rows(), 1);
        inputs[0]->setGrad<T>(g);
    }
    else if (G.cols() == 1)
    {
        g = G.replicate(1, A->cols());
        inputs[0]->setGrad<T>(g);
    }
    else
    {
        inputs[0]->setGrad<T>(G);
    }
}

/// --- Minimizaer Operation ----

template <typename T>
Minimizer<T>::Minimizer(GradientDescentOptimizer *grd, BaseNode *loss)
{
    grdOpt_ = Locking_ptr<GradientDescentOptimizer>(grd, &gMtx_);
    loss_ = Locking_ptr<BaseNode>(loss);
    learningRate_ = grd->learningRate_;
}

template <typename T>
Minimizer<T>::Minimizer(Minimizer<T> &&other)
{
    //std::cout << " Minimizer move contructor..." << std::endl;
    // move members
    grdOpt_ = std::move(other.grdOpt_);
    loss_ = std::move(other.loss_);
    learningRate_ = std::move(other.learningRate_);
}

template <typename T>
Minimizer<T> &Minimizer<T>::operator=(Minimizer<T> &&other)
{
    //std::cout << " Minimizer move assignment contructor..." << std::endl;
    if (this != &other)
    {
        // move members
        grdOpt_ = std::move(other.grdOpt_);
        loss_ = std::move(other.loss_);
        learningRate_ = std::move(other.learningRate_);
    }
    return *this;
}

// Compute updates the variable gradients based on learning rate
template <typename T>
void Minimizer<T>::compute()
{
    //std::cout << "Compute Minimization operation ..." << std::endl;
    // compute grdients 
    grdOpt_->computeGradients(loss_.get());
    // iterate through nodes and update variable values
    auto list = grdOpt_->getNodeQueue(loss_.get());
    for (auto &n : list)
    {
        if (n->getNodeType() == nodeType::variable)
        {
            auto v = static_cast<Variable<T> *>(n.get());
            v->updateValue(learningRate_);
        }
    }
}

template <typename T>
void Minimizer<T>::gradient() { return; }


template <typename T>
void Session::Run(BaseNode *n, std::unordered_map<std::string, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> *> feed)
{
    Locking_ptr<BaseNode> nptr(n);
    // empty node list
    _nodesList.clear();
    // obtain inputs for node n in post-order, to resolve inputs befor computation of an operation
    updateNodesList(nptr.get());

    for (auto m : _nodesList)
    {
        // if it's a placeholder then feed the data
        if (m->getNodeType() == nodeType::placeholder)
        {
            // set the output value
            Locking_ptr<Placeholder<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>> plcptr(static_cast<Placeholder<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> *>(m.get()));
            plcptr->setValue(std::move(*feed[plcptr->getName()]));
        } // if it's a operation then compute the value
        else if (m->getNodeType() == nodeType::operation)
        {
            // compute add and set the output value
            m->compute();
        }
    }
}

// get post order list of nodes
void Session::updateNodesList(BaseNode *n)
{
    Locking_ptr<BaseNode> nptr(n);
    // only operations have input nodes
    if (nptr->getNodeType() == nodeType::operation)
    {
        // go through input nodes recoursively
        auto list = nptr->getInputs();
        for (auto &m : list)
        {
            updateNodesList(m.get());
        }
    }
    // add node to the list
    _nodesList.push_back(Locking_ptr<BaseNode>(nptr.get()));
}

// Return  nodes list
std::vector<Locking_ptr<BaseNode>> Session::getNodesList()
{
    return _nodesList;
}