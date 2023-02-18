#pragma once

#include <memory>
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <Eigen/Dense>

#define BARWIDTH 30


namespace ts {

    // Enum for channel splitting directions in CNN
    // (declared outside for now because scoped enum declarationb seems
    // impossible)
    enum class ChannelSplit : int {
        NOSPLIT,
        SPLIT_HOR,	// Splits lines
        SPLIT_VERT	// Splits columns
    };

	std::vector<std::string> split(std::string str, char delimeter);

	std::string serializeUnsignedVec2D(
		std::vector<std::vector<unsigned>> &vec2d
	);

	std::vector<std::vector<unsigned>> parseUnsignedVec2D(
		std::ifstream &in
	);

	void progressBar(unsigned current, unsigned max);

	template <typename T> class Node;
	template <typename T> class InputNode;
	template <typename T> class ElementWiseNode;
	template <typename T> class MatProdNode;
	template <typename T> class ScalarNode;

	template <typename T> class WengertList;
	template <typename T> class Tensor;
	template <typename T> class Gradient;


	// This helper function allows us to create Tensor instances without
	// template syntax. This way, the type will be the same as its parent
	// WengertList.

	template <typename T>
	ts::Tensor<T> NewTensor(
		Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> newValue,
		ts::WengertList<T> * newWList
	);

	template <typename T>
	ts::Tensor<T> operator+(const ts::Tensor<T> &x, const ts::Tensor<T> &y);
	template <typename T>
	ts::Tensor<T> operator-(const ts::Tensor<T> &x, const ts::Tensor<T> &y);
	template <typename T>
	ts::Tensor<T> operator*(const ts::Tensor<T> &x, const ts::Tensor<T> &y);
	template <typename T>
	ts::Tensor<T> operator/(const ts::Tensor<T> &x, const ts::Tensor<T> &y);

	template <typename T>
	ts::Tensor<T> matProd(const ts::Tensor<T> &x, const ts::Tensor<T> &y);
	template <typename T>
	ts::Tensor<T> sigmoid(const ts::Tensor<T> &x);
	template <typename T>
	ts::Tensor<T> relu(const ts::Tensor<T> &x);
	template <typename T>
	ts::Tensor<T> leakyRelu(const ts::Tensor<T> &x);
	template <typename T>
	ts::Tensor<T> rescale(const ts::Tensor<T> &x);
	template <typename T>
	ts::Tensor<T> squaredNorm(const ts::Tensor<T> &x);


	// Forward declaration of friends
	// (grad accumulators and other autodiff operations)
	template <typename T> class GaElement;
	template <typename T> class GradientAccumulator;
	template <typename T> class AdamOptimizer;

	enum class ChannelSplit : int;

	template <typename T>
	ts::Tensor<T> convolution(const ts::Tensor<T> &mat, const ts::Tensor<T> &ker);

	template <typename T>
	ts::Tensor<T> maxPooling(const ts::Tensor<T> &x, std::vector<unsigned> pool);

	template <typename T>
	std::vector<ts::Tensor<T>> split(
		const ts::Tensor<T> &x,
		ChannelSplit channelSplit,
		unsigned nInputChannels
	);

	template <typename T>
	ts::Tensor<T> vertCat(const std::vector<ts::Tensor<T>> &x);

	template <typename T>
	ts::Tensor<T> flattening(const ts::Tensor<T> &x);

	template <typename T>
	ts::Tensor<T> im2col(
		const std::vector<ts::Tensor<T>> &x,
		std::vector<unsigned> kernelDim
	);

	template <typename T>
	std::vector<ts::Tensor<T>> col2im(
		const ts::Tensor<T> &x,
		std::vector<unsigned> outputDim
	);


	// ts::Node
    template <typename T>
    class ts::Node {
    protected:

        Node() {}

        // Represents an input variable
        Node(std::vector<long> shape);

        // Represents a unary operator
        Node(std::vector<long> shape,
            Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> xVal, int xDep
        );

        // Represents a binary operator
        Node(
            std::vector<long> shape,
            Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> xVal, int xDep,
            Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> yVal, int yDep
        );


        std::vector<int> dependencies{};

        virtual Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> incrementGradient(
                Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> &childDerivative,
                unsigned &j
        ) = 0;

        std::vector< Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> > values{};

        // Shape of the corresponding tensor
        long rows, cols;

    public:

        friend ts::Tensor<T>;
        friend ts::WengertList<T>;
        friend ts::GradientAccumulator<T>;
        friend ts::AdamOptimizer<T>;	// Needed to initialize moment estimates

        friend ts::Tensor<T> operator+<>(const ts::Tensor<T> &x, const ts::Tensor<T> &y);
        friend ts::Tensor<T> operator-<>(const ts::Tensor<T> &x, const ts::Tensor<T> &y);
        friend ts::Tensor<T> operator*<>(const ts::Tensor<T> &x, const ts::Tensor<T> &y);
        friend ts::Tensor<T> operator/<>(const ts::Tensor<T> &x, const ts::Tensor<T> &y);

        friend ts::Tensor<T> matProd<>(const ts::Tensor<T> &x, const ts::Tensor<T> &y);
        friend ts::Tensor<T> sigmoid<>(const ts::Tensor<T> &x);
        friend ts::Tensor<T> relu<>(const ts::Tensor<T> &x);
        friend ts::Tensor<T> leakyRelu<>(const ts::Tensor<T> &x);
        friend ts::Tensor<T> rescale<>(const ts::Tensor<T> &x);
        friend ts::Tensor<T> squaredNorm<>(const ts::Tensor<T> &x);

        friend ts::Tensor<T> convolution<>(const ts::Tensor<T> &mat, const ts::Tensor<T> &ker);
        friend ts::Tensor<T> maxPooling<>(const ts::Tensor<T> &x, std::vector<unsigned> pool);
        friend std::vector<ts::Tensor<T>> split<>(
            const ts::Tensor<T> &x,
            ChannelSplit channelSplit,
            unsigned nInputChannels
        );
        friend ts::Tensor<T> vertCat<>(const std::vector<ts::Tensor<T>> &x);
        friend ts::Tensor<T> flattening<>(const ts::Tensor<T> &x);
        friend ts::Tensor<T> im2col<>(
            const std::vector<ts::Tensor<T>> &x,
            std::vector<unsigned> kernelDim
        );
        friend std::vector<ts::Tensor<T>> col2im<>(
            const ts::Tensor<T> &x,
            std::vector<unsigned> outputDim
        );

    };



    template <typename T>
    class ts::InputNode : public ts::Node<T> {
    private:
        using ts::Node<T>::Node;
        InputNode(std::vector<long> shape, bool model);

        Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> incrementGradient(
                Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> &childDerivative,
                unsigned &j
        );

        // We will need this to optimize the tensor value in a ts::Model
        ts::Tensor<T> * optimizedTensor = NULL;

        // If true, node won't be removed on wList reset
        bool isModel = false;

    public:

        friend ts::WengertList<T>;
        friend ts::Tensor<T>;
        friend ts::GradientAccumulator<T>;
    };



    template <typename T>
    class ts::ElementWiseNode : public ts::Node<T> {
    private:
        using ts::Node<T>::Node;

        Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> incrementGradient(
                Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> &childDerivative,
                unsigned &j
        );
    };



    template <typename T>
    class ts::MatProdNode : public ts::Node<T> {
    private:
        using ts::Node<T>::Node;

        MatProdNode(
            std::vector<long> shape,
            Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> xVal, int xDep,
            Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> yVal, int yDep,
            std::vector<long int> newXSize, std::vector<long int> newYSize
        );

        // Size of the operands to figure out how to increment their partial
        // derivatives
        std::vector<long int> xSize;
        std::vector<long int> ySize;

        Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> incrementGradient(
                Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> &childDerivative,
                unsigned &j
        );

        friend ts::Tensor<T> matProd<>(const ts::Tensor<T> &x, const ts::Tensor<T> &y);
    };



    template <typename T>
    class ts::ScalarNode : public ts::Node<T> {
    private:
        using ts::Node<T>::Node;

        Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> incrementGradient(
                Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> &childDerivative,
                unsigned &j
        );
    };



        // ts::WengertList

    template <typename T>
    class ts::WengertList {
    private:
        bool elementWiseOnly = true;
        std::vector< std::shared_ptr<ts::Node<T>> > nodes{};

    public:
        int size();
        int reset();

        // Make a tensor optimizable
        void toggleOptimize(ts::Tensor<T> * tensor, bool enable);

        friend class ts::Tensor<T>;
        friend class ts::GradientAccumulator<T>;
        friend class ts::AdamOptimizer<T>;	// Needed to initialize moment estimates

        // Other non-element wise operations (to change elementWiseOnly)
        friend ts::Tensor<T> matProd<>(const ts::Tensor<T> &x, const ts::Tensor<T> &y);
        friend ts::Tensor<T> sigmoid<>(const ts::Tensor<T> &x);
        friend ts::Tensor<T> relu<>(const ts::Tensor<T> &x);
        friend ts::Tensor<T> leakyRelu<>(const ts::Tensor<T> &x);
        friend ts::Tensor<T> rescale<>(const ts::Tensor<T> &x);
        friend ts::Tensor<T> squaredNorm<>(const ts::Tensor<T> &x);

        friend ts::Tensor<T> convolution<>(const ts::Tensor<T> &mat, const ts::Tensor<T> &ker);
        friend ts::Tensor<T> maxPooling<>(const ts::Tensor<T> &x, std::vector<unsigned> pool);
        friend std::vector<ts::Tensor<T>> split<>(
            const ts::Tensor<T> &x,
            ChannelSplit channelSplit,
            unsigned nInputChannels
        );
        friend ts::Tensor<T> vertCat<>(const std::vector<ts::Tensor<T>> &x);
        friend ts::Tensor<T> flattening<>(const ts::Tensor<T> &x);
        friend ts::Tensor<T> im2col<>(
            const std::vector<ts::Tensor<T>> &x,
            std::vector<unsigned> kernelDim
        );
        friend std::vector<ts::Tensor<T>> col2im<>(
            const ts::Tensor<T> &x,
            std::vector<unsigned> outputDim
        );
    };



        // ts::Tensor

    template <typename T>
    class ts::Tensor {
    private:
        Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> value;
        ts::WengertList<T> * wList = NULL;
        int index;

        // We want this constructor to be private as it is supposed to be called by
        // our friends overloaded operators and functions only. This constructor
        // thus allows us to create a Tensor with dependencies in the Wengert list.
        Tensor(
            Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> newValue,
            ts::WengertList<T> * newWList, std::shared_ptr<ts::Node<T>> node
        );

    public:

        Tensor() {};

        // Input tensor, part of model
        Tensor(
            Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> newValue,
            ts::WengertList<T> * newWList
        );

        // Non part of model input tensor
        // (equivalent to calling previous constructor with model = false)
        Tensor(
            Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> newValue,
            ts::WengertList<T> * newWList, bool model
        );

        Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> getValue();
        ts::Gradient<T> grad();


        friend ts::WengertList<T>;

        friend ts::Gradient<T>;
        friend ts::GaElement<T>;
        friend ts::GradientAccumulator<T>;

        friend ts::Tensor<T> operator+<>(const ts::Tensor<T> &x, const ts::Tensor<T> &y);
        friend ts::Tensor<T> operator-<>(const ts::Tensor<T> &x, const ts::Tensor<T> &y);
        friend ts::Tensor<T> operator*<>(const ts::Tensor<T> &x, const ts::Tensor<T> &y);
        friend ts::Tensor<T> operator/<>(const ts::Tensor<T> &x, const ts::Tensor<T> &y);

        friend ts::Tensor<T> matProd<>(const ts::Tensor<T> &x, const ts::Tensor<T> &y);
        friend ts::Tensor<T> sigmoid<>(const ts::Tensor<T> &x);
        friend ts::Tensor<T> relu<>(const ts::Tensor<T> &x);
        friend ts::Tensor<T> leakyRelu<>(const ts::Tensor<T> &x);
        friend ts::Tensor<T> rescale<>(const ts::Tensor<T> &x);
        friend ts::Tensor<T> squaredNorm<>(const ts::Tensor<T> &x);

        friend ts::Tensor<T> convolution<>(const ts::Tensor<T> &mat, const ts::Tensor<T> &ker);
        friend ts::Tensor<T> maxPooling<>(const ts::Tensor<T> &x, std::vector<unsigned> pool);
        friend std::vector<ts::Tensor<T>> split<>(
            const ts::Tensor<T> &x,
            ChannelSplit channelSplit,
            unsigned nInputChannels
        );
        friend ts::Tensor<T> vertCat<>(const std::vector<ts::Tensor<T>> &x);
        friend ts::Tensor<T> flattening<>(const ts::Tensor<T> &x);
        friend ts::Tensor<T> im2col<>(
            const std::vector<ts::Tensor<T>> &x,
            std::vector<unsigned> kernelDim
        );
        friend std::vector<ts::Tensor<T>> col2im<>(
            const ts::Tensor<T> &x,
            std::vector<unsigned> outputDim
        );
    };

    template <typename T> std::string serializeTensor(ts::Tensor<T> &tensor);

	template <typename T> ts::Tensor<T> parseTensor(
		std::ifstream &in, ts::WengertList<T> * wList
	);

	template <typename T> std::string serializeTensorsVector(
		std::vector<ts::Tensor<T>> &tensorsVector
	);

	template <typename T> std::vector<ts::Tensor<T>> parseTensorsVector(
		std::ifstream &in, ts::WengertList<T> * wList
	);


    template <typename T>
    Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> convArray(
        const Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> &mat,
        const Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> &ker
    );

    template <typename T> class ConvolutionNode;
    template <typename T>
    ts::Tensor<T> convolution(const ts::Tensor<T> &mat, const ts::Tensor<T> &ker);

    template <typename T> class SplitNode;
    template <typename T>
    std::vector<ts::Tensor<T>> split(
        const ts::Tensor<T> &x,
        ChannelSplit channelSplit,
        unsigned nInputChannels
    );

    template <typename T> class PoolingNode;
    template <typename T>
    ts::Tensor<T> maxPooling(const ts::Tensor<T> &x, std::vector<unsigned> pool);

    template <typename T> class VertCatNode;
    template <typename T>
    ts::Tensor<T> vertCat(const std::vector<ts::Tensor<T>> &x);

    template <typename T> class FlatteningNode;
    template <typename T>
    ts::Tensor<T> flattening(const ts::Tensor<T> &x);

    template <typename T> class Im2ColNode;
    template <typename T>
    ts::Tensor<T> im2col(
        const std::vector<ts::Tensor<T>> &x,
        std::vector<unsigned> kernelDim
    );

    template <typename T> class Col2ImNode;
    template <typename T>
    std::vector<ts::Tensor<T>> col2im(
        const ts::Tensor<T> &x,
        std::vector<unsigned> outputDim
    );
}



    // ts::ConvolutionNode

template <typename T>
class ts::ConvolutionNode : public ts::Node<T> {
private:
    using ts::Node<T>::Node;

    Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> incrementGradient(
            Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> &childDerivative,
            unsigned &j
    );
};



    // ts::PoolingNode

template <typename T>
class ts::PoolingNode : public ts::Node<T> {
private:
    using ts::Node<T>::Node;

    PoolingNode(
        std::vector<long> shape,
        Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> xVal, int xDep,
        std::vector<unsigned> newPool
    );

    Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> incrementGradient(
            Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> &childDerivative,
            unsigned &j
    );

    std::vector<unsigned> pool = {};

    friend ts::Tensor<T> ts::maxPooling<>(
        const ts::Tensor<T> &x, std::vector<unsigned> pool
    );
};



    // ts::SplitNode

template <typename T>
class ts::SplitNode : public ts::Node<T> {
private:
    using ts::Node<T>::Node;

    SplitNode(
        std::vector<long> shape,
        int xDep,
        std::vector<long> originalShape,
        ChannelSplit newSplitDirection,
        unsigned newPosition
    );

    Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> incrementGradient(
            Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> &childDerivative,
            unsigned &j
    );

    long originalRows, originalCols;
    ChannelSplit splitDirection;
    unsigned position;

    friend std::vector<ts::Tensor<T>> ts::split<>(
        const ts::Tensor<T> &x,
        ChannelSplit channelSplit,
        unsigned nInputChannels
    );
};



    // ts::VertCatNode

template <typename T>
class ts::VertCatNode : public ts::Node<T> {
private:
    using ts::Node<T>::Node;

    // This node can have n parents !
    VertCatNode(
        std::vector<long> shape,
        std::vector<int> newDependencies,
        std::vector<long> newHeights
    );

    Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> incrementGradient(
            Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> &childDerivative,
            unsigned &j
    );

    std::vector<long> heights = {};

    friend ts::Tensor<T> ts::vertCat<>(const std::vector<ts::Tensor<T>> &x);
};



    // ts::FlatteningNode

template <typename T>
class ts::FlatteningNode : public ts::Node<T> {
private:
    using ts::Node<T>::Node;

    FlatteningNode(
        std::vector<long> shape,
        Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> xVal, int xDep,
        std::vector<long> newSize
    );

    Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> incrementGradient(
            Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> &childDerivative,
            unsigned &j
    );

    std::vector<long> size = {};

    friend ts::Tensor<T> ts::flattening<>(const ts::Tensor<T> &x);
};



    // ts::Im2ColNode

template <typename T>
class ts::Im2ColNode : public ts::Node<T> {
private:
    using ts::Node<T>::Node;

    // This node can have n parents !
    Im2ColNode(
        std::vector<long> shape,
        std::vector<int> newDependencies,
        std::vector<long> newKernelDim,
        std::vector<long> newMatrixDim,
        unsigned newNChannels
    );

    Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> incrementGradient(
            Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> &childDerivative,
            unsigned &j
    );

    std::vector<long> kernelDim = {};
    std::vector<long> matrixDim = {};	// Size of one channel
    unsigned nChannels;	// Input nChannels

    friend ts::Tensor<T> ts::im2col<>(
        const std::vector<ts::Tensor<T>> &x,
        std::vector<unsigned> kernelDim
    );
};



    // ts::Col2ImNode

template <typename T>
class ts::Col2ImNode : public ts::Node<T> {
private:
    using ts::Node<T>::Node;

    // This node can have n parents !
    Col2ImNode(
        std::vector<long> shape,
        int xDep,
        unsigned newPosition,
        long newNChannels
    );

    Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> incrementGradient(
            Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> &childDerivative,
            unsigned &j
    );

    unsigned position;
    unsigned nChannels;

    friend std::vector<ts::Tensor<T>> ts::col2im<>(
        const ts::Tensor<T> &x,
        std::vector<unsigned> outputDim
    ) ;
};    

template <typename T>
ts::Tensor<T> NewTensor(
    Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> newValue,
    ts::WengertList<T> * newWList
);



    // ts::Gradient

template <typename T>
class ts::Gradient {
private:
    // Constructor is private since we want instances of this class to be
    // generated by the Tensor::grad() method only
    Gradient(
        std::vector< Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> > newDerivatives
    );

    std::vector< Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> > derivatives;

public:
    Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> getValue(ts::Tensor<T> a);
    bool isEmpty();

    friend class ts::Tensor<T>;
    friend class ts::GradientAccumulator<T>;
    friend class ts::AdamOptimizer<T>;
};
