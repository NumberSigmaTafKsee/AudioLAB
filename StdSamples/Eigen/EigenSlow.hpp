#pragma once

#include <vector>
#include <memory>
#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>

#define BARWIDTH 30

namespace ts {
	std::vector<std::string> split(std::string str, char delimeter);

	std::string serializeUnsignedVec2D(
		std::vector<std::vector<unsigned>> &vec2d
	);

	std::vector<std::vector<unsigned>> parseUnsignedVec2D(
		std::ifstream &in
	);

	void progressBar(unsigned current, unsigned max);
}


/*
* General automatic differentiation engine based on a Wengert list
* implementation. Reverse mode only.
*/



namespace ts {
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
}



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



namespace ts {
	template <typename T> class Model;

	template <typename T> class Polynom;
	template <typename T> class MultiLayerPerceptron;
	template <typename T> class ConvolutionalNetwork;

	// Friends forward declaration
	template <typename T> class GradientAccumulator;
}

namespace ts {
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
}

// Enum for channel splitting directions in CNN
// (declared outside for now because scoped enum declarationb seems
// impossible)
enum class ts::ChannelSplit : int {
	NOSPLIT,
	SPLIT_HOR,	// Splits lines
	SPLIT_VERT	// Splits columns
};


namespace ts {

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



	// ts::Model

template <typename T>
class ts::Model {
private:

public:
	ts::WengertList<T> wList;

	// Call the WengertList toggleOptimize method
	void toggleOptimize(ts::Tensor<T> * tensor, bool enable);

	// Helper function to optimize the whole model
	virtual void toggleGlobalOptimize(bool enable) = 0;

	// General method for computing the model forward pass
	virtual ts::Tensor<T> compute(ts::Tensor<T> input) = 0;

	// Serializes / parses model into / from a file
	virtual void save(std::string filePath) = 0;
	virtual void load(std::string filePath) = 0;

	friend ts::GradientAccumulator<T>;
};



	// ts::Polynom
	// (element-wise polynom for nxn tensors)

template <typename T>
class ts::Polynom : public ts::Model<T> {
private:
	long nRows = 0;
	long nCols = 0;

public:
	Polynom(unsigned order, std::vector<long> size);

	std::vector<ts::Tensor<T>> coefficients = {};

	void toggleGlobalOptimize(bool enable);

	ts::Tensor<T> compute(ts::Tensor<T> input);

	void save(std::string filePath);
	void load(std::string filePath);

	long rows();
	long cols();
};



	// ts::MultiLayerPerceptron

template <typename T>
class ts::MultiLayerPerceptron : public ts::Model<T> {
private:

public:
	MultiLayerPerceptron(unsigned inputSize, std::vector<unsigned> layers);

	ts::Tensor<T> (*activationFunction)(const ts::Tensor<T>&) = &(ts::relu);
	ts::Tensor<T> (*finalActivation)(const ts::Tensor<T>&) = &(ts::sigmoid);

	std::vector<ts::Tensor<T>> weights = {};
	std::vector<ts::Tensor<T>> biases = {};

	void toggleGlobalOptimize(bool enable);

	ts::Tensor<T> compute(ts::Tensor<T> input);

	void save(std::string filePath);
	void load(std::string filePath);
};



	// ts::ConvolutionalNetwork

template <typename T>
class ts::ConvolutionalNetwork : public ts::Model<T> {
private:

public:
	ConvolutionalNetwork(
		std::vector<unsigned> inputSize,
		ChannelSplit splitDirection, unsigned inputChannels,
		std::vector<std::vector<unsigned>> convLayers,
		std::vector<std::vector<unsigned>> poolingLayers,
		std::vector<unsigned> denseLayers
	);

	ts::Tensor<T> (*convActivation)(const ts::Tensor<T>&) = &(ts::leakyRelu);
	ts::Tensor<T> (*denseActivation)(const ts::Tensor<T>&) = &(ts::relu);
	ts::Tensor<T> (*finalActivation)(const ts::Tensor<T>&) = &(ts::sigmoid);



	// Convolution section
	std::vector<ts::Tensor<T>> convKernels = {};
	std::vector<ts::Tensor<T>> convBiases = {};
	std::vector<std::vector<unsigned>> pooling;
	std::vector<std::vector<unsigned>> kernelDims;
	std::vector<std::vector<unsigned>> outputDims;	// Outputs right after convs

	// Dense section
	std::vector<ts::Tensor<T>> weights = {};
	std::vector<ts::Tensor<T>> fullBiases = {};

	ChannelSplit channelSplit = ChannelSplit::NOSPLIT;
	unsigned nInputChannels = 1;

	void toggleGlobalOptimize(bool enable);

	ts::Tensor<T> compute(ts::Tensor<T> input);

	void save(std::string filePath);
	void load(std::string filePath);
};

namespace ts {
	template <typename T> class TrainingData;

	template <typename T> class GaElement;
	template <typename T> class GradientAccumulator;

	template <typename T> class Optimizer;
	template <typename T> class GradientDescentOptimizer;
	template <typename T> class AdamOptimizer;
};



	// ts::TrainingData
	// (helper class containing both input data and its expected result)

template <typename T>
class ts::TrainingData {
private:

public:
	TrainingData(
		Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> newInput,
		Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> newExpected
	);

	Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> input;
	Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> expected;
};



	// ts::GaElement
	// (Gradient accumulator element that keeps summed derivatives and index of
	// an optimized node/tensor)

template <typename T>
class ts::GaElement {
private:
	GaElement(ts::Tensor<T> * inputTensor);

	Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> gradSum;
	unsigned index;	// in the ts::WengertList / ts::Gradient (see WARNING below)

	void reset();

public:

	friend ts::GradientAccumulator<T>;
	friend ts::Optimizer<T>;
	friend ts::GradientDescentOptimizer<T>;
	friend ts::AdamOptimizer<T>;
};



	// ts::GradientAccumulator
	// A collection of accumulated gradient elements for all optimizable tensors
	// of a model.
	// WARNING The indices in this array are completely independent of those in
	// the Wengert List, since some nodes may not be input or optimizable.

template <typename T>
class ts::GradientAccumulator {
private:
	GradientAccumulator();
	GradientAccumulator(ts::Model<T> &model);

	std::vector<ts::GaElement<T>> elements = {};

	void reset();
	void increment(ts::Gradient<T> &gradient);
	void updateTensor(
		ts::Model<T> &model, unsigned i,
		Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> value
	);
	void clear();


public:

	friend ts::Optimizer<T>;
	friend ts::GradientDescentOptimizer<T>;
	friend ts::AdamOptimizer<T>;
};



	// ts::Optimizer

template <typename T>
class ts::Optimizer {
private:

protected:
	ts::GradientAccumulator<T> gradAccumulator;

	void resetGradAccumulator();	// Set values to 0
	void setupGradAccumulator(ts::Model<T> &model);	// Generate 0-filled elements

	// Dependent of optimizer type. Applies and the accumulated gradient.
	virtual void updateModel(ts::Model<T> &model, unsigned batchSize) = 0;

public:
	Optimizer();

	ts::Tensor<T> (*normFunction)(const ts::Tensor<T>&) = &(ts::squaredNorm);

	unsigned epochs = 1;

	// Optimizes the model by running its compute() method on the batches data
	virtual std::vector<std::vector<std::vector< T >>> run(
		ts::Model<T> &model, std::vector<std::vector< ts::TrainingData<T> >> &batches
	) = 0;

};



	// ts::GradientDescentOptimizer

template <typename T>
class ts::GradientDescentOptimizer : public ts::Optimizer<T> {
private:
	void updateModel(ts::Model<T> &model, unsigned batchSize);

public:
	using ts::Optimizer<T>::Optimizer;

	GradientDescentOptimizer(T newLearningRate);

	T learningRate = 0.1;

	std::vector<std::vector<std::vector< T >>> run(
		ts::Model<T> &model, std::vector<std::vector< ts::TrainingData<T> >> &batches
	);
};



	// ts::AdamOptimizer

template <typename T>
class ts::AdamOptimizer : public ts::Optimizer<T> {
private:
	void updateModel(ts::Model<T> &model, unsigned batchSize);

	// Moment estimates, initialized at run time
	std::vector<Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic>> m = {};
	std::vector<Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic>> v = {};

	// NOTE For now, all mHat and vHat are stored in a vector. It would be
	// possible to only use one Eigen::Array, and to resize it every time
	// instead of storing all values.
	// -> This is a tradeoff between time and memory, it would be interesting
	// to benchmark both methods in the future.
	std::vector<Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic>> mHat = {};
	std::vector<Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic>> vHat = {};

	void initMomentEstimates(
		std::vector< std::shared_ptr<ts::Node<T>> > nodes
	);

	void computeIncrement(
		std::vector< Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> >& derivatives,
		std::vector<ts::GaElement<T>>& elements
	);

	T decayedBeta1;
	T decayedBeta2;


public :
	using ts::Optimizer<T>::Optimizer;

	AdamOptimizer(
		T newAlpha,
		T newBeta1, T newBeta2,
		T newEpsilon
	);

	// Default values set according to original paper
	T alpha = 0.001;
	T beta1 = 0.9;
	T beta2 = 0.999;
	T epsilon = 0.00000001;

	std::vector<std::vector<std::vector< T >>> run(
		ts::Model<T> &model, std::vector<std::vector< ts::TrainingData<T> >> &batches
	);
};


/*
* General automatic differentiation engine based on a Wengert list
* implementation. Reverse mode only.
*/

#
	// ts::Node

template <typename T>
ts::Node<T>::Node(std::vector<long> shape) {
	rows = shape[0];
	cols = shape[1];
};



template <typename T>
ts::Node<T>::Node(
	std::vector<long> shape,
	Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> xVal, int xDep
) {
	rows = shape[0];
	cols = shape[1];

	values =  {xVal};	// [da/dx]
	dependencies =  {xDep};
}



template <typename T>
ts::Node<T>::Node(
	std::vector<long> shape,
	Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> xVal, int xDep,
	Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> yVal, int yDep
) {
	rows = shape[0];
	cols = shape[1];

	values =  {xVal, yVal};	// [da/dx, da/dy]
	dependencies =  {xDep, yDep};
}



template <typename T>
Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> ts::InputNode<T>::incrementGradient(
		Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> &childDerivative,
		unsigned &j
) {

	// Used in the  ts::Tensor::grad() method. For an input node, this
	// function should never be called

	Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> increment;
	return increment;
}



template <typename T>
Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> ts::ElementWiseNode<T>::incrementGradient(
		Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> &childDerivative,
		unsigned &j
) {

	// Used in the  ts::Tensor::grad() method. Computes the increment of a derivative
	// for an element-wise operation

	Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> increment;
	increment = this->values[j] * childDerivative;

	return increment;
}



template <typename T>
ts::InputNode<T>::InputNode(std::vector<long> shape, bool model) {
	this->rows = shape[0];
	this->cols = shape[1];

	isModel = model;
};



template <typename T>
ts::MatProdNode<T>::MatProdNode(
	std::vector<long> shape,
	Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> xVal, int xDep,
	Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> yVal, int yDep,
	std::vector<long int> newXSize, std::vector<long int> newYSize
) {

	// MatProdNode specific constructor to store the size of the operands

	this->rows = shape[0];
	this->cols = shape[1];

	this->values =  {xVal, yVal};	// [da/dx, da/dy]
	this->dependencies =  {xDep, yDep};

	xSize = newXSize;
	ySize = newYSize;
}



template <typename T>
Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> ts::MatProdNode<T>::incrementGradient(
		Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> &childDerivative,
		unsigned &j
) {

	// Used in the  ts::Tensor::grad() method. Computes the increment of a derivative
	// for a matrix-matrix product.

	Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> increment;

	// Make sure operands are at the correct position for product x * y

	// Incremening x
	if(
		xSize[1] == this->values[j].rows() &&
		xSize[0] == this->values[j].cols()
	) {
		increment = (this->values[j].matrix() * childDerivative.matrix()).array();
	}

	// Incrementing y
	else if(
		ySize[1] == this->values[j].rows() &&
		ySize[0] == this->values[j].cols()
	) {
		increment = (childDerivative.matrix() * this->values[j].matrix() ).array();
	}

	else {
		exit(-1);
	}

	return increment;
}



template <typename T>
Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> ts::ScalarNode<T>::incrementGradient(
		Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> &childDerivative,
		unsigned &j
) {

	// Used in the ts::Tensor::grad() method. Computes the increment of a derivative
	// for a tensor to scalar operation.

	Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> increment;
	increment = this->values[j] * childDerivative(0, 0);

	return increment;
}



	// ts::WengertList

template <typename T>
int ts::WengertList<T>::size() {
	return nodes.size();
}



template <typename T>
int ts::WengertList<T>::reset() {
	// Used to remove all nodes but the optimizable input nodes, so the input tensors can
	// be reused in new computations. Returns the new size of the list.

	// First pass : remove non optimizable variables
	for(unsigned i = nodes.size(); i-- > 0; ) {

		// If the node is not an input (has dependencies)
		if(nodes[i]->dependencies.size() != 0) {
			nodes.erase(nodes.begin() + i);
		}

		// Input node
		else {
			std::shared_ptr<ts::InputNode<T>> inputPtr =
			std::static_pointer_cast<ts::InputNode<T>>(nodes[i]);


			// If the node is not part of model (probably model input)
			if(!(inputPtr->isModel)) {
				nodes.erase(nodes.begin() + i);
			}
		}
	}

	// Second pass : update tensors indices
	for(unsigned i = nodes.size(); i-- > 0; ) {
		std::shared_ptr<ts::InputNode<T>> inputPtr =
		std::static_pointer_cast<ts::InputNode<T>>(nodes[i]);

		if(inputPtr->optimizedTensor != NULL) {
			inputPtr->optimizedTensor->index = i;
		}
	}


	return nodes.size();
}



template <typename T>
void ts::WengertList<T>::toggleOptimize(ts::Tensor<T> * tensor, bool enable) {

	std::shared_ptr<ts::InputNode<T>> inputPtr =
	std::static_pointer_cast<ts::InputNode<T>>(nodes[tensor->index]);

	if(enable) {
		inputPtr->optimizedTensor = tensor;
	}
	else {
		inputPtr->optimizedTensor = NULL;
	}
}



	// ts::Tensor

// Input and not part of model
template <typename T>
ts::Tensor<T>::Tensor(
	Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> newValue,
	ts::WengertList<T> * newWList
) {
	value = newValue;
	wList = newWList;

	if(wList != NULL) {
		// Add new Tensor to the Wengert list
		index = wList->nodes.size();

		// Node without dependencies (input var,)
		std::shared_ptr<ts::Node<T>> nodePtr (
			new ts::InputNode<T>({newValue.rows(), newValue.cols()}, false)
		);

		wList->nodes.push_back(nodePtr);
	} else {
		index = -1;
	}
};



// Input and part of model
template <typename T>
ts::Tensor<T>::Tensor(
	Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> newValue,
	ts::WengertList<T> * newWList,
	bool model
) {
	value = newValue;
	wList = newWList;

	if(wList != NULL) {
		// Add new Tensor to the Wengert list
		index = wList->nodes.size();

		// Node without dependencies (input var,)
		std::shared_ptr<ts::Node<T>> nodePtr (
			new ts::InputNode<T>({newValue.rows(), newValue.cols()}, model)
		);

		wList->nodes.push_back(nodePtr);
	} else {
		index = -1;
	}
};



// Tensor with dependencies
template <typename T>
ts::Tensor<T>::Tensor(
	Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> newValue,
	ts::WengertList<T> * newWList, std::shared_ptr<ts::Node<T>> node
) {
	value = newValue;
	wList = newWList;

	if(wList != NULL) {
		// Add new Tensor to the Wengert list
		index = wList->nodes.size();
		wList->nodes.push_back(node);	// This node can contain dependencies & values
	} else {
		index = -1;
	}
}



// Helper function to create new instances without syntax template
template <typename T>
ts::Tensor<T> ts::NewTensor(
	Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> newValue,
	ts::WengertList<T> * newWList
) {
	return ts::Tensor<T>(newValue, newWList);
}



template <typename T>
Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> ts::Tensor<T>::getValue() {
	return value;
}



template <typename T>
ts::Gradient<T> ts::Tensor<T>::grad() {
	// Computes the gradient of this variable with respect to all the Wengert
	// list's nodes. Derivatives are stored in a vector wich size equals the
	// Wengert list's.

	// 2 possibilities :
	// - All operations are element wise, so we allow this tensor not to be a scalar
	// - Some operations change shapes of tensors, we only allow this tensor to be scalar

	// Making sure that we're not in case 2 with a non-scalar tensor
	if(!wList->elementWiseOnly && value.rows() != 1 && value.cols() != 1) {
		return ts::Gradient<T>({});
	}


	std::vector< Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> > derivatives(
		wList->nodes.size(),
		Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic>()
	);

	// Initialize all gradients with correct size zero-filled arrays
	for(unsigned i = 0; i < derivatives.size(); i++) {
		derivatives[i].setZero(wList->nodes[i]->rows, wList->nodes[i]->cols);
	}

	// Initialize gradient of self with respect to itself
	derivatives[index].fill(1.0);


	// Iterate over the Wengert list backwards
	for (unsigned i = wList->nodes.size(); i-- > 0; ) {

		std::shared_ptr<ts::Node<T>> node = wList->nodes[i];

		// Increment parent nodes
		for(unsigned j = 0; j < node->dependencies.size(); j++) {
			derivatives[node->dependencies[j]] += node->incrementGradient(
				derivatives[i], j
			);
		}
	}

	return ts::Gradient<T>(derivatives);
}



	// ts::Gradient

template <typename T>
ts::Gradient<T>::Gradient(
	std::vector< Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> > newDerivatives
) {
	derivatives = newDerivatives;
}



template <typename T>
Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> ts::Gradient<T>::getValue(ts::Tensor<T> a) {
	// Prevents segfault if tensor is out of bound for some reason
	if((unsigned) a.index >= derivatives.size()) {
		return Eigen::Array<T, 0, 0>();
	}

	return derivatives[a.index];
}



template <typename T>
bool ts::Gradient<T>::isEmpty() {
	// Used to look for errors after computing a gradient
	return derivatives.size() == 0 ? true : false;
}


	// Overloaded arithmetic operators

template <typename T>
ts::Tensor<T> ts::operator+(const ts::Tensor<T> &x, const ts::Tensor<T> &y){
	// Element-wise sum operation

	if(
		x.wList != y.wList ||
		x.value.rows() != y.value.rows() ||
		x.value.cols() != y.value.cols()
	) {
		return ts::Tensor<T>(Eigen::Array<T, 0, 0>(), NULL);
	}


	// a = x + y
	// da / dx = 1
	// da / dy = 1

	Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> grad;
	grad.setOnes(x.value.rows(), x.value.cols());

	std::shared_ptr<ts::Node<T>> nodePtr (
		new ts::ElementWiseNode<T>(
			{x.value.rows(), x.value.cols()},
			grad, x.index,
			grad, y.index
		)
	);

	return ts::Tensor<T>(x.value + y.value, x.wList, nodePtr);
}



template <typename T>
ts::Tensor<T> ts::operator-(const ts::Tensor<T> &x, const ts::Tensor<T> &y){
	// Element-wise difference operation

	if(
		x.wList != y.wList ||
		x.value.rows() != y.value.rows() ||
		x.value.cols() != y.value.cols()
	) {
		return ts::Tensor<T>(Eigen::Array<T, 0, 0>(), NULL);
	}


	// a = x - y
	// da / dx = 1
	// da / dy = -1

	Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> grad;
	grad.setOnes(x.value.rows(), x.value.cols());

	std::shared_ptr<ts::Node<T>> nodePtr (
		new ts::ElementWiseNode<T>(
			{x.value.rows(), x.value.cols()},
			grad, x.index,
			-1 * grad, y.index
		)
	);

	return ts::Tensor<T>(x.value - y.value,x.wList, nodePtr);
}



template <typename T>
ts::Tensor<T> ts::operator*(const ts::Tensor<T> &x, const ts::Tensor<T> &y){
	// Element-wise (Hadamard) product operation

	if(
		x.wList != y.wList ||
		x.value.rows() != y.value.rows() ||
		x.value.cols() != y.value.cols()
	) {
		return ts::Tensor<T>(Eigen::Array<T, 0, 0>(), NULL);
	}


	// a = x * y
	// da / dx = y
	// da / dy = x

	std::shared_ptr<ts::Node<T>> nodePtr (
		new ts::ElementWiseNode<T>(
			{x.value.rows(), x.value.cols()},
			y.value, x.index,
			x.value, y.index
		)
	);

	return ts::Tensor<T>(x.value * y.value,x.wList, nodePtr);
}



template <typename T>
ts::Tensor<T> ts::operator/(const ts::Tensor<T> &x, const ts::Tensor<T> &y){
	// Element-wise quotient operation

	if(
		x.wList != y.wList ||
		x.value.rows() != y.value.rows() ||
		x.value.cols() != y.value.cols()
	) {
		return ts::Tensor<T>(Eigen::Array<T, 0, 0>(), NULL);
	}


	// a = x / y
	// da / dx = 1 / y
	// da / dy = -x / y^2

	std::shared_ptr<ts::Node<T>> nodePtr (
		new ts::ElementWiseNode<T>(
			{x.value.rows(), x.value.cols()},
			1.0 / y.value, x.index,
			-x.value / (y.value * y.value), y.index
		)
	);

	return ts::Tensor<T>(x.value + y.value, x.wList, nodePtr);
}



	// Matrix product

template <typename T>
ts::Tensor<T> ts::matProd(const ts::Tensor<T> &x, const ts::Tensor<T> &y) {
	// Classic matrix-matrix product

	if(x.wList != y.wList || x.value.cols() != y.value.rows()) {
		return ts::Tensor<T>(Eigen::Array<T, 0, 0>(), NULL);
	}

	// The gradient will have to be computed for a scalar
	x.wList->elementWiseOnly = false;

	// a = x.y
	// dx = y^T	(transposed)
	// dy = x^T
	// (will be used in matrix product when computing gradient)

	std::shared_ptr<ts::Node<T>> nodePtr (
		new ts::MatProdNode<T>(
			{x.value.rows(), y.value.cols()},
			y.value.matrix().transpose(), x.index,
			x.value.matrix().transpose(), y.index,
			{x.value.rows(), x.value.cols()}, {y.value.rows(), y.value.cols()}
		)
	);

	return ts::Tensor<T>( x.value.matrix() * y.value.matrix(), x.wList, nodePtr);
}



	// Activation functions

template <typename T>
ts::Tensor<T> ts::sigmoid(const ts::Tensor<T> &x) {
	// Element-wise sigmoid function

	// a = e^x / (e^x + 1) = 1 / (1 + e^-x)
	// da / dx = e^x / (e^x + 1)^2

	std::shared_ptr<ts::Node<T>> nodePtr (
		new ts::ElementWiseNode<T>(
			{x.value.rows(), x.value.cols()},
			x.value.exp() / (x.value.exp() + 1).pow(2), x.index
		)
	);

	return ts::Tensor<T>(x.value.exp() / (x.value.exp() + 1), x.wList, nodePtr);
}



template <typename T>
ts::Tensor<T> ts::relu(const ts::Tensor<T> &x) {
	// Element-wise ReLU function
	// a = max(0, x)
	// da / dx = 0 if x<= 0 ; 1 if x > 0
	// Output is then rescaled between 0 and 1

	// Apply cwise max function
	Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> res;
	res.resize(x.value.rows(), x.value.cols());
	Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> dx;
	dx.resize(x.value.rows(), x.value.cols());


	// #pragma omp parallel for schedule(auto)
	for(unsigned i=0; i<res.size(); i++) {
		res(i) = (x.value(i) <= 0) ? 0 : x.value(i);
		dx(i) = (x.value(i) <= 0) ? 0 : 1;
	}


	// Return value
	std::shared_ptr<ts::Node<T>> nodePtr (
		new ts::ElementWiseNode<T>(
			{x.value.rows(), x.value.cols()},
			dx, x.index
		)
	);

	return ts::Tensor<T>(res, x.wList, nodePtr);

}



template <typename T>
ts::Tensor<T> ts::leakyRelu(const ts::Tensor<T> &x) {
	// Element-wise ReLU function
	// a = max(0, x)
	// da / dx = 0 if x<= 0 ; 1 if x > 0
	// Output is then rescaled between 0 and 1

	// Apply cwise max function
	Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> res;
	res.resize(x.value.rows(), x.value.cols());
	Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> dx;
	dx.resize(x.value.rows(), x.value.cols());


	// #pragma omp parallel for schedule(auto)
	for(unsigned i=0; i<res.size(); i++) {
			res(i) = (x.value(i) <= 0) ? 0.1 * x.value(i) : x.value(i);
			dx(i) = (x.value(i) <= 0) ? 0.1 : 1;
	}


	// Return value
	std::shared_ptr<ts::Node<T>> nodePtr (
		new ts::ElementWiseNode<T>(
			{x.value.rows(), x.value.cols()},
			dx, x.index
		)
	);

	return ts::Tensor<T>(res, x.wList, nodePtr);

}



template <typename T>
ts::Tensor<T> ts::rescale(const ts::Tensor<T> &x) {
	// Rescales tensor to 1
	// a = a / max(a)
	// da / dx = 1 / max(a)
	// Output is then rescaled between 0 and 1


	T max = x.value.maxCoeff();

	Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> res = x.value;
	if(max != 0) {
		res = res / max;
	}

	Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> dx;
	dx.setZero(x.value.rows(), x.value.cols());
	dx = dx + max;

	// Return value
	std::shared_ptr<ts::Node<T>> nodePtr (
		new ts::ElementWiseNode<T>(
			{x.value.rows(), x.value.cols()},
			dx, x.index
		)
	);

	return ts::Tensor<T>(res, x.wList, nodePtr);

}



	// Norm functions

template <typename T>
ts::Tensor<T> ts::squaredNorm(const ts::Tensor<T> &x) {
	// Returns the square of the 2-norm / euclidean norm of a vector

	// The gradient will have to be computed for a scalar
	x.wList->elementWiseOnly = false;

	// a = norm(x)^2
	// da / dx = 2x

	std::shared_ptr<ts::Node<T>> nodePtr (
		new ts::ScalarNode<T>(
			{1, 1}, 2 * x.value.matrix(), x.index
		)
	);

	Eigen::Array<T, 1, 1> res;
	res << (T) x.value.matrix().squaredNorm();

	return ts::Tensor<T>(res, x.wList, nodePtr);
}


	// Convolution operation on Eigen arrays
	// (LEGACY, for benchmarking purpose only)

template <typename T>
Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> ts::convArray(
	const Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> &mat,
	const Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> &ker
) {

	// Make sure kernel is smaller
	if(
		mat.rows() < ker.rows() ||
		mat.cols() < ker.cols()
	) {
		return Eigen::Array<T, 0, 0>();
	}

	unsigned newRows = mat.rows() - ker.rows() + 1;
	unsigned newCols = mat.cols() - ker.cols() + 1;


	Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> res;
	res.resize(newRows, newCols);

	for(unsigned i=0; i<newCols; i++) {
		for(unsigned j=0; j<newRows; j++) {
			// Compute one element of feature map
			res(j, i) =
			(mat.block(j, i, ker.rows(), ker.cols()) * ker).sum();
		}
	}

	return res;
}



	// Convolution
	// (LEGACY, for benchmarking purpose only)

template <typename T>
Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> ts::ConvolutionNode<T>::incrementGradient(
		Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> &childDerivative,
		unsigned &j
) {

	// Used in the  ts::Tensor::grad() method. Computes the increment of a derivative
	// for a convolution operation.

	Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> increment;

	// Matrices are already prepared at this stage, so we only need to put the
	// operands in the correct order for convolution.

	if(
		childDerivative.rows() > this->values[j].rows() &&
		childDerivative.rows() > this->values[j].cols()
	) {
		increment = ts::convArray(childDerivative, this->values[j]);
	} else {
		increment = ts::convArray(this->values[j], childDerivative);
	}

	return increment;
}



template <typename T>
ts::Tensor<T> ts::convolution(const ts::Tensor<T> &mat, const ts::Tensor<T> &ker) {
	// Convolution operation
	// Resulting matrix is of size : (mat.x - ker.x + 1, mat.y - ker.y + 1)
	// (where mat.x >= ker.x and mat.y >= mat.y)

	if(
		mat.wList != ker.wList ||
		mat.value.rows() < ker.value.rows() ||
		mat.value.cols() < ker.value.cols()
	) {
		return ts::Tensor<T>(Eigen::Array<T, 0, 0>(), NULL);
	}

	// The gradient will have to be computed for a scalar
	mat.wList->elementWiseOnly = false;


	// Compute res
	Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> res = ts::convArray(
		mat.value, ker.value
	);

	// Init dMat matrix (for matrix partial derivative)
	Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> dMat;
	dMat.setZero(
		2 * res.rows() + ker.value.rows() - 2,
		2 * res.cols() + ker.value.cols() - 2
	);

	dMat.block(
		res.rows() - 1,
		res.cols() - 1,
		ker.value.rows(), ker.value.cols()
	) = ker.value.rowwise().reverse().colwise().reverse();

	std::shared_ptr<ts::Node<T>> nodePtr (
		new ts::ConvolutionNode<T>(
			{res.rows(), res.cols()},
			dMat, mat.index,
			mat.value, ker.index
		)
	);

	return ts::Tensor<T>(res, mat.wList, nodePtr);
}



	// Max pooling

template <typename T>
ts::PoolingNode<T>::PoolingNode(
	std::vector<long> shape,
	Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> xVal, int xDep,
	std::vector<unsigned> newPool
) {

	// PoolingNode specific constructor to store the size of pools
	// (this allows us to easily upscale the matrix in grad computation)

	this->rows = shape[0];
	this->cols = shape[1];

	this->values =  {xVal};	// [da/dx]
	this->dependencies =  {xDep};

	pool = newPool;
}



template <typename T>
Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> ts::PoolingNode<T>::incrementGradient(
	Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> &childDerivative,
	unsigned &j
) {

	// Used in the  ts::Tensor::grad() method. Computes the increment of a derivative
	// for a max pooling / downsample operation.


	// Upsample matrix of child derivative to match original size

	Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> upsample;
	upsample.setZero(this->values[j].rows(), this->values[j].cols());


	// Affect coefficients of childDerivative to upsample pools by filling each
	// pool with the corresponding value

	for(unsigned i=0; i<childDerivative.cols(); i++) {
		for(unsigned j=0; j<childDerivative.rows(); j++) {

			// Fill one pool with one value
			for(unsigned k=0; k<pool[1]; k++) {
				for(unsigned l=0; l<pool[0]; l++) {
					upsample(j * pool[0] + l, i * pool[1] + k) =
					childDerivative(j, i);
				}
			}
		}
	}


	// Compute & return element-wise product with the this->values
	// (since this->values is 0/1-flled, we will only get the coefficients in
	// the desired positions, and 0 anywhere else)
	return upsample * this->values[j];
}



template <typename T>
ts::Tensor<T> ts::maxPooling(const ts::Tensor<T> &x, std::vector<unsigned> pool) {
	// Max pooling operation : we keep only the biggest element in each pool
	// in order to reduce the size of a matrix
	// Resulting matrix is of size : (mat.x / pool.x, mat.y / pool.y)
	// (where mat.x >= pool.x and mat.y >= pool.y)

	if(pool.size() != 2) {
		return ts::Tensor<T>(Eigen::Array<T, 0, 0>(), NULL);
	}

	if(
		x.value.rows() % pool[0] != 0 ||
		x.value.cols() % pool[1] != 0
	) {
		return ts::Tensor<T>(Eigen::Array<T, 0, 0>(), NULL);
	}

	// The gradient will have to be computed for a scalar
	x.wList->elementWiseOnly = false;


	// Init result
	Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> res;
	res.setZero(x.value.rows() / pool[0], x.value.cols() / pool[1]);


	// Init dx
	// (dx is 1 for each max element, 0 elsewhere)
	Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> dx;
	dx.setZero(x.value.rows(), x.value.cols());


	unsigned xMax, yMax;
	T maxVal;


	// Compute both pooled matrix (res) and dx
	// #pragma omp parallel for collapse(2) schedule(auto)
	for(unsigned i=0; i<res.cols(); i++) {
		for(unsigned j=0; j<res.rows(); j++) {

			// Get index of pool's max element
			// (for now it seems the best way is to manually iterate over
			// elements)

			xMax = j * pool[0];
			yMax = i * pool[1];
			maxVal = x.value(j * pool[0], i * pool[1]);

			for(unsigned k=0; k<pool[1]; k++) {
				for(unsigned l=0; l<pool[0]; l++) {

					if(x.value(j * pool[1] + l, i * pool[0] + k) > maxVal) {
						maxVal = x.value(j * pool[1] + l, i * pool[0] + k);
						xMax = j * pool[1] + l;
						yMax = i * pool[0] + k;
					}

				}
			}

			// Assigning values for result and derivative
			res(j, i) = maxVal;
			dx(xMax, yMax) = 1.0;

		}
	}


	std::shared_ptr<ts::Node<T>> nodePtr (
		new ts::PoolingNode<T>(
			{res.rows(), res.cols()},
			dx, x.index,
			pool
		)
	);

	return ts::Tensor<T>(res, x.wList, nodePtr);
}



	// Splitting

template <typename T>
ts::SplitNode<T>::SplitNode(
	std::vector<long> shape,
	int xDep,
	std::vector<long> originalShape,
	ChannelSplit newSplitDirection,
	unsigned newPosition
) {
	// SplitNode specific constructor to store the split direction

	this->dependencies =  {xDep};

	// New tensor shape (dimension of split matrices)
	this->rows = shape[0];
	this->cols = shape[1];

	// Original matrix shape
	originalRows = originalShape[0];
	originalCols = originalShape[1];


	splitDirection = newSplitDirection;
	position = newPosition;
}



template <typename T>
Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> ts::SplitNode<T>::incrementGradient(
	Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> &childDerivative,
	unsigned &j
) {

	// Used in the  ts::Tensor::grad() method. Computes the increment of a derivative
	// for a matrix split.

	// childDerivative is one of the resulting matrices. We will reconstruct
	// the partial derivative with regard to this considered matrix in order to
	// compute the increment. Index of the corresponding matrix is given by j.

	// Shape of base matrix derivative (initally zero filled)
	Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> increment;
	increment.setZero(originalRows, originalCols);

	// Affect childDerivative values to correct positions, according to
	// split direction & matrix index (j)
	if(splitDirection == ChannelSplit::SPLIT_VERT) {
		increment.block(0, position * this->cols, this->rows, this->cols) =
		childDerivative;
	}

	else if(splitDirection == ChannelSplit::SPLIT_HOR) {
		increment.block(position * this->rows, 0, this->rows, this->cols) =
		childDerivative;
	}

	return increment;
}



template <typename T>
std::vector<ts::Tensor<T>> ts::split(
	const ts::Tensor<T> &x,
	ChannelSplit channelSplit,
	unsigned nInputChannels
) {

	// The gradient will have to be computed for a scalar
	x.wList->elementWiseOnly = false;

	std::vector<ts::Tensor<T>> matrices = {};

	if(channelSplit == ChannelSplit::SPLIT_HOR) {
		unsigned channelSize = x.value.rows() / nInputChannels;

		for(unsigned i=0; i<nInputChannels; i++) {

			// Get matrix form of block
			Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> tmp =
			x.value.block(
				i * channelSize, 0, channelSize, x.value.cols()
			);

			// Create associated Tensor
			std::shared_ptr<ts::Node<T>> nodePtr (
				new ts::SplitNode<T>(
					{channelSize, x.value.cols()},
					x.index,
					{x.value.rows(), x.value.cols()},
					channelSplit,
					i
				)
			);

			matrices.push_back(ts::Tensor<T>(tmp, x.wList, nodePtr));
		}
	}

	if(channelSplit == ChannelSplit::SPLIT_VERT) {
		unsigned channelSize = x.value.cols() / nInputChannels;

		for(unsigned i=0; i<nInputChannels; i++) {

			// Get matrix form of block
			Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> tmp =
			x.value.block(
				0, i * channelSize, x.value.rows(), channelSize
			);

			// Create associated Tensor
			std::shared_ptr<ts::Node<T>> nodePtr (
				new ts::SplitNode<T>(
					{x.value.rows(), channelSize},
					x.index,
					{x.value.rows(), x.value.cols()},
					channelSplit,
					i
				)
			);

			matrices.push_back(ts::Tensor<T>(tmp, x.wList, nodePtr));

		}
	}

	if(channelSplit == ChannelSplit::NOSPLIT) {
		matrices.push_back(x);
	}

	return matrices;
}



	// Vertical concatenation

template <typename T>
ts::VertCatNode<T>::VertCatNode(
	std::vector<long> shape,
	std::vector<int> newDependencies,
	std::vector<long> newHeights
) {

	// VertCatNode specific constructor to store the height of first matrix
	// This way we can copy correct elements in inrementGradient

	// New tensor shape (vector)
	this->rows = shape[0];
	this->cols = shape[1];

	this->dependencies =  newDependencies;

	// Height of the first matrix
	heights = newHeights;
}



template <typename T>
Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> ts::VertCatNode<T>::incrementGradient(
	Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> &childDerivative,
	unsigned &j
) {

	// Used in the  ts::Tensor::grad() method. Computes the increment of a derivative
	// for a matrix flattening.

	// childDerivative is a flattened vector. We need to convert it back to a
	// matrix with the dimensions of the original matrix.

	Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> mat;
	mat.resize(heights[j+1] - heights[j], childDerivative.cols());

	mat = childDerivative.block(
		heights[j], 0,
		heights[j+1] - heights[j], childDerivative.cols()
	);


	return mat;
}



template <typename T>
ts::Tensor<T> ts::vertCat(const std::vector<ts::Tensor<T>> &x) {
	// Vertical concatenation operation
	// x[i] will be under x[i-1]

	if(x.size() == 0) {
		return ts::Tensor<T>(Eigen::Array<T, 0, 0>(), NULL);
	}


	// The gradient will have to be computed for a scalar
	x[0].wList->elementWiseOnly = false;


	// Compute size of resulting matrix, and storing each input matrix position
	// We will also make sure that all matrices have the same width / wList

	std::vector<long> heights = {0}; // Values are cumulative starting heights
	long height = 0;
	long width = x[0].value.cols();
	std::vector<int> dependencies = {};

	for(unsigned i=0; i<x.size(); i++) {

		if(x[i].value.cols() != width || x[i].wList != x[0].wList) {
			return ts::Tensor<T>(Eigen::Array<T, 0, 0>(), NULL);
		}

		height += x[i].value.rows();
		heights.push_back(height);
		dependencies.push_back(x[i].index);

	}


	// Set res vector
	Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> res;
	res.resize(height, width);

	for(unsigned i=0; i<x.size(); i++) {
		res.block(heights[i], 0, heights[i+1] - heights[i], width) = x[i].value;
	}


	// Return
	std::shared_ptr<ts::Node<T>> nodePtr (
		new ts::VertCatNode<T>(
			{res.rows(), res.cols()},
			dependencies,
			heights
		)
	);

	return ts::Tensor<T>(res, x[0].wList, nodePtr);
}



	// Flattening

template <typename T>
ts::FlatteningNode<T>::FlatteningNode(
	std::vector<long> shape,
	Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> xVal, int xDep,
	std::vector<long> newSize
) {

	// FlatteningNode specific constructor to store the size of original matrix
	// (this allows us to easily rescale the flattened vector in grad
	// computation)

	// New tensor shape (vector)
	this->rows = shape[0];
	this->cols = shape[1];

	this->values =  {xVal};	// [da/dx]
	this->dependencies =  {xDep};

	// Original matrix size
	size = newSize;
}



template <typename T>
Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> ts::FlatteningNode<T>::incrementGradient(
	Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> &childDerivative,
	unsigned &j
) {

	// Used in the  ts::Tensor::grad() method. Computes the increment of a derivative
	// for a matrix flattening.

	// childDerivative is a flattened vector. We need to convert it back to a
	// matrix with the dimensions of the original matrix.

	Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> mat = childDerivative;
	mat = Eigen::Map<Eigen::Array<T, -1, 1>>(
		mat.data(), mat.cols() * mat.rows()
	);
	mat.resize(size[0], size[1]);

	return mat;
}



template <typename T>
ts::Tensor<T> ts::flattening(const ts::Tensor<T> &x) {
	// Flattening operation to convert matrix to vector
	// A x matrix of size m*n becomes a (m * n, 1) vector
	// Conversion is column major by default (because of Eigen storing order)

	// The gradient will have to be computed for a scalar
	x.wList->elementWiseOnly = false;


	// Set res vector
	Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> res = x.value;
	res = Eigen::Map<Eigen::Array<T, -1, 1>>(
		res.data(), res.cols() * res.rows()
	);


	// Set dx matrix
	// It should be 1-filled since we're keeping all values of x in res, but
	// storing the full matrix would not be memory-efficient
	Eigen::Array<T, 0, 0> dx;


	// Return
	std::shared_ptr<ts::Node<T>> nodePtr (
		new ts::FlatteningNode<T>(
			{res.rows(), res.cols()},
			dx, x.index,
			{x.value.rows(), x.value.cols()}
		)
	);

	return ts::Tensor<T>(res, x.wList, nodePtr);
}



	// Im2col

template <typename T>
ts::Im2ColNode<T>::Im2ColNode(
	std::vector<long> shape,
	std::vector<int> newDependencies,
	std::vector<long> newKernelDim,
	std::vector<long> newMatrixDim,
	unsigned newNChannels
) {
	// New tensor shape (vector)
	this->rows = shape[0];
	this->cols = shape[1];

	this->dependencies =  newDependencies;

	// Original matrix size
	kernelDim = newKernelDim;
	matrixDim = newMatrixDim;
	nChannels = newNChannels;
}



template <typename T>
Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> ts::Im2ColNode<T>::incrementGradient(
		Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> &childDerivative,
		unsigned &j
) {

	// Used in the  ts::Tensor::grad() method. Computes the increment of a derivative
	// for a im2col operation.

	// childDerivative has the shape of the final matrix.
	// The increment will have the shape of one input matrix (this method will
	// be called once for each channel)

	Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> mat;
	mat.setZero(matrixDim[0], matrixDim[1]);

	// This matrix will be converted back to "normal" shape
	Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> im2colMat = childDerivative.block(
		j * kernelDim[0] * kernelDim[1], 0,
		kernelDim[0] * kernelDim[1], childDerivative.cols()
	);

	for(unsigned i=0; i<im2colMat.cols(); i++) {
		// Each column is a col-major flattened submatrix
		for(unsigned j=0; j<im2colMat.rows(); j++) {
			// Get top left coords of submatrix
			int submatTopX = i / (matrixDim[0] - kernelDim[0] + 1);
			int submatTopY = i % (matrixDim[0] - kernelDim[0] + 1);

			// Get coords in submatrix
			int submatX = j / kernelDim[1];
			int submatY = j % kernelDim[1];

			// Add derivative to coords in original matrix
			mat(submatTopX + submatX, submatTopY + submatY) =
			mat(submatTopX + submatX, submatTopY + submatY) + im2colMat(j, i);

		}
	}

	return mat;
}



template <typename T>
ts::Tensor<T> ts::im2col(
	const std::vector<ts::Tensor<T>> &x,
	std::vector<unsigned> kernelDim
) {
	// Turns a tensor vector into a single im2col matrix
	// Using a kernels matrix, one entire conv layer could be computed in
	// only one matrix product

	std::vector<int> dependencies = {};

	Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> res;
	res.resize(
		kernelDim[0] * kernelDim[1] * x.size(),
		(x[0].value.rows() - kernelDim[0] + 1) * (x[0].value.cols() - kernelDim[1] + 1)
	);

	for(unsigned i=0; i<x.size(); i++) {
		// #pragma omp parallel for collapse(2) schedule(auto)
		for(unsigned j=0; j<x[i].value.cols() - kernelDim[0] + 1; j++) {
			for(unsigned k=0; k<x[i].value.rows() - kernelDim[1] + 1; k++) {

				Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> tmp =
				x[i].value.block(k, j, kernelDim[0], kernelDim[1]);

				Eigen::Map<Eigen::Array<T, -1, 1>> map =
				Eigen::Map<Eigen::Array<T, -1, 1>>(
					tmp.data(), tmp.cols() * tmp.rows()
				);

				res.block(i * kernelDim[0] * kernelDim[1], k * (x[i].value.rows() - kernelDim[0] + 1)  + j, kernelDim[0] * kernelDim[1], 1)
				= map;

			}
		}

		dependencies.push_back(x[i].index);
	}


	// Return
	std::shared_ptr<ts::Node<T>> nodePtr (
		new ts::Im2ColNode<T>(
			{res.rows(), res.cols()},
			dependencies,
			{kernelDim[0], kernelDim[1]},
			{x[0].value.rows(), x[0].value.cols()},
			x.size()
		)
	);

	return ts::Tensor<T>(res, x[0].wList, nodePtr);
}



	// Col2im

template <typename T>
ts::Col2ImNode<T>::Col2ImNode(
	std::vector<long> shape,
	int xDep,
	unsigned newPosition,
	long newNChannels
) {
	// New tensor shape (vector)
	this->rows = shape[0];
	this->cols = shape[1];

	this->dependencies =  {xDep};

	// Original matrix size
	position = newPosition;
	nChannels = newNChannels;
}



template <typename T>
Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> ts::Col2ImNode<T>::incrementGradient(
		Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> &childDerivative,
		unsigned &j
) {

	Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> flat = childDerivative;
	flat = Eigen::Map<Eigen::Array<T, 1, -1>>(
		flat.data(), flat.cols() * flat.rows()
	);

	Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> res;
	res.setZero(nChannels, flat.cols());
	res.block(position, 0, 1, flat.cols()) = flat;

	return res;
}



template <typename T>
std::vector<ts::Tensor<T>> ts::col2im(
	const ts::Tensor<T> &x,
	std::vector<unsigned> outputDim
) {
	// Turns an im2col matrix into a channels vector
	// The output can be reused in another im2col, or
	// flattened before dense layers.

	std::vector<ts::Tensor<T>> res = {};

	// Each line contains some channel's coefficients in row-major order
	for(unsigned i=0; i<x.value.rows(); i++) {
		Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> tmp =
		x.value.block(i, 0, 1, x.value.cols());

		tmp.resize(outputDim[0], outputDim[1]);


		// Convert it back to matrix form
		std::shared_ptr<ts::Node<T>> nodePtr (
			new ts::Col2ImNode<T>(
				{tmp.rows(), tmp.cols()},
				x.index,
				i,
				x.value.rows()
			)
		);

		res.push_back(ts::Tensor<T>(tmp, x.wList, nodePtr));
	}

	return res;
}


// ts::Model

template <typename T>
void ts::Model<T>::toggleOptimize(ts::Tensor<T> * tensor, bool enable) {
	wList.toggleOptimize(tensor, enable);
}



	// ts::Polynom

template <typename T>
ts::Polynom<T>::Polynom(unsigned order, std::vector<long> size) {

	// +1 for deg 0 coefficient
	for(unsigned i=0; i<order+1; i++) {
		coefficients.push_back(ts::Tensor<T>(
			Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic>()
			.setRandom(size[0], size[1]),
			&(this->wList), true)
		);
	}

	// Size of tensors
	nRows = size[0];
	nCols = size[1];

}



template <typename T>
void ts::Polynom<T>::toggleGlobalOptimize(bool enable) {
	for(unsigned i=0; i<coefficients.size(); i++) {
		this->wList.toggleOptimize(&(coefficients[i]), enable);
	}
}



template <typename T>
ts::Tensor<T> ts::Polynom<T>::compute(ts::Tensor<T> input) {

	// Assert input and coefficients have the same size
	for(unsigned i=0; i<coefficients.size(); i++) {
		if(
			input.getValue().cols() != coefficients[i].getValue().cols() ||
			input.getValue().rows() != coefficients[i].getValue().rows()
		) {
			return ts::Tensor<T>(Eigen::Array<T, 0, 0>(), NULL);
		}
	}


	ts::Tensor<T> result = ts::Tensor<T>(
		Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(
			coefficients[0].getValue().cols(),
			coefficients[0].getValue().rows()
		),
		&(this->wList)
	);

	ts::Tensor<T> element;


	// Begin computation loop
	for(unsigned i=0; i<coefficients.size(); i++) {
		// Reset element
		element = ts::Tensor<T>(
			Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(
				coefficients[i].getValue().cols(),
				coefficients[i].getValue().rows()
			),
			&(this->wList)
		);

		element = element + coefficients[i];

		// Compute element
		for(unsigned j=0; j<i; j++) {
			element = element * input;
		}

		// Increment result
		result = result + element;
	}

	return result;
}



template <typename T>
void ts::Polynom<T>::save(std::string filePath) {
	std::ofstream out(filePath);
	out << ts::serializeTensorsVector(coefficients);
	out.close();
}



template <typename T>
void ts::Polynom<T>::load(std::string filePath) {
	// Delete current tensors and reset wList
	coefficients = {};
	this->wList.reset();

	// Load new tensors
	std::ifstream in(filePath);
	coefficients = ts::parseTensorsVector(in, &(this->wList));
	in.close();

	// Set number of wors and cols
	if(coefficients.size() > 0) {
		nRows = coefficients[0].getValue().rows();
		nCols = coefficients[0].getValue().cols();
	} else {
		nRows = 0;
		nCols = 0;
	}
}



template <typename T>
long ts::Polynom<T>::rows() {
	return nRows;
}



template <typename T>
long ts::Polynom<T>::cols() {
	return nCols;
}



	// ts::MultiLayerPerceptron

template <typename T>
ts::MultiLayerPerceptron<T>::MultiLayerPerceptron(
	unsigned inputSize, std::vector<unsigned> layers
) {
	// Each element of the layers vector is a new layer, its value represents
	// the layer size. Values are randomly initialized between 0 and 1.

	// Make sure layers/input are not of size 0
	if(inputSize == 0) {
		return;
	}

	for(unsigned i=1; i<layers.size(); i++) {
		if(layers[i] == 0) {
			return;
		}
	}


	layers.insert(layers.begin(), inputSize);

	for(unsigned i=1; i<layers.size(); i++) {
		// Initializing values according to He Initialization
		T variance = sqrt(2.0 / layers[i-1]);

		// Add layer i-1 -> layer i weights
		weights.push_back(ts::Tensor<T>(
			Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic>()
			.setRandom(layers[i], layers[i-1]) * variance,
			&(this->wList), true)
		);

		// Add layer i biases
		biases.push_back(ts::Tensor<T>(
			Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic>()
			.setRandom(layers[i], 1) * variance,
			&(this->wList), true)
		);
	}

}



template <typename T>
void ts::MultiLayerPerceptron<T>::toggleGlobalOptimize(bool enable) {
	for(unsigned i=0; i<weights.size(); i++) {
		this->toggleOptimize(&(weights[i]), enable);
		this->toggleOptimize(&(biases[i]), enable);
	}
}



template <typename T>
ts::Tensor<T> ts::MultiLayerPerceptron<T>::compute(ts::Tensor<T> input) {

	// Assert expected size
	if(input.getValue().rows() != weights[0].getValue().cols() ||
	input.getValue().cols() != 1) {
		return ts::Tensor<T>(Eigen::Array<T, 0, 0>(), NULL);
	}

	// Assert weights and biases vectors have the same size
	if(weights.size() != biases.size()) {
		return ts::Tensor<T>(Eigen::Array<T, 0, 0>(), NULL);
	}

	// Begin computation loop
	for(unsigned i=0; i<weights.size(); i++) {
		// Hidden layer
		if(i < weights.size() - 1) {
			input = (*activationFunction)(matProd(weights[i], input) + biases[i]);
		}
		// Final layer (we might want another activation function)
		else {
			input = (*finalActivation)(matProd(weights[i], input) + biases[i]);
		}
	}

	return input;
}



template <typename T>
void ts::MultiLayerPerceptron<T>::save(std::string filePath) {
	std::ofstream out(filePath);

	out << ts::serializeTensorsVector(weights);
	out << ts::serializeTensorsVector(biases);

	out.close();
}



template <typename T>
void ts::MultiLayerPerceptron<T>::load(std::string filePath) {
	// Delete current tensors and reset wList
	weights = {};
	biases = {};
	this->wList.reset();

	// Load new tensors
	std::ifstream in(filePath);

	weights = ts::parseTensorsVector(in, &(this->wList));
	biases = ts::parseTensorsVector(in, &(this->wList));

	in.close();
}



	// ts::ConvolutionalNetwork

template <typename T>
ts::ConvolutionalNetwork<T>::ConvolutionalNetwork(
	std::vector<unsigned> inputSize,
	ChannelSplit splitDirection, unsigned inputChannels,
	std::vector<std::vector<unsigned>> convLayers,
	std::vector<std::vector<unsigned>> poolingLayers,
	std::vector<unsigned> denseLayers
) {
	// inputSize : std::vector of size 3 for dimensions of 2D image / matrix
	//	+ number of channels (number of conv kernels for each layer)
	// convLayers : sizes of convolution kernels (std::vector of dimension 2)
	// fullLayers: sizes of fully connected layers


		// Validate dimensions of network

	if(inputSize.size() != 2) {
		std::cout << "ERROR: Input is not of dimension 2" << std::endl;
		return;
	}
	if(inputSize[0] == 0 || inputSize[1] == 0) {
		std::cout << "ERROR: Input is of size 0" << std::endl;
		return;
	}

	// Do we have an equal number of convolution and pooling layers ?
	if(convLayers.size() != poolingLayers.size()) {
		std::cout << "ERROR: Different numbers for convolution and pooling layers"
		<< std::endl;
		return;
	}


	std::vector<int> intermediarySize = {(int) inputSize[0], (int) inputSize[1]};

	// Make sure channel splitting is possible

	// Splitting rows
	if(splitDirection == ChannelSplit::SPLIT_HOR) {
		if(
			inputSize[0] % inputChannels != 0 ||
			inputSize[0] < inputChannels
		) {
			std::cout << "ERROR: Impossible to split horizontally"
			<< std::endl;
			return;
		}

		intermediarySize = {(int) inputSize[0] / (int) inputChannels, (int) inputSize[1]};
	}

	// Splitting cols
	else if(splitDirection == ChannelSplit::SPLIT_VERT) {
		if(
			inputSize[1] % inputChannels != 0 ||
			inputSize[1] < inputChannels
		) {
			std::cout << "ERROR: Impossible to split vertically"
			<< std::endl;
			return;
		}

		intermediarySize = {(int) inputSize[0], (int) inputSize[1] / (int) inputChannels};
	}

	// No split
	else {
		intermediarySize = {(int) inputSize[0], (int) inputSize[1]};
	}


	// Make sure convolutions / poolings are possible
	for(unsigned i=0; i<convLayers.size(); i++) {
		// Is size of kernel correctly described
		if(convLayers[i].size() != 3) {
			std::cout << "ERROR: Convolution layer " << i <<
			" is not of dimension 3" << std::endl;
			return;
		}
		// Are the different numbers of channels > 0 ?
		// (the numbers must be in growing order, and multipliers between each others)
		if(i != 0) {
			if(
				convLayers[i][2] == 0
			) {
				std::cout << "ERROR: Number of channels for " << i <<
				" is 0" << std::endl;
				return;
			}
		}

		// Is size of pooling correctly described
		if(poolingLayers[i].size() != 2) {
			std::cout << "ERROR: Pooling layer " << i <<
			" is not of dimension 2" << std::endl;
			return;
		}


		// Compute size of matrix after convolution
		intermediarySize[0] = intermediarySize[0] - convLayers[i][0] + 1;
		intermediarySize[1] = intermediarySize[1] - convLayers[i][1] + 1;


		if(intermediarySize[0] <= 0 || intermediarySize[1] <= 0) {
			std::cout << "ERROR: Convolution layer " << i <<
			" is impossible" << std::endl;
			return;
		}

		// Compute size of matrix after pooling
		if(poolingLayers[i][0] != 0 && poolingLayers[i][1] != 0) {
			if(
				intermediarySize[0] % poolingLayers[i][0] != 0 ||
				intermediarySize[1] % poolingLayers[i][1] != 0
			) {
				std::cout << "ERROR: Pooling layer " << i <<
				" is impossible" << std::endl;
				return;
			}

			intermediarySize[0] = intermediarySize[0] / poolingLayers[i][0];
			intermediarySize[1] = intermediarySize[1] / poolingLayers[i][1];
		}

	}


		// Randomly init kernels, weights and biases

	// Convolution layers

	// Splitting rows
	if(splitDirection == ChannelSplit::SPLIT_HOR) {
		intermediarySize = {(int) inputSize[0] / (int) inputChannels, (int) inputSize[1]};
	}

	// Splitting cols
	else if(splitDirection == ChannelSplit::SPLIT_VERT) {
		intermediarySize = {(int) inputSize[0], (int) inputSize[1] / (int) inputChannels};
	}

	// No split
	else {
		intermediarySize = {(int) inputSize[0], (int) inputSize[1]};
	}

	channelSplit = splitDirection;
	nInputChannels = inputChannels;

	// This is the input of the network, of dimension 1
	convLayers.insert(convLayers.begin(), {0, 0, inputChannels});

	convKernels = {};
	convBiases = {};

	for(unsigned i=1; i<convLayers.size(); i++) {

		// Initializing values according to He Initialization
		T variance = sqrt(2.0 / (convLayers[i][0] * convLayers[i][1] * convLayers[i-1][2]));

		convKernels.push_back(ts::Tensor<T>(
			Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic>()
			.setRandom(
				convLayers[i][2],
				convLayers[i][0] * convLayers[i][1] * convLayers[i-1][2]
			) * variance,
			&(this->wList), true)
		);

		intermediarySize[0] = intermediarySize[0] - convLayers[i][0] + 1;
		intermediarySize[1] = intermediarySize[1] - convLayers[i][1] + 1;

		outputDims.push_back(
			{(unsigned) intermediarySize[0], (unsigned) intermediarySize[1]}
		);

		convBiases.push_back(ts::Tensor<T>(
			Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic>()
			.setZero(
				convLayers[i][2],
				intermediarySize[0] * intermediarySize[1]
			),
			&(this->wList), true)
		);

		if(poolingLayers[i-1][0] != 0 && poolingLayers[i-1][1] != 0) {
			intermediarySize[0] = intermediarySize[0] / poolingLayers[i-1][0];
			intermediarySize[1] = intermediarySize[1] / poolingLayers[i-1][1];
		}
	}

	// Fully connected layers
	denseLayers.insert(
		// First dense layer input will be the flattened convolution output
		denseLayers.begin(),
		intermediarySize[0] * intermediarySize[1] * convLayers[convLayers.size() - 1][2]
	);

	for(unsigned i=1; i<denseLayers.size(); i++) {
		// Initializing values according to He Initialization
		T variance = sqrt(2.0 / denseLayers[i-1]);

		// Add layer i-1 -> layer i weights
		weights.push_back(ts::Tensor<T>(
			Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic>()
			.setRandom(denseLayers[i], denseLayers[i-1]) * variance,
			&(this->wList), true)
		);

		// Add layer i biases
		fullBiases.push_back(ts::Tensor<T>(
			Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic>()
			.setRandom(denseLayers[i], 1) * variance,
			&(this->wList), true)
		);
	}

	// Set up data fields
	pooling = poolingLayers;
	convLayers.erase(convLayers.begin());
	kernelDims = convLayers;
}



template <typename T>
void ts::ConvolutionalNetwork<T>::toggleGlobalOptimize(bool enable) {
	if(convKernels.size() != convBiases.size()) {
		return;
	}

	if(weights.size() != fullBiases.size()) {
		return;
	}

	for(unsigned i=0; i<convKernels.size(); i++) {
		this->toggleOptimize(&(convKernels[i]), enable);
		this->toggleOptimize(&(convBiases[i]), enable);
	}

	for(unsigned i=0; i<weights.size(); i++) {
		this->toggleOptimize(&(weights[i]), enable);
		this->toggleOptimize(&(fullBiases[i]), enable);
	}
}



template <typename T>
ts::Tensor<T> ts::ConvolutionalNetwork<T>::compute(ts::Tensor<T> input) {

	// NOTE It might be a good idea to add an entire function to make sure that
	// all parameters are compatible (in terms of size), and that output is
	// computable


	// Convert input to 2D vector (for number of channels) for use with the
	// im2col method. This should be a faster way to compute convolutions.

	std::vector<ts::Tensor<T>> inputVec = {};

	if(channelSplit != ChannelSplit::NOSPLIT) {
		inputVec = ts::split(input, channelSplit, nInputChannels);
	}
	else {
		inputVec.push_back(input);
	}


	// 1) Convolution / pooling computation loop
	for(unsigned i=0; i<convKernels.size(); i++) {
		// Compute the im2col multichannel convolution
		input = ts::im2col(inputVec, kernelDims[i]);
		input = (*convActivation)(matProd(convKernels[i], input) + convBiases[i]);
		inputVec = ts::col2im(input,  outputDims[i]);

		// A pooling layer of size 0 means we want to skip it
		if(pooling[i][0] != 0 || pooling[i][1] != 0) {
			for(unsigned j=0; j<inputVec.size(); j++) {
				inputVec[j] = ts::maxPooling(inputVec[j], pooling[i]);
			}
		}
	}


	// 2) Gather all channels back to input tensor,
	// and flatten convolution outputs
	input = vertCat(inputVec);
	input = flattening(input);


	// 3) Dense layers computation loop
	for(unsigned i=0; i<weights.size(); i++) {
		if(i < weights.size() - 1) {
			input = (*denseActivation)(matProd(weights[i], input) + fullBiases[i]);
		}
		// Final layer (we might want another activation function)
		else {
			input = (*finalActivation)(matProd(weights[i], input) + fullBiases[i]);
		}
	}

	return input;
}



template <typename T>
void ts::ConvolutionalNetwork<T>::save(std::string filePath) {
	std::ofstream out(filePath);

	out << static_cast<std::underlying_type<ts::ChannelSplit>::type>(channelSplit) << std::endl;
	out << nInputChannels << std::endl;

	out << ts::serializeUnsignedVec2D(pooling);
	out << ts::serializeUnsignedVec2D(kernelDims);
	out << ts::serializeUnsignedVec2D(outputDims);

	out << ts::serializeTensorsVector(convKernels);
	out << ts::serializeTensorsVector(convBiases);

	out << ts::serializeTensorsVector(weights);
	out << ts::serializeTensorsVector(fullBiases);

	out.close();
}



template <typename T>
void ts::ConvolutionalNetwork<T>::load(std::string filePath) {
	// Delete current model, reset wList
	convKernels = {};
	convBiases = {};
	weights = {};
	fullBiases = {};
	this->wList.reset();

	pooling = {};
	kernelDims = {};
	outputDims = {};


	// Load new model
	std::string line;
	std::ifstream in(filePath);

	std::getline(in, line);
	channelSplit = static_cast<ts::ChannelSplit>(std::stoi(line));

	std::getline(in, line);
	nInputChannels = std::stoi(line);

	pooling = ts::parseUnsignedVec2D(in);
	kernelDims = ts::parseUnsignedVec2D(in);
	outputDims = ts::parseUnsignedVec2D(in);

	convKernels = ts::parseTensorsVector(in, &(this->wList));
	convBiases = ts::parseTensorsVector(in, &(this->wList));

	weights = ts::parseTensorsVector(in, &(this->wList));
	fullBiases = ts::parseTensorsVector(in, &(this->wList));

	in.close();
}


	// ts::TrainingData

template <typename T>
ts::TrainingData<T>::TrainingData(
		Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> newInput,
		Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> newExpected
) {
	input = newInput;
	expected = newExpected;
}



	// ts::GaElement

template <typename T>
ts::GaElement<T>::GaElement(ts::Tensor<T> * inputTensor) {

	unsigned rows = inputTensor->value.rows();
	unsigned cols = inputTensor->value.cols();

	gradSum = Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(rows, cols);

	index = inputTensor->index;

}



template <typename T>
void ts::GaElement<T>::reset()  {
	// Reset the value of accumulated gradient
	gradSum.setZero();
}



	// ts::GradientAccumulator

template <typename T>
ts::GradientAccumulator<T>::GradientAccumulator() {
}



template <typename T>
ts::GradientAccumulator<T>::GradientAccumulator(ts::Model<T> &model) {

	// Reset wengertList in case it has been used before
	model.wList.reset();

	for(unsigned i=0; i<model.wList.nodes.size(); i++) {

		std::shared_ptr<ts::InputNode<T>> inputPtr =
		std::static_pointer_cast<ts::InputNode<T>>(model.wList.nodes[i]);

		// Check if it is associated with a tensor (== optimizable)
		if(inputPtr->optimizedTensor != NULL) {

			// Then append it to the gradient accumulator
			elements.push_back(ts::GaElement<T>(inputPtr->optimizedTensor));
		}
	}
}



template <typename T>
void ts::GradientAccumulator<T>::reset() {
	// #pragma omp parallel for
	for(unsigned i=0; i<elements.size(); i++) {
		elements[i].reset();
	}
}



template <typename T>
void ts::GradientAccumulator<T>::increment(ts::Gradient<T> &gradient) {
	// Increment all elements of gradAccumulator according to gradient

	// #pragma omp parallel for
	for(unsigned i=0; i<elements.size(); i++) {
		// We use two different indices systems here
		// (one for the wList/grad and one for the gradient accumulator)
		elements[i].gradSum += gradient.derivatives[elements[i].index];
	}
}



template <typename T>
void ts::GradientAccumulator<T>::updateTensor(
	ts::Model<T> &model, unsigned i,
	Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> value
) {
	// Update a tensor via the gradient accumulator
	std::shared_ptr<ts::InputNode<T>> inputPtr =
	std::static_pointer_cast<ts::InputNode<T>>(model.wList.nodes[i]);

	inputPtr->optimizedTensor->value -= value;
}



template <typename T>
void ts::GradientAccumulator<T>::clear() {
	// Empty elements
	elements.clear();
}



	// ts::Optimizer

template <typename T>
ts::Optimizer<T>::Optimizer() {
}



	// ts::GradientDescentOptimizer

template <typename T>
ts::GradientDescentOptimizer<T>::GradientDescentOptimizer(T newLearningRate) {
	learningRate = newLearningRate;
}



template <typename T>
void ts::GradientDescentOptimizer<T>::updateModel(
	ts::Model<T> &model, unsigned batchSize
) {
	// #pragma omp parallel for
	for(unsigned i=0; i<this->gradAccumulator.elements.size(); i++) {
		this->gradAccumulator.updateTensor(
			model, i,
			learningRate * this->gradAccumulator.elements[i].gradSum / batchSize
		);
	}
}



template <typename T>
std::vector<std::vector<std::vector< T >>> ts::GradientDescentOptimizer<T>::run(
	ts::Model<T> &model, std::vector<std::vector< ts::TrainingData<T> >> &batches
) {

		// Set up gradient accumulator (this also resets wList)

	this->gradAccumulator = ts::GradientAccumulator<T>(model);


		// Start running and training the model

		std::vector<std::vector<std::vector< T >>> losses(this->epochs, (std::vector<std::vector<T>>) {});

	// Epochs
	for(unsigned i=0; i<this->epochs; i++) {

		std::cout << "Epoch " << i + 1 << "/" <<  this->epochs << ": " << std::endl;

		losses[i] = std::vector<std::vector<T>>(batches.size(), (std::vector<T>) {});

		// Batches
		for(unsigned j=0; j<batches.size(); j++) {

			losses[i][j] = std::vector<T>(batches[j].size(), 0);

			// Data instance
			for(unsigned k=0; k<batches[j].size(); k++) {

				ts::Tensor<T> input = ts::Tensor<T>(
					batches[j][k].input, &(model.wList)
				);
				ts::Tensor<T> expected = ts::Tensor<T>(
					batches[j][k].expected, &(model.wList)
				);

				// Compute model and norm
				ts::Tensor<T> output = model.compute(input);
				ts::Tensor<T> norm = (*this->normFunction)(output - expected);

				// Get gradient and increment gradient accumulator
				ts::Gradient<T> gradient = norm.grad();
				this->gradAccumulator.increment(gradient);

				model.wList.reset();

				losses[i][j][k] = norm.getValue()(0, 0);
			}

			updateModel(model, batches[j].size());
			this->gradAccumulator.reset();

			ts::progressBar(j + 1, batches.size());
		}
		std::cout << std::endl << std::endl;
	}


		// Clean

	this->gradAccumulator.clear();
	model.wList.reset();

	return losses;
}



	// ts::AdamOptimizer

template <typename T>
ts::AdamOptimizer<T>::AdamOptimizer(
	T newAlpha,
	T newBeta1, T newBeta2,
	T newEpsilon
) {
	alpha = newAlpha;
	beta1 = newBeta1;
	beta2 = newBeta2;
	epsilon = newEpsilon;
}



template <typename T>
void ts::AdamOptimizer<T>::updateModel(
	ts::Model<T> &model, unsigned batchSize
) {
	// #pragma omp parallel for
	for(unsigned i=0; i<this->gradAccumulator.elements.size(); i++) {
		this->gradAccumulator.updateTensor(
			model, i,
			alpha * this->gradAccumulator.elements[i].gradSum / batchSize
		);
	}
}



template <typename T>
void ts::AdamOptimizer<T>::initMomentEstimates(
	std::vector< std::shared_ptr<ts::Node<T>> > nodes
) {
	// Initialize shape of moment estimates
	// (same size as reset wList, zero filled)

	Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> tmp;
	for(unsigned i = 0; i<nodes.size(); i++) {

		tmp.setZero(nodes[i]->rows, nodes[i]->cols);
		m.push_back(tmp);
		v.push_back(tmp);

		mHat.push_back(Eigen::Array<T, 0, 0>());
		vHat.push_back(Eigen::Array<T, 0, 0>());
	}
}



template <typename T>
void ts::AdamOptimizer<T>::computeIncrement(
	std::vector< Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> >& derivatives,
	std::vector<ts::GaElement<T>>& elements
) {
	// gradAccumulator will only be used to convert the indices of optimizable
	// tensors from the gradient / wList indices system (iGrad) to the
	// gradAccumulator indices system (iAcc)

	unsigned iGrad;
	for(unsigned iAcc=0; iAcc<elements.size(); iAcc++) {

		// Get index in the gradAccumulator system
		iGrad = elements[iAcc].index;

		// Compute biased moment estimates
		m[iAcc] = beta1 * m[iAcc] + (1-beta1) * derivatives[iGrad];
		v[iAcc] = beta2 * v[iAcc] + (1-beta2) * derivatives[iGrad] * derivatives[iGrad];

		// Compute bias-corrected moment estimates
		mHat[iAcc] = m[iAcc] / (1 - decayedBeta1);
		vHat[iAcc] = v[iAcc] / (1 - decayedBeta2);

		// Replace gradient with its corrected value
		// (since gradient is used in the gradAccumulator increment method)
		derivatives[iGrad] = mHat[iAcc] / (vHat[iAcc].sqrt() + epsilon);
	}
}



template <typename T>
std::vector<std::vector<std::vector< T >>> ts::AdamOptimizer<T>::run(
	ts::Model<T> &model, std::vector<std::vector< ts::TrainingData<T> >> &batches
) {

		// Set up gradient accumulator (this also resets wList)
	this->gradAccumulator = ts::GradientAccumulator<T>(model);


		//Init parameters

	initMomentEstimates(model.wList.nodes);
	decayedBeta1 = beta1;
	decayedBeta2 = beta2;


		// Start running and training the model

		std::vector<std::vector<std::vector< T >>> losses(this->epochs, (std::vector<std::vector<T>>) {});

	// Epochs
	for(unsigned i=0; i<this->epochs; i++) {

		std::cout << "Epoch " << i + 1 << "/" <<  this->epochs << ": " << std::endl;

		losses[i] = std::vector<std::vector<T>>(batches.size(), (std::vector<T>) {});

		// Batches
		for(unsigned j=0; j<batches.size(); j++) {

			losses[i][j] = std::vector<T>(batches[j].size(), 0);

			// Data instance
			for(unsigned k=0; k<batches[j].size(); k++) {

				ts::Tensor<T> input = ts::Tensor<T>(
					batches[j][k].input, &(model.wList)
				);
				ts::Tensor<T> expected = ts::Tensor<T>(
					batches[j][k].expected, &(model.wList)
				);

				// Compute model and norm
				ts::Tensor<T> output = model.compute(input);
				ts::Tensor<T> norm = (*this->normFunction)(output - expected);

				// Get & correct gradient, then increment gradient accumulator
				ts::Gradient<T> gradient = norm.grad();
				computeIncrement(
					gradient.derivatives,
					this->gradAccumulator.elements
				);
				this->gradAccumulator.increment(gradient);

				model.wList.reset();

				losses[i][j][k] = norm.getValue()(0, 0);
			}

			updateModel(model, batches[j].size());
			this->gradAccumulator.reset();

			// Decay betas
			decayedBeta1 = decayedBeta1 * beta1;
			decayedBeta2 = decayedBeta2 * beta2;

			ts::progressBar(j + 1, batches.size());
		}
		std::cout << std::endl << std::endl;
	}


		// Clean

	this->gradAccumulator.clear();
	model.wList.reset();

	m = {};
	v = {};
	mHat = {};
	vHat = {};


	return losses;
}

// Serializes a ts::Tensor in a string.
// The output string has the following format:
// *ROWS*
// *COLS*
// *VAL*,*VAL*, *VAL*, ..., *VAL*
template <typename T>
std::string ts::serializeTensor(ts::Tensor<T> &tensor) {

	// Save array to string
	std::ostringstream stringStream;
	stringStream << tensor.getValue() << std::endl;

	std::string arrayString =  stringStream.str();


	// Because some extra spaces can be added by Eigen, we remove consecutive
	// spaces/linebreaks
	for(int i=arrayString.size()-2; i>=0; i--) {
		if(
			(arrayString[i] == ' ' || arrayString[i] == '\n') &&
			(arrayString[i+1] == ' ' || arrayString[i+1] == '\n')
		) {
			arrayString.erase(i, 1);
		}
	}

	// We might have a trailing space remaining at index 0 ...
	if(arrayString[0] == ' ') {
		arrayString.erase(0, 1);
	}

	// ... as well as a \n at last position
	if(arrayString[arrayString.size()-1] == '\n') {
		arrayString.erase(arrayString.size()-1, 1);
	}


	// Finally, remove spaces and linebreaks
	std::replace(arrayString.begin(), arrayString.end(), ' ', ',');
	std::replace(arrayString.begin(), arrayString.end(), '\n', ',');


	// Create out stream, return it as string
	std::ostringstream outStream;

	outStream << tensor.getValue().rows() << std::endl;
	outStream << tensor.getValue().cols() << std::endl;
	outStream << arrayString;

	return outStream.str();
}



// Reads an ifstream starting at a serialized tensor, and parses it to a
// ts::Tensor
template <typename T>
ts::Tensor<T> ts::parseTensor(
	std::ifstream &in, ts::WengertList<T> * wList
) {

	std::string line;

	// Get rows
	std::getline(in, line);
	unsigned rows = std::stoi(line);

	// Get cols
	std::getline(in, line);
	unsigned cols = std::stoi(line);

	// Get elements vector
	std::getline(in, line);
	std::vector<std::string> stringElements = ts::split(line, ',');


	// Initialize Eigen::Array
	Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> array;
	array.resize(rows, cols);

	for(unsigned i=0; i<rows; i++) {
		for(unsigned j=0; j<cols; j++) {
			array(i, j) = std::stof(stringElements[i * cols + j]);
		}
	}


	return ts::Tensor<T>(array, wList);
}



// Serializes a vector of ts::Tensor using the ts::serializeTensor
// method. The outputted string has the following format :
// *N == VECTOR SIZE*
// *TENSOR 1 (3 lines)*
// *TENSOR 2 (3 lines)*
// ...
// *TENSOR N (3 lines)*
template <typename T>
std::string ts::serializeTensorsVector(
	std::vector<ts::Tensor<T>> &tensorsVector
) {

	std::ostringstream stringStream;

	stringStream << tensorsVector.size() << std::endl;

	for(unsigned i=0; i<tensorsVector.size(); i++) {
		stringStream << ts::serializeTensor(tensorsVector[i]) << std::endl;
	}

	return stringStream.str();
}



// Reads an ifstream starting at a serialized tensors vector, and parses it to
// a std::vector<ts::Tensor>
template <typename T>
std::vector<ts::Tensor<T>> ts::parseTensorsVector(
	std::ifstream &in, ts::WengertList<T> * wList
) {
	std::vector<ts::Tensor<T>> vector = {};

	std::string line;

	// Get tensor size
	std::getline(in, line);
	unsigned size = std::stoi(line);

	for(unsigned i=0; i<size; i++) {
		vector.push_back(ts::parseTensor(in, wList));
	}

	return vector;
}

// Splits a string by a char delimiter and returns substrings in a vector
std::vector<std::string> ts::split(std::string str, char delimeter){

	std::stringstream ss(str);
	std::string item;
	std::vector<std::string> splitted;

	while (std::getline(ss, item, delimeter)) {
		splitted.push_back(item);
	}

	return splitted;
}



// Read / write 2D unsigned vector

std::string ts::serializeUnsignedVec2D(
	std::vector<std::vector<unsigned>> &vec2d
) {
	std::ostringstream stringStream;

	stringStream << vec2d.size() << std::endl;

	for(unsigned i=0; i<vec2d.size(); i++) {
		stringStream << vec2d[i].size() << std::endl;
		for(unsigned j=0; j<vec2d[i].size(); j++) {
			stringStream << vec2d[i][j] << std::endl;
		}
	}

	return stringStream.str();
}



std::vector<std::vector<unsigned>> ts::parseUnsignedVec2D(
	std::ifstream &in
) {

	std::vector<std::vector<unsigned>> vec2d = {};

	std::string line;

	// Get vec2dsize size
	std::getline(in, line);
	unsigned size2d = std::stoi(line);

	for(unsigned i=0; i<size2d; i++) {
		vec2d.push_back({});

		std::getline(in, line);
		unsigned size1d = std::stoi(line);

		for(unsigned j=0; j<size1d; j++) {
			std::getline(in, line);
			unsigned val = std::stoi(line);
			vec2d[i].push_back(val);
		}
	}

	return vec2d;
}



void ts::progressBar(unsigned current, unsigned max) {
	float progress = (float) current / (float) max;

	std::cout << "\r[";
	int pos = BARWIDTH * progress;

	for (int i = 0; i < BARWIDTH; ++i) {
	    if (i < pos) {
			std::cout << "=";
		}
	    else if (i == pos) {
			std::cout << ">";
		}
	    else {
			std::cout << " ";
		}
	}
	std::cout << "] " << (int) (progress * 100) << "% ";

	std::cout << "(" << current << "/" << max << " batches)";

	std::cout.flush();

}




// float

template class ts::Node<float>;
template class ts::InputNode<float>;
template class ts::ElementWiseNode<float>;
template class ts::MatProdNode<float>;
template class ts::ScalarNode<float>;

template class ts::WengertList<float>;
template class ts::Tensor<float>;
template class ts::Gradient<float>;

template ts::Tensor<float> ts::NewTensor(
	Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic> newValue,
	ts::WengertList<float> * newWList
);

template ts::Tensor<float> ts::operator+(const ts::Tensor<float> &x, const ts::Tensor<float> &y);
template ts::Tensor<float> ts::operator-(const ts::Tensor<float> &x, const ts::Tensor<float> &y);
template ts::Tensor<float> ts::operator*(const ts::Tensor<float> &x, const ts::Tensor<float> &y);
template ts::Tensor<float> ts::operator/(const ts::Tensor<float> &x, const ts::Tensor<float> &y);

template ts::Tensor<float> ts::matProd(const ts::Tensor<float> &x, const ts::Tensor<float> &y);
template ts::Tensor<float> ts::sigmoid(const ts::Tensor<float> &x);
template ts::Tensor<float> ts::relu(const ts::Tensor<float> &x);
template ts::Tensor<float> ts::leakyRelu(const ts::Tensor<float> &x);
template ts::Tensor<float> ts::rescale(const ts::Tensor<float> &x);
template ts::Tensor<float> ts::squaredNorm(const ts::Tensor<float> &x);


template class ts::Model<float>;
template class ts::Polynom<float>;
template class ts::MultiLayerPerceptron<float>;
template class ts::ConvolutionalNetwork<float>;

template class ts::TrainingData<float>;
template class ts::GaElement<float>;
template class ts::GradientAccumulator<float>;
template class ts::Optimizer<float>;
template class ts::GradientDescentOptimizer<float>;
template class ts::AdamOptimizer<float>;

template std::string ts::serializeTensor(ts::Tensor<float> &tensor);
template ts::Tensor<float> ts::parseTensor(
	std::ifstream &in, ts::WengertList<float> * wList
);
template std::string ts::serializeTensorsVector(
	std::vector<ts::Tensor<float>> &tensorsVector
);
template std::vector<ts::Tensor<float>> ts::parseTensorsVector(
	std::ifstream &in, ts::WengertList<float> * wList
);

template Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic> ts::convArray(
	const Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic> &mat,
	const Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic> &ker
);
template class ts::ConvolutionNode<float>;
template ts::Tensor<float> ts::convolution(
	const ts::Tensor<float> &mat, const ts::Tensor<float> &ker
);
template class ts::PoolingNode<float>;
template ts::Tensor<float> ts::maxPooling(
	const ts::Tensor<float> &x, std::vector<unsigned> pool
);
template class ts::SplitNode<float>;
template std::vector<ts::Tensor<float>> ts::split(
	const ts::Tensor<float> &x, ChannelSplit channelSplit, unsigned nInputChannels
);
template class ts::VertCatNode<float>;
template ts::Tensor<float> ts::vertCat<float>(
	const std::vector<ts::Tensor<float>> &x
);
template class ts::FlatteningNode<float>;
template ts::Tensor<float> ts::flattening<float>(const ts::Tensor<float> &x);
template class ts::Im2ColNode<float>;
template ts::Tensor<float> ts::im2col<float>(
	const std::vector<ts::Tensor<float>> &x,
	std::vector<unsigned> kernelDim
);
template class ts::Col2ImNode<float>;
template std::vector<ts::Tensor<float>> ts::col2im<float>(
	const ts::Tensor<float> &x,
	std::vector<unsigned> outputDim
);



	// double

template class ts::Node<double>;
template class ts::InputNode<double>;
template class ts::ElementWiseNode<double>;
template class ts::MatProdNode<double>;
template class ts::ScalarNode<double>;

template class ts::WengertList<double>;
template class ts::Tensor<double>;
template class ts::Gradient<double>;

template ts::Tensor<double> ts::NewTensor(
	Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic> newValue,
	ts::WengertList<double> * newWList
);

template ts::Tensor<double> ts::operator+(const ts::Tensor<double> &x, const ts::Tensor<double> &y);
template ts::Tensor<double> ts::operator-(const ts::Tensor<double> &x, const ts::Tensor<double> &y);
template ts::Tensor<double> ts::operator*(const ts::Tensor<double> &x, const ts::Tensor<double> &y);
template ts::Tensor<double> ts::operator/(const ts::Tensor<double> &x, const ts::Tensor<double> &y);

template ts::Tensor<double> ts::matProd(const ts::Tensor<double> &x, const ts::Tensor<double> &y);
template ts::Tensor<double> ts::sigmoid(const ts::Tensor<double> &x);
template ts::Tensor<double> ts::relu(const ts::Tensor<double> &x);
template ts::Tensor<double> ts::leakyRelu(const ts::Tensor<double> &x);
template ts::Tensor<double> ts::rescale(const ts::Tensor<double> &x);
template ts::Tensor<double> ts::squaredNorm(const ts::Tensor<double> &x);

template class ts::Model<double>;
template class ts::Polynom<double>;
template class ts::MultiLayerPerceptron<double>;
template class ts::ConvolutionalNetwork<double>;

template class ts::TrainingData<double>;
template class ts::GaElement<double>;
template class ts::GradientAccumulator<double>;
template class ts::Optimizer<double>;
template class ts::GradientDescentOptimizer<double>;
template class ts::AdamOptimizer<double>;

template std::string ts::serializeTensor(ts::Tensor<double> &tensor);
template ts::Tensor<double> ts::parseTensor(
	std::ifstream &in, ts::WengertList<double> * wList
);
template std::string ts::serializeTensorsVector(
	std::vector<ts::Tensor<double>> &tensorsVector
);
template std::vector<ts::Tensor<double>> ts::parseTensorsVector(
	std::ifstream &in, ts::WengertList<double> * wList
);

template Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic> ts::convArray(
	const Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic> &mat,
	const Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic> &ker
);
template class ts::ConvolutionNode<double>;
template ts::Tensor<double> ts::convolution(
	const ts::Tensor<double> &mat, const ts::Tensor<double> &ker
);
template class ts::PoolingNode<double>;
template ts::Tensor<double> ts::maxPooling(
	const ts::Tensor<double> &x, std::vector<unsigned> pool
);
template class ts::SplitNode<double>;
template std::vector<ts::Tensor<double>> ts::split(
	const ts::Tensor<double> &x, ChannelSplit channelSplit, unsigned nInputChannels
);
template class ts::VertCatNode<double>;
template ts::Tensor<double> ts::vertCat<double>(
	const std::vector<ts::Tensor<double>> &x
);
template class ts::FlatteningNode<double>;
template ts::Tensor<double> ts::flattening<double>(const ts::Tensor<double> &x);
template class ts::Im2ColNode<double>;
template ts::Tensor<double> ts::im2col<double>(
	const std::vector<ts::Tensor<double>> &x,
	std::vector<unsigned> kernelDim
);
template class ts::Col2ImNode<double>;
template std::vector<ts::Tensor<double>> ts::col2im<double>(
	const ts::Tensor<double> &x,
	std::vector<unsigned> outputDim
);