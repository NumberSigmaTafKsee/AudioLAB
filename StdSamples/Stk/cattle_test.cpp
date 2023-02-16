#include "Cattle.hpp"
#include <gtest/gtest.h>

namespace cattle {

/**
 * A namespace for C-ATTL3's test functions.
 */
namespace test {

bool verbose;

/**
 * A trait struct for the name of a scalar type and the default numeric constants used for gradient
 * verification depending on the scalar type.
 */
template<typename Scalar>
struct ScalarTraits {
	static constexpr Scalar step_size = 1e-5;
	static constexpr Scalar abs_epsilon = 1e-2;
	static constexpr Scalar rel_epsilon = 1e-2;
	inline static std::string name() {
		return "double";
	}
};

/**
 * Template specialization for single precision floating point scalars.
 */
template<>
struct ScalarTraits<float> {
	static constexpr float step_size = 5e-4;
	static constexpr float abs_epsilon = 1.5e-1;
	static constexpr float rel_epsilon = 1e-1;
	inline static std::string name() {
		return "float";
	}
};

/**
 * An alias for a unique pointer to an optimizer.
 */
template<typename Scalar, std::size_t Rank, bool Sequential>
using OptimizerPtr = std::unique_ptr<Optimizer<Scalar,Rank,Sequential>>;

/**
 * @param dims The dimensions of the prospective tensor.
 * @return The number of samples.
 */
template<typename Scalar, std::size_t Rank, bool Sequential>
inline std::size_t get_rows(const typename std::enable_if<!Sequential,
		std::array<std::size_t,Rank>>::type& dims) {
	return dims[0];
}

/**
 * @param dims The dimensions of the prospective tensor.
 * @return The number of samples multiplied by the number of time steps.
 */
template<typename Scalar, std::size_t Rank, bool Sequential>
inline std::size_t get_rows(const typename std::enable_if<Sequential,
		std::array<std::size_t,Rank>>::type& dims) {
	return dims[0] * dims[1];
}

/**
 * @param dims The dimensions of the random tensor to create.
 * @return A tensor of the specified dimensions filled with random values in the range of
 * -1 to 1.
 */
template<typename Scalar, std::size_t Rank>
inline TensorPtr<Scalar,Rank> random_tensor(const std::array<std::size_t,Rank>& dims) {
	TensorPtr<Scalar,Rank> tensor_ptr(new Tensor<Scalar,Rank>(dims));
	tensor_ptr->setRandom();
	return tensor_ptr;
}

/**
 * @param dims The dimensions of the random tensor to create.
 * @param non_label The value of the non-label elements of the objective tensor (0 for cross
 * entropy, -1 for hinge loss).
 * @return A one-hot tensor of the specified dimensions in which only one element per sample
 * and time step is 1 while all others equal non_label.
 */
template<typename Scalar, std::size_t Rank, bool Sequential>
inline TensorPtr<Scalar,Rank> random_one_hot_tensor(const std::array<std::size_t,Rank>& dims,
		Scalar non_label = 0) {
	TensorPtr<Scalar,Rank> tensor_ptr(new Tensor<Scalar,Rank>(dims));
	tensor_ptr->setConstant(non_label);
	int rows = get_rows<Scalar,Rank,Sequential>(dims);
	MatrixMap<Scalar> mat_map(tensor_ptr->data(), rows, tensor_ptr->size() / rows);
	std::random_device rd;
	std::mt19937 rng(rd());
	std::uniform_int_distribution<int> dist(0, mat_map.cols() - 1);
	for (int i = 0; i < mat_map.rows(); ++i)
		mat_map(i,dist(rng)) = 1;
	return tensor_ptr;
}

/**
 * @param dims The dimensions of the random tensor to create.
 * @param non_label The value of the non-label elements of the objective tensor (0 for cross
 * entropy, -1 for hinge loss).
 * @return A multi-label objective tensor whose elements equal either 1 or non_label.
 */
template<typename Scalar, std::size_t Rank>
inline TensorPtr<Scalar,Rank> random_multi_hot_tensor(const std::array<std::size_t,Rank>& dims,
		Scalar non_label = 0) {
	auto tensor_ptr = random_tensor<Scalar,Rank>(dims);
	Tensor<bool,Rank> if_tensor = (*tensor_ptr) > tensor_ptr->constant((Scalar) 0);
	Tensor<Scalar,Rank> then_tensor = tensor_ptr->constant((Scalar) 1);
	Tensor<Scalar,Rank> else_tensor = tensor_ptr->constant(non_label);
	*tensor_ptr = if_tensor.select(then_tensor, else_tensor);
	return tensor_ptr;
}

/**
 * @param input_dims The input dimensions of the layer.
 * @return A fully connected kernel layer.
 */
template<typename Scalar, std::size_t Rank>
inline KernelPtr<Scalar,Rank> kernel_layer(const typename std::enable_if<Rank != 1,
		Dimensions<std::size_t,Rank>>::type& input_dims) {
	return KernelPtr<Scalar,Rank>(new ConvKernelLayer<Scalar,Rank>(
			input_dims, input_dims.template extend<3 - Rank>()(2),
			std::make_shared<HeParameterInitialization<Scalar>>()));
}

/**
 * @param input_dims The input dimensions of the layer.
 * @return A convolutional kernel layer.
 */
template<typename Scalar, std::size_t Rank>
inline KernelPtr<Scalar,Rank> kernel_layer(const typename std::enable_if<Rank == 1,
		Dimensions<std::size_t,Rank>>::type& input_dims) {
	return KernelPtr<Scalar,Rank>(new DenseKernelLayer<Scalar,Rank>(
			input_dims, input_dims.get_volume(),
			std::make_shared<GlorotParameterInitialization<Scalar>>()));
}

/**
 * @param input_dims The input dimensions of the neural network.
 * @return A simple feed-forward neural network with unrestricted output values.
 */
template<typename Scalar, std::size_t Rank>
inline NeuralNetPtr<Scalar,Rank,false> neural_net(const Dimensions<std::size_t,Rank>& input_dims) {
	std::vector<LayerPtr<Scalar,Rank>> layers(1);
	layers[0] = kernel_layer<Scalar,Rank>(input_dims);
	return NeuralNetPtr<Scalar,Rank,false>(new FeedforwardNeuralNetwork<Scalar,Rank>(std::move(layers)));
}

/**
 * @param input_dims The input dimensions of the neural network.
 * @return A simple feed-forward neural network with a sigmoid activation function and hence output values
 * between 0 and 1.
 */
template<typename Scalar, std::size_t Rank>
inline NeuralNetPtr<Scalar,Rank,false> sigmoid_neural_net(const Dimensions<std::size_t,Rank>& input_dims) {
	std::vector<LayerPtr<Scalar,Rank>> layers(2);
	layers[0] = kernel_layer<Scalar,Rank>(input_dims);
	layers[1] = LayerPtr<Scalar,Rank>(new SigmoidActivationLayer<Scalar,Rank>(layers[0]->get_output_dims()));
	return NeuralNetPtr<Scalar,Rank,false>(new FeedforwardNeuralNetwork<Scalar,Rank>(std::move(layers)));
}

/**
 * @param input_dims The input dimensions of the neural network.
 * @return A simple feed-forward neural network with a tanh activation function and hence output values
 * between -1 and 1.
 */
template<typename Scalar, std::size_t Rank>
inline NeuralNetPtr<Scalar,Rank,false> tanh_neural_net(const Dimensions<std::size_t,Rank>& input_dims) {
	std::vector<LayerPtr<Scalar,Rank>> layers(2);
	layers[0] = kernel_layer<Scalar,Rank>(input_dims);
	layers[1] = LayerPtr<Scalar,Rank>(new TanhActivationLayer<Scalar,Rank>(layers[0]->get_output_dims()));
	return NeuralNetPtr<Scalar,Rank,false>(new FeedforwardNeuralNetwork<Scalar,Rank>(std::move(layers)));
}

/**
 * @param input_dims The input dimensions of the neural network.
 * @return A simple feed-forward neural network with a softmax activation function and hence output values
 * between 0 and 1.
 */
template<typename Scalar, std::size_t Rank>
inline NeuralNetPtr<Scalar,Rank,false> softmax_neural_net(const Dimensions<std::size_t,Rank>& input_dims) {
	std::vector<LayerPtr<Scalar,Rank>> layers(2);
	layers[0] = kernel_layer<Scalar,Rank>(input_dims);
	layers[1] = LayerPtr<Scalar,Rank>(new SoftmaxActivationLayer<Scalar,Rank>(layers[0]->get_output_dims()));
	return NeuralNetPtr<Scalar,Rank,false>(new FeedforwardNeuralNetwork<Scalar,Rank>(std::move(layers)));
}

/**
 * @param input_dims The input dimensions of the neural network.
 * @return A simple recurrent neural network without an identity output activation function and with a
 * single output time step.
 */
template<typename Scalar, std::size_t Rank>
inline NeuralNetPtr<Scalar,Rank,true> recurrent_neural_net(const Dimensions<std::size_t,Rank>& input_dims,
		std::function<std::pair<std::size_t,std::size_t>(std::size_t)> seq_schedule_func) {
	return NeuralNetPtr<Scalar,Rank,true>(new RecurrentNeuralNetwork<Scalar,Rank>(kernel_layer<Scalar,Rank>(input_dims),
			kernel_layer<Scalar,Rank>(input_dims), kernel_layer<Scalar,Rank>(input_dims),
			ActivationPtr<Scalar,Rank>(new TanhActivationLayer<Scalar,Rank>(input_dims)),
			ActivationPtr<Scalar,Rank>(new IdentityActivationLayer<Scalar,Rank>(input_dims)), seq_schedule_func));
}

/**
 * @param input_dims The input dimensions of the neural network.
 * @return A simple single-output recurrent neural network without an identity output activation function
 * and with a single output time step.
 */
template<typename Scalar, std::size_t Rank>
inline NeuralNetPtr<Scalar,Rank,true> single_output_recurrent_neural_net(const Dimensions<std::size_t,Rank>& input_dims,
		std::function<std::pair<std::size_t,std::size_t>(std::size_t)> seq_schedule_func) {
	auto init = std::make_shared<GlorotParameterInitialization<Scalar>>();
	return NeuralNetPtr<Scalar,Rank,true>(new RecurrentNeuralNetwork<Scalar,Rank>(kernel_layer<Scalar,Rank>(input_dims),
			kernel_layer<Scalar,Rank>(input_dims), KernelPtr<Scalar,Rank>(new DenseKernelLayer<Scalar,Rank>(input_dims, 1, init)),
			ActivationPtr<Scalar,Rank>(new TanhActivationLayer<Scalar,Rank>(input_dims)),
			ActivationPtr<Scalar,Rank>(new IdentityActivationLayer<Scalar,Rank>(Dimensions<std::size_t,Rank>())),
			seq_schedule_func));
}

/**
 * @param net The first module of the composite single-output network.
 * @return A stacked network with an output dimensionality of one.
 */
template<typename Scalar, std::size_t Rank>
inline NeuralNetPtr<Scalar,Rank,false> single_output_net(NeuralNetPtr<Scalar,Rank,false> net) {
	auto init = std::make_shared<GlorotParameterInitialization<Scalar>>();
	std::vector<NeuralNetPtr<Scalar,Rank,false>> modules(2);
	modules[0] = std::move(net);
	modules[1] = NeuralNetPtr<Scalar,Rank,false>(new FeedforwardNeuralNetwork<Scalar,Rank>(
			LayerPtr<Scalar,Rank>(new DenseKernelLayer<Scalar,Rank>(modules[0]->get_output_dims(), 1, init))));
	return NeuralNetPtr<Scalar,Rank,false>(new StackedNeuralNetwork<Scalar,Rank,false>(std::move(modules)));
}

/**
 * @param test_case_name The name of the test case/suite.
 * @param test_name The name of the test.
 */
template<typename Scalar, std::size_t Rank, bool Sequential>
inline void print_test_header(std::string test_case_name, std::string test_name) {
	std::transform(test_case_name.begin(), test_case_name.end(), test_case_name.begin(), ::toupper);
	std::transform(test_name.begin(), test_name.end(), test_name.begin(), ::toupper);
	std::string header = "|   " + test_case_name + ": " + test_name + "; SCALAR TYPE: " +
			ScalarTraits<Scalar>::name() + "; RANK: " + std::to_string(Rank) +
			"; SEQ: " + std::to_string(Sequential) + "   |";
	std::size_t header_length = header.length();
	int padding_content_length = header_length - 2;
	std::string header_border = " " + std::string(padding_content_length, '-') + " ";
	std::string upper_header_padding = "/" + std::string(padding_content_length, ' ') + "\\";
	std::string lower_header_padding = "\\" + std::string(padding_content_length, ' ') + "/";
	std::cout << std::endl << header_border << std::endl << upper_header_padding << std::endl <<
			header << std::endl << lower_header_padding << std::endl << header_border << std::endl;
}

} /* namespace test */

} /* namespace cattle */


namespace cattle {
namespace test {

/**
 * Determines the verbosity of the training tests.
 */
extern bool verbose;

/**
 * @param name The name of the training test.
 * @param train_prov The training data provider.
 * @param net The neural network to train.
 * @param opt The optimizer to use to train the network.
 * @param epochs The number of epochs for whicht the network is to be trained.
 * @param epsilon The maximum acceptable absolute difference between 0 and
 * the final training loss.
 */
template<typename Scalar, std::size_t Rank, bool Sequential>
inline void train_test(std::string name, DataProvider<Scalar,Rank,Sequential>& train_prov,
		NeuralNetwork<Scalar,Rank,Sequential>& net, Optimizer<Scalar,Rank,Sequential>& opt,
		unsigned epochs, Scalar epsilon) {
	print_test_header<Scalar,Rank,Sequential>("training test", name);
	net.init();
	opt.fit(net);
	Scalar loss = opt.train(net, train_prov, epochs, 0, epsilon, verbose);
	EXPECT_TRUE(NumericUtils<Scalar>::almost_equal((Scalar) 0, loss, epsilon));
}

/**
 * @param name The name of the training test.
 * @param opt The optimizer to use to train the feed-forward network.
 * @param input_dims The input dimensions of the network.
 * @param samples The total number of data samples to include in the training data.
 * @param epochs The number of epochs for which the network is to be trained.
 * @param epsilon The maximum acceptable absolute difference between 0 and
 * the final training loss.
 */
template<typename Scalar, std::size_t Rank>
inline void ffnn_train_test(std::string name, OptimizerPtr<Scalar,Rank,false> opt,
		const Dimensions<std::size_t,Rank> input_dims, unsigned samples,
		unsigned epochs, Scalar epsilon = ScalarTraits<Scalar>::abs_epsilon) {
	NeuralNetPtr<Scalar,Rank,false> net = single_output_net<Scalar,Rank>(tanh_neural_net<Scalar,Rank>(input_dims));
	std::array<std::size_t,Rank + 1> input_batch_dims = input_dims.template promote<>();
	input_batch_dims[0] = samples;
	TensorPtr<Scalar,Rank + 1> input_tensor(new Tensor<Scalar,Rank + 1>(input_batch_dims));
	input_tensor->setRandom();
	std::array<std::size_t,Rank + 1> output_batch_dims = net->get_output_dims().template promote<>();
	output_batch_dims[0] = samples;
	TensorPtr<Scalar,Rank + 1> output_tensor(new Tensor<Scalar,Rank + 1>(output_batch_dims));
	output_tensor->setRandom();
	MemoryDataProvider<Scalar,Rank,false> prov(std::move(input_tensor), std::move(output_tensor));
	train_test(name, prov, *net, *opt, epochs, epsilon);
}

/**
 * @param name The name of the training test.
 * @param opt The optimizer to use to train the network.
 * @param input_dims The input dimensions of the network.
 * @param samples The total number of data samples to include in the training data.
 * @param time_steps The number of time steps each sample is to contain.
 * @param epochs The number of epochs for which the network is to be trained.
 * @param epsilon The maximum acceptable absolute difference between 0 and
 * the final training loss.
 */
template<typename Scalar, std::size_t Rank>
inline void rnn_train_test(std::string name, OptimizerPtr<Scalar,Rank,true> opt,
		const Dimensions<std::size_t,Rank> input_dims, unsigned samples, unsigned time_steps,
		unsigned epochs, Scalar epsilon = ScalarTraits<Scalar>::abs_epsilon) {
	NeuralNetPtr<Scalar,Rank,true> net = single_output_recurrent_neural_net<Scalar,Rank>(input_dims,
			[](int input_seq_length) { return std::make_pair(1, input_seq_length - 1); });
	std::array<std::size_t,Rank + 2> input_batch_dims = input_dims.template promote<2>();
	input_batch_dims[0] = samples;
	input_batch_dims[1] = time_steps;
	TensorPtr<Scalar,Rank + 2> input_tensor(new Tensor<Scalar,Rank + 2>(input_batch_dims));
	input_tensor->setRandom();
	std::array<std::size_t,Rank + 2> output_batch_dims = net->get_output_dims().template promote<2>();
	output_batch_dims[0] = samples;
	output_batch_dims[1] = 1;
	TensorPtr<Scalar,Rank + 2> output_tensor(new Tensor<Scalar,Rank + 2>(output_batch_dims));
	output_tensor->setRandom();
	MemoryDataProvider<Scalar,Rank,true> prov(std::move(input_tensor), std::move(output_tensor));
	train_test(name, prov, *net, *opt, epochs, epsilon);
}

/**
 * Performs training tests using a vanilla stochastic gradient descent optimizer.
 */
template<typename Scalar>
inline void vanilla_sgd_train_test() {
	const unsigned samples = 20;
	const unsigned time_steps = 3;
	const unsigned epochs = 500;
	const Scalar epsilon = .2;
	const Scalar seq_epsilon = .25;
	const Dimensions<std::size_t,3> dims({ 6u, 6u, 2u });
	auto loss = std::make_shared<SquaredLoss<Scalar,3,false>>();
	ffnn_train_test<Scalar,3>("vanilla sgd batch",
			OptimizerPtr<Scalar,3,false>(new VanillaSGDOptimizer<Scalar,3,false>(loss, samples, 2e-3)),
			dims, samples, epochs, epsilon);
	ffnn_train_test<Scalar,3>("vanilla sgd mini-batch",
			OptimizerPtr<Scalar,3,false>(new VanillaSGDOptimizer<Scalar,3,false>(loss, 10, 5e-4)),
			dims, samples, epochs, epsilon);
	ffnn_train_test<Scalar,3>("vanilla sgd online",
			OptimizerPtr<Scalar,3,false>(new VanillaSGDOptimizer<Scalar,3,false>(loss, 1, 1e-4)),
			dims, samples, epochs, epsilon);
	auto seq_loss = std::make_shared<SquaredLoss<Scalar,3,true>>();
	rnn_train_test<Scalar,3>("vanilla sgd batch",
			OptimizerPtr<Scalar,3,true>(new VanillaSGDOptimizer<Scalar,3,true>(seq_loss, samples, 2e-3)),
			dims, samples, time_steps, epochs, seq_epsilon);
	rnn_train_test<Scalar,3>("vanilla sgd mini-batch",
			OptimizerPtr<Scalar,3,true>(new VanillaSGDOptimizer<Scalar,3,true>(seq_loss, 10, 5e-4)),
			dims, samples, time_steps, epochs, seq_epsilon);
	rnn_train_test<Scalar,3>("vanilla sgd online",
			OptimizerPtr<Scalar,3,true>(new VanillaSGDOptimizer<Scalar,3,true>(seq_loss, 1, 1e-4)),
			dims, samples, time_steps, epochs, seq_epsilon);
}

TEST(TrainingTest, VanillaSGD) {
	vanilla_sgd_train_test<float>();
	vanilla_sgd_train_test<double>();
}

/**
 * Performs training tests using a momentum accelerated stochastic gradient descent optimizer.
 */
template<typename Scalar>
inline void momentum_sgd_train_test() {
	const unsigned samples = 20;
	const unsigned time_steps = 3;
	const unsigned epochs = 500;
	const Scalar epsilon = .18;
	const Scalar seq_epsilon = .22;
	const Dimensions<std::size_t,3> dims({ 6u, 6u, 2u });
	auto loss = std::make_shared<SquaredLoss<Scalar,3,false>>();
	ffnn_train_test<Scalar,3>("momentum accelerated sgd batch",
			OptimizerPtr<Scalar,3,false>(new MomentumSGDOptimizer<Scalar,3,false>(loss, samples, 2e-3)),
			dims, samples, epochs, epsilon);
	ffnn_train_test<Scalar,3>("momentum accelerated sgd mini-batch",
			OptimizerPtr<Scalar,3,false>(new MomentumSGDOptimizer<Scalar,3,false>(loss, 10, 1e-3, 1e-2, .99)),
			dims, samples, epochs, epsilon);
	ffnn_train_test<Scalar,3>("momentum accelerated sgd online",
			OptimizerPtr<Scalar,3,false>(new MomentumSGDOptimizer<Scalar,3,false>(loss, 1, 5e-4)),
			dims, samples, epochs, epsilon);
	auto seq_loss = std::make_shared<SquaredLoss<Scalar,3,true>>();
	rnn_train_test<Scalar,3>("momentum accelerated sgd batch",
			OptimizerPtr<Scalar,3,true>(new MomentumSGDOptimizer<Scalar,3,true>(seq_loss, samples, 2e-3)),
			dims, samples, time_steps, epochs, seq_epsilon);
	rnn_train_test<Scalar,3>("momentum accelerated sgd mini-batch",
			OptimizerPtr<Scalar,3,true>(new MomentumSGDOptimizer<Scalar,3,true>(seq_loss, 10, 1e-3, 1e-2, .99)),
			dims, samples, time_steps, epochs, seq_epsilon);
	rnn_train_test<Scalar,3>("momentum accelerated sgd online",
			OptimizerPtr<Scalar,3,true>(new MomentumSGDOptimizer<Scalar,3,true>(seq_loss, 1, 5e-4)),
			dims, samples, time_steps, epochs, seq_epsilon);
}

TEST(TrainingTest, MomentumSGD) {
	momentum_sgd_train_test<float>();
	momentum_sgd_train_test<double>();
}

/**
 * Performs training tests using a nesterov momentum accelerated stochastic gradient descent
 * optimizer.
 */
template<typename Scalar>
inline void nesterov_momentum_sgd_train_test() {
	const unsigned samples = 20;
	const unsigned time_steps = 3;
	const unsigned epochs = 500;
	const Scalar epsilon = .18;
	const Scalar seq_epsilon = .22;
	const Dimensions<std::size_t,3> dims({ 6u, 6u, 2u });
	auto loss = std::make_shared<SquaredLoss<Scalar,3,false>>();
	ffnn_train_test<Scalar,3>("nesterov momentum accelerated sgd batch",
			OptimizerPtr<Scalar,3,false>(new NesterovMomentumSGDOptimizer<Scalar,3,false>(loss, samples)),
			dims, samples, epochs, epsilon);
	ffnn_train_test<Scalar,3>("nesterov momentum accelerated sgd mini-batch",
			OptimizerPtr<Scalar,3,false>(new NesterovMomentumSGDOptimizer<Scalar,3,false>(loss, 10, 1e-3, 1e-2, .99)),
			dims, samples, epochs, epsilon);
	ffnn_train_test<Scalar,3>("nesterov momentum accelerated sgd online",
			OptimizerPtr<Scalar,3,false>(new NesterovMomentumSGDOptimizer<Scalar,3,false>(loss, 1, 5e-4)),
			dims, samples, epochs, epsilon);
	auto seq_loss = std::make_shared<SquaredLoss<Scalar,3,true>>();
	rnn_train_test<Scalar,3>("nesterov momentum accelerated sgd batch",
			OptimizerPtr<Scalar,3,true>(new NesterovMomentumSGDOptimizer<Scalar,3,true>(seq_loss, samples)),
			dims, samples, time_steps, epochs, seq_epsilon);
	rnn_train_test<Scalar,3>("nesterov momentum accelerated sgd mini-batch",
			OptimizerPtr<Scalar,3,true>(new NesterovMomentumSGDOptimizer<Scalar,3,true>(seq_loss, 10, 1e-3, 1e-2, .99)),
			dims, samples, time_steps, epochs, seq_epsilon);
	rnn_train_test<Scalar,3>("nesterov momentum accelerated sgd online",
			OptimizerPtr<Scalar,3,true>(new NesterovMomentumSGDOptimizer<Scalar,3,true>(seq_loss, 1, 5e-4)),
			dims, samples, time_steps, epochs, seq_epsilon);
}

TEST(TrainingTest, NesterovMomentumSGD) {
	nesterov_momentum_sgd_train_test<float>();
	nesterov_momentum_sgd_train_test<double>();
}

/**
 * Performs training tests using an AdaGrad optimizer.
 */
template<typename Scalar>
inline void adagrad_train_test() {
	const unsigned samples = 20;
	const unsigned time_steps = 3;
	const unsigned epochs = 500;
	const Scalar epsilon = .03;
	const Scalar seq_epsilon = .08;
	const Dimensions<std::size_t,3> dims({ 6u, 6u, 2u });
	auto loss = std::make_shared<SquaredLoss<Scalar,3,false>>();
	ffnn_train_test<Scalar,3>("adagrad batch",
			OptimizerPtr<Scalar,3,false>(new AdaGradOptimizer<Scalar,3,false>(loss, samples, 2.5e-2)),
			dims, samples, epochs, epsilon);
	ffnn_train_test<Scalar,3>("adagrad mini-batch",
			OptimizerPtr<Scalar,3,false>(new AdaGradOptimizer<Scalar,3,false>(loss, 10, 1e-2)),
			dims, samples, epochs, epsilon);
	ffnn_train_test<Scalar,3>("adagrad online",
			OptimizerPtr<Scalar,3,false>(new AdaGradOptimizer<Scalar,3,false>(loss)),
			dims, samples, epochs, epsilon);
	auto seq_loss = std::make_shared<SquaredLoss<Scalar,3,true>>();
	rnn_train_test<Scalar,3>("adagrad batch",
			OptimizerPtr<Scalar,3,true>(new AdaGradOptimizer<Scalar,3,true>(seq_loss, samples, 1e-2)),
			dims, samples, time_steps, epochs, seq_epsilon);
	rnn_train_test<Scalar,3>("adagrad mini-batch",
			OptimizerPtr<Scalar,3,true>(new AdaGradOptimizer<Scalar,3,true>(seq_loss, 10, 5e-3)),
			dims, samples, time_steps, epochs, seq_epsilon);
	rnn_train_test<Scalar,3>("adagrad online",
			OptimizerPtr<Scalar,3,true>(new AdaGradOptimizer<Scalar,3,true>(seq_loss)),
			dims, samples, time_steps, epochs, seq_epsilon);
}

TEST(TrainingTest, AdaGrad) {
	adagrad_train_test<float>();
	adagrad_train_test<double>();
}

/**
 * Performs training tests using an RMSprop optimizer.
 */
template<typename Scalar>
inline void rmsprop_train_test() {
	const unsigned samples = 20;
	const unsigned time_steps = 3;
	const unsigned epochs = 500;
	const Scalar epsilon = .02;
	const Scalar seq_epsilon = .08;
	const Dimensions<std::size_t,3> dims({ 6u, 6u, 2u });
	auto loss = std::make_shared<SquaredLoss<Scalar,3,false>>();
	ffnn_train_test<Scalar,3>("rmsprop batch",
			OptimizerPtr<Scalar,3,false>(new RMSPropOptimizer<Scalar,3,false>(loss, samples, 1e-2)),
			dims, samples, epochs, epsilon);
	ffnn_train_test<Scalar,3>("rmsprop mini-batch",
			OptimizerPtr<Scalar,3,false>(new RMSPropOptimizer<Scalar,3,false>(loss, 10, 5e-3, 1e-2)),
			dims, samples, epochs, epsilon);
	ffnn_train_test<Scalar,3>("rmsprop online",
			OptimizerPtr<Scalar,3,false>(new RMSPropOptimizer<Scalar,3,false>(loss)),
			dims, samples, epochs, epsilon);
	auto seq_loss = std::make_shared<SquaredLoss<Scalar,3,true>>();
	rnn_train_test<Scalar,3>("rmsprop batch",
			OptimizerPtr<Scalar,3,true>(new RMSPropOptimizer<Scalar,3,true>(seq_loss, samples, 1e-2)),
			dims, samples, time_steps, epochs, seq_epsilon);
	rnn_train_test<Scalar,3>("rmsprop mini-batch",
			OptimizerPtr<Scalar,3,true>(new RMSPropOptimizer<Scalar,3,true>(seq_loss, 10, 5e-3, 1e-2)),
			dims, samples, time_steps, epochs, seq_epsilon);
	rnn_train_test<Scalar,3>("rmsprop online",
			OptimizerPtr<Scalar,3,true>(new RMSPropOptimizer<Scalar,3,true>(seq_loss)),
			dims, samples, time_steps, epochs, seq_epsilon);
}

TEST(TrainingTest, RMSProp) {
	rmsprop_train_test<float>();
	rmsprop_train_test<double>();
}

/**
 * Performs training tests using an AdaDelta optimizer.
 */
template<typename Scalar>
inline void adadelta_train_test() {
	const unsigned samples = 20;
	const unsigned time_steps = 3;
	const unsigned epochs = 500;
	const Scalar epsilon = .02;
	const Scalar seq_epsilon = .1;
	const Dimensions<std::size_t,3> dims({ 6u, 6u, 2u });
	auto loss = std::make_shared<SquaredLoss<Scalar,3,false>>();
	ffnn_train_test<Scalar,3>("adadelta batch",
			OptimizerPtr<Scalar,3,false>(new AdaDeltaOptimizer<Scalar,3,false>(loss, samples)),
			dims, samples, epochs, epsilon);
	ffnn_train_test<Scalar,3>("adadelta mini-batch",
			OptimizerPtr<Scalar,3,false>(new AdaDeltaOptimizer<Scalar,3,false>(loss, 10, 1e-2)),
			dims, samples, epochs, epsilon);
	ffnn_train_test<Scalar,3>("adadelta online",
			OptimizerPtr<Scalar,3,false>(new AdaDeltaOptimizer<Scalar,3,false>(loss, 1, 1e-1)),
			dims, samples, epochs, epsilon);
	auto seq_loss = std::make_shared<SquaredLoss<Scalar,3,true>>();
	rnn_train_test<Scalar,3>("adadelta batch",
			OptimizerPtr<Scalar,3,true>(new AdaDeltaOptimizer<Scalar,3,true>(seq_loss, samples, 5e-2)),
			dims, samples, time_steps, epochs, seq_epsilon);
	rnn_train_test<Scalar,3>("adadelta mini-batch",
			OptimizerPtr<Scalar,3,true>(new AdaDeltaOptimizer<Scalar,3,true>(seq_loss, 10, 1e-2)),
			dims, samples, time_steps, epochs, seq_epsilon);
	rnn_train_test<Scalar,3>("adadelta online",
			OptimizerPtr<Scalar,3,true>(new AdaDeltaOptimizer<Scalar,3,true>(seq_loss, 1, 1e-1)),
			dims, samples, time_steps, epochs, seq_epsilon);
}

TEST(TrainingTest, AdaDelta) {
	adadelta_train_test<float>();
	adadelta_train_test<double>();
}

/**
 * Performs training tests using an Adam optimizer.
 */
template<typename Scalar>
inline void adam_train_test() {
	const unsigned samples = 20;
	const unsigned time_steps = 3;
	const unsigned epochs = 500;
	const Scalar epsilon = .015;
	const Scalar seq_epsilon = .1;
	const Dimensions<std::size_t,3> dims({ 6u, 6u, 2u });
	auto loss = std::make_shared<SquaredLoss<Scalar,3,false>>();
	ffnn_train_test<Scalar,3>("adam batch",
			OptimizerPtr<Scalar,3,false>(new AdamOptimizer<Scalar,3,false>(loss, samples, 1e-2)),
			dims, samples, epochs, epsilon);
	ffnn_train_test<Scalar,3>("adam mini-batch",
			OptimizerPtr<Scalar,3,false>(new AdamOptimizer<Scalar,3,false>(loss, 10, 1e-3, 5e-2, 5e-3)),
			dims, samples, epochs, epsilon);
	ffnn_train_test<Scalar,3>("adam online",
			OptimizerPtr<Scalar,3,false>(new AdamOptimizer<Scalar,3,false>(loss)),
			dims, samples, epochs, epsilon);
	auto seq_loss = std::make_shared<SquaredLoss<Scalar,3,true>>();
	rnn_train_test<Scalar,3>("adam batch",
			OptimizerPtr<Scalar,3,true>(new AdamOptimizer<Scalar,3,true>(seq_loss, samples, 5e-3)),
			dims, samples, time_steps, epochs, seq_epsilon);
	rnn_train_test<Scalar,3>("adam mini-batch",
			OptimizerPtr<Scalar,3,true>(new AdamOptimizer<Scalar,3,true>(seq_loss, 10, 1e-3, 5e-2, 5e-3)),
			dims, samples, time_steps, epochs, seq_epsilon);
	rnn_train_test<Scalar,3>("adam online",
			OptimizerPtr<Scalar,3,true>(new AdamOptimizer<Scalar,3,true>(seq_loss)),
			dims, samples, time_steps, epochs, seq_epsilon);
}

TEST(TrainingTest, Adam) {
	adam_train_test<float>();
	adam_train_test<double>();
}

/**
 * Performs training tests using an AdaMax optimizer.
 */
template<typename Scalar>
inline void adamax_train_test() {
	const unsigned samples = 20;
	const unsigned time_steps = 3;
	const unsigned epochs = 500;
	const Scalar epsilon = .015;
	const Scalar seq_epsilon = .1;
	const Dimensions<std::size_t,3> dims({ 6u, 6u, 2u });
	auto loss = std::make_shared<SquaredLoss<Scalar,3,false>>();
	ffnn_train_test<Scalar,3>("adamax batch",
			OptimizerPtr<Scalar,3,false>(new AdaMaxOptimizer<Scalar,3,false>(loss, samples, 1e-2)),
			dims, samples, epochs, epsilon);
	ffnn_train_test<Scalar,3>("adamax mini-batch",
			OptimizerPtr<Scalar,3,false>(new AdaMaxOptimizer<Scalar,3,false>(loss, 10, 1e-3, 5e-2, 5e-3)),
			dims, samples, epochs, epsilon);
	ffnn_train_test<Scalar,3>("adamax online",
			OptimizerPtr<Scalar,3,false>(new AdaMaxOptimizer<Scalar,3,false>(loss)),
			dims, samples, epochs, epsilon);
	auto seq_loss = std::make_shared<SquaredLoss<Scalar,3,true>>();
	rnn_train_test<Scalar,3>("adamax batch",
			OptimizerPtr<Scalar,3,true>(new AdaMaxOptimizer<Scalar,3,true>(seq_loss, samples, 5e-3)),
			dims, samples, time_steps, epochs, seq_epsilon);
	rnn_train_test<Scalar,3>("adamax mini-batch",
			OptimizerPtr<Scalar,3,true>(new AdaMaxOptimizer<Scalar,3,true>(seq_loss, 10, 1e-3, 5e-2, 5e-3)),
			dims, samples, time_steps, epochs, seq_epsilon);
	rnn_train_test<Scalar,3>("adamax online",
			OptimizerPtr<Scalar,3,true>(new AdaMaxOptimizer<Scalar,3,true>(seq_loss)),
			dims, samples, time_steps, epochs, seq_epsilon);
}

TEST(TrainingTest, AdaMax) {
	adamax_train_test<float>();
	adamax_train_test<double>();
}

/**
 * Performs training tests using an Nadam optimizer.
 */
template<typename Scalar>
inline void nadam_train_test() {
	const unsigned samples = 20;
	const unsigned time_steps = 3;
	const unsigned epochs = 500;
	const Scalar epsilon = .01;
	const Scalar seq_epsilon = .06;
	const Dimensions<std::size_t,3> dims({ 6u, 6u, 2u });
	auto loss = std::make_shared<SquaredLoss<Scalar,3,false>>();
	ffnn_train_test<Scalar,3>("nadam batch",
			OptimizerPtr<Scalar,3,false>(new NadamOptimizer<Scalar,3,false>(loss, samples, 1e-2)),
			dims, samples, epochs, epsilon);
	ffnn_train_test<Scalar,3>("nadam mini-batch",
			OptimizerPtr<Scalar,3,false>(new NadamOptimizer<Scalar,3,false>(loss, 10, 1e-3, 5e-2, 5e-3)),
			dims, samples, epochs, epsilon);
	ffnn_train_test<Scalar,3>("nadam online",
			OptimizerPtr<Scalar,3,false>(new NadamOptimizer<Scalar,3,false>(loss)),
			dims, samples, epochs, epsilon);
	auto seq_loss = std::make_shared<SquaredLoss<Scalar,3,true>>();
	rnn_train_test<Scalar,3>("nadam batch",
			OptimizerPtr<Scalar,3,true>(new NadamOptimizer<Scalar,3,true>(seq_loss, samples, 5e-3)),
			dims, samples, time_steps, epochs, seq_epsilon);
	rnn_train_test<Scalar,3>("nadam mini-batch",
			OptimizerPtr<Scalar,3,true>(new NadamOptimizer<Scalar,3,true>(seq_loss, 10, 1e-3, 5e-2, 5e-3)),
			dims, samples, time_steps, epochs, seq_epsilon);
	rnn_train_test<Scalar,3>("nadam online",
			OptimizerPtr<Scalar,3,true>(new NadamOptimizer<Scalar,3,true>(seq_loss)),
			dims, samples, time_steps, epochs, seq_epsilon);
}

TEST(TrainingTest, Nadam) {
	nadam_train_test<float>();
	nadam_train_test<double>();
}

/**
 * Performs training tests using an AMSGrad optimizer.
 */
template<typename Scalar>
inline void amsgrad_train_test() {
	const unsigned samples = 20;
	const unsigned time_steps = 3;
	const unsigned epochs = 500;
	const Scalar epsilon = .015;
	const Scalar seq_epsilon = .12;
	const Dimensions<std::size_t,3> dims({ 6u, 6u, 2u });
	auto loss = std::make_shared<SquaredLoss<Scalar,3,false>>();
	ffnn_train_test<Scalar,3>("amsgrad batch",
			OptimizerPtr<Scalar,3,false>(new AMSGradOptimizer<Scalar,3,false>(loss, samples, 1e-2)),
			dims, samples, epochs, epsilon);
	ffnn_train_test<Scalar,3>("amsgrad mini-batch",
			OptimizerPtr<Scalar,3,false>(new AMSGradOptimizer<Scalar,3,false>(loss, 10, 1e-3)),
			dims, samples, epochs, epsilon);
	ffnn_train_test<Scalar,3>("amsgrad online",
			OptimizerPtr<Scalar,3,false>(new AMSGradOptimizer<Scalar,3,false>(loss)),
			dims, samples, epochs, epsilon);
	auto seq_loss = std::make_shared<SquaredLoss<Scalar,3,true>>();
	rnn_train_test<Scalar,3>("amsgrad batch",
			OptimizerPtr<Scalar,3,true>(new AMSGradOptimizer<Scalar,3,true>(seq_loss, samples, 5e-3)),
			dims, samples, time_steps, epochs, seq_epsilon);
	rnn_train_test<Scalar,3>("amsgrad mini-batch",
			OptimizerPtr<Scalar,3,true>(new AMSGradOptimizer<Scalar,3,true>(seq_loss, 10, 1e-3)),
			dims, samples, time_steps, epochs, seq_epsilon);
	rnn_train_test<Scalar,3>("amsgrad online",
			OptimizerPtr<Scalar,3,true>(new AMSGradOptimizer<Scalar,3,true>(seq_loss)),
			dims, samples, time_steps, epochs, seq_epsilon);
}

TEST(TrainingTest, AMSGrad) {
	amsgrad_train_test<float>();
	amsgrad_train_test<double>();
}

} /* namespace test */
} /* namespace cattle */


int main(int argc, char** argv) {
	static const char* verbose_flag = "-verbose";
	for (int i = 1; i < argc; ++i) {
		if (!strcmp(argv[i], verbose_flag)) {
			cattle::test::verbose = true;
			break;
		}
	}
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}