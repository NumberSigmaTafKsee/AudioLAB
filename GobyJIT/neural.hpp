#pragma once


#include <vector>
#include <memory>
#include <iostream>
#include <unsupported/Eigen/CXX11/Tensor>

namespace nn {
    /**
     * @brief How to initialize the dense weights
     */
    enum class InitializationScheme {
        GlorotUniform,
        GlorotNormal
    };

    template <typename Dtype>
    class WeightDistribution {
    public:
        /**
         * @brief Create a weight distribution to draw from
         * @param scheme [in]: The scheme to initialize with
         * @param fanIn [in]: The fan in of the layer
         * @param fanOut [in]: The fan out of the layer
         */
        explicit WeightDistribution(InitializationScheme scheme, int fanIn, int fanOut):
                m_scheme(scheme),
                m_randomNumberGenerator(std::random_device()())
        {
            if (m_scheme == InitializationScheme::GlorotUniform) {
                Dtype limit = std::sqrt(6.0 / (fanIn + fanOut));
                m_uniformDist.reset(new std::uniform_real_distribution<Dtype>(-limit, limit));
            } else if (m_scheme == InitializationScheme::GlorotNormal) {
                Dtype std = std::sqrt(2.0 / (fanIn + fanOut));
                m_normalDist.reset(new std::normal_distribution<Dtype>(0, std));
            }
        }

        /**
         * @brief Get a value from the distribution
         * @return
         */
        Dtype get() {
            if (m_scheme == InitializationScheme::GlorotUniform) {
                return (*m_uniformDist)(m_randomNumberGenerator);
            } else if (m_scheme == InitializationScheme::GlorotNormal) {
                return (*m_normalDist)(m_randomNumberGenerator);
            } else {
                std::cerr << "Tried to draw from distribution that is uninitialized" << std::endl;
                exit(-1);
            }
        }

    private:
        InitializationScheme m_scheme;                                        ///< Our init scheme
        std::mt19937 m_randomNumberGenerator;                                 ///< Our random number generator
        std::unique_ptr<std::uniform_real_distribution<Dtype>> m_uniformDist; ///< Our uniform distribution
        std::unique_ptr<std::normal_distribution<Dtype>> m_normalDist;        ///< Our normal distribution
    };


    /**
     * @brief Initialize a tensor of dimension (input x output) with a specified scheme
     * @tparam Dtype [in]: Datatype of the tensor (float/double)
     * @param inputDimensions [in]: The input dimensions of the layer
     * @param outputDimensions [in]: The output dimensions of the layer
     * @param scheme [in]: Initialization Scheme
     * @return A randomly initialized tensor
     *
     * @note This function only exists because I can't seem to get Tensor.setRandom<Generator> to work
     *       with their builtins. This is way, way less efficient, but is only called on creation of a new layer
     */
    template <typename Dtype>
    Eigen::Tensor<Dtype, 2> getRandomWeights(int inputDimensions, int outputDimensions,
                                             InitializationScheme scheme = InitializationScheme::GlorotUniform) {
        Eigen::Tensor<Dtype, 2> weights(inputDimensions, outputDimensions);
        weights.setZero();

        auto distribution = WeightDistribution<Dtype>(scheme, inputDimensions, outputDimensions);
        for (unsigned int ii = 0; ii < inputDimensions; ++ii) {
            for (unsigned int jj = 0; jj < outputDimensions; ++jj) {
                weights(ii, jj) = distribution.get();
            }
        }
        return weights;
    };
}

namespace nn {
    template <typename Dtype, int Dims>
    class OptimizerImpl {
    public:
        virtual Eigen::Tensor<Dtype, Dims> weightUpdate(const Eigen::Tensor<Dtype, Dims> &weights) = 0;
    };
}

namespace nn {
    namespace internal {
        template<typename Dtype, int Dims>
        class StochasticGradientDescentImpl : public OptimizerImpl<Dtype, Dims> {
        public:

            // TODO: Add momentum
            /**
             * @brief Initialize our SGD Solver
             * @param learningRate [in]: The learning rate of SGD
             */
            explicit StochasticGradientDescentImpl(Dtype learningRate):
                    m_learningRate(learningRate) {}

            /**
             * @brief Get the update to apply to the weights
             * @param gradWeights [in]: Weights to update
             * @return The factor to update the weights by
             */
            Eigen::Tensor<Dtype, Dims> weightUpdate(const Eigen::Tensor<Dtype, Dims> &gradWeights) {
                return gradWeights * gradWeights.constant(m_learningRate);
            };

        private:
            Dtype m_learningRate; ///< Our current learning rate
        };
    }
}

namespace nn {
    namespace internal {
        template<typename Dtype, int Dims>
        class AdamImpl : public OptimizerImpl<Dtype, Dims> {
        public:

            /**
             * @brief Initialize our Adam Solver
             * @param learningRate [in]: Base learning rate
             * @param beta1 [in]: The first moment decay factor (default = 0.9)
             * @param beta2 [in]: The second moment decay factor (default = 0.999)
             * @param epsilon [in]: A stabilizing factor for division (default = 1e-8)
             */
            explicit AdamImpl(Dtype learningRate, Dtype beta1, Dtype beta2, Dtype epsilon):
                    m_learningRate(learningRate), m_beta1(beta1), m_beta2(beta2), m_epsilon(epsilon),
                    m_isInitialized(false), m_currentTimestep(1)
            {}

            /**
             * @brief Get the update to apply to the weights
             * @param gradWeights [in]: Weights to update
             * @return The factor to update the weights by
             */
            Eigen::Tensor<Dtype, Dims> weightUpdate(const Eigen::Tensor<Dtype, Dims> &gradWeights) {
                if (!m_isInitialized) {
                    m_firstMoment = Eigen::Tensor<Dtype, Dims>(gradWeights.dimensions());
                    m_firstMoment.setZero();

                    m_secondMoment = Eigen::Tensor<Dtype, Dims>(gradWeights.dimensions());
                    m_secondMoment.setZero();
                    m_isInitialized = true;
                }

                // m_t = B_1 * m_(t-1) + (1 - B_1) * g_t
                m_firstMoment = m_firstMoment.constant(m_beta1) * m_firstMoment +
                        gradWeights.constant(1 - m_beta1) * gradWeights;

                // v_t = B_2 * v_(t-1) + (1 - B_2) * g_t^2
                m_secondMoment = m_secondMoment.constant(m_beta2) * m_secondMoment +
                        gradWeights.constant(1 - m_beta2) * gradWeights.square();
//
//                std::cout << "First moment: " << m_firstMoment << std::endl;
//                std::cout << "Second moment: " << m_secondMoment << std::endl;
//                std::cout << std::endl << std::endl << std::endl;

                auto biasCorrectedFirstMoment = m_firstMoment / m_firstMoment.constant(1 - pow(m_beta1, m_currentTimestep));
                auto biasCorrectedSecondMoment = m_secondMoment / m_secondMoment.constant(1 - pow(m_beta2, m_currentTimestep));
//
//                std::cout << "Bias corrected first: " << biasCorrectedFirstMoment << std::endl;
//                std::cout << "Bias corrected second: " << biasCorrectedSecondMoment << std::endl;
//                std::cout << std::endl << std::endl << std::endl;


                m_currentTimestep ++;
                // Return firstMoment  * (learning_rate) / (sqrt(secondMoment) + epsilon)
                return biasCorrectedFirstMoment * (
                              (gradWeights.constant(m_learningRate) /
                               (biasCorrectedSecondMoment.sqrt() + gradWeights.constant(m_epsilon))
                ));
            };

        private:
            Dtype m_learningRate; ///< The learning rate of our optimizer
            Dtype m_beta1;        ///< Our B1 parameter (first moment decay)
            Dtype m_beta2;        ///< Our B2 parameter (second moment decay)
            Dtype m_epsilon;      ///< Stability factor

            bool m_isInitialized;      ///< On our first iteration, set the first and second order gradients to zero
            size_t m_currentTimestep;  ///< Our current timestep (iteration)

            // Our exponentially decaying average of past gradients
            Eigen::Tensor<Dtype, Dims> m_firstMoment;  ///< Our m_t term that represents the first order gradient decay
            Eigen::Tensor<Dtype, Dims> m_secondMoment; ///< Our v_t term that represents the second order gradient decay
        };
    }
}

namespace nn {
    /**
     * The current design is that you declare your optimizer in your main training area
     * and the Net class propagates this to all layers, which create their own Impls
     * with one Impl per weight. The design is geared towards more complex optimizers
     */


    /**
     * @brief Factory method of SGD
     *
     * @tparam Dtype : The floating point type of the optimizer
     */
    template <typename Dtype>
    class StochasticGradientDescent {
    public:
        /**
         * @brief Create a SGD factory w/ learning rate
         * @param learningRate [in]: Learning rate of SGD optimizer
         */
        explicit StochasticGradientDescent(Dtype learningRate):
            m_learningRate(learningRate) {}

        /**
         * @brief Create an optimizer Impl for our given type
         * @tparam Dims [in]: The dimensionality of the tensor the optimizer will update
         * @return An optimizer impl that can update weights and keep track of state
         */
        template <int Dims>
        std::unique_ptr<OptimizerImpl<Dtype, Dims>> createOptimizer() const {
            return std::unique_ptr<OptimizerImpl<Dtype, Dims>>(new internal::StochasticGradientDescentImpl<Dtype, Dims>(m_learningRate));
        }

    private:
        Dtype m_learningRate; ///< The learning rate
    };

    template <typename Dtype>
    class Adam {
    public:
        /**
         * @brief Create an Adam optimizer
         * @param learningRate [in]: Base learning rate
         * @param beta1 [in]: The first moment decay factor (default = 0.9)
         * @param beta2 [in]: The second moment decay factor (default = 0.999)
         * @param epsilon [in]: A stabilizing factor for division (default = 1e-8)
         */
        explicit Adam(Dtype learningRate, Dtype beta1 = 0.9, Dtype beta2 = 0.999, Dtype epsilon = 1e-8):
                m_learningRate(learningRate), m_beta1(beta1), m_beta2(beta2), m_epsilon(epsilon)
        {}

        /**
         * Create an optimizer Impl for our given type
         * @tparam Dims [in]: The dimensionality of the tensor the optimizer will update
         * @return An optimizer impl that can update weights and keep track of state
         */
        template <int Dims>
        std::unique_ptr<OptimizerImpl<Dtype, Dims>> createOptimizer() const {
            return std::unique_ptr<OptimizerImpl<Dtype, Dims>>(new internal::AdamImpl<Dtype, Dims>(m_learningRate, m_beta1, m_beta2, m_epsilon));
        };

    private:
        Dtype m_learningRate; ///< The learning rate of our optimizer
        Dtype m_beta1;        ///< Our B1 parameter (first moment decay)
        Dtype m_beta2;        ///< Our B2 parameter (second moment decay)
        Dtype m_epsilon;      ///< Stability factor
    };


}

namespace nn {
    template <typename Dtype = float, int Dims = 2>
    class Layer {
    public:
        /**
         * @brief Return the name of the layer
         */
        virtual const std::string& getName() = 0;

        /**
         * @brief Take an input tensor, perform an operation on it, and return a new tensor
         * @param input [in]: The input tensor (from the previous layer)
         * @return An output tensor, which is fed into the next layer
         */
        virtual Eigen::Tensor<Dtype, Dims> forward(const Eigen::Tensor<Dtype, Dims> &input) = 0;

        /**
         * @brief Perform the backwards operation on the layer.
         * @param input [in]: The input tensor (from next layer)
         * @return The output tensor, which is fed into the previous layer
         */
        virtual Eigen::Tensor<Dtype, Dims> backward(const Eigen::Tensor<Dtype, Dims> &output) = 0;

        /**
         * @brief Update the weights after a backwards pass
         */
        virtual void step() = 0;

        // TODO: Need to find a clean way to inherit optimizers to reduce repetition
        // TODO: If anyone is reading this, and has ideas... please let me know virtual methods w/ templates is impossible :(
        /**
         * @brief Registers the optimizer with the layer
         * @param optimizer [in]: The optimizer to register
         */
        virtual void registerOptimizer(std::shared_ptr<StochasticGradientDescent<Dtype>> optimizer) = 0;

        /**
         * @brief Registers the optimizer with the layer
         * @param optimizer [in]: The optimizer to register
         */
        virtual void registerOptimizer(std::shared_ptr<Adam<Dtype>> optimizer) = 0;
    };
}

namespace nn {

    template <typename Dtype = float, int Dims = 2>
    class Dense : public Layer<Dtype, Dims> {
    public:
        /**
         * @brief Create a Dense layer
         * @param batchSize [in]: The batch size going through the network
         * @param inputDimension [in]: Expected input dimension
         * @param outputDimension [in]: The output dimensionality (number of neurons)
         * @param useBias [in]: Whether to use a bias term
         * @param weightInitializer [in]: The weight initializer scheme to use. Defaults to GlorotUniform
         */
        explicit Dense(int batchSize, int inputDimension, int outputDimension, bool useBias,
                       InitializationScheme weightInitializer = InitializationScheme::GlorotUniform);

        /**
         * @brief Return the name of the layer
         * @return The layer name
         */
        const std::string& getName() {
            const static std::string name = "Dense";
            return name;
        }

        /**
         * @brief Forward through the layer (compute the output)
         * @param input [in]: The input to the layer (either data or previous layer output)
         * @return The output of this layer
         */
        Eigen::Tensor<Dtype, Dims> forward(const Eigen::Tensor<Dtype, Dims> &input);

        /**
         * @brief Compute the gradient (backward pass) of the layer
         * @param accumulatedGrad [in]: The input to the backwards pass. (from next layer)
         * @return The output of the backwards pass (sent to previous layer)
         */
        Eigen::Tensor<Dtype, Dims> backward(const Eigen::Tensor<Dtype, Dims> &accumulatedGrad);

        /**
         * @brief Get the input shape
         * @return The input shape
         */

        Eigen::array<Eigen::Index, Dims> getOutputShape() {
            return m_outputShape;
        };

        /**
         * @brief Update weights of the layer w.r.t. gradient
         */
        void step();

        // TODO: Find a nicer way then duplication for subtype optimizer factories
        /**
         * @brief Set up the optimizer for our weights
         */
        void registerOptimizer(std::shared_ptr<StochasticGradientDescent<Dtype>> optimizer);

        /**
         * @brief Set up the optimizer for our weights
         */
        void registerOptimizer(std::shared_ptr<Adam<Dtype>> optimizer);

    private:
        Eigen::array<Eigen::Index, Dims> m_outputShape;                ///< The output shape of this layer
        Eigen::Tensor<Dtype, Dims> m_inputCache;                       ///< Cache the input to calculate gradient
        Eigen::Tensor<Dtype, Dims> m_weights;                          ///< Our weights of the layer
        Eigen::Tensor<Dtype, Dims> m_bias;                             ///< The bias weights if specified

        // Gradients
        Eigen::Tensor<Dtype, Dims> m_weightsGrad;                      ///< The gradient of the weights
        Eigen::Tensor<Dtype, Dims> m_biasGrad;                         ///< The gradient of the bias
        std::unique_ptr<OptimizerImpl<Dtype, Dims>> m_weightOptimizer; ///< The optimizer of our weights
        std::unique_ptr<OptimizerImpl<Dtype, Dims>> m_biasOptimizer;   ///< The optimizer of our bias

        bool m_useBias;                                                ///< Whether we use the bias
    };

    template <typename Dtype, int Dims>
    Dense<Dtype, Dims>::Dense(int batchSize, int inputDimension, int outputDimension, bool useBias,
                              InitializationScheme weightInitializer):
            m_outputShape({batchSize, outputDimension}),
            m_useBias(useBias)
    {
        m_weights = getRandomWeights<Dtype>(inputDimension, outputDimension, weightInitializer);

        m_weightsGrad = Eigen::Tensor<Dtype, Dims>(inputDimension, outputDimension);
        m_weightsGrad.setZero();

        if (useBias) {
            m_bias = getRandomWeights<Dtype>(1, outputDimension, weightInitializer);

            m_biasGrad = Eigen::Tensor<Dtype, Dims>(1, outputDimension);
            m_biasGrad.setZero();
        }
    };

    template <typename Dtype, int Dims>
    Eigen::Tensor<Dtype, Dims> Dense<Dtype, Dims>::forward(const Eigen::Tensor<Dtype, Dims> &input) {
        assert(input.dimensions()[1] == m_weights.dimensions()[0] &&
                            "Dense::forward dimensions of input and weights do not match");
        m_inputCache = input;

        Eigen::array<Eigen::IndexPair<int>, 1> productDims = { Eigen::IndexPair<int>(1, 0) };
        auto result = input.contract(m_weights, productDims);

        if (m_useBias) {
            // Copy the bias from (1, outputSize) to (inputBatchSize, outputDimension)
            return result + m_bias.broadcast(Eigen::array<Eigen::Index, 2>{input.dimensions()[0], 1});
        } else {
            return result;
        }
    }

    template <typename Dtype, int Dims>
    Eigen::Tensor<Dtype, Dims> Dense<Dtype, Dims>::backward(const Eigen::Tensor<Dtype, Dims> &accumulatedGrad) {
        assert(accumulatedGrad.dimensions()[0] == m_inputCache.dimensions()[0] &&
                       "Dense::backward dimensions of accumulatedGrad and inputCache do not match");
        // m_inputCache is of shape (batchSize, inputDimension)
        // accumulatedGrad is of shape (batchSize, outputDimension)
        // So we want to contract along dimensions (0, 0), aka m_inputCache.T * accumulatedGrad
        // Where dimensions would be (inputDimension, batchSize) * (batchSize, outputDimension)
        static const Eigen::array<Eigen::IndexPair<int>, 1> transposeInput = { Eigen::IndexPair<int>(0, 0) };

        m_weightsGrad = m_inputCache.contract(accumulatedGrad, transposeInput);
        if (m_useBias) {
            m_biasGrad = accumulatedGrad.sum(Eigen::array<int, 1>{0}).eval().reshape(Eigen::array<Eigen::Index, 2>{1, m_outputShape[1]});
        }

        // accumulatedGrad is of shape (batchSize, outputDimensions)
        // m_weights is of shape (inputDimensions, outputDimensions)
        // So we want to contract along dimensions (1, 1), which would be accumulatedGrad * m_weights.T
        // Where dimensions would be (batchSize, outputDimension) * (outputDimension, inputDimension)
        static const Eigen::array<Eigen::IndexPair<int>, 1> transposeWeights = { Eigen::IndexPair<int>(1, 1)};
        return accumulatedGrad.contract(m_weights, transposeWeights);
    }

    template <typename Dtype, int Dims>
    void Dense<Dtype, Dims>::step() {
        m_weights -= m_weightOptimizer->weightUpdate(m_weightsGrad);

        if (m_useBias) {
            m_bias -= m_biasOptimizer->weightUpdate(m_biasGrad);
        }
    }

    // TODO: Find a nicer way then duplication for subtype optimizer factories
    template <typename Dtype, int Dims>
    void Dense<Dtype, Dims>::registerOptimizer(std::shared_ptr<StochasticGradientDescent<Dtype>> optimizer) {
        m_weightOptimizer = std::move(optimizer->template createOptimizer<Dims>());

        if (m_useBias) {
            m_biasOptimizer = std::move(optimizer->template createOptimizer<Dims>());
        }
    }

    template <typename Dtype, int Dims>
    void Dense<Dtype, Dims>::registerOptimizer(std::shared_ptr<Adam<Dtype>> optimizer) {
        m_weightOptimizer = std::move(optimizer->template createOptimizer<Dims>());

        if (m_useBias) {
            m_biasOptimizer = std::move(optimizer->template createOptimizer<Dims>());
        }

    }
}

namespace nn {
    template <typename Dtype = float, int Dims = 2>
    class Relu : public Layer<Dtype, Dims> {
    public:

        /**
         * @brief initialize Relu
         */
        Relu() = default;

        /**
         * @brief Return the name of the layer
         * @return The layer name
         */
        const std::string& getName() {
            const static std::string name = "Relu";
            return name;
        }

        /**
         * @brief Forward through the layer (compute the output)
         * @param input [in]: The input tensor to apply Relu to
         * @return max(0, input)
         */
        Eigen::Tensor<Dtype, Dims> forward(const Eigen::Tensor<Dtype, Dims> &input);

        /**
         * @brief Compute the gradient (backward pass) of the layer
         * @param accumulatedGrad [in]: The input to the backwards pass (from the next layer)
         * @return The output of the backwards pass (sent to the previous layer)
         */
        Eigen::Tensor<Dtype, Dims> backward(const Eigen::Tensor<Dtype, Dims> &accumulatedGrad);

        /**
         * @brief Void function in relu
         */
        void step() {}

        /**
         * @brief Void function in relu
         */
        void registerOptimizer(std::shared_ptr<StochasticGradientDescent<Dtype>> optimizer) {}

        /**
         * @brief Void function in relu
         */
        void registerOptimizer(std::shared_ptr<Adam<Dtype>> optimizer) {}

    private:
        Eigen::Tensor<Dtype, Dims> m_output; ///< The output of the forward pass
    };

    template <typename Dtype, int Dims>
    Eigen::Tensor<Dtype, Dims> Relu<Dtype, Dims>::forward(const Eigen::Tensor<Dtype, Dims> &input) {
        m_output = input.cwiseMax(static_cast<Dtype>(0));
        return m_output;
    };

    template <typename Dtype, int Dims>
    Eigen::Tensor<Dtype, Dims> Relu<Dtype, Dims>::backward(const Eigen::Tensor<Dtype, Dims> &accumulatedGrad) {
        // Could also check a cached input to Relu::forward, but since
        // the output is simply (x, 0), we can just check our already cached output.
        auto inputPositive = m_output > static_cast<Dtype>(0);
        return inputPositive.select(accumulatedGrad, accumulatedGrad.constant(0.0));
    }
}

namespace nn {
    template <typename Dtype = float, int Dims = 2>
    class Softmax : public Layer<Dtype, Dims> {
    public:

        /**
         * @brief initialize Softmax
         */
        Softmax() = default;

        /**
         * @brief Return the name of the layer
         * @return The layer name
         */
        const std::string& getName() {
            const static std::string name = "Softmax";
            return name;
        }

        /**
         * @brief Forward through the layer (compute the output)
         * @param input [in]: The input tensor to apply softmax to
         * @return
         */
        Eigen::Tensor<Dtype, Dims> forward(const Eigen::Tensor<Dtype, Dims> &input);

        /**
         * @brief Compute the gradient (backwards pass) of the layer
         * @param accumulatedGrad [in]: The input tensor to the backwards pass (from the next layer). This should be one hot encoded labels
         * @return The output of the backwards pass (sent ot the previous layer)
         */
        Eigen::Tensor<Dtype, Dims> backward(const Eigen::Tensor<Dtype, Dims> &accumulatedGrad);

        /**
         * @brief Update Weights (doesn't do anything w/ softmax)
         */
        void step() {}

        /**
         * @brief Void function in softmax
         */
        void registerOptimizer(std::shared_ptr<StochasticGradientDescent<Dtype>> optimizer) {}

        /**
         * @brief Void function in softmax
         */
        void registerOptimizer(std::shared_ptr<Adam<Dtype>> optimizer) {}

    private:
        Eigen::Tensor<Dtype, Dims> m_output; ///< The output of the forward pass
    };

    template <typename Dtype, int Dims>
    Eigen::Tensor<Dtype, Dims> Softmax<Dtype, Dims>::forward(const Eigen::Tensor<Dtype, Dims> &input) {
        int batchSize = input.dimensions()[0];
        int classDims = input.dimensions()[1];
        auto shiftedInput = input - input.maximum(Eigen::array<int, 1>{1})
                                    .eval().reshape(Eigen::array<int, 2>{batchSize, 1})
                                    .broadcast(Eigen::array<int, 2>{1, classDims});

        auto exponentiated = shiftedInput.exp();
        m_output = exponentiated * exponentiated.sum(Eigen::array<int, 1>{1})
                                   .inverse().eval()
                                   .reshape(Eigen::array<int, 2>({batchSize, 1}))
                                   .broadcast(Eigen::array<int, 2>({1, classDims}));
        return m_output;
    }

    template <typename Dtype, int Dims>
    Eigen::Tensor<Dtype, Dims> Softmax<Dtype, Dims>::backward(const Eigen::Tensor<Dtype, Dims> &accumulatedGrad) {
        const int batchSize = accumulatedGrad.dimensions()[0];
        assert(batchSize == m_output.dimensions()[0] && "Dimensions of number of batches does not match");
        return accumulatedGrad / accumulatedGrad.constant(batchSize);
    }
}

namespace nn {
    template<typename Dtype, int Dims>
    class CrossEntropyLoss {
    public:
        /**
         * @brief Create a cross entropy loss layer
         */
        CrossEntropyLoss() = default;

        /**
         * @brief Calculate the cross entropy loss
         * @param probabilities [in]: "Probabilities" as in 0-1 values output by a layer like Softmax
         * @param labels [in]: One hot encoded labels
         * @return The loss
         */
        Dtype loss(const Eigen::Tensor<Dtype, Dims> &probabilities, const Eigen::Tensor<Dtype, Dims> &labels);

        /**
         * @brief Calculate the accuracy of our labels
         * @param probabilities [in]: "Probabilities" as in 0-1 values output by a layer like Softmax
         * @param labels [in]: One hot encoded labels
         * @return The total accuracy (num_correct / total)
         */
        Dtype accuracy(const Eigen::Tensor<Dtype, Dims> &probabilities, const Eigen::Tensor<Dtype, Dims> &labels);

        /**
         * @brief Compute the gradient for Cross Entropy Loss
         * @param probabilities [in]: "Probabilities" as in 0-1 values output by a layer like Softmax
         * @param labels [in]: One hot encoded labels
         * @return The gradient of this loss layer
         */
        Eigen::Tensor<Dtype, Dims>
        backward(const Eigen::Tensor<Dtype, Dims> &probabilities, const Eigen::Tensor<Dtype, Dims> &labels);
    };

    template<typename Dtype, int Dims>
    Dtype CrossEntropyLoss<Dtype, Dims>::loss(const Eigen::Tensor<Dtype, Dims> &probabilities,
                                              const Eigen::Tensor<Dtype, Dims> &labels) {
        int batchSize = probabilities.dimensions()[0];

        // TODO: Do I need a stabilizing const here?
        static const Dtype stabilizingVal = 0.0001;
        Eigen::Tensor<Dtype, 0> summedLoss = (labels *
                                              (probabilities.constant(stabilizingVal) + probabilities).log()).sum();
        return (-1.0 / batchSize) * summedLoss(0);
    }

    template<typename Dtype, int Dims>
    Dtype CrossEntropyLoss<Dtype, Dims>::accuracy(const Eigen::Tensor<Dtype, Dims> &probabilities,
                                                  const Eigen::Tensor<Dtype, Dims> &labels) {
        assert(probabilities.dimensions()[0] == labels.dimensions()[0] &&
               "CrossEntropy::accuracy dimensions did not match");
        assert(probabilities.dimensions()[1] == labels.dimensions()[1] &&
               "CrossEntropy::accuracy dimensions did not match");

        auto batchSize = static_cast<Dtype>(labels.dimensions()[0]);

        // Argmax across dimension = 1 (so we get a column vector)
        Eigen::Tensor<bool, 1> ifTensor = probabilities.argmax(1) == labels.argmax(1);
        Eigen::Tensor<Dtype, 1> thenTensor(batchSize);
        auto result = ifTensor.select(thenTensor.constant(1.0), thenTensor.constant(0));
        Eigen::Tensor<Dtype, 0> count = result.sum();
        return static_cast<Dtype>(count(0)) / batchSize;
    }

    template<typename Dtype, int Dims>
    Eigen::Tensor<Dtype, Dims> CrossEntropyLoss<Dtype, Dims>::backward(const Eigen::Tensor<Dtype, Dims> &probabilities,
                                                                       const Eigen::Tensor<Dtype, Dims> &labels) {
        return probabilities - labels;
    }
}

namespace nn {
    template<typename Dtype, int Dims>
    class HuberLoss {
    public:

        /**
         * @brief Initialize a SmoothL1Loss loss function
         */
        explicit HuberLoss(Dtype threshold = 1.0): m_threshold(threshold) {}

        /**
         * @brief Compute the loss
         * @param predictions [in]: Predictions from the network
         * @param labels [in]: Labels to compute loss with
         * @return The loss as a scalar
         */
        Dtype loss(const Eigen::Tensor<Dtype, Dims> &predictions, const Eigen::Tensor<Dtype, Dims> &labels);

        /**
         * @brief Compute the gradient of the loss given this data
         * @param predictions [in]: Predictions from the network
         * @param labels [in]: Labels from dataset
         * @return The gradient of the loss layer
         */
        Eigen::Tensor<Dtype, Dims>
        backward(const Eigen::Tensor<Dtype, Dims> &predictions, const Eigen::Tensor<Dtype, Dims> &labels);

    private:
        Dtype m_threshold;                                ///< The threshold used to determine which part of the piecewise loss
        Eigen::Tensor<bool, Dims> m_cachedSwitchResults;  ///< Whether abs(y - y_hat) <= m_threshold
    };

    template<typename Dtype, int Dims>
    Dtype HuberLoss<Dtype, Dims>::loss(const Eigen::Tensor<Dtype, Dims> &predictions,
                                              const Eigen::Tensor<Dtype, Dims> &labels) {
        assert(predictions.dimensions()[0] == labels.dimensions()[0] &&
               "HuberLoss::loss dimensions don't match");
        assert(predictions.dimensions()[1] == labels.dimensions()[1] &&
               "HuberLoss::loss dimensions don't match");
        int batchSize = predictions.dimensions()[0];
        // Definition taken from: https://en.wikipedia.org/wiki/Huber_loss

        // Precalculate y_hat - y
        auto error = predictions - labels;
        auto absoluteError = error.abs();

        // Set up our switch statement and cache it
        m_cachedSwitchResults = absoluteError <= m_threshold;

        // Calculate both terms for the huber loss
        auto lessThanThreshold = error.constant(0.5) * error.square();
        auto moreThanThreshold = error.constant(m_threshold) * absoluteError - error.constant(0.5 * pow(m_threshold, 2));

        // If abs(y_hat - y) <= threshold
        auto perItemLoss = m_cachedSwitchResults.select(
                lessThanThreshold, // Then use 0.5 * (y_hat - y)^2
                moreThanThreshold); // Else use thresh * |y_hat - y| - (0.5 * threshold^2)

        Eigen::Tensor<Dtype, 0> sum = perItemLoss.sum();
        // Sum and divide by N
        return sum(0) / batchSize;
    }

    template<typename Dtype, int Dims>
    Eigen::Tensor<Dtype, Dims> HuberLoss<Dtype, Dims>::backward(const Eigen::Tensor<Dtype, Dims> &predictions,
                                                                       const Eigen::Tensor<Dtype, Dims> &labels) {

        auto error = predictions - labels;

        // Note: Grad of linear part of error is threshold * (error / abs(error)), which
        // simplifies to threshold * sign(error)
        auto errorPositiveOrZero = error >= static_cast<Dtype>(0);
        auto absoluteErrorGrad = errorPositiveOrZero.select(error.constant(m_threshold), error.constant(-m_threshold));
        return m_cachedSwitchResults.select(error, absoluteErrorGrad);
    }
}

namespace nn {
    template<typename Dtype, int Dims>
    class MeanSquaredError {
    public:

        /**
         * @brief Initialize a Mean Squared Error loss function
         */
        MeanSquaredError() = default;

        /**
         * @brief Compute the MSE loss
         * @param predictions [in]: Predictions from the network
         * @param labels [in]: Labels to compute loss with
         * @return The loss as a scalar
         */
        Dtype loss(const Eigen::Tensor<Dtype, Dims> &predictions, const Eigen::Tensor<Dtype, Dims> &labels);

        /**
         * @brief Compute the gradient of Mean Squared Error given this data
         * @param predictions [in]: Predictions from the network
         * @param labels [in]: Labels from dataset
         * @return The gradient of the loss layer
         */
        Eigen::Tensor<Dtype, Dims>
        backward(const Eigen::Tensor<Dtype, Dims> &predictions, const Eigen::Tensor<Dtype, Dims> &labels);
    };

    template<typename Dtype, int Dims>
    Dtype MeanSquaredError<Dtype, Dims>::loss(const Eigen::Tensor<Dtype, Dims> &predictions,
                                              const Eigen::Tensor<Dtype, Dims> &labels) {
        assert(predictions.dimensions()[0] == labels.dimensions()[0] &&
               "MeanSquaredError::loss dimensions don't match");
        assert(predictions.dimensions()[1] == labels.dimensions()[1] &&
               "MeanSquaredError::loss dimensions don't match");

        int batchSize = predictions.dimensions()[0];

        Eigen::Tensor<Dtype, 0> squaredSum = (predictions - labels).square().sum();
        return squaredSum(0) / batchSize;
    }

    template<typename Dtype, int Dims>
    Eigen::Tensor<Dtype, Dims> MeanSquaredError<Dtype, Dims>::backward(const Eigen::Tensor<Dtype, Dims> &predictions,
                                                                       const Eigen::Tensor<Dtype, Dims> &labels) {
        return predictions - labels;
    }
}

namespace nn {

    /**
     * @brief A neural network class
     */
    template <typename Dtype = float>
    class Net {
    public:
        /**
         * @brief Init a neural network wrapper
         */
        Net() = default;

        template <int inputDim, int outputDim>
        Eigen::Tensor<Dtype, outputDim> forward(Eigen::Tensor<Dtype, inputDim> input) {
            if (m_layers.empty()) {
                std::cerr << "No layers specified" << std::endl;
                return {};
            }

            // TODO: How to ensure each forward call returns a lazily evaluated expression instead of a Tensor
            // That way we can use this to autogenerate the evaluation chain for efficiency.
            // Right now it seems to evaluate each layer individually.
            auto currentInput = input;
            for (const auto &layer : m_layers) {
                currentInput = layer->forward(currentInput);
            }
            return currentInput;
        }

        template <int labelDims>
        void backward(Eigen::Tensor<Dtype, labelDims> input) {
            if (!m_hasOptimizer) {
                std::cerr << "No registered optimizer" << std::endl;
                return;
            }

            if (m_layers.empty()) {
                std::cerr << "No layers specified" << std::endl;
                return;
            }

            auto accumulatedGrad = input;
            for (auto rit = m_layers.rbegin(); rit != m_layers.rend(); ++rit) {
                accumulatedGrad = (*rit)->backward(accumulatedGrad);
            }
        }

        void registerOptimizer(nn::StochasticGradientDescent<Dtype> *optimizer) {
            m_hasOptimizer = true;
            // TODO: Pulled this out of a private member var cause I can't supertype
            std::shared_ptr<nn::StochasticGradientDescent<Dtype>> optimizerPtr(optimizer);
            for (auto &layer : m_layers) {
                layer->registerOptimizer(optimizerPtr);
            }
        }

        void registerOptimizer(nn::Adam<Dtype> *optimizer) {
            m_hasOptimizer = true;
            // TODO: Pulled this out of a private member var cause I can't supertype
            std::shared_ptr<nn::Adam<Dtype>> optimizerPtr(optimizer);
            for (auto &layer : m_layers) {
                layer->registerOptimizer(optimizerPtr);
            }
        }

        /**
         * @brief Update weights for each layer
         */
        void step() {
            for (auto &layer : m_layers) {
                layer->step();
            }
        }

        /**
         * @brief Add a layer to the neural network
         * @param layer [in]: A layer to add
         * @return A reference to *this for method chaining
         */
        template <int Dims>
        Net<Dtype>& add(std::unique_ptr<Layer<Dtype, Dims>> layer) {
            m_layers.push_back(layer);
            return *this;
        }

        /**
         * Add a dense layer
         * @param denseLayer [in]: The dense layer to add
         * @return A reference to *this for method chaining
         */
        template <int Dims>
        Net<Dtype>& add(Dense<Dtype, Dims> *denseLayer) {
            // Do shape checks here
            m_layers.push_back(std::unique_ptr<Layer<Dtype, Dims>>(denseLayer));
            return *this;
        }

        /**
         * Add a relu layer
         * @param reluLayer [in]: The relu layer to add
         * @return A reference to *this for method chaining
         */
        template <int Dims>
        Net<Dtype>& add(Relu<Dtype, Dims> *reluLayer) {
            m_layers.push_back(std::unique_ptr<Layer<Dtype, Dims>>(reluLayer));
            return *this;
        }

        /**
         * Add a softmax layer
         * @param softmaxLayer [in]: The softmax layer to add
         * @return A reference to *this for method chaining
         */
        template <int Dims>
        Net<Dtype>& add(Softmax<Dtype, Dims> *softmaxLayer) {
            m_layers.push_back(std::unique_ptr<Layer<Dtype, Dims>>(softmaxLayer));
            return *this;
        }


    private:
        std::vector<std::unique_ptr<Layer<Dtype>>> m_layers; ///< A vector of all our layers
        bool m_hasOptimizer;                                 ///< An optimizer has been added to the net
    };
}