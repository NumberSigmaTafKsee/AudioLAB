#pragma once

#include <stdlib.h>
#include <math.h>
#include <numeric>
#include <chrono>
#include <string>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <iterator>
#include <unordered_map>
#include <vector>
#include <cmath>
#include <functional>
#include <typeinfo>
#include <typeindex>
#include <cassert>
#include "microunit.hpp"
#include "easylogging++.hpp"

//#include "Chrono.h"
#ifdef _WIN32
#include <time.h>
#else
#include <sys/time.h>
#endif

namespace utils {
//Typical sigmoid function created from input x
//Returns the sigmoided value
inline double sigmoid(double x) {
  return 1 / (1 + exp(-x));
}

// Derivative of sigmoid function
inline double deriv_sigmoid(double x) {
  return sigmoid(x)*(1 - sigmoid(x));
};

//Compute hyperbolic tangent (tanh)
//Returns the hyperbolic tangent of x.
inline double hyperbolic_tan(double x) {
  return (tanh)(x);
}

// Derivative of hyperbolic tangent function
inline double deriv_hyperbolic_tan(double x) {
  return 1 - (std::pow)(hyperbolic_tan(x), 2);
};

inline double linear(double x) {
  return x;
}

// Derivative of linear function
inline double deriv_linear(double x) {
  return 1;
};

struct ActivationFunctionsManager {
  bool GetActivationFunctionPair(const std::string & activation_name,
                                    std::pair<std::function<double(double)>,
                                    std::function<double(double)> > **pair) {
    auto iter = activation_functions_map.find(activation_name);
    if (iter != activation_functions_map.end())
      *pair = &(iter->second);
    else
      return false;
    return true;
  }

  static ActivationFunctionsManager & Singleton() {
    static ActivationFunctionsManager instance;
    return instance;
  }
private:
  void AddNewPair(std::string function_name,
                  std::function<double(double)> function,
                  std::function<double(double)> deriv_function) {
    activation_functions_map.insert(std::make_pair(function_name,
                                                   std::make_pair(function,
                                                                  deriv_function)));
  };

  ActivationFunctionsManager() {
    AddNewPair("sigmoid", sigmoid, deriv_sigmoid);
    AddNewPair("tanh", hyperbolic_tan, deriv_hyperbolic_tan);
    AddNewPair("linear", linear, deriv_linear);
  };

  std::unordered_map<std::string,
    std::pair<std::function<double(double)>, std::function<double(double)> > >
    activation_functions_map;
};

struct gen_rand {
  double factor;
  double offset;
public:
  gen_rand(double r = 2.0) : factor(r / RAND_MAX), offset(r / 2) {}
  double operator()() {
    return rand() * factor - offset;
  }
};

inline void Softmax(std::vector<double> *output) {
  size_t num_elements = output->size();
  std::vector<double> exp_output(num_elements);
  double exp_total = 0.0;
  for (size_t i = 0; i < num_elements; i++) {
    exp_output[i] = exp((*output)[i]);
    exp_total += exp_output[i];
  }
  for (size_t i = 0; i < num_elements; i++) {
    (*output)[i] = exp_output[i] / exp_total;
  }
}

inline void  GetIdMaxElement(const std::vector<double> &output, size_t * class_id) {
  *class_id = std::distance(output.begin(),
                            std::max_element(output.begin(),
                                             output.end()));
}
}

#define CONSTANT_WEIGHT_INITIALIZATION 0

class Node {
public:
  Node() {
    m_num_inputs = 0;
    m_bias = 0.0;
    m_weights.clear();
  };
  Node(int num_inputs,
       bool use_constant_weight_init = true,
       double constant_weight_init = 0.5) {
    m_num_inputs = num_inputs;
    m_bias = 0.0;
    m_weights.clear();
    //initialize weight vector
    WeightInitialization(m_num_inputs,
                         use_constant_weight_init,
                         constant_weight_init);
  };

  ~Node() {
    m_num_inputs = 0;
    m_bias = 0.0;
    m_weights.clear();
  };

  void WeightInitialization(int num_inputs,
                            bool use_constant_weight_init = true,
                            double constant_weight_init = 0.5) {
    m_num_inputs = num_inputs;
    //initialize weight vector
    if (use_constant_weight_init) {
      m_weights.resize(m_num_inputs, constant_weight_init);
    } else {
      m_weights.resize(m_num_inputs);
      std::generate_n(m_weights.begin(),
                      m_num_inputs,
                      utils::gen_rand());
    }
  }

  int GetInputSize() const {
    return m_num_inputs;
  }

  void SetInputSize(int num_inputs) {
    m_num_inputs = num_inputs;
  }

  double GetBias() const {
    return m_bias;
  }
  void SetBias(double bias) {
    m_bias = bias;
  }

  std::vector<double> & GetWeights() {
    return m_weights;
  }

  const std::vector<double> & GetWeights() const {
    return m_weights;
  }

  void SetWeights( std::vector<double> & weights ){
      // check size of the weights vector
      if( weights.size() == m_num_inputs )
          m_weights = weights;
      else
          throw new std::logic_error("Incorrect weight size in SetWeights call");
  }

  size_t GetWeightsVectorSize() const {
    return m_weights.size();
  }

  void GetInputInnerProdWithWeights(const std::vector<double> &input,
                                    double * output) const {
    assert(input.size() == m_weights.size());
    double inner_prod = std::inner_product(begin(input),
                                           end(input),
                                           begin(m_weights),
                                           0.0);
    *output = inner_prod;
  }

  void GetOutputAfterActivationFunction(const std::vector<double> &input,
                                        std::function<double(double)> activation_function,
                                        double * output) const {
    double inner_prod = 0.0;
    GetInputInnerProdWithWeights(input, &inner_prod);
    *output = activation_function(inner_prod);
  }

  void GetBooleanOutput(const std::vector<double> &input,
                        std::function<double(double)> activation_function,
                        bool * bool_output,
                        double threshold = 0.5) const {
    double value;
    GetOutputAfterActivationFunction(input, activation_function, &value);
    *bool_output = (value > threshold) ? true : false;
  };

  void UpdateWeights(const std::vector<double> &x,
                     double error,
                     double learning_rate) {
    assert(x.size() == m_weights.size());
    for (size_t i = 0; i < m_weights.size(); i++)
      m_weights[i] += x[i] * learning_rate *  error;
  };

  void UpdateWeight(int weight_id,
                    double increment,
                    double learning_rate) {
    m_weights[weight_id] += learning_rate*increment;
  }

  void SaveNode(FILE * file) const {
    fwrite(&m_num_inputs, sizeof(m_num_inputs), 1, file);
    fwrite(&m_bias, sizeof(m_bias), 1, file);
    fwrite(&m_weights[0], sizeof(m_weights[0]), m_weights.size(), file);
  };
  void LoadNode(FILE * file) {
    m_weights.clear();

    fread(&m_num_inputs, sizeof(m_num_inputs), 1, file);
    fread(&m_bias, sizeof(m_bias), 1, file);
    m_weights.resize(m_num_inputs);
    fread(&m_weights[0], sizeof(m_weights[0]), m_weights.size(), file);
  };

protected:
  size_t m_num_inputs{ 0 };
  double m_bias{ 0.0 };
  std::vector<double> m_weights;
};

class Layer {
public:
  Layer() {
    m_num_nodes = 0;
    m_nodes.clear();
  };

  Layer(int num_inputs_per_node,
        int num_nodes,
        const std::string & activation_function,
        bool use_constant_weight_init = true,
        double constant_weight_init = 0.5) {
    m_num_inputs_per_node = num_inputs_per_node;
    m_num_nodes = num_nodes;
    m_nodes.resize(num_nodes);

    for (int i = 0; i < num_nodes; i++) {
      m_nodes[i].WeightInitialization(num_inputs_per_node,
                                      use_constant_weight_init,
                                      constant_weight_init);
    }

    std::pair<std::function<double(double)>,
      std::function<double(double)> > *pair;
    bool ret_val = utils::ActivationFunctionsManager::Singleton().
      GetActivationFunctionPair(activation_function,
                                &pair);
    assert(ret_val);
    m_activation_function = (*pair).first;
    m_deriv_activation_function = (*pair).second;
    m_activation_function_str = activation_function;
  };

  ~Layer() {
    m_num_inputs_per_node = 0;
    m_num_nodes = 0;
    m_nodes.clear();
  };

  int GetInputSize() const {
    return m_num_inputs_per_node;
  };

  int GetOutputSize() const {
    return m_num_nodes;
  };

  const std::vector<Node> & GetNodes() const {
    return m_nodes;
  }

  /**
   * Return the internal list of nodes, but modifiable.
   */
  std::vector<Node> & GetNodesChangeable() {
    return m_nodes;
  }


  void GetOutputAfterActivationFunction(const std::vector<double> &input,
                                        std::vector<double> * output) const {
    assert(input.size() == m_num_inputs_per_node);

    output->resize(m_num_nodes);

    for (size_t i = 0; i < m_num_nodes; ++i) {
      m_nodes[i].GetOutputAfterActivationFunction(input,
                                                  m_activation_function,
                                                  &((*output)[i]));
    }
  }

  void UpdateWeights(const std::vector<double> &input_layer_activation,
                     const std::vector<double> &deriv_error,
                     double m_learning_rate,
                     std::vector<double> * deltas) {
    assert(input_layer_activation.size() == m_num_inputs_per_node);
    assert(deriv_error.size() == m_nodes.size());

    deltas->resize(m_num_inputs_per_node, 0);

    for (size_t i = 0; i < m_nodes.size(); i++) {
      double net_sum;
      m_nodes[i].GetInputInnerProdWithWeights(input_layer_activation,
                                              &net_sum);

      //dE/dwij = dE/doj . doj/dnetj . dnetj/dwij
      double dE_doj = 0.0;
      double doj_dnetj = 0.0;
      double dnetj_dwij = 0.0;

      dE_doj = deriv_error[i];
      doj_dnetj = m_deriv_activation_function(net_sum);

      for (size_t j = 0; j < m_num_inputs_per_node; j++) {
        (*deltas)[j] += dE_doj * doj_dnetj * m_nodes[i].GetWeights()[j];

        dnetj_dwij = input_layer_activation[j];

        m_nodes[i].UpdateWeight(j,
                                -(dE_doj * doj_dnetj * dnetj_dwij),
                                m_learning_rate);
      }
    }
  };


  void SetWeights( std::vector<std::vector<double>> & weights )
  {
      if( 0 <= weights.size() && weights.size() <= m_num_nodes )
      {
          // traverse the list of nodes
          size_t node_i = 0;
          for( Node & node : m_nodes )
          {
              node.SetWeights( weights[node_i] );
              node_i++;
          }
      }
      else
          throw new std::logic_error("Incorrect layer number in SetWeights call");
  };

  void SaveLayer(FILE * file) const {
    fwrite(&m_num_nodes, sizeof(m_num_nodes), 1, file);
    fwrite(&m_num_inputs_per_node, sizeof(m_num_inputs_per_node), 1, file);

    size_t str_size = m_activation_function_str.size();
    fwrite(&str_size, sizeof(size_t), 1, file);
    fwrite(m_activation_function_str.c_str(), sizeof(char), str_size, file);

    for (size_t i = 0; i < m_nodes.size(); i++) {
      m_nodes[i].SaveNode(file);
    }
  };
  void LoadLayer(FILE * file) {
    m_nodes.clear();

    fread(&m_num_nodes, sizeof(m_num_nodes), 1, file);
    fread(&m_num_inputs_per_node, sizeof(m_num_inputs_per_node), 1, file);

    size_t str_size = 0;
    fread(&str_size, sizeof(size_t), 1, file);
    m_activation_function_str.resize(str_size);
    fread(&(m_activation_function_str[0]), sizeof(char), str_size, file);

    std::pair<std::function<double(double)>,
      std::function<double(double)> > *pair;
    bool ret_val = utils::ActivationFunctionsManager::Singleton().
      GetActivationFunctionPair(m_activation_function_str,
                                &pair);
    assert(ret_val);
    m_activation_function = (*pair).first;
    m_deriv_activation_function = (*pair).second;
    
    m_nodes.resize(m_num_nodes);
    for (size_t i = 0; i < m_nodes.size(); i++) {
      m_nodes[i].LoadNode(file);
    }

  };

protected:
  size_t m_num_inputs_per_node{ 0 };
  size_t m_num_nodes{ 0 };
  std::vector<Node> m_nodes;

  std::string m_activation_function_str;
  std::function<double(double)>  m_activation_function;
  std::function<double(double)>  m_deriv_activation_function;
};

#include <iostream>
#include <stdlib.h>
#include <vector>

class Sample {
public:
  Sample(const std::vector<double> & input_vector) {

    m_input_vector = input_vector;
  }
  const std::vector<double> & input_vector() const {
    return m_input_vector;
  }
  size_t GetInputVectorSize() const {
    return m_input_vector.size();
  }
  void AddBiasValue(double bias_value) {
    m_input_vector.insert(m_input_vector.begin(), bias_value);
  }
  friend std::ostream & operator<<(std::ostream &stream, Sample const & obj) {
    obj.PrintMyself(stream);
    return stream;
  };
protected:
  virtual void PrintMyself(std::ostream& stream) const {
    stream << "Input vector: [";
    for (size_t i = 0; i < m_input_vector.size(); i++) {
      if (i != 0)
        stream << ", ";
      stream << m_input_vector[i];
    }
    stream << "]";
  }

  std::vector<double> m_input_vector;
};


class TrainingSample : public Sample {
public:
  TrainingSample(const std::vector<double> & input_vector,
                 const std::vector<double> & output_vector) :
    Sample(input_vector) {
    m_output_vector = output_vector;
  }
  const std::vector<double> & output_vector() const {
    return m_output_vector;
  }
  size_t GetOutputVectorSize() const {
    return m_output_vector.size();
  }

protected:
  virtual void PrintMyself(std::ostream& stream) const {
    stream << "Input vector: [";
    for (size_t i = 0; i < m_input_vector.size(); i++) {
      if (i != 0)
        stream << ", ";
      stream << m_input_vector[i];
    }
    stream << "]";

    stream << "; ";

    stream << "Output vector: [";
    for (size_t i = 0; i < m_output_vector.size(); i++) {
      if (i != 0)
        stream << ", ";
      stream << m_output_vector[i];
    }
    stream << "]";
  }

  std::vector<double> m_output_vector;
};

class MLP {
public:
  //desired call syntax :  MLP({64*64,20,4}, {"sigmoid", "linear"},
  MLP(const std::vector<uint64_t> & layers_nodes,
      const std::vector<std::string> & layers_activfuncs,
      bool use_constant_weight_init = false,
      double constant_weight_init = 0.5);
  MLP(const std::string & filename);
  ~MLP();

  void SaveMLPNetwork(const std::string & filename)const;
  void LoadMLPNetwork(const std::string & filename);

  void GetOutput(const std::vector<double> &input,
                 std::vector<double> * output,
                 std::vector<std::vector<double>> * all_layers_activations = nullptr) const;
  void GetOutputClass(const std::vector<double> &output, size_t * class_id) const;

  void Train(const std::vector<TrainingSample> &training_sample_set_with_bias,
                       double learning_rate,
                       int max_iterations = 5000,
                       double min_error_cost = 0.001,
                       bool output_log = true);
  size_t GetNumLayers();
  std::vector<std::vector<double>> GetLayerWeights( size_t layer_i );
  void SetLayerWeights( size_t layer_i, std::vector<std::vector<double>> & weights );

protected:
  void UpdateWeights(const std::vector<std::vector<double>> & all_layers_activations,
                     const std::vector<double> &error,
                     double learning_rate);
private:
  void CreateMLP(const std::vector<uint64_t> & layers_nodes,
                 const std::vector<std::string> & layers_activfuncs,
                 bool use_constant_weight_init,
                 double constant_weight_init = 0.5);
  size_t m_num_inputs{ 0 };
  int m_num_outputs{ 0 };
  int m_num_hidden_layers{ 0 };
  std::vector<uint64_t> m_layers_nodes;
  std::vector<Layer> m_layers;
};


//desired call sintax :  MLP({64*64,20,4}, {"sigmoid", "linear"},
MLP::MLP(const std::vector<uint64_t> & layers_nodes,
         const std::vector<std::string> & layers_activfuncs,
         bool use_constant_weight_init,
         double constant_weight_init) {
  assert(layers_nodes.size() >= 2);
  assert(layers_activfuncs.size() + 1 == layers_nodes.size());

  CreateMLP(layers_nodes,
            layers_activfuncs,
            use_constant_weight_init,
            constant_weight_init);
};

MLP::MLP(const std::string & filename) {
  LoadMLPNetwork(filename);
}

MLP::~MLP() {
  m_num_inputs = 0;
  m_num_outputs = 0;
  m_num_hidden_layers = 0;
  m_layers_nodes.clear();
  m_layers.clear();
};

void MLP::CreateMLP(const std::vector<uint64_t> & layers_nodes,
                    const std::vector<std::string> & layers_activfuncs,
                    bool use_constant_weight_init,
                    double constant_weight_init) {
  m_layers_nodes = layers_nodes;
  m_num_inputs = m_layers_nodes[0];
  m_num_outputs = m_layers_nodes[m_layers_nodes.size() - 1];
  m_num_hidden_layers = m_layers_nodes.size() - 2;

  for (size_t i = 0; i < m_layers_nodes.size() - 1; i++) {
    m_layers.emplace_back(Layer(m_layers_nodes[i],
                                m_layers_nodes[i + 1],
                                layers_activfuncs[i],
                                use_constant_weight_init,
                                constant_weight_init));
  }
};

void MLP::SaveMLPNetwork(const std::string & filename)const {
  FILE * file;
  file = fopen(filename.c_str(), "wb");
  fwrite(&m_num_inputs, sizeof(m_num_inputs), 1, file);
  fwrite(&m_num_outputs, sizeof(m_num_outputs), 1, file);
  fwrite(&m_num_hidden_layers, sizeof(m_num_hidden_layers), 1, file);
  if (!m_layers_nodes.empty())
    fwrite(&m_layers_nodes[0], sizeof(m_layers_nodes[0]), m_layers_nodes.size(), file);
  for (size_t i = 0; i < m_layers.size(); i++) {
    m_layers[i].SaveLayer(file);
  }
  fclose(file);
};
void MLP::LoadMLPNetwork(const std::string & filename) {
  m_layers_nodes.clear();
  m_layers.clear();

  FILE * file;
  file = fopen(filename.c_str(), "rb");
  fread(&m_num_inputs, sizeof(m_num_inputs), 1, file);
  fread(&m_num_outputs, sizeof(m_num_outputs), 1, file);
  fread(&m_num_hidden_layers, sizeof(m_num_hidden_layers), 1, file);
  m_layers_nodes.resize(m_num_hidden_layers + 2);
  if (!m_layers_nodes.empty())
    fread(&m_layers_nodes[0], sizeof(m_layers_nodes[0]), m_layers_nodes.size(), file);
  m_layers.resize(m_layers_nodes.size() - 1);
  for (size_t i = 0; i < m_layers.size(); i++) {
    m_layers[i].LoadLayer(file);
  }
  fclose(file);
};

void MLP::GetOutput(const std::vector<double> &input,
                    std::vector<double> * output,
                    std::vector<std::vector<double>> * all_layers_activations) const {
  assert(input.size() == m_num_inputs);
  int temp_size;
  if (m_num_hidden_layers == 0)
    temp_size = m_num_outputs;
  else
    temp_size = m_layers_nodes[1];

  std::vector<double> temp_in(m_num_inputs, 0.0);
  std::vector<double> temp_out(temp_size, 0.0);
  temp_in = input;

  for (size_t i = 0; i < m_layers.size(); ++i) {
    if (i > 0) {
      //Store this layer activation
      if (all_layers_activations != nullptr)
        all_layers_activations->emplace_back(std::move(temp_in));

      temp_in.clear();
      temp_in = temp_out;
      temp_out.clear();
      temp_out.resize(m_layers[i].GetOutputSize());
    }
    m_layers[i].GetOutputAfterActivationFunction(temp_in, &temp_out);
  }

  if (temp_out.size() > 1)
    utils::Softmax(&temp_out);
  *output = temp_out;

  //Add last layer activation
  if (all_layers_activations != nullptr)
    all_layers_activations->emplace_back(std::move(temp_in));
}

void MLP::GetOutputClass(const std::vector<double> &output, size_t * class_id) const {
  utils::GetIdMaxElement(output, class_id);
}

void MLP::UpdateWeights(const std::vector<std::vector<double>> & all_layers_activations,
                        const std::vector<double> &deriv_error,
                        double learning_rate) {

  std::vector<double> temp_deriv_error = deriv_error;
  std::vector<double> deltas{};
  //m_layers.size() equals (m_num_hidden_layers + 1)
  for (int i = m_num_hidden_layers; i >= 0; --i) {
    m_layers[i].UpdateWeights(all_layers_activations[i], temp_deriv_error, learning_rate, &deltas);
    if (i > 0) {
      temp_deriv_error.clear();
      temp_deriv_error = std::move(deltas);
      deltas.clear();
    }
  }
};

void MLP::Train(const std::vector<TrainingSample> &training_sample_set_with_bias,
                          double learning_rate,
                          int max_iterations,
                          double min_error_cost,
                          bool output_log) {
  //rlunaro.03/01/2019. the compiler says that these variables are unused
  //int num_examples = training_sample_set_with_bias.size();
  //int num_features = training_sample_set_with_bias[0].GetInputVectorSize();

  //{
  //  int layer_i = -1;
  //  int node_i = -1;
  //  std::cout << "Starting weights:" << std::endl;
  //  for (const auto & layer : m_layers) {
  //    layer_i++;
  //    node_i = -1;
  //    std::cout << "Layer " << layer_i << " :" << std::endl;
  //    for (const auto & node : layer.GetNodes()) {
  //      node_i++;
  //      std::cout << "\tNode " << node_i << " :\t";
  //      for (auto m_weightselement : node.GetWeights()) {
  //        std::cout << m_weightselement << "\t";
  //      }
  //      std::cout << std::endl;
  //    }
  //  }
  //}

  int i = 0;
  double current_iteration_cost_function = 0.0;

  for (i = 0; i < max_iterations; i++) {
    current_iteration_cost_function = 0.0;
    for (auto & training_sample_with_bias : training_sample_set_with_bias) {

      std::vector<double> predicted_output;
      std::vector< std::vector<double> > all_layers_activations;

      GetOutput(training_sample_with_bias.input_vector(),
                &predicted_output,
                &all_layers_activations);

      const std::vector<double> &  correct_output =
        training_sample_with_bias.output_vector();

      assert(correct_output.size() == predicted_output.size());
      std::vector<double> deriv_error_output(predicted_output.size());

      if (output_log && ((i % (max_iterations / 10)) == 0)) {
        std::stringstream temp_training;
        temp_training << training_sample_with_bias << "\t\t";

        temp_training << "Predicted output: [";
        for (size_t i = 0; i < predicted_output.size(); i++) {
          if (i != 0)
            temp_training << ", ";
          temp_training << predicted_output[i];
        }
        temp_training << "]";
	
        LOG(INFO) << temp_training.str();
	
      }

      for (size_t j = 0; j < predicted_output.size(); j++) {
        current_iteration_cost_function +=
          (std::pow)((correct_output[j] - predicted_output[j]), 2);
        deriv_error_output[j] =
          -2 * (correct_output[j] - predicted_output[j]);
      }

      UpdateWeights(all_layers_activations,
                    deriv_error_output,
                    learning_rate);
    }


    if (output_log && ((i % (max_iterations / 10)) == 0))
      LOG(INFO) << "Iteration " << i << " cost function f(error): "
      << current_iteration_cost_function;

    if (current_iteration_cost_function < min_error_cost)
      break;
  }
  LOG(INFO) << "Iteration " << i << " cost function f(error): "
    << current_iteration_cost_function;

  LOG(INFO) << "******************************";
  LOG(INFO) << "******* TRAINING ENDED *******";
  LOG(INFO) << "******* " << i << " iters *******";
  LOG(INFO) << "******************************";

  //{
  //  int layer_i = -1;
  //  int node_i = -1;
  //  std::cout << "Final weights:" << std::endl;
  //  for (const auto & layer : m_layers) {
  //    layer_i++;
  //    node_i = -1;
  //    std::cout << "Layer " << layer_i << " :" << std::endl;
  //    for (const auto & node : layer.GetNodes()) {
  //      node_i++;
  //      std::cout << "\tNode " << node_i << " :\t";
  //      for (auto m_weightselement : node.GetWeights()) {
  //        std::cout << m_weightselement << "\t";
  //      }
  //      std::cout << std::endl;
  //    }
  //  }
  //}
};


size_t MLP::GetNumLayers()
{
    return m_layers.size();
}

std::vector<std::vector<double>> MLP::GetLayerWeights( size_t layer_i )
{
    std::vector<std::vector<double>> ret_val;
    // check parameters
    if( 0 <= layer_i && layer_i < m_layers.size() )
    {
        Layer current_layer = m_layers[layer_i];
        for( Node & node : current_layer.GetNodesChangeable() )
        {
            ret_val.push_back( node.GetWeights() );
        }
        return ret_val;
    }
    else
        throw new std::logic_error("Incorrect layer number in GetLayerWeights call");

}

void MLP::SetLayerWeights( size_t layer_i, std::vector<std::vector<double>> & weights )
{
    // check parameters
    if( 0 <= layer_i && layer_i < m_layers.size() )
    {
        m_layers[layer_i].SetWeights( weights );
    }
    else
        throw new std::logic_error("Incorrect layer number in SetLayerWeights call");
}
