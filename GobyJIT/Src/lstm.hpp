/**
 * MIT License
 *
 * Copyright (c) 2018 Prabhsimran Singh
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <map>
#include <set>
#include <vector>

#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/string_file.hpp>
#include <boost/timer/timer.hpp>

#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
using namespace Eigen;


namespace F {

    inline float glorot_uniform(const int &n_in, const int &n_out) {
        return 2.0 / (n_in + n_out);
    }

    inline MatrixXf sigmoid(const MatrixXf &x) {
        return (1.0 + (-x).array().exp()).inverse().matrix();
    }

    inline MatrixXf sigmoid_derivative(const MatrixXf &x) {
        return x * (1 - x.array()).matrix();
    }

    inline MatrixXf relu(const MatrixXf &x) {
        return x.cwiseMax(0.0);
    }

    MatrixXf tanh(const MatrixXf &x) {
        return sigmoid(x);
    }

    MatrixXf log_softmax(const MatrixXf &x) {
        auto e = x.array().exp();
        return (e / e.sum()).log().matrix();
    }
}

namespace util {

/**
 * Load Dataset.
 * 
 * Loads a dataset by reading file into memory.
 * 
 * @param filename the name of the file containing data.
 * @returns the list of lists of tokens by line.
 */
void load_dataset(std::vector<std::vector<std::string>> &contents, const std::string &filename) {
    std::string line;

    boost::filesystem::ifstream file(filename);
    while (std::getline(file, line)) {
        std::vector<std::string> strs;
        boost::split(strs, line, boost::is_any_of(" "));
        contents.push_back(strs);
    }
}

/**
 * Generate Vocabs.
 * 
 * Generates token-to-id and id-to-token vocabularies using training corpus.
 * 
 * @param corpus the corpus containing sequences of tokens.
 * @returns the pair of token-to-id and id-to-token maps.
 */
void generate_vocabs(const std::vector<std::vector<std::string>> &corpus,
                     std::map<std::string, int> &word_to_idx,
                     std::map<int, std::string> &idx_to_word) {
    std::set<std::string> tokens;

    for (const auto &line : corpus) {
        for (const auto &word : line) {
            tokens.insert(word);
        }
    }

    int index = 0;
    for (const auto &token : tokens) {
        word_to_idx[token] = index;
        idx_to_word[index] = token;
        index++;
    }
}

class Trainer {
};
} // namespace util

namespace nn {

std::string SOS_TOKEN = "~";
std::string EOS_TOKEN = "#";

class Dense {
  private:
    MatrixXf W;
    RowVectorXf b;

  public:
    explicit Dense(const int &, const int &);

    MatrixXf operator()(const MatrixXf &);

    MatrixXf forward(const MatrixXf &);

    // MatrixXf backward(const MatrixXf &, const MatrixXf &);
};

Dense::Dense(const int &num_inputs, const int &num_outputs) {
    W = MatrixXf(num_inputs, num_outputs).setRandom() * F::glorot_uniform(num_inputs, num_outputs);
    b = RowVectorXf(num_outputs).setZero();
}

MatrixXf Dense::operator()(const MatrixXf &inputs) {
    return forward(inputs);
}

MatrixXf Dense::forward(const MatrixXf &inputs) {
    MatrixXf Z = inputs * W;
    Z += b.replicate(inputs.rows(), 1);
    return Z;
}

// MatrixXf Dense::backward(const MatrixXf &inputs, const MatrixXf &gradients) {
// }

class Embedding {
  private:
    MatrixXf embeddings;
    int n_tokens;
    int embedding_dim;

  public:
    explicit Embedding(const int &, const int &);

    inline MatrixXf operator()(const MatrixXf &);

    inline MatrixXf forward(const MatrixXf &);

    // MatrixXf backward(const MatrixXf &, const MatrixXf &);
};

Embedding::Embedding(const int &n_tokens, const int &embedding_dim) : n_tokens(n_tokens), embedding_dim(embedding_dim) {
    this->embeddings = MatrixXf(n_tokens, embedding_dim).setRandom();
}

inline MatrixXf Embedding::operator()(const MatrixXf &inputs) {
    return forward(inputs);
}

inline MatrixXf Embedding::forward(const MatrixXf &inputs) {
    // here inputs are one-hot encodings
    return inputs * embeddings;
}

// MatrixXf Embedding::backward(const MatrixXf &inputs, const MatrixXf &gradients) {

// }


struct LSTMState {
    // hidden state
    MatrixXf h;
    // cell state
    MatrixXf c;
};

class LSTMCell {
  private:
    int batch_size;
    int hidden_size;
    int embedding_dim;

    Dense i2h;
    Dense h2h;

  protected:
    LSTMState state;

    friend class LSTMNetwork;

  public:
    explicit LSTMCell(const int &, const int &, const int &);

    MatrixXf &operator()(const MatrixXf &);

    MatrixXf &forward(const MatrixXf &);

    // MatrixXf backward(const MatrixXf &, const MatrixXf &);
};

/**
 * LSTMCell Constructor.
 * 
 * @param hidden_size the size of the hidden state and cell state.
 * @param batch_size the size of batch used during training (for vectorization purposes).
 */
LSTMCell::LSTMCell(const int &hidden_size, const int &embedding_dim, const int &batch_size)
    : i2h(Dense(embedding_dim, 4 * hidden_size)), h2h(Dense(hidden_size, 4 * hidden_size)) {

    this->hidden_size = hidden_size;
    this->embedding_dim = embedding_dim;
    this->batch_size = batch_size;

    this->state = LSTMState{MatrixXf(batch_size, hidden_size).setRandom() * F::glorot_uniform(batch_size, hidden_size),
                            MatrixXf(batch_size, hidden_size).setRandom() * F::glorot_uniform(batch_size, hidden_size)};
}

MatrixXf &LSTMCell::operator()(const MatrixXf &xt) {
    return forward(xt);
}

/**
 * LSTMCell Forward Pass.
 * 
 * @param xt the input vector at time-step t.
 * @returns the next hidden state for input into next lstm layer.
 */
MatrixXf &LSTMCell::forward(const MatrixXf &xt) {
    // i2h + h2h = [it_pre, ft_pre, ot_pre, x_pre]
    MatrixXf preactivations = i2h(xt) + h2h(state.h);
    // all pre sigmoid gates chunk
    MatrixXf pre_sigmoid_chunk = preactivations.block(0, 0, batch_size, 3 * hidden_size);
    // compute sigmoid on gates chunk
    MatrixXf all_gates = F::sigmoid(pre_sigmoid_chunk);
    // compute c_in (x_transform) i.e. information vector
    MatrixXf x_pre = preactivations.block(0, 3 * hidden_size, batch_size, hidden_size);
    MatrixXf x_transform = F::tanh(x_pre);
    // single out all the gates
    MatrixXf it = all_gates.block(0, 0, batch_size, hidden_size);
    MatrixXf ft = all_gates.block(0, hidden_size, batch_size, hidden_size);
    MatrixXf ot = all_gates.block(0, 2 * hidden_size, batch_size, hidden_size);
    // update cell state
    MatrixXf c_forget = ft.cwiseProduct(state.c);
    MatrixXf c_input = it.cwiseProduct(x_transform);
    state.c = c_forget + c_input;
    // compute next hidden state
    MatrixXf c_transform = F::tanh(state.c);
    state.h = ot.cwiseProduct(c_transform);
    return state.h;
}

/**
 * LSTMCell Backward Pass.
 * 
 * @param inputs the input vector given at time-step t.
 * @param gradients the gradients from upper layers computed using chain rule.
 * @returns the output i.e. the hidden state at time-step t.
 */
// MatrixXf LSTMCell::backward(const MatrixXf &inputs, const MatrixXf &gradients) {
// }

class LSTMNetwork {
  private:
    size_t n_layers;
    int batch_size;
    int embedding_dim;
    int n_tokens;
    int hidden_size;

    Embedding embedding;
    std::vector<LSTMCell> layers;
    Dense out;

    std::vector<std::vector<LSTMState>> states;

    friend class util::Trainer;

  public:
    explicit LSTMNetwork(const size_t &, const int &, const int &, const int &, const int &);

    MatrixXf operator()(const MatrixXf &);

    MatrixXf forward(const MatrixXf &);

    void backward(const MatrixXf &);
};

LSTMNetwork::LSTMNetwork(const size_t &n_layers, const int &hidden_size, const int &n_tokens, const int &embedding_dim, const int &batch_size)
    : embedding(Embedding(n_tokens, embedding_dim)), out(Dense(hidden_size, n_tokens)) {

    this->n_layers = n_layers;
    this->hidden_size = hidden_size;
    this->n_tokens = n_tokens;
    this->embedding_dim = embedding_dim;
    this->batch_size = batch_size;

    // initial layer (embedding -> hidden)
    layers.push_back(LSTMCell(hidden_size, embedding_dim, batch_size));
    for (size_t i = 1; i < n_layers; i++) {
        // rest of the layers (hidden -> hidden)
        layers.push_back(LSTMCell(hidden_size, hidden_size, batch_size));
    }
}

MatrixXf LSTMNetwork::operator()(const MatrixXf &inputs) {
    return forward(inputs);
}

MatrixXf LSTMNetwork::forward(const MatrixXf &inputs) {
    std::vector<LSTMState> t_states;
    MatrixXf output = embedding(inputs);
    for (size_t i = 0; i < layers.size(); i++) {
        output = layers[i](output);
        t_states.push_back(layers[i].state);
    }
    states.push_back(t_states);
    MatrixXf logits = out(output);
    return F::log_softmax(logits);
}

// void LSTMNetwork::backward(const MatrixXf &loss_grad) {

//     states.clear();
// }

} // namespace nn