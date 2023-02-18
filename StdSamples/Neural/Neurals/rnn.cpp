#pragma once

#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <random>

using Matrix = Eigen::MatrixXd;

class RNN
{
public:
	RNN(unsigned lexicSize, unsigned hiddenSize, unsigned seqLength);
	void forward(std::vector<unsigned> &inputs);
	void backProp(std::vector<unsigned> &targets);
	void update();
    std::vector<unsigned> generate(unsigned seed, unsigned int iter);

	static double tanh(double x)
	{
		return std::tanh(x);
	}

	static double dtanh(double x)
	{
		double t = std::tanh(x);
		return 1 - t * t;
	}

	static double exp(double x)
	{
		return std::exp(x);
	}

    static double adagradInv(double x)
    {
        return 1 / std::sqrt(x + 1e-8);
    }
	static double clip(double x)
	{
		if (x > 5)
			return 5;
		else if (x < -5)
			return -5;
		else
			return x;
	}

private:
	void adagrad(Matrix& param, Matrix& dparam, Matrix& mem);

	// hyperparameters
	unsigned lexicSize = 26;
	unsigned hiddenSize = 100;  // size of hidden layer of neuron
	unsigned seqLength = 30;  // number of steps to unroll the RNN for
	double learningRate = 1e-1;

	// parameters
	Matrix Wxh;
	Matrix Whh;
	Matrix Why;
	Matrix bh;
	Matrix by;

	Matrix hprev;

	// gradients
	Matrix dWxh;
	Matrix dWhh;
	Matrix dWhy;
	Matrix dbh;
	Matrix dby;

	// rnn pipeline
	std::vector<Matrix> hs; // hiddenstate
	std::vector<Matrix> xs; // one hot vectors
	std::vector<Matrix> ys; // unormalized log prob for next
	std::vector<Matrix> ps; // prob for next

	// memory for adagrad
	Matrix mWxh;
	Matrix mWhh;
	Matrix mWhy;
	Matrix mbh;
	Matrix mby;

  std::random_device rd;
  std::mt19937 gen;

};
RNN::RNN(unsigned lexicSize, unsigned hiddenSize, unsigned seqLength)
	:lexicSize(lexicSize),
		hiddenSize(hiddenSize),
		seqLength(seqLength),
		learningRate(1e-1),
		Wxh(Matrix::Random(hiddenSize,lexicSize) * 0.01),
		Whh(Matrix::Random(hiddenSize, hiddenSize) * 0.01),
		Why(Matrix::Random(lexicSize, hiddenSize) * 0.01),
		bh(Matrix::Zero(hiddenSize, 1)),
		by(Matrix::Zero(lexicSize, 1)),
        hprev(Matrix::Zero(hiddenSize, 1)),
		dWxh(Matrix::Zero(hiddenSize, lexicSize)),
		dWhh(Matrix::Zero(hiddenSize, hiddenSize)),
		dWhy(Matrix::Zero(lexicSize, hiddenSize)),
		dbh(Matrix::Zero(hiddenSize, 1)),
		dby(Matrix::Zero(lexicSize, 1)),
		mWxh(Matrix::Zero(hiddenSize,lexicSize)),
		mWhh(Matrix::Zero(hiddenSize, hiddenSize)),
		mWhy(Matrix::Zero(lexicSize, hiddenSize)),
		mbh(Matrix::Zero(hiddenSize, 1)),
		mby(Matrix::Zero(lexicSize, 1)),
    gen(rd())
{
	std::cout << Wxh.unaryExpr(std::ptr_fun(RNN::tanh)) << std::endl;
	std::cout << std::endl;
	std::cout << bh << std::endl;
}

void RNN::forward(std::vector<unsigned> &inputs)
{
	hs.push_back(hprev);
	for (int i = 0; i < static_cast<int>(inputs.size()); ++i) {
        //Matrix x = Matrix::Zero(lexicSize, 1);
		xs.push_back(Matrix::Zero(lexicSize, 1));
		xs.back()(inputs[i], 0) = 1;

		hs.push_back((Wxh * xs.back() + Whh * hs.back() + bh).unaryExpr(std::ptr_fun(RNN::tanh)));

		ys.push_back(Why * hs.back() + by);

		Matrix ysExp = ys.back().unaryExpr(std::ptr_fun(RNN::exp));
		ps.push_back(ysExp / ysExp.sum());
	}

	// save hprev for next forward pass
	hprev = hs.back();
}

void RNN::backProp(std::vector<unsigned> &targets)
{
	dWxh.setZero();
	dWhh.setZero();
	dWhy.setZero();
	dbh.setZero();
	dby.setZero();

	Matrix dhnext = Matrix::Zero(hiddenSize, 1);
	for (int i = targets.size() - 1; i >= 0; --i) {
		Matrix dy = ps[i];
		dy(targets[i], 0) -= 1;
		dWhy += dy * hs[i + 1].transpose();
		dby += dy;
		Matrix dh = Why.transpose() * dy + dhnext;
		Matrix dhraw = hs[i + 1].unaryExpr(std::ptr_fun(RNN::dtanh)).cwiseProduct(dh);
		dbh += dhraw;
		dWxh += dhraw * xs[i].transpose();
		dWhh += dhraw * hs[i].transpose();
		dhnext = Whh.transpose() * dhraw;
	}

	dWxh.unaryExpr(std::ptr_fun(RNN::clip));
	dWhh.unaryExpr(std::ptr_fun(RNN::clip));
	dWhy.unaryExpr(std::ptr_fun(RNN::clip));
	dbh.unaryExpr(std::ptr_fun(RNN::clip));
	dby.unaryExpr(std::ptr_fun(RNN::clip));
}

void RNN::adagrad(Matrix& param, Matrix& dparam, Matrix& mem)
{
	mem += dparam.cwiseProduct(dparam);
  //std::cout << mem << std::endl;

  // alpha * dparam / sqrt(mem + 1e-8)
	param += (-learningRate * dparam).cwiseProduct(mem.unaryExpr(std::ptr_fun(RNN::adagradInv)));
  //std::cout << param << std::endl;
}

void RNN::update()
{
	adagrad(Wxh, dWxh, mWxh);
	adagrad(Whh, dWhh, mWhh);
	adagrad(Why, dWhy, mWhy);
	adagrad(bh, dbh, mbh);
	adagrad(by, dby, mby);

	xs.clear();
	hs.clear();
	ys.clear();
	ps.clear();
}

std::vector<unsigned> RNN::generate(unsigned seed, unsigned iter)
{
    std::vector<unsigned> out;
    Matrix x = Matrix::Zero(lexicSize, 1);
    x(seed, 0) = 1;
    unsigned iprev = seed;
    Matrix h = hprev;
    for (unsigned i = 0; i < iter; ++i) {
        h = Wxh * x + Whh * h + bh;
        h = h.unaryExpr(std::ptr_fun(RNN::tanh));
        Matrix y = Why * h + by;
        Matrix yexp = y.unaryExpr(std::ptr_fun(RNN::exp));
        Matrix p = yexp / yexp.sum();
        //if (i == 0)
        //    std::cout << p;
        //std::cout << std::endl;

        double *data = p.data();
        if (i == 0) {
          for (unsigned j = 0; j < lexicSize; ++j)
              std::cout << data[j] << std::endl;
          std::cout << std::endl;
        }
        std::discrete_distribution<> dd (data, data + lexicSize);
        unsigned ix = dd(gen);
        out.push_back(ix);
        //ix = np.random.choice(range(vocab_size), p=p.ravel())
        //x = np.zeros((vocab_size, 1))
        //x[ix] = 1
        //ixes.append(ix)
    }

    return out;
}

#include <iostream>
#include <fstream>
#include <set>
#include <map>


std::vector<char> load(std::ifstream &datafs, unsigned count, bool training)
{
    char *cinput = new char[count];
	datafs.read(cinput, count);
    //std::cout << datafs.gcount() << " chars read :\n";
    //std::cout << cinput << '\n';
    //std::cin.get();

    if (training && datafs.gcount() < count) {
  		datafs.clear();
		datafs.seekg(0, std::ios::beg);
    }

	std::vector<char> input;
    for (int i = 0; i < datafs.gcount(); ++i) {
       input.push_back(cinput[i]);
    }
    delete[] cinput;
	return input;
}

void addchars(std::set<char> &vocab, std::vector<char> &input)
{
	for (auto &c : input) {
		vocab.insert(c);
	}
}

void buildDicts(std::set<char> &vocab,
				std::vector<char> &intToChar,
				std::map<char, unsigned> &charToInt)
{
	intToChar.clear();
	charToInt.clear();
	unsigned i = 0;
	for (auto &c : vocab) {
		intToChar.push_back(c);
		charToInt.insert(std::pair<char,unsigned>(c, i));
		++i;
	}
}

int main(int argc, const char ** argv)
{
    const unsigned seqLen = 31;

	std::set<char> vocab;
	std::vector<char> intToChar;
	std::map<char, unsigned> charToInt;

	std::ifstream datafs("input.txt");
	while (!datafs.eof()) {
		std::vector<char> input = load(datafs, seqLen, false);
		addchars(vocab, input);
	}
	buildDicts(vocab, intToChar, charToInt);
	for (unsigned i = 0; i < intToChar.size(); ++i) {
		std::cout << i << " => " << intToChar[i] << std::endl;
		std::cout << intToChar[i] << " => " << charToInt[intToChar[i]] << std::endl;
	}

	RNN rnn(vocab.size(), 40, seqLen - 1);

    unsigned long long iter = 0;
	datafs.seekg(0, std::ios::beg);
	while (!datafs.eof()) {
		std::vector<char> seq = load(datafs, seqLen, true);
        if (seq.size() < seqLen)
            continue;

        if (iter % 500 == 0) {
            std::cout << "Iteration " << iter << std::endl;
            std::cout << std::string(seq.begin(), seq.end()) << std::endl;
            std::vector<unsigned> results = rnn.generate(0, 30);
            std::stringstream sstream;
            for (auto &c : results) {
            sstream << intToChar[c];
            }
            std::cout << sstream.str() << std::endl;
        }

        // Building inputs and targets
            std::vector<unsigned> inputs;
            for (auto &c: seq){
                inputs.push_back(charToInt[c]);
            }
        std::vector<unsigned> targets = inputs;
        inputs.pop_back();
        targets.erase(targets.begin());

        rnn.forward(inputs);
        rnn.backProp(targets);
        rnn.update();
        iter++;
	}

	return 0;
}