// Copyright 2020-present pytorch-cpp Authors

#include <torch/torch.h>
#include <iostream>
#include <iomanip>
#include <cassert>
#include <algorithm>
#include <functional>

typedef float floatType;

#include "Plot.hpp"
//#include "ml_helpers.hpp"

#include "Analog/VAAnalogSVF.hpp"
#include "Analog/VABlitOscillators.hpp"

#define MATRIXRC(w) std::cout << w.rows() << "," << w.cols() << std::endl
#define PRINT(X) std::cout << (X) << std::endl
#define LOOP(x,start,end) for(size_t x = start; x < end; x++)


Plot<floatType> plt;

class NeuralNetImpl : public torch::nn::Module {
 public:
    NeuralNetImpl(int64_t input_size, int64_t hidden_size, int64_t num_classes);

    torch::Tensor forward(torch::Tensor x);

 private:
    torch::nn::Linear fc1;
    torch::nn::Linear fc2;
    torch::nn::Linear fc3;
    
};

TORCH_MODULE(NeuralNet);



NeuralNetImpl::NeuralNetImpl(int64_t input_size, int64_t hidden_size, int64_t num_classes)
    : fc1(input_size, hidden_size), fc2(hidden_size, hidden_size), fc3(hidden_size,num_classes) {
    register_module("fc1", fc1);
    register_module("fc2", fc2);
    register_module("fc3", fc3);
    register_parameter("cutoff",torch::zeros({1,1}),true);
    register_parameter("q",torch::zeros({1,1}),true);
    
}

torch::Tensor NeuralNetImpl::forward(torch::Tensor x) {    
    x = fc1->forward(x);
    x = (torch::tanh(fc2->forward(x)));
    return fc3->forward(x);
}


struct DataSet : public torch::data::datasets::Dataset<DataSet>
{
    std::vector<floatType> data,train;
    Analog::Oscillators::Blit::BlitSaw osc;
	Analog::Filters::AnalogSVF::AnalogSVF filter;
	
    
    DataSet() : filter(44100,250,9.5) {
		for(int i = 0; i < 256; i++) {
			floatType s = osc.Tick();
			floatType x = filter.Tick(s);
			data.push_back(s);
			train.push_back(x);
		}        		

		plt.plot_x(data.data(),data.size(),"examples");
		plt.plot_x(train.data(),train.size(),"training");	
    }
    using Example = torch::data::Example<>;
    Example get(size_t i) {
        auto tdata = torch::tensor(torch::ArrayRef<floatType>(data.data(),data.size()));
        auto tlabel= torch::tensor(torch::ArrayRef<floatType>(train.data(),train.size()));
        return {tdata,tlabel};
    }
    torch::optional<size_t> size() const { return 1; }
};

/*
int main()
{
    auto t = XORDataSet();
    std::cout << t.get(0).data << std::endl;
}
*/

int main() {
    std::cout << "FeedForward Neural Network\n\n";

    // Device
    auto cuda_available = torch::cuda::is_available();
    torch::Device device(cuda_available ? torch::kCUDA : torch::kCPU);
    std::cout << (cuda_available ? "CUDA available. Training on GPU." : "Training on CPU.") << '\n';

    // Hyper parameters
    const int64_t input_size  = 256;
    const int64_t hidden_size = 256;
    const int64_t num_classes = 256;
    const int64_t batch_size = 1;
    const size_t num_epochs = 1000;
    const floatType learning_rate = 0.3;

	plt.setstyle("lines");
    
    //auto train_dataset = XORDataSet().map(torch::data::transforms::Stack<>());;
    auto train_dataset = DataSet().map(torch::data::transforms::Stack<>());;

    // Number of samples in the training set
    auto num_train_samples = train_dataset.size().value();

    auto test_dataset = DataSet().map(torch::data::transforms::Stack<>());

    // Number of samples in the testset
    auto num_test_samples = test_dataset.size().value();

		
    // Data loaders
    auto train_loader = torch::data::make_data_loader<torch::data::samplers::RandomSampler>(
        std::move(train_dataset), batch_size);

    auto test_loader = torch::data::make_data_loader<torch::data::samplers::SequentialSampler>(
        std::move(test_dataset), batch_size);

    // Neural Network model
    NeuralNet model(input_size, hidden_size, num_classes);
    model->to(device);

    // Optimizer
    torch::optim::SGD optimizer(model->parameters(), torch::optim::SGDOptions(learning_rate));

    // Set floating point output precision
    std::cout << std::fixed << std::setprecision(4);

    std::cout << "Training...\n";

    // Train the model
    for (size_t epoch = 0; epoch != num_epochs; ++epoch) {
        // Initialize running metrics
        floatType running_loss = 0.0;
        size_t num_correct = 0;

        for (auto& batch : *train_loader) {
            auto data = batch.data.to(device);
            auto target = batch.target.to(device);

            // Forward pass
            auto output = model->forward(data);
            filter.reset();
            filter.setCutoff(model->cutoff({0,0}).item<float>();
            filter.setQ(model->q({0,i}).item<float>();
            auto loss = torch::nn::functional::mse_loss(output, target);

            // Update running loss
            running_loss += loss.item<floatType>() * data.size(0);
            
            // Update number of correctly classified samples
            //num_correct += prediction.eq(target).sum().item<int64_t>();

            // Backward and optimize
            optimizer.zero_grad();
            loss.backward();
            optimizer.step();
        }

        auto sample_mean_loss = running_loss / num_train_samples;
        auto accuracy = static_cast<floatType>(num_correct) / num_train_samples;

        std::cout << "Epoch [" << (epoch + 1) << "/" << num_epochs << "], Trainset - Loss: "
            << sample_mean_loss << ", Accuracy: " << accuracy << '\n';
    }

    std::cout << "Training finished!\n\n";
    std::cout << "Testing...\n";

    // Test the model
    model->eval();
    torch::NoGradGuard no_grad;

    floatType running_loss = 0.0;
    size_t num_correct = 0;

    for (const auto& batch : *train_loader) {
        auto data = batch.data.view({batch_size, -1}).to(device);
        auto target = batch.target.to(device);

        auto output = model->forward(data);
        //std::cout << output << std::endl;
        auto loss = torch::nn::functional::mse_loss(output, target);

        running_loss += loss.item<floatType>() * data.size(0);

        std::vector<floatType> out(output.size(1));
        // Calculate prediction
        //auto prediction = output.argmax(1);                    
        for(int i = 0; i < output.size(1); i++)
            out[i] = output.index({0,i}).item<float>();
        plt.plot_x(out.data(),out.size(),"neural");
                    
    }

    std::cout << "Testing finished!\n";

    auto test_accuracy = static_cast<floatType>(num_correct) / num_test_samples;
    auto test_sample_mean_loss = running_loss / num_test_samples;

    std::cout << "Testset - Loss: " << test_sample_mean_loss << ", Accuracy: " << test_accuracy << '\n';
    sleep(10);
}

