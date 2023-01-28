#pragma once

#include <memory>
#include <vector>
#include <string>
#include <tuple>

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>
#include <string>
#include <functional>
#include <numeric>

namespace tf_core
{
    enum class DiscretizationMethod
    {
        Forward = 0u,
        Backward = 1u,
        Tustin = 2u
    };

    class Complex
    {
        public:
            Complex(void) = default;
            Complex(const Complex & obj) = default;
            //Complex(const Complex && obj) = default;
            ~Complex(void) = default;
            Complex(const float real, const float img);

            //Complex operator=(const Complex & rhs) = default;

            Complex operator+(const Complex & rhs) const;
            Complex & operator+=(const Complex & rhs);
            Complex operator-(const Complex & rhs) const;
            Complex & operator-=(const Complex & rhs);
            Complex operator*(const Complex & rhs) const;
            Complex & operator*=(const Complex & rhs);
            Complex operator*(const float rhs) const;
            Complex & operator*=(const float rhs);
            Complex operator/(const Complex & rhs) const;
            Complex & operator/=(const Complex & rhs);
            Complex operator/(const float rhs) const;
            Complex & operator/=(const float rhs);

            Complex Power(const int power) const;
            Complex operator^(const int power) const;

            bool operator==(const Complex & rhs) const;
            bool operator!=(const Complex & rhs) const;

            Complex Conjugate(void) const;

            float Real(void) const;
            float Img(void) const;

            float Abs(void) const;
            float Phase(void) const;

        private:
            float real_ = 0.0f;
            float img_ = 0.0f;
    };


    class Polynomial
    {
        public:
            using Coefficients = std::vector<float>;

            Polynomial(void);
            Polynomial(const Polynomial & polynomial) = default;
            Polynomial(Polynomial && polynomial) = default;
            ~Polynomial(void) = default;
            Polynomial & operator=(const Polynomial & polynomial) = default;
            Polynomial & operator=(Polynomial && polynomial) = delete;

            explicit Polynomial(const Coefficients & coefficients);

            Polynomial operator+(const Polynomial & rhs) const;
            Polynomial & operator+=(const Polynomial & rhs);
            Polynomial operator-(const Polynomial & rhs) const;
            Polynomial & operator-=(const Polynomial & rhs);
            Polynomial operator*(const Polynomial & rhs) const;
            Polynomial & operator*=(const Polynomial & rhs);
            Polynomial operator*(const float & rhs) const;
            Polynomial & operator*=(const float rhs);
            Polynomial operator/(const float rhs) const;
            Polynomial & operator/=(const float rhs);

            bool operator==(const Polynomial & rhs) const;

            friend Polynomial operator*=(const float lhs, const Polynomial & rhs) {
                auto mult = rhs * lhs;
                return mult;
            }

            Polynomial Power(const size_t n) const;

            const Coefficients & GetCoefficients(void) const;
            size_t GetPolynomialOrder(void) const;

            bool IsStable(void) const;

            std::string ToString(void) const;

        private:
            void ReduceUnnecessaryElements(void);
            Eigen::MatrixXf GetSubHurwitzianMatrixWithOrder(const size_t order) const;

            Coefficients coefficients_;
    };

     class PolynomialFraction
    {
        public:
            PolynomialFraction(const Polynomial & num, const Polynomial & den);
            ~PolynomialFraction(void) = default;

            PolynomialFraction operator+(const PolynomialFraction & rhs) const;
            PolynomialFraction & operator+=(const PolynomialFraction & rhs);
            PolynomialFraction operator-(const PolynomialFraction & rhs) const;
            PolynomialFraction & operator-=(const PolynomialFraction & rhs);
            PolynomialFraction operator*(const float & rhs);
            bool operator==(const PolynomialFraction & rhs) const;
            bool operator!=(const PolynomialFraction & rhs) const;

            PolynomialFraction Power(const size_t n);
            
            const Polynomial & GetNum(void) const;
            const Polynomial & GetDen(void) const;

        private:
            Polynomial num_;
            Polynomial den_;
    };
    class PolynomialDiscretizer
    {
        public:
            PolynomialDiscretizer(void) = default;
            ~PolynomialDiscretizer(void) = default;

            using DiscretizedPolynomial = std::tuple<Polynomial, Polynomial>;

            static DiscretizedPolynomial Discretize(const Polynomial & continuous,
                const float discretization_time,
                const DiscretizationMethod discretization_method);

        private:
            static PolynomialFraction GetApproximation(const float discretization_time,
                const DiscretizationMethod discretization_method);
    };
    class CoreTransferFunction
    {
        public:
            CoreTransferFunction(void) = delete;
            CoreTransferFunction(const CoreTransferFunction & obj) = default;
            CoreTransferFunction(const Polynomial & num, const Polynomial & den);
            ~CoreTransferFunction(void) = default;

            CoreTransferFunction Series(const CoreTransferFunction & rhs) const;
            CoreTransferFunction Parallel(const CoreTransferFunction & rhs) const;
            CoreTransferFunction Feedback(const CoreTransferFunction & feedback_branch, bool positive = false) const;

            CoreTransferFunction operator+(const CoreTransferFunction & rhs) const;
            CoreTransferFunction operator*(const CoreTransferFunction & rhs) const;

            bool operator==(const CoreTransferFunction & rhs) const;
            bool operator!=(const CoreTransferFunction & rhs) const;

            const Polynomial & GetNum(void) const;
            const Polynomial & GetDen(void) const;
            
        private:
            Polynomial num_;
            Polynomial den_;            
    };    
    
    class ComplexTransferFunction;
    class CoreTransferFunction;
    class Polynomial;

    using Signal = std::vector<float>;
    using FrequencyCharacteristic = std::vector<std::pair<float, float>>;

    class TransferFunction
    {
        public:
            using CoefficientsVector = std::vector<float>;

            TransferFunction(const CoefficientsVector & num, const CoefficientsVector & den);
            TransferFunction(void) = delete;
            ~TransferFunction(void);

            TransferFunction Series(const TransferFunction & rhs) const;
            TransferFunction Parallel(const TransferFunction & rhs) const;
            TransferFunction Feedback(const TransferFunction & feedback_branch, bool positive = false) const;

            TransferFunction operator+(const TransferFunction & rhs) const;
            TransferFunction operator*(const TransferFunction & rhs) const;

            bool operator==(const TransferFunction & rhs) const;
            bool operator!=(const TransferFunction & rhs) const;

            const CoefficientsVector & GetNum(void) const;
            const CoefficientsVector & GetDen(void) const;

            TransferFunction Discretize(const float discretization_time,
                const DiscretizationMethod discretization_method = DiscretizationMethod::Tustin) const;
            
            Signal SimulateDiscrete(const Signal & input_signal, const float sampling_time) const;
            Signal Step(float simulation_time) const;
            Signal Impulse(float simulation_time) const;

            bool IsStable(void) const;

            FrequencyCharacteristic Bode(const std::vector<float> & omeha) const;
            FrequencyCharacteristic Nyquist(const std::vector<float> & omeha) const;

            std::string ToString(void) const;

        private:
            explicit TransferFunction(const CoreTransferFunction & core_tf);

            std::unique_ptr<Polynomial> num_;
            std::unique_ptr<Polynomial> den_;
            std::unique_ptr<CoreTransferFunction> tf_;
            std::unique_ptr<ComplexTransferFunction> ctf_;
    };

     using Frequencies = std::vector<float>;
    using Characteristic = std::vector<std::pair<float, float>>;
    using ComplexCharacteristic = std::vector<Complex>;

    class ComplexTransferFunction
    {
        public:
            explicit ComplexTransferFunction(const CoreTransferFunction & tf);
            ComplexTransferFunction(void) = delete;
            ~ComplexTransferFunction(void) = default;

            Characteristic CalculateBode(const Frequencies & omega) const;
            Characteristic CalculateNyquist(const Frequencies & omega) const;

        private:
            ComplexCharacteristic CalculateCharacteristics(const Frequencies & omega) const;
            Complex CalculateValue(const float w) const;

            CoreTransferFunction tf_;
    };

    using Signal = std::vector<float>;

    class DiscreteSimulator
    {
        public:
            DiscreteSimulator(void) = delete;
            DiscreteSimulator(const DiscreteSimulator & obbj) = delete;
            DiscreteSimulator(const DiscreteSimulator && obj) = delete;
            ~DiscreteSimulator(void) = delete;

            static Signal Simulate(const CoreTransferFunction & tf, const Signal & input_signal);        
    };

    class TransferFunctionDiscretizer
    {
        public:
            TransferFunctionDiscretizer() = delete;
            TransferFunctionDiscretizer(const TransferFunctionDiscretizer & obj) = delete;
            TransferFunctionDiscretizer(const TransferFunctionDiscretizer && obj) = delete;
            ~TransferFunctionDiscretizer(void) = default;

            static CoreTransferFunction Discretize(const CoreTransferFunction & tf,
                const float discretization_time,
                const DiscretizationMethod discretization_method);
    };

    Polynomial::Polynomial(void)
        : coefficients_{{0.0f}} {
    }

    Polynomial::Polynomial(const Polynomial::Coefficients & coefficients)
        : coefficients_{coefficients} {
        ReduceUnnecessaryElements();
    }

    Polynomial Polynomial::operator+(const Polynomial & rhs) const {
        auto output_size = std::max(this->GetPolynomialOrder(), rhs.GetPolynomialOrder());
        Coefficients new_coefficients(output_size);
        std::copy(this->coefficients_.begin(), this->coefficients_.end(), new_coefficients.begin());

        for (auto idx = 0u; idx < rhs.GetPolynomialOrder(); idx++)
            new_coefficients.at(idx) += rhs.coefficients_.at(idx);

        return Polynomial(new_coefficients);
    }

    Polynomial & Polynomial::operator+=(const Polynomial & rhs) {
        auto sum = (*this) + rhs;
        this->coefficients_ = sum.GetCoefficients();
        ReduceUnnecessaryElements();
        return *this;
    }

    Polynomial Polynomial::operator-(const Polynomial & rhs) const {
        auto output_size = std::max(this->GetPolynomialOrder(), rhs.GetPolynomialOrder());
        Coefficients new_coefficients(output_size);
        std::copy(this->coefficients_.begin(), this->coefficients_.end(), new_coefficients.begin());

        for (auto idx = 0u; idx < rhs.GetPolynomialOrder(); idx++)
            new_coefficients.at(idx) -= rhs.coefficients_.at(idx);

        return Polynomial(new_coefficients);
    }

    Polynomial & Polynomial::operator-=(const Polynomial & rhs) {
        auto sum = (*this) - rhs;
        this->coefficients_ = sum.GetCoefficients();
        ReduceUnnecessaryElements();
        return *this;
    }

    Polynomial Polynomial::operator*(const Polynomial & rhs) const {
        auto new_polynomial_size = this->GetPolynomialOrder() + rhs.GetPolynomialOrder() - 1u;
        Coefficients new_coeficients(new_polynomial_size, 0.0f);

        for (auto a_idx = 0u; a_idx < this->GetPolynomialOrder(); a_idx++) {
            for (auto b_idx = 0u; b_idx < rhs.GetPolynomialOrder(); b_idx++) {
                new_coeficients.at(a_idx + b_idx) += this->coefficients_.at(a_idx) * rhs.coefficients_.at(b_idx);
            }
        }

        return Polynomial(new_coeficients);
    }

    Polynomial & Polynomial::operator*=(const Polynomial & rhs) {
        auto mult = (*this) * rhs;
        this->coefficients_ = mult.coefficients_;
        ReduceUnnecessaryElements();
        return *this;
    }

    Polynomial Polynomial::operator*(const float & rhs) const {
        Coefficients new_coeficients = this->coefficients_;
        std::transform(new_coeficients.begin(), new_coeficients.end(), new_coeficients.begin(),
            [=](auto coefficient) {
                return coefficient * rhs;
            }
        );

        return Polynomial(new_coeficients);
    }

    Polynomial & Polynomial::operator*=(const float rhs) {
        auto mult = (*this) * rhs;
        this->coefficients_ = mult.coefficients_;
        ReduceUnnecessaryElements();

        return (*this);
    }

    Polynomial Polynomial::operator/(const float rhs) const {
        if (std::abs(rhs) < 1e-9f)
            throw std::invalid_argument("Polynomial::operator/(const float & rhs) ; division by something close to 0.0f");
        
        Coefficients new_coeficients = this->coefficients_;
        std::transform(new_coeficients.begin(), new_coeficients.end(), new_coeficients.begin(),
            [=](auto coefficient) {
                return coefficient / rhs;
            }
        );

        return Polynomial(new_coeficients);
    }

    Polynomial & Polynomial::operator/=(const float rhs) {
        auto div = (*this) / rhs;
        this->coefficients_ = div.coefficients_;
        ReduceUnnecessaryElements();

        return (*this);
    }

    bool Polynomial::operator==(const Polynomial & rhs) const {
        return this->coefficients_ == rhs.coefficients_;
    }

    Polynomial Polynomial::Power(const size_t n) const {
        auto output = Polynomial({1.0f});

        for (auto i = 0u; i < n; i++)
            output *= (*this);

        return output;
    }

    const Polynomial::Coefficients & Polynomial::GetCoefficients(void) const {
        return coefficients_;
    }

    size_t Polynomial::GetPolynomialOrder(void) const {
        return coefficients_.size();
    }

    void Polynomial::ReduceUnnecessaryElements(void) {
        auto unnecessery_elements_number = 0u;

        for (auto it = coefficients_.rbegin(); it != coefficients_.rend(); it++) {
            if (*it == 0.0f)
                unnecessery_elements_number++;
            else
                break;
        }

        if ((unnecessery_elements_number > 0u) && (unnecessery_elements_number < coefficients_.size()))
        {
            coefficients_.resize(coefficients_.size() - unnecessery_elements_number);
        }
        else if (unnecessery_elements_number == coefficients_.size())
        {
            coefficients_ = {0.0f};
        }
        else
        {
            // Do nothing
        }
    }

    bool Polynomial::IsStable(void) const {
        auto first_negative = std::find_if(coefficients_.begin(), coefficients_.end(), [](float coefficient){ return coefficient <= 0.0f; });
        if (first_negative != coefficients_.end()) {
            return false;
        } else {
            for (auto index = 0u; index < GetPolynomialOrder(); index++) {
                auto sub_hurwitz = GetSubHurwitzianMatrixWithOrder(index);
                auto determinant = sub_hurwitz.determinant();
                if (determinant <= 0.0f)
                    return false;
            }
            return true;
        }
    }

    Eigen::MatrixXf Polynomial::GetSubHurwitzianMatrixWithOrder(const size_t order) const {
        auto sub_hurwitz = Eigen::MatrixXf(order, order);
        
        auto get_element = [this](int n) { return ((n >= 0) && (n < static_cast<int>(GetPolynomialOrder()))) ? coefficients_.at(n) : 0.0f; };

        for (auto row = 0; row < static_cast<int>(order); row++) {
            auto index = static_cast<int>(GetPolynomialOrder()) - 1 - (2 * row);
            for (auto col = 0; col < static_cast<int>(order); col++) {
                sub_hurwitz(col, row) = get_element(index + col);
            }
        }

        return sub_hurwitz;
    }

    std::string Polynomial::ToString(void) const {
        std::string output = "";
        auto idx = 0u;
        std::for_each(coefficients_.rbegin(), coefficients_.rend(),
            [&](float coeficient) {
                auto sign = (coeficient < 0.0f) ? std::string(" - ") : std::string(" + ");
                auto abs_coefficient = std::to_string(std::abs(coeficient));
                auto power = coefficients_.size() - 1u - idx;

                if (idx != 0u)
                {
                    output.append(sign);
                    output.append(abs_coefficient);
                }
                else
                {
                    output.append(std::to_string(coeficient));
                }                

                switch (power)
                {
                    case 0u:
                    {
                        break;
                    }
                    case 1:
                    {
                        output.append(std::string(" s"));
                        break;
                    }
                    default:
                    {
                        output.append(std::string(" s^")).append(std::to_string(power));
                    }
                }
                idx++;
            }
        );

        return output;
    }
    
        PolynomialDiscretizer::DiscretizedPolynomial PolynomialDiscretizer::Discretize(const Polynomial & continuous,
        const float discretization_time,
        const DiscretizationMethod discretization_method) {
        auto discrete_polynomial = PolynomialFraction(Polynomial({0.0f}), Polynomial({1.0f}));
        auto s = GetApproximation(discretization_time, discretization_method);

        for (auto idx = 0u; idx < continuous.GetPolynomialOrder(); idx++)
            discrete_polynomial += s.Power(idx) * continuous.GetCoefficients().at(idx);
        
        return std::make_tuple(discrete_polynomial.GetNum(), discrete_polynomial.GetDen());
    }

    PolynomialFraction PolynomialDiscretizer::GetApproximation(const float discretization_time,
        const DiscretizationMethod discretization_method) {
        switch (discretization_method) {
            case DiscretizationMethod::Forward:
                return PolynomialFraction(Polynomial({-1.0f, 1.0f}), Polynomial({discretization_time}));
            case DiscretizationMethod::Backward:
                return PolynomialFraction(Polynomial({-1.0f, 1.0f}), Polynomial({0.0f, discretization_time}));
            case DiscretizationMethod::Tustin:
                return PolynomialFraction(Polynomial({-1.0f, 1.0f}), Polynomial({1.0f, 1.0f})) * (2.0f / discretization_time);
            default:
                throw std::invalid_argument("PolynomialDiscretizer::GetApproximation - Invali discretization method");
        }
    }

     PolynomialFraction::PolynomialFraction(const Polynomial & num, const Polynomial & den)
        : num_{num}
        , den_{den} {
    }

    PolynomialFraction PolynomialFraction::operator+(const PolynomialFraction & rhs) const {
        auto new_den = this->den_ * rhs.den_;
        auto new_num = this->num_ * rhs.den_ + this->den_ * rhs.num_;

        return PolynomialFraction(new_num, new_den);
    }

    PolynomialFraction & PolynomialFraction::operator+=(const PolynomialFraction & rhs) {
        auto sum = (*this) + rhs;
        
        this->num_ = sum.GetNum();
        this->den_ = sum.GetDen();

        return (*this);
    }

    PolynomialFraction PolynomialFraction::operator-(const PolynomialFraction & rhs) const {
        auto new_den = this->den_ * rhs.den_;
        auto new_num = this->num_ * rhs.den_ - this->den_ * rhs.num_;

        return PolynomialFraction(new_num, new_den);
    }

    PolynomialFraction & PolynomialFraction::operator-=(const PolynomialFraction & rhs) {
        auto diff = (*this) - rhs;

        this->num_ = diff.GetNum();
        this->den_ = diff.GetDen();

        return (*this);
    }

    PolynomialFraction PolynomialFraction::operator*(const float & rhs) {
        return PolynomialFraction(num_ * rhs, den_);
    }

    bool PolynomialFraction::operator==(const PolynomialFraction & rhs) const {
        return (this->num_ == rhs.num_) && (this->den_ == rhs.den_);
    }
    bool PolynomialFraction::operator!=(const PolynomialFraction & rhs) const {
        return !(*this == rhs);
    }

    PolynomialFraction PolynomialFraction::Power(const size_t n) {
        return PolynomialFraction(num_.Power(n), den_.Power(n));
    }

    const Polynomial & PolynomialFraction::GetNum(void) const {
        return num_;
    }

    const Polynomial & PolynomialFraction::GetDen(void) const {
        return den_;
    }

    CoreTransferFunction TransferFunctionDiscretizer::Discretize(const CoreTransferFunction & tf,
        const float discretization_time,
        const DiscretizationMethod discretization_method) {
        auto [num_num, num_den] = tf_core::PolynomialDiscretizer::Discretize(tf.GetNum(), discretization_time, discretization_method);
        auto [den_num, den_den] = tf_core::PolynomialDiscretizer::Discretize(tf.GetDen(), discretization_time, discretization_method);

        auto num = CoreTransferFunction(num_num, num_den);
        auto inverted_den = CoreTransferFunction(den_den, den_num);

        auto result = num * inverted_den;

        return CoreTransferFunction(result);
    }

    TransferFunction::TransferFunction(const std::vector<float> & num, const std::vector<float> & den)
        : num_{std::make_unique<Polynomial>(num)}
        , den_{std::make_unique<Polynomial>(den)}
        , tf_{std::make_unique<CoreTransferFunction>(*num_, *den_)}
        , ctf_{std::make_unique<ComplexTransferFunction>(*tf_)} {
        std::ignore = num;
        std::ignore = den;
    }

    TransferFunction::TransferFunction(const CoreTransferFunction & core_tf)
        : num_{std::make_unique<Polynomial>(core_tf.GetNum())}
        , den_{std::make_unique<Polynomial>(core_tf.GetDen())}
        , tf_{std::make_unique<CoreTransferFunction>(core_tf.GetNum(), core_tf.GetDen())}  { 
    }

    TransferFunction::~TransferFunction(void) = default;

    TransferFunction TransferFunction::Series(const TransferFunction & rhs) const {
        return TransferFunction((*tf_) * (*rhs.tf_));
    }

    TransferFunction TransferFunction::Parallel(const TransferFunction & rhs) const {
        return TransferFunction((*tf_) + (*rhs.tf_));
    }

    TransferFunction TransferFunction::Feedback(const TransferFunction & feedback_branch, bool positive) const {
        return TransferFunction(tf_->Feedback(*feedback_branch.tf_, positive));
    }

    TransferFunction TransferFunction::operator+(const TransferFunction & rhs) const {
        return TransferFunction((*tf_) + (*rhs.tf_));
    }

    TransferFunction TransferFunction::operator*(const TransferFunction & rhs) const {
        return TransferFunction((*tf_) * (*rhs.tf_));
    }

    bool TransferFunction::operator==(const TransferFunction & rhs) const {
        return (*tf_) == (*rhs.tf_);
    }

    bool TransferFunction::operator!=(const TransferFunction & rhs) const {
        return (*tf_) != (*rhs.tf_);
    }

    const TransferFunction::CoefficientsVector & TransferFunction::GetNum(void) const {
        return num_->GetCoefficients();
    }

    const TransferFunction::CoefficientsVector & TransferFunction::GetDen(void) const {
        return den_->GetCoefficients();
    }

    TransferFunction TransferFunction::Discretize(const float discretization_time, const DiscretizationMethod discretization_method) const {
        auto discrete_tf = tf_core::TransferFunctionDiscretizer::Discretize(*tf_, discretization_time, discretization_method);
        return TransferFunction(discrete_tf);
    }

    Signal TransferFunction::SimulateDiscrete(const Signal & input_signal, const float sampling_time) const {
        auto discrete_tf = Discretize(sampling_time);
        return DiscreteSimulator::Simulate(*discrete_tf.tf_, input_signal);
    }

    Signal TransferFunction::Step(float simulation_time) const {
        float characteristics_samples_number_ = 100.0f;

        auto sampling_time = simulation_time / 100.0f;
        auto discrete_tf = Discretize(sampling_time);
        
        auto input_signal = Signal(static_cast<size_t>(characteristics_samples_number_), 1.0f);

        return DiscreteSimulator::Simulate(*discrete_tf.tf_, input_signal);
    }

    Signal TransferFunction::Impulse(float simulation_time) const {
        float characteristics_samples_number_ = 100.0f;

        auto sampling_time = simulation_time / 100.0f;
        auto discrete_tf = Discretize(sampling_time);
        
        auto input_signal = Signal(static_cast<size_t>(characteristics_samples_number_), 0.0f);
        input_signal.at(0u) = 100.0f;
        
        return DiscreteSimulator::Simulate(*discrete_tf.tf_, input_signal);
    }

    FrequencyCharacteristic TransferFunction::Bode(const std::vector<float> & omega) const {
        return ctf_->CalculateBode(omega);
    }

    FrequencyCharacteristic TransferFunction::Nyquist(const std::vector<float> & omega) const {
        return ctf_->CalculateNyquist(omega);
    }

    bool TransferFunction::IsStable(void) const {
        return den_->IsStable();
    }

    std::string TransferFunction::ToString(void) const {
        auto num_string = num_->ToString();
        auto num_string_length = static_cast<int>(num_string.size());

        auto den_string = den_->ToString();
        auto den_string_length = static_cast<int>(den_string.size());
        
        auto length_difference = std::abs(den_string_length - num_string_length);
        auto offset_size = length_difference / 2u;
        auto line_size = std::max(num_string_length, den_string_length);

        std::string output("");
        auto offset = std::string(offset_size, ' ');
        auto line = std::string(line_size, '-');

        if (num_string_length < den_string_length) {
            output.append(offset + num_string + offset);
            output.append("\n");
            output.append(line);
            output.append("\n");
            output.append(den_string);
        } else {
            output.append(num_string);
            output.append("\n");
            output.append(line);
            output.append("\n");
            output.append(offset + den_string + offset);
        }

        return output;
    }

     Complex::Complex(const float real, const float img)
        : real_{real}
        , img_{img} {
    }

    Complex Complex::operator+(const Complex & rhs) const {
        return Complex(this->Real() + rhs.Real(), this->Img() + rhs.Img());
    }

    Complex & Complex::operator+=(const Complex & rhs) {
        *this = (*this) + rhs;
        return *this;
    }

    Complex Complex::operator-(const Complex & rhs) const {
        return Complex(this->Real() - rhs.Real(), this->Img() - rhs.Img());
    }

    Complex & Complex::operator-=(const Complex & rhs) {
        *this = (*this) - rhs;
        return *this;
    }

    Complex Complex::operator*(const Complex & rhs) const {
        auto new_real = this->Real() * rhs.Real() - this->Img() * rhs.Img();
        auto new_img = this->Real() * rhs.Img() + this->Img() * rhs.Real();
        return Complex(new_real, new_img);
    }

    Complex & Complex::operator*=(const Complex & rhs) {
        *this = (*this) * rhs;
        return *this;
    }

    Complex Complex::operator*(const float rhs) const {
        return Complex(this->Real() * rhs, this->Img() * rhs);
    }

    Complex & Complex::operator*=(const float rhs) {
        *this = (*this) * rhs;
        return *this;
    }

    Complex Complex::operator/(const Complex & rhs) const {
        return ((*this) * rhs.Conjugate()) / std::pow(rhs.Abs(), 2.0f);
    }

    Complex & Complex::operator/=(const Complex & rhs) {
        *this = (*this) / rhs;
        return *this;
    }

    Complex Complex::operator/(const float rhs) const {
        return Complex(this->Real() / rhs, this->Img() / rhs);
    }

    Complex & Complex::operator/=(const float rhs) {
        *this = (*this) / rhs;
        return *this;
    }

    Complex Complex::Power(const int power) const {
        auto output = Complex(1.0f, 0.0f);
        for (auto idx = 0; idx < power; idx++)
            output *= (*this);
        
        return output;
    }

    Complex Complex::operator^(const int power) const {
        return this->Power(power);
    }

    bool Complex::operator==(const Complex & rhs) const {
        return ((this->Real() == rhs.Real()) && (this->Img() == rhs.Img()));
    }

    bool Complex::operator!=(const Complex & rhs) const {
        return !((*this) == rhs);
    }

    Complex Complex::Conjugate(void) const {
        return Complex(Real(), -Img());
    }

    float Complex::Real(void) const {
        return real_;
    }

    float Complex::Img(void) const {
        return img_;
    }

    float Complex::Abs(void) const {
        return std::hypot(Real(), Img());
    }

    float Complex::Phase(void) const {
        return std::atan2(Img(), Real());
    }

    ComplexTransferFunction::ComplexTransferFunction(const CoreTransferFunction & tf)
        : tf_{tf} {
    }

    Characteristic ComplexTransferFunction::CalculateBode(const Frequencies & omega) const {
        auto complex_characteristic = CalculateCharacteristics(omega);

        Characteristic bode_characteristic(omega.size());
        std::transform(complex_characteristic.begin(), complex_characteristic.end(), bode_characteristic.begin(),
            [](const Complex v) {
                return std::pair<float, float>(v.Abs(), v.Phase());
            }
        );

        return bode_characteristic;
    }

    Characteristic ComplexTransferFunction::CalculateNyquist(const Frequencies & omega) const {
        auto complex_characteristic = CalculateCharacteristics(omega);

        Characteristic nyquist_characteristic(omega.size());
        std::transform(complex_characteristic.begin(), complex_characteristic.end(), nyquist_characteristic.begin(),
            [](const Complex v) {
                return std::pair<float, float>(v.Real(), v.Img());
            }
        );

        return nyquist_characteristic;
    }

    ComplexCharacteristic ComplexTransferFunction::CalculateCharacteristics(const Frequencies & omega) const {
        ComplexCharacteristic characteristic(omega.size());

        std::transform(omega.begin(), omega.end(), characteristic.begin(),
            [this](const float w) {
                return CalculateValue(w);
            }
        );

        return characteristic;
    }

    Complex  ComplexTransferFunction::CalculateValue(const float w) const {
        const auto omega = Complex(0.0f, w);

        auto init_num = Complex(0.0f, 0.0f);
        int pow_num = 0;
        auto num = std::accumulate(tf_.GetNum().GetCoefficients().begin(), tf_.GetNum().GetCoefficients().end(), init_num,
            [=,&pow_num](Complex sum, const float coefficient) {
                auto w_power = omega^pow_num;
                pow_num++;

                return (sum + w_power * coefficient);
            }
        );

        auto init_den = Complex(0.0f, 0.0f);
        int pow_den = 0;
        auto den = std::accumulate(tf_.GetDen().GetCoefficients().begin(), tf_.GetDen().GetCoefficients().end(), init_den,
            [=,&pow_den](Complex sum, const float coefficient) {
                auto w_power = omega^pow_den;
                pow_den++;

                return (sum + w_power * coefficient);
            }
        );

        return (num / den);
    }

    Signal DiscreteSimulator::Simulate(const CoreTransferFunction & tf, const Signal & input_signal) {
        auto input_history_size = tf.GetNum().GetPolynomialOrder();
        auto output_history_size = tf.GetDen().GetPolynomialOrder() - 1u;

        Signal input_history(input_history_size, 0.0f);
        Signal output_history(output_history_size, 0.0f);

        Signal output_signal(input_signal.size(), 0.0f);

        std::transform(input_signal.begin(), input_signal.end(), output_signal.begin(),
            [&](float input) {
                std::rotate(input_history.rbegin(), input_history.rbegin() + 1u, input_history.rend());
                *input_history.rbegin() = input;

                auto index = 0u;
                auto input_sum = std::accumulate(input_history.begin(), input_history.end(), 0.0f,
                    [=, &index] (float result, float input) {
                        return result + input * tf.GetNum().GetCoefficients().at(index++);
                    }
                );

                index = 0u;
                auto output_sum = std::accumulate(output_history.begin(), output_history.end(), 0.0f,
                    [=, &index] (float result, float output) {
                        return result + output * tf.GetDen().GetCoefficients().at(index++);
                    }
                );

                auto output = (input_sum - output_sum) / *(tf.GetDen().GetCoefficients().end() - 1u);

                std::rotate(output_history.rbegin(), output_history.rbegin() + 1u, output_history.rend());
                *output_history.rbegin() = output;

                return output;
            }
        );

        return output_signal;
    }

    CoreTransferFunction::CoreTransferFunction(const Polynomial & num, const Polynomial & den)
        : num_{num}
        , den_{den} {
    }

    CoreTransferFunction CoreTransferFunction::Series(const CoreTransferFunction & rhs) const {
        auto new_num = this->num_ * rhs.num_;
        auto new_den = this->den_ * rhs.den_;

        return CoreTransferFunction(new_num, new_den);
    }

    CoreTransferFunction CoreTransferFunction::Parallel(const CoreTransferFunction & rhs) const {
        auto new_num = this->num_ * rhs.den_ + rhs.num_ * this->den_;
        auto new_den = this->den_ * rhs.den_;
        
        return CoreTransferFunction(new_num, new_den);
    }

    CoreTransferFunction CoreTransferFunction::Feedback(const CoreTransferFunction & feedback_branch, bool positive) const {
        auto feedback_sign = positive ? -1.0f : 1.0f;

        auto new_num = this->num_ * feedback_branch.den_;
        auto new_den = this->den_ * feedback_branch.den_ + this->num_ * feedback_branch.num_ * feedback_sign;

        return CoreTransferFunction(new_num, new_den);
    }

    CoreTransferFunction CoreTransferFunction::operator+(const CoreTransferFunction & rhs) const {
        return this->Parallel(rhs);
    }

    CoreTransferFunction CoreTransferFunction::operator*(const CoreTransferFunction & rhs) const {
        return this->Series(rhs);
    }

    bool CoreTransferFunction::operator==(const CoreTransferFunction & rhs) const {
        return (this->num_ == rhs.num_) && (this->den_ == rhs.den_);
    }

    bool CoreTransferFunction::operator!=(const CoreTransferFunction & rhs) const {
        return !((*this) == rhs);
    }

    const Polynomial & CoreTransferFunction::GetNum(void) const {
        return num_;
    }

    const Polynomial & CoreTransferFunction::GetDen(void) const {
        return den_;
    }

}    