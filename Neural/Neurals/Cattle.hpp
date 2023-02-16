#pragma once


#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <cstddef>
#include <string>
#include <cassert>
#include <Eigen/Dense>
#include <string>
#include <thread>
#include <unsupported/Eigen/CXX11/Tensor>
#include <memory>
#include <type_traits>
#include <utility>
#include <array>
#include <cassert>
#include <initializer_list>
#include <sstream>
#include <stdexcept>
#include <iomanip>
#include <set>
#include <array>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <sys/types.h>
#include <regex>
#include <cctype>
#include <dirent.h>
#include <map>
#include <pthread.h>
#include <functional>

namespace cattle {

    /**
    * An alias for a single row matrix of an arbitrary scalar type.
    */
    template<typename Scalar>
    using RowVector = Eigen::Matrix<Scalar,1,Eigen::Dynamic,Eigen::RowMajor, 1,Eigen::Dynamic>;

    /**
    * An alias for a single column matrix of an arbitrary scalar type.
    */
    template <typename Scalar>
    using ColVector = Eigen::Matrix<Scalar,Eigen::Dynamic,1,Eigen::ColMajor,Eigen::Dynamic,1>;

    /**
    * An alias for a dynamically sized matrix of an arbitrary scalar type.
    */
    template<typename Scalar>
    using Matrix = Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor,
            Eigen::Dynamic,Eigen::Dynamic>;

    /**
    * An alias for a class that can be used to map raw pointer data to a dynamically
    * sized Matrix of an arbitrary scalar type.
    */
    template<typename Scalar>
    using MatrixMap = Eigen::Map<Matrix<Scalar>>;

    /**
    * An alias for a tensor of arbitrary rank and scalar type with dynamic dimensionality.
    */
    template<typename Scalar, std::size_t Rank>
    using Tensor = Eigen::Tensor<Scalar,Rank,Eigen::ColMajor,std::size_t>;

    /**
    * An for a class that can be used to map raw pointer data to a tensor of arbitrary
    * rank and scalar type with dynamic dimensionality.
    */
    template<typename Scalar, std::size_t Rank>
    using TensorMap = Eigen::TensorMap<Tensor<Scalar,Rank>>;

    /**
    * An alias for permutation matrices.
    */
    using PermMatrix = Eigen::PermutationMatrix<Eigen::Dynamic,Eigen::Dynamic>;

    /**
    * An alias for self-adjoint eigen solvers.
    */
    template<typename Scalar>
    using EigenSolver = Eigen::SelfAdjointEigenSolver<Matrix<Scalar>>;

    /**
    * An alias for Eigen's bi-diagonal divide and conquer singular-value decomposition.
    */
    template<typename Scalar>
    using SVD = Eigen::BDCSVD<Matrix<Scalar>>;

    /**
    * An alias for Eigen's singular-value decomposition options.
    */
    using SVDOptions = Eigen::DecompositionOptions;

    /**
    * @return The number of threads used by Eigen to accelerate operations
    * supporting multithreading.
    */
    inline int num_of_eval_threads() {
        return Eigen::nbThreads();
    }

    /**
    * @param num_of_threads The number of threads Eigen should use to accelerate
    * operations supporting multithreading. The lower bound of the actual value
    * applied is 1 while the upper bound is the maximum of 1 and the level of
    * hardware concurrency detected.
    */
    inline void set_num_of_eval_threads(int num_of_threads) {
        int max = std::max(1, (int) std::thread::hardware_concurrency());
        Eigen::setNbThreads(std::max(1, std::min(num_of_threads, max)));
    }

    /**
    * It serializes the matrix in a format such that the first two numbers denote the
    * matrix's number of rows and columns respectively and the remaining numbers represent
    * the coefficients of the matrix in column-major order.
    *
    * @param matrix The matrix to serialize.
    * @param out_stream The non-binary stream to serialize the matrix to.
    */
    template<typename Scalar>
    inline void serialize(const Matrix<Scalar>& matrix, std::ostream& out_stream) {
        out_stream << sizeof(Scalar);
        out_stream << " " << matrix.rows();
        out_stream << " " << matrix.cols();
        for (std::size_t i = 0; i < matrix.size(); ++i)
            out_stream << " " << *(matrix.data() + i);
        out_stream << std::flush;
    }

    /**
    * It serializes the matrix into a file at the specified file path.
    *
    * @param matrix The matrix to serialize.
    * @param file_path The path to the file to which the matrix is to be serialized.
    */
    template<typename Scalar>
    inline void serialize(const Matrix<Scalar>& matrix, const std::string& file_path) {
        std::ofstream out_stream(file_path);
        assert(out_stream.is_open());
        serialize<Scalar>(matrix, out_stream);
    }

    /**
    * It serializes the matrix in a format such that the first 2 bytes denote the size of
    * a single coefficient of the matrix in bytes, the second and third 4 bytes denote the
    * matrix's number of rows and columns respectively, and the remaining bytes contain
    * the coefficients of the matrix in column-major order.
    *
    * @param matrix The matrix to serialize.
    * @param out_stream The binary stream to serialize the matrix to.
    */
    template<typename Scalar>
    inline void serialize_binary(const Matrix<Scalar>& matrix, std::ostream& out_stream) {
        unsigned short scalar_size = static_cast<unsigned short>(sizeof(Scalar));
        out_stream.write(reinterpret_cast<const char*>(&scalar_size),
                std::streamsize(sizeof(unsigned short)));
        unsigned rows = static_cast<unsigned>(matrix.rows());
        unsigned cols = static_cast<unsigned>(matrix.cols());
        out_stream.write(reinterpret_cast<const char*>(&rows), std::streamsize(sizeof(unsigned)));
        out_stream.write(reinterpret_cast<const char*>(&cols), std::streamsize(sizeof(unsigned)));
        out_stream.write(reinterpret_cast<const char*>(matrix.data()),
                std::streamsize(matrix.size() * sizeof(Scalar)));
        out_stream << std::flush;
    }

    /**
    * It serializes the matrix into a binary file at the specified file path.
    *
    * @param matrix The matrix to serialize.
    * @param file_path The path to the binary file to which the matrix is to be serialized.
    */
    template<typename Scalar>
    inline void serialize_binary(const Matrix<Scalar>& matrix, const std::string& file_path) {
        std::ofstream out_stream(file_path, std::ios::binary);
        assert(out_stream.is_open());
        serialize_binary<Scalar>(matrix, out_stream);
    }

    /**
    * It deserializes a matrix assuming the serialized format matches that used by the
    * serialize() method.
    *
    * @param in_stream The stream to the serialized matrix.
    * @return The unserialized matrix.
    */
    template<typename Scalar>
    inline Matrix<Scalar> deserialize(std::istream& in_stream) {
        unsigned rows, cols;
        in_stream >> rows;
        in_stream >> cols;
        Matrix<Scalar> matrix(rows, cols);
        for (std::size_t i = 0; i < matrix.size(); ++i)
            in_stream >> *(matrix.data() + i);
        return matrix;
    }

    /**
    * It deserializes a matrix from the file at the provided file path.
    *
    * @param file_path The path to the file containing the serialized matrix.
    * @return The deserialized matrix.
    */
    template<typename Scalar>
    inline Matrix<Scalar> deserialize(const std::string& file_path) {
        std::ifstream in_stream(file_path);
        assert(in_stream.is_open());
        return deserialize<Scalar>(in_stream);
    }

    /**
    * It deserializes a matrix assuming the serialized format matches that used by the
    * serialize_binary() method.
    *
    * @param in_stream The binary stream to the serialized matrix.
    * @return The unserialized matrix.
    */
    template<typename Scalar>
    inline Matrix<Scalar> deserialize_binary(std::istream& in_stream) {
        unsigned short scalar_size;
        in_stream.read(reinterpret_cast<char*>(&scalar_size), std::streamsize(sizeof(unsigned short)));
        assert(scalar_size == sizeof(Scalar));
        unsigned rows, cols;
        in_stream.read(reinterpret_cast<char*>(&rows), std::streamsize(sizeof(unsigned)));
        in_stream.read(reinterpret_cast<char*>(&cols), std::streamsize(sizeof(unsigned)));
        Matrix<Scalar> matrix(rows, cols);
        in_stream.read(reinterpret_cast<char*>(matrix.data()),
                std::streamsize(matrix.size() * sizeof(Scalar)));
        return matrix;
    }

    /**
    * It deserializes a matrix from the binary file at the provided file path.
    *
    * @param file_path The path to the binary file containing the serialized matrix.
    * @return The deserialized matrix.
    */
    template<typename Scalar>
    inline Matrix<Scalar> deserialize_binary(const std::string& file_path) {
        std::ifstream in_stream(file_path, std::ios::binary);
        assert(in_stream.is_open());
        return deserialize_binary<Scalar>(in_stream);
    }

    /**
    * An interface template for coder-decoders.
    */
    template<typename Scalar, std::size_t Rank>
    class Codec {
    public:
        virtual ~Codec() = default;
        /**
        * Encodes the tensor and writes it to a file.
        *
        * @param data The tensor whose contents are to be encoded.
        * @param file_path The path to the file to which the encoded data should be written.
        * If it does not exist, it will be created; if it does, the encoded data is appended
        * to the contents of the file.
        */
        virtual void encode(const Tensor<Scalar,Rank>& data, const std::string& file_path) const = 0;
        /**
        * Decodes the contents of a file into a tensor.
        *
        * @param file_path The path to the file containing the encoded data.
        * @return The decoded data in the form of a tensor.
        */
        virtual Tensor<Scalar,Rank> decode(const std::string& file_path) const = 0;
    };

    /**
    * An enumeration for PPM format types.
    */
    enum PPMFormatType {
        P2, P3, P5, P6
    };

    /**
    * A PPM image encoder-decoder.
    */
    template<typename Scalar, PPMFormatType Type = P6>
    class PPMCodec : public Codec<Scalar,3> {
        static_assert(Type >= P2 && Type <= P6, "illegal ppm format type argument");
        static constexpr int MAX_SINGLE_BYTE_VAL = 255;
        static constexpr int MAX_DOUBLE_BYTE_VAL = 65535;
        static constexpr int MAX_LINE_LENGTH = 70;
        static constexpr int MAX_VAL_STRING_LENGTH = 5;
        static constexpr int BUFFER_SIZE = 3072;
        static constexpr bool GRAY_SCALE = Type == P2 || Type == P5;
        static constexpr bool BINARY = Type == P5 || Type == P6;
    public:
        PPMCodec() : type(resolve_type_string()) { }
        inline void encode(const Tensor<Scalar,3>& data, const std::string& file_path) const {
            assert(data.dimension(0) > 0 && data.dimension(1) > 0 &&
                    data.dimension(2) == (GRAY_SCALE ? 1 : 3));
            std::ofstream file_stream(file_path, BINARY ? std::ios::binary : std::ios::out);
            assert(file_stream.is_open());
            Tensor<Scalar,0> max_tensor = data.maximum();
            const int max_val = (int) max_tensor(0u);
            assert(max_val >= 0 && max_val <= MAX_DOUBLE_BYTE_VAL);
            const bool single_byte = max_val <= MAX_SINGLE_BYTE_VAL;
            std::string header = type + "\n" + std::to_string(data.dimension(1)) + " " +
                    std::to_string(data.dimension(0)) + "\n" + std::to_string(max_val) + "\n";
            file_stream.write(header.c_str(), header.length());
            int ind = 0;
            // For non-binary formats.
            int last_line_break = 0;
            unsigned char buffer[+BUFFER_SIZE];
            for (std::size_t i = 0; i < data.dimension(0); ++i) {
                for (std::size_t j = 0; j < data.dimension(1); ++j) {
                    for (std::size_t k = 0; k < data.dimension(2); ++k) {
                        int val = data(i,j,k);
                        assert(val >= 0);
                        if (BINARY) {
                            if (ind == +BUFFER_SIZE) {
                                file_stream.write(reinterpret_cast<char*>(buffer), ind);
                                ind = 0;
                            }
                            // The buffer size is divisible by 2; no need to worry about buffer overflow.
                            if (!single_byte)
                                buffer[ind++] = (unsigned char) (val >> 8);
                            buffer[ind++] = (unsigned char) val;
                        } else {
                            if (ind >= +BUFFER_SIZE - (MAX_VAL_STRING_LENGTH + 1)) {
                                file_stream.write(reinterpret_cast<char*>(buffer), ind);
                                last_line_break -= ind;
                                ind = 0;
                            }
                            std::string val_string = std::to_string(val);
                            for (int l = 0; l < val_string.length(); ++l)
                                buffer[ind++] = *(val_string.c_str() + l);
                            if ((ind + 1 - last_line_break >= (MAX_LINE_LENGTH - MAX_VAL_STRING_LENGTH))) {
                                buffer[ind++] =  '\n';
                                last_line_break = ind;
                            } else
                                buffer[ind++] =  ' ';
                        }
                    }
                }
            }
            if (ind != 0)
                file_stream.write(reinterpret_cast<char*>(buffer), ind);
        }
        inline Tensor<Scalar,3> decode(const std::string& file_path) const {
            std::ifstream file_stream(file_path, BINARY ? std::ios::binary : std::ios::in);
            assert(file_stream.is_open());
            std::string format_type, dims, max_val_string;
            int width, height, max_val;
            std::getline(file_stream, format_type);
            assert(type == format_type);
            std::getline(file_stream, dims);
            std::istringstream dims_stream(dims);
            dims_stream >> width;
            dims_stream >> height;
            assert(width > 0 && height > 0);
            std::getline(file_stream, max_val_string);
            std::istringstream max_val_stream(max_val_string);
            max_val_stream >> max_val;
            assert(max_val >= 0 && max_val <= MAX_DOUBLE_BYTE_VAL);
            const std::size_t depth = GRAY_SCALE ? 1u : 3u;
            const int total_values = height * width * depth;
            Tensor<Scalar,3> data((std::size_t) height, (std::size_t) width, depth);
            unsigned char buffer[+BUFFER_SIZE];
            int ind = 0;
            if (BINARY) {
                const bool single_byte = max_val <= MAX_SINGLE_BYTE_VAL;
                int values_in_buffer = std::min(+BUFFER_SIZE, (2 - single_byte) * total_values);
                int read_values = values_in_buffer;
                file_stream.read(reinterpret_cast<char*>(&buffer), values_in_buffer);
                assert(file_stream.gcount() == values_in_buffer);
                for (std::size_t i = 0; i < height; ++i) {
                    for (std::size_t j = 0; j < width; ++j) {
                        for (std::size_t k = 0; k < depth; ++k) {
                            if (ind == values_in_buffer) {
                                values_in_buffer = std::min(+BUFFER_SIZE,
                                        (2 - single_byte) * total_values - read_values);
                                file_stream.read(reinterpret_cast<char*>(&buffer), values_in_buffer);
                                assert(file_stream.gcount() == values_in_buffer);
                                read_values += values_in_buffer;
                                ind = 0;
                            }
                            unsigned val;
                            if (single_byte)
                                val = (unsigned) buffer[ind++];
                            else { // No buffer overflow possible due to the even buffer size.
                                val = (unsigned) buffer[ind++];
                                val |= ((unsigned) buffer[ind++]) << 8;
                            }
                            data(i,j,k) = (Scalar) val;
                        }
                    }
                }
            } else {
                file_stream.read(reinterpret_cast<char*>(&buffer), +BUFFER_SIZE);
                for (std::size_t i = 0; i < height; ++i) {
                    for (std::size_t j = 0; j < width; ++j) {
                        for (std::size_t k = 0; k < depth; ++k) {
                            std::vector<unsigned char> chars;
                            bool found_num = false;
                            for (;;) {
                                if (ind == +BUFFER_SIZE) {
                                    file_stream.read(reinterpret_cast<char*>(&buffer), +BUFFER_SIZE);
                                    ind = 0;
                                } else if (ind == file_stream.gcount())
                                    break;
                                unsigned char curr_char = buffer[ind++];
                                if (curr_char >= '0' && curr_char <= '9') {
                                    chars.push_back(curr_char);
                                    found_num = true;
                                } else if (found_num)
                                    break;
                            }
                            assert(found_num);
                            std::string val_string(chars.begin(), chars.end());
                            data(i,j,k) = (Scalar) std::stoi(val_string);
                        }
                    }
                }
            }
            return data;
        }
    private:
        /**
        * @return The string representation of the PPM format type.
        */
        inline static std::string resolve_type_string() {
            switch (Type) {
                case P2: return "P2";
                case P3: return "P3";
                case P5: return "P5";
                case P6: return "P6";
                default: return "";
            }
        }
        const std::string type;
    };

/**
 * NOTE: Expression templates are rather unwarranted for this class as its objects are unlikely to
 * be used in complex expressions. They have only been implemented for the sake of practical learning.
 */

template<typename Derived, typename IndexType, std::size_t Rank>
class DimExpression;

/**
 * A class representing dimensions along one or more ranks. It can describe the dimensionality of
 * tensors of arbitrary ranks.
 */
template<typename IndexType, std::size_t Rank>
class Dimensions : public DimExpression<Dimensions<IndexType,Rank>,IndexType,Rank> {
	template<typename OtherIndexType, std::size_t OtherRank>
	friend class Dimensions;
public:
	inline Dimensions() :
			values(Rank, 1) { }
	inline Dimensions(const std::initializer_list<IndexType>& values) :
			Dimensions() {
		assert(values.size() <= Rank);
		std::copy(values.begin(), values.end(), this->values.begin());
	}
	inline Dimensions(const std::array<IndexType,Rank>& array) :
			Dimensions() {
		std::copy(array.begin(), array.end(), values.begin());
	}
	template<typename OtherDerived>
	inline Dimensions(const DimExpression<OtherDerived,IndexType,Rank>& dims) :
			Dimensions() {
		for (std::size_t i = 0; i < Rank; ++i)
			values[i] = dims(i);
	}
	/**
	 * A constant method that returns a copy of the instance with n additional ranks prepended to it.
	 *
	 * @return A new Dimensions instance with additional n ranks of size 1 prepended to it.
	 */
	template<std::size_t Ranks = 1>
	inline Dimensions<IndexType,Rank + Ranks> promote() const {
		Dimensions<IndexType,Rank + Ranks> promoted;
		std::copy(values.begin(), values.end(), promoted.values.begin() + Ranks);
		return promoted;
	}
	/**
	 * A constant method that returns a copy of the instance without the n most-significant ranks.
	 *
	 * @return A new Dimensions instance with the first n ranks removed.
	 */
	template<std::size_t Ranks = 1>
	inline Dimensions<IndexType,Rank - Ranks> demote() const {
		static_assert(Rank > Ranks, "rank must be greater than the number of ranks to demote by");
		Dimensions<IndexType,Rank - Ranks> demoted;
		std::copy(values.begin() + Ranks, values.end(), demoted.values.begin());
		return demoted;
	}
	/**
	 * A constant method that returns a copy of the instance with n ranks appended to it.
	 *
	 * @return A new Dimensions instance with additional n ranks with size 1 appended to it.
	 */
	template<std::size_t Ranks = 1>
	inline Dimensions<IndexType,Rank + Ranks> extend() const {
		Dimensions<IndexType,Rank + Ranks> extended;
		std::copy(values.begin(), values.end(), extended.values.begin());
		return extended;
	}
	/**
	 * A constant method that returns a copy of the instance without the n least-significant ranks.
	 *
	 * @return A new Dimensions instance with the last n ranks removed.
	 */
	template<std::size_t Ranks = 1>
	inline Dimensions<IndexType,Rank - Ranks> contract() const {
		static_assert(Rank > Ranks, "rank must be greater than the number of ranks to contract by");
		Dimensions<IndexType,Rank - Ranks> contracted;
		std::copy(values.begin(), values.end() - Ranks, contracted.values.begin());
		return contracted;
	}
	/**
	 * A simple constant method that returns a copy of the numeral representing the
	 * number of dimensions along a given rank.
	 *
	 * @param i The index of the rank whose dimensionality is to be returned.
	 * @return The dimensionality of the i-th rank.
	 */
	inline IndexType operator()(std::size_t i) const {
		if (i < 0 || i >= Rank)
			throw std::invalid_argument("illegal index value: " + std::to_string(i));
		return values[i];
	}
	/**
	 * A that returns a non-constant reference to the numeral representing the
	 * number of dimensions along a given rank.
	 *
	 * @param i The index of the rank whose dimensionality is to be returned.
	 * @return The dimensionality of the i-th rank.
	 */
	inline IndexType& operator()(std::size_t i) {
		if (i < 0 || i >= Rank)
			throw std::invalid_argument("illegal index value: " + std::to_string(i));
		return values[i];
	}
	/**
	 * A constant conversion operator that returns an array with the contents
	 * of the instance.
	 */
	inline operator std::array<IndexType,Rank>() const {
		std::array<IndexType,Rank> array;
		std::copy(values.begin(), values.end(), array.begin());
		return array;
	}
    private:
        std::vector<IndexType> values;
    };

    /**
    * An expression representing an operation between a numeral and all ranks of a dimension expression.
    */
    template<typename IndexType, std::size_t Rank, typename LhsExpr, typename OpType>
    class UnaryDimExpression :
            public DimExpression<UnaryDimExpression<IndexType,Rank,LhsExpr,OpType>,IndexType,Rank> {
    public:
        inline UnaryDimExpression(const LhsExpr& lhs, const IndexType& rhs) :
                lhs(lhs),
                rhs(rhs) { };
        inline IndexType operator()(std::size_t i) const {
            if (i < 0 || i >= Rank)
                throw std::invalid_argument("illegal index value: " + std::to_string(i));
            return OpType::apply(lhs(i), rhs);
        }
    private:
        const LhsExpr& lhs;
        IndexType rhs;
    };

    /**
    * An expression representing an operation between a numeral and a single rank of a dimension expression.
    */
    template<typename IndexType, std::size_t Rank, typename LhsExpr, typename OpType>
    class UnaryRankWiseDimExpression :
            public DimExpression<UnaryRankWiseDimExpression<IndexType,Rank,LhsExpr,OpType>,IndexType,Rank> {
    public:
        inline UnaryRankWiseDimExpression(const LhsExpr& lhs, const IndexType& rhs, std::size_t rank) :
                lhs(lhs),
                rhs(rhs),
                rank(rank) {
            assert(rank < Rank);
        };
        inline IndexType operator()(std::size_t i) const {
            if (i < 0 || i >= Rank)
                throw std::invalid_argument("illegal index value: " + std::to_string(i));
            return i == rank ? OpType::apply(lhs(i), rhs) : lhs(i);
        }
    private:
        const LhsExpr& lhs;
        IndexType rhs;
        std::size_t rank;
    };

    /**
    * An expression representing an operation between two dimension expressions of the same rank.
    */
    template<typename IndexType, std::size_t Rank, typename LhsExpr, typename RhsExpr, typename OpType>
    class BinaryDimExpression :
            public DimExpression<BinaryDimExpression<IndexType,Rank,LhsExpr,RhsExpr,OpType>,IndexType,Rank> {
    public:
        inline BinaryDimExpression(const LhsExpr& lhs, const RhsExpr& rhs) :
                lhs(lhs),
                rhs(rhs) { };
        inline IndexType operator()(std::size_t i) const {
            if (i < 0 || i >= Rank)
                throw std::invalid_argument("illegal index value: " + std::to_string(i));
            return OpType::apply(lhs(i), rhs(i));
        }
    protected:
        const LhsExpr& lhs;
        const RhsExpr& rhs;
    };

    /**
    * An expression representing an operation along a single rank of two dimension expressions of the
    * same rank.
    */
    template<typename IndexType, std::size_t Rank, typename LhsExpr, typename RhsExpr, typename OpType>
    class BinaryRankWiseDimExpression :
            public DimExpression<BinaryRankWiseDimExpression<IndexType,Rank,LhsExpr,RhsExpr,OpType>,IndexType,Rank> {
    public:
        inline BinaryRankWiseDimExpression(const LhsExpr& lhs, const RhsExpr& rhs, std::size_t rank) :
                lhs(lhs),
                rhs(rhs),
                rank(rank) {
            assert(rank < Rank);
        }
        inline IndexType operator()(std::size_t i) const {
            if (i < 0 || i >= Rank)
                throw std::invalid_argument("illegal index value: " + std::to_string(i));
            return i == rank ? OpType::apply(lhs(i), rhs(i)) : lhs(i);
        }
    protected:
        const LhsExpr& lhs;
        const RhsExpr& rhs;
        std::size_t rank;
    };

    // Arithmetic operations.
    template<typename Operand> class SumOp {
    public: inline static Operand apply(Operand lhs, Operand rhs) { return lhs + rhs; }
    };
    template<typename Operand> class SubOp {
    public: inline static Operand apply(Operand lhs, Operand rhs) { return lhs - rhs; }
    };
    template<typename Operand> class MulOp {
    public: inline static Operand apply(Operand lhs, Operand rhs) { return lhs * rhs; }
    };
    template<typename Operand> class DivOp {
    public: inline static Operand apply(Operand lhs, Operand rhs) { return lhs / rhs; }
    };

    /**
    * The base class of all dimension expressions.
    */
    template<typename Derived, typename IndexType, std::size_t Rank>
    class DimExpression {
        static_assert(Rank > 0, "illegal rank");
        typedef DimExpression<Derived,IndexType,Rank> Self;
        template<typename OtherDerived> using Other = DimExpression<OtherDerived,IndexType,Rank>;
    public:
        inline operator Derived&() {
            return static_cast<Derived&>(*this);
        }
        inline operator const Derived&() const {
            return static_cast<const Derived&>(*this);
        }
        /**
        * Evaluates the expression along the given rank.
        *
        * @param i The index of the single rank along which the expression is to be evaluated.
        * @return The result of the expression evaluated along the i-th rank.
        */
        inline IndexType operator()(std::size_t i) const {
            return static_cast<const Derived&>(*this)(i);
        }
        /**
        * A constant method the returns the volume of the expression.
        *
        * @return The product of the dimensions of each rank of the instance.
        */
        inline IndexType get_volume() const {
            int volume = 1;
            for (std::size_t i = 0; i < Rank; ++i)
                volume *= (*this)(i);
            return volume;
        }
        /**
        * A method for forcing the evaluation of the expression.
        *
        * @return A Dimensions instance containing the results of the evaluated
        * expression.
        */
        inline Dimensions<IndexType,Rank> eval() {
            return Dimensions<IndexType,Rank>(*this);
        }
        /**
        * It evaluates the expression and returns a string containing the results.
        *
        * @return A string representation of the evaluated expression.
        */
        inline std::string to_string() const {
            std::stringstream strm;
            strm << "[" + std::to_string((*this)(0));
            for (std::size_t i = 1; i < Rank; ++i)
                strm << "," << std::to_string((*this)(i));
            strm << "]";
            return strm.str();
        }
        /**
        * @param n The value to add.
        * @param rank The rank to which the value is to be added.
        * @return An expression representing the addition of the specified value to the
        * specified rank of the dimension expression.
        */
        inline UnaryRankWiseDimExpression<IndexType,Rank,Self,SumOp<IndexType>>
        add_along_rank(const IndexType& n, std::size_t rank) const {
            return UnaryRankWiseDimExpression<IndexType,Rank,Self,SumOp<IndexType>>(*this, n, rank);
        }
        /**
        * @param n The value to subtract.
        * @param rank The rank from which the value is to be subtracted.
        * @return An expression representing the subtraction of the specified value from the
        * specified rank of the dimension expression.
        */
        inline UnaryRankWiseDimExpression<IndexType,Rank,Self,SubOp<IndexType>>
        subtract_along_rank(const IndexType& n, std::size_t rank) const {
            return UnaryRankWiseDimExpression<IndexType,Rank,Self,SubOp<IndexType>>(*this, n, rank);
        }
        /**
        * @param n The value to multiply by.
        * @param rank The rank which is to be multiplied by the value.
        * @return An expression representing the multiplication of the specified rank of the
        * dimension expression by the specified value.
        */
        inline UnaryRankWiseDimExpression<IndexType,Rank,Self,MulOp<IndexType>>
        multiply_along_rank(const IndexType& n, std::size_t rank) const {
            return UnaryRankWiseDimExpression<IndexType,Rank,Self,MulOp<IndexType>>(*this, n, rank);
        }
        /**
        * @param n The value to divide by.
        * @param rank The rank which is to be divide by the value.
        * @return An expression representing the division of the specified rank of the
        * dimension expression by the specified value.
        */
        inline UnaryRankWiseDimExpression<IndexType,Rank,Self,DivOp<IndexType>>
        divide_along_rank(const IndexType& n, std::size_t rank) const {
            return UnaryRankWiseDimExpression<IndexType,Rank,Self,DivOp<IndexType>>(*this, n, rank);
        }
        /**
        * @param dims The dimension expression to add.
        * @param rank The rank along which the expressions are to be added.
        * @return An expression representing the addition of the specified expression
        * dimension to the instance on which the method is invoked along the specified rank.
        */
        template<typename OtherDerived>
        inline BinaryRankWiseDimExpression<IndexType,Rank,Self,Other<OtherDerived>,SumOp<IndexType>>
        add_along_rank(const Other<OtherDerived>& dims, std::size_t rank) const {
            return BinaryRankWiseDimExpression<IndexType,Rank,Self,
                    Other<OtherDerived>,SumOp<IndexType>>(*this, dims, rank);
        }
        /**
        * @param dims The dimension expression to subtract.
        * @param rank The rank along which the dimension expression is to be subtracted from
        * the instance on which the method is invoked.
        * @return An expression representing the subtraction of the specified expression
        * dimension from the instance on which the method is invoked along the specified rank.
        */
        template<typename OtherDerived>
        inline BinaryRankWiseDimExpression<IndexType,Rank,Self,Other<OtherDerived>,SubOp<IndexType>>
        subtract_along_rank(const Other<OtherDerived>& dims, std::size_t rank) const {
            return BinaryRankWiseDimExpression<IndexType,Rank,Self,
                    Other<OtherDerived>,SubOp<IndexType>>(*this, dims, rank);
        }
        template<typename OtherDerived>
        inline BinaryDimExpression<IndexType,Rank,Self,Other<OtherDerived>,SumOp<IndexType>>
        operator+(const Other<OtherDerived>& dims) const {
            return BinaryDimExpression<IndexType,Rank,Self,Other<OtherDerived>,SumOp<IndexType>>(*this, dims);
        }
        template<typename OtherDerived>
        inline BinaryDimExpression<IndexType,Rank,Self,Other<OtherDerived>,SubOp<IndexType>>
        operator-(const Other<OtherDerived>& dims) const {
            return BinaryDimExpression<IndexType,Rank,Self,Other<OtherDerived>,SubOp<IndexType>>(*this, dims);
        }
        template<typename OtherDerived>
        inline BinaryDimExpression<IndexType,Rank,Self,Other<OtherDerived>,MulOp<IndexType>>
        operator*(const Other<OtherDerived>& dims) const {
            return BinaryDimExpression<IndexType,Rank,Self,Other<OtherDerived>,MulOp<IndexType>>(*this, dims);
        }
        inline UnaryDimExpression<IndexType,Rank,Self,SumOp<IndexType>>
        operator+(const IndexType& n) const {
            return UnaryDimExpression<IndexType,Rank,Self,SumOp<IndexType>>(*this, n);
        };
        inline UnaryDimExpression<IndexType,Rank,Self,SubOp<IndexType>>
        operator-(const IndexType& n) const {
            return UnaryDimExpression<IndexType,Rank,Self,SubOp<IndexType>>(*this, n);
        }
        inline UnaryDimExpression<IndexType,Rank,Self,DivOp<IndexType>>
        operator/(const IndexType& n) const {
            return UnaryDimExpression<IndexType,Rank,Self,DivOp<IndexType>>(*this, n);
        }
        inline UnaryDimExpression<IndexType,Rank,Self,MulOp<IndexType>>
        operator*(const IndexType& n) const {
            return UnaryDimExpression<IndexType,Rank,Self,MulOp<IndexType>>(*this, n);
        }
        inline friend UnaryDimExpression<IndexType,Rank,Self,MulOp<IndexType>>
        operator*(const IndexType& n, const Self& dims) {
            return UnaryDimExpression<IndexType,Rank,Self,MulOp<IndexType>>(dims, n);
        };
        template<typename OtherDerived, typename OtherIndexType, std::size_t OtherRank>
        inline bool operator==(const DimExpression<OtherDerived,OtherIndexType,OtherRank>& dims) const {
            return false;
        }
        template<typename OtherDerived>
        inline bool operator==(const Other<OtherDerived>& dims) const {
            for (std::size_t i = 0; i < Rank; ++i) {
                if ((*this)(i) != dims(i))
                    return false;
            }
            return true;
        }
        inline bool operator==(const std::array<IndexType,Rank>& dims) const {
            for (std::size_t i = 0; i < Rank; ++i) {
                if ((*this)(i) != dims[i])
                    return false;
            }
            return true;
        }
        inline friend std::ostream& operator<<(std::ostream& os, const Self& dims) {
            return os << dims.to_string();
        }
    };

    /**
    * An alias for a pair of two tensors of the same rank. It represents observation-objective pairs.
    */
    template<typename Scalar, std::size_t Rank, bool Sequential>
    using DataPair = std::pair<Tensor<Scalar,Rank + Sequential + 1>,Tensor<Scalar,Rank + Sequential + 1>>;

    /**
    * A class template for fetching data from memory or disk.
    */
    template<typename Scalar, std::size_t Rank, bool Sequential>
    class DataProvider {
        static_assert(std::is_floating_point<Scalar>::value, "non floating-point scalar type");
        static_assert(Rank > 0 && Rank < 4, "illegal data provider rank");
    protected:
        static constexpr std::size_t DATA_RANK = Rank + Sequential + 1;
        typedef Tensor<Scalar,DATA_RANK> Data;
        typedef Dimensions<std::size_t,Rank> Dims;
    public:
        virtual ~DataProvider() = default;
        /**
        * A simple constant getter method for the dimensions of the observations.
        *
        * @return A constant reference to the dimensions of the observations.
        */
        virtual const Dims& get_obs_dims() const = 0;
        /**
        * A simple constant getter method for the dimensions of the objectives.
        *
        * @return A constant reference to the dimensions of the objectives.
        */
        virtual const Dims& get_obj_dims() const = 0;
        /**
        * A method that returns whether the data provider instance has more data to provide.
        * It should always be called before calling get_data(std::size_t).
        *
        * @return Whether there are more observation-objective pairs to read from the
        * instance.
        */
        virtual bool has_more() = 0;
        /**
        * Reads and returns the specified number of observation-objective pairs. It also
        * offsets the reader by the specified number. If has_more() returns false, the
        * invocation of this method results in the throwing of a std::out_of_range exception.
        *
        * @param batch_size The maximum number of observation-objective pairs to read and
        * return.
        * @return At most batch_size number of observation-objective pairs. If the number
        * of unread pairs is less than batch_size, the number of returned pairs is that
        * of the unread ones.
        */
        virtual DataPair<Scalar,Rank,Sequential> get_data(std::size_t batch_size) = 0;
        /**
        * It resets the reader head to the beginning of the data storage.
        */
        virtual void reset() = 0;
        /**
        * It skips the specified number of data points. If has_more() returns false,
        * the invocation of the method has no effect.
        *
        * @param instances The number of instances to skip.
        */
        virtual void skip(std::size_t instances) = 0;
    };

    template<typename Scalar>
    class Parameters {
    public:
        virtual ~Parameters() = default;
        /**
        * @return A pointer to copy of the instance on which the method is
        * invoked.
        */
        virtual Parameters<Scalar>* clone() const = 0;
        /**
        * Determines whether the parameters are optimizable. Non-optimizable
        * parameters are ignored by optimizers and do not have to maintain
        * gradients or worry about regularization.
        *
        * @return Whether the parameters are learnable via gradient descent
        * or any other optimization method.
        */
        virtual bool are_optimizable() const = 0;
        /**
        * @return The number of rows of the parameter matrix.
        */
        virtual std::size_t get_rows() const = 0;
        /**
        * @return The number of columns of the parameter matrix.
        */
        virtual std::size_t get_cols() const = 0;
        /**
        * It initializes the values of the parameters.
        */
        virtual void init_values() = 0;
        /**
        * It initializes the gradient of the parameters.
        */
        virtual void init_grad() = 0;
        /**
        * @return A constant reference to the values of the parameters.
        */
        virtual const Matrix<Scalar>& get_values() const = 0;
        /**
        * @param values The new values of the parameters. The matrix is expected to havey
        * the dimensions specified by the get_rows() and get_cols() methods.
        */
        virtual void set_values(Matrix<Scalar> values) = 0;
        /**
        * @return A constant reference to the gradient of the parameters.
        */
        virtual const Matrix<Scalar>& get_grad() const = 0;
        /**
        * @param grad The values to add to the current gradient of the parameters. The matrix
        * is expected to have the dimensions specified by the get_rows() and get_cols() methods.
        */
        virtual void accumulate_grad(const Matrix<Scalar>& grad) = 0;
        /**
        * It resets the gradient to all zeroes.
        */
        virtual void reset_grad() = 0;
        /**
        * @return The regularization penalty imposed on the parameters.
        */
        virtual Scalar get_regularization_penalty() const = 0;
        /**
        * It adds the derivative of the regularization function w.r.t.t values of the
        * parameters to the gradient of the parameters.
        */
        virtual void regularize() = 0;
        /**
        * @return Whether the parameters should not be updated.
        */
        virtual bool are_frozen() const = 0;
        /**
        * @param frozen Whether the parameters are to be frozen, i.e. not to be updated.
        */
        virtual void set_frozen(bool frozen) = 0;
        /**
        * It initializes both the values and the gradient of the parameters.
        */
        inline virtual void init() {
            init_values();
            init_grad();
        }
    };



    /**
    * An alias for a shared pointer to a Parameters instance.
    */
    template<typename Scalar>
    using ParamsSharedPtr = std::shared_ptr<Parameters<Scalar>>;

    /**
    * An abstract class template representing layers in a neural network.
    */
    template<typename Scalar, std::size_t Rank>
    class Layer {
        static_assert(std::is_floating_point<Scalar>::value, "non floating-point scalar type");
        static_assert(Rank > 0 && Rank < 4, "illegal rank");
    protected:
        // Rank is increased by one to allow for batch training.
        static constexpr std::size_t DATA_RANK = Rank + 1;
        typedef Tensor<Scalar,DATA_RANK> Data;
        typedef Dimensions<std::size_t,Rank> Dims;
    public:
        virtual ~Layer() = default;
        /**
        * It returns a clone of the layer instance.
        *
        * @return A pointer to a copy of the instance. The instance does not take ownership of
        * the returned pointer (i.e. the caller is responsible for deleting it).
        */
        virtual Layer<Scalar,Rank>* clone() const = 0;
        /**
        * It returns a clone of the layer instance using a reference to the original's parameters.
        * Non-parametric layers do not need to support parameter sharing and thus are just expected
        * to return a normal clone.
        *
        * @return A clone of the original layer instance sharing the same parameters with the
        * original.
        */
        virtual Layer<Scalar,Rank>* clone_with_shared_params() = 0;
        /**
        * It returns a reference to the layer owning the parameters used. If this owner goes out
        * of scope (in case this one is a clone with shared parameters), the behaviour of the clone
        * is undefined.
        *
        * @return A reference to the layer owning the parameters. If this layer is not using
        * shared parameters, it returns a reference to itself.
        */
        virtual const Layer<Scalar,Rank>& get_params_owner() const = 0;
        /**
        * A simple constant getter method for the input dimensionality of the layer.
        *
        * @return A constant reference to the member variable denoting the dimensions of the
        * tensors accepted by the layer as its input (except for the first rank which denotes
        * the variable sample size).
        */
        virtual const Dims& get_input_dims() const = 0;
        /**
        * A simple constant getter method for the output dimensionality of the layer.
        *
        * @return A constant reference to the member variable denoting the dimensions of the
        * tensors output by the layer along all ranks except the first one.
        */
        virtual const Dims& get_output_dims() const = 0;
        /**
        * A constant method that returns whether this layer functions as an input layer. An input
        * layer does not need to propagate the gradients all the way during the backward pass as
        * it is assumed that no other layer needs them derive the gradient on its parameters. It
        * is therefore possible for an input layer to simply return a null tensor as the output of
        * its backward pass.
        *
        * @return Whether this layer is the input layer of the neural network that contains it.
        */
        virtual bool is_input_layer() const = 0;

        /**
        * Sets this instance's input layer status to the given value.
        *
        * @param input_layer Whether this layer is to be an input layer or not.
        */
        virtual void set_input_layer(bool input_layer) = 0;
        /**
        * It empties the layer's caches such as those required for the derivation of the function
        * represented by the layer.
        */
        virtual void empty_cache() = 0;
        /**
        * It returns a vector of constant non-owning pointers to the parameters of the layer.
        *
        * @return A vector of constant pointers to the parameters of the layer.
        */
        virtual std::vector<const Parameters<Scalar>*> get_params() const = 0;
        /**
        * It returns a vector of non-owning pointers to the parameters of the layer.
        *
        * @return A vector of pointers to the parameters of the layer.
        */
        virtual std::vector<Parameters<Scalar>*> get_params() = 0;
        /**
        * It has the function represented by the layer applied to the input tensor.
        *
        * @param in A tensor representing a batch of observations. The observations are of
        * the rank specified by the layer's template parameter and the input tensors rank is
        * one greater.
        * @param training Whether the input is to be processed in training or inference mode.
        * If the forward pass is performed in inference mode, the backward pass is not
        * guaranteed to work.
        * @return The output of the function represented by the layer applied to the input
        * tensor.
        */
        virtual Data pass_forward(Data in, bool training) = 0;
        /**
        * It back-propagates the derivative of the error function w.r.t. the output of the
        * layer updating the gradient of its learnable parameters along the way if there are
        * any.
        *
        * @param out_grad The derivative of the loss function w.r.t. the output of the
        * layer
        * @return The derivative of the loss function w.r.t. the output of the previous layer
        * or a null tensor if the layer is an input layer.
        */
        virtual Data pass_back(Data out_grad) = 0;
        /**
        * It determines whether the layer instance is a clone using the shared parameters of
        * another instance.
        *
        * @return Whether the layer instance is a shared-parameter clone.
        */
        inline bool is_shared_params_clone() const {
            return this != &get_params_owner();
        }
        /**
        * A method that returns whether the layer has parameters.
        *
        * @return Whether the layer uses parameters.
        */
        inline bool is_parametric() const {
            return get_params().size() > 0;
        }
    };
    /**
    * An abstract neural network class template. It allows for inference and training via
    * back-propagation.
    */
    template<typename Scalar, std::size_t Rank, bool Sequential>
    class NeuralNetwork {
        static_assert(std::is_floating_point<Scalar>::value, "non floating-point scalar type");
        static_assert(Rank > 0 && Rank < 4, "illegal neural network rank");
    protected:
        static constexpr std::size_t DATA_RANK = Rank + Sequential + 1;
        typedef Tensor<Scalar,DATA_RANK> Data;
        typedef Dimensions<std::size_t,Rank> Dims;
        static std::string PARAM_SERIAL_PREFIX;
    public:
        virtual ~NeuralNetwork() = default;
        /**
        * A constant method implementing the clone pattern.
        *
        * @return A pointer to a copy of the instance. The instance does not take ownership of
        * the returned pointer (i.e. the caller is responsible for deleting it).
        */
        virtual NeuralNetwork<Scalar,Rank,Sequential>* clone() const = 0;
        /**
        * @return A constant reference to the member variable denoting the dimensions of the
        * tensors accepted by the network as its input along (except for the first rank which
        * denotes the variable sample size and in case of sequential networks the second rank
        * which denotes the variable time steps).
        */
        virtual const Dims& get_input_dims() const = 0;
        /**
        * @return A constant reference to the member variable denoting the dimensions of the
        * tensors output by the network (except for the first rank which denotes the variable
        * sample size and in case of sequential networks the second rank which denotes the
        * variable time steps).
        */
        virtual const Dims& get_output_dims() const = 0;
        /**
        * @return A vector of pointers to constant layers constituting the network. The ownership
        * of the layers remains with the network.
        */
        virtual std::vector<const Layer<Scalar,Rank>*> get_layers() const = 0;
        /**
        * @return A vector of pointers to the layers of the network. The ownership of the
        * layers remains with the network.
        */
        virtual std::vector<Layer<Scalar,Rank>*> get_layers() = 0;
        /**
        * @return Whether the instance is a foremost network. If the instance is not a stand-alone
        * network and it is not the first module of a complex network, it is not a foremost
        * network. Foremost networks do not need to back-propagate the gradients all the way
        * given that no other network is expected to depend on them.
        */
        virtual bool is_foremost() const = 0;
        /**
        * Sets the foremost status of the network.
        *
        * @param foremost Whether the network is to function as a foremost network.
        */
        virtual void set_foremost(bool foremost) = 0;
        /**
        * Empties the caches of every layer of the network.
        */
        virtual void empty_caches() = 0;
        /**
        * It propagates the input tensor through the network and outputs its prediction.
        *
        * @param input The input tensor to propagate through.
        * @param training Whether the input is to be propagated in training mode or not.
        * Propagating the input in training mode may be more time and memory consuming, but
        * is a prerequisite of back-propagation.
        * @return The output tensor of the network in response to the input.
        */
        virtual Data propagate(Data input, bool training) = 0;
        /**
        * It back-propagates the derivative of the loss function w.r.t. the output of the
        * network through its layers updating the gradients on their parameters.
        *
        * @param out_grad The derivative of the loss function w.r.t. the output of the
        * network.
        * @return The derivative of the loss function w.r.t. the input of the network or
        * a null tensor if the network is a foremost network.
        */
        virtual Data backpropagate(Data out_grad) = 0;
        /**
        * @return A vector of pointers to the constant parameters of the network. The pointers
        * are not necessarily unique.
        */
        inline std::vector<const Parameters<Scalar>*> get_all_params() const {
            std::vector<const Parameters<Scalar>*> params_vec;
            for (auto layer_ptr : get_layers()) {
                if (!layer_ptr)
                    continue;
                for (auto params_ptr : layer_ptr->get_params()) {
                    if (params_ptr)
                        params_vec.push_back(params_ptr);
                }
            }
            return params_vec;
        }
        /**
        * @return A vector of pointers to the parameters of the network. The pointers are
        * not necessarily unique.
        */
        inline std::vector<Parameters<Scalar>*> get_all_params() {
            std::vector<Parameters<Scalar>*> params_vec;
            for (auto layer_ptr : get_layers()) {
                if (!layer_ptr)
                    continue;
                for (auto params_ptr : layer_ptr->get_params()) {
                    if (params_ptr)
                        params_vec.push_back(params_ptr);
                }
            }
            return params_vec;
        }
        /**
        * @return A vector of pointers to the constant, unique parameters of the network.
        */
        inline std::vector<const Parameters<Scalar>*> get_all_unique_params() const {
            std::vector<const Parameters<Scalar>*> params_vec;
            std::set<const Parameters<Scalar>*> params_set;
            for (auto params_ptr : get_all_params()) {
                if (params_set.find(params_ptr) == params_set.end()) {
                    params_set.insert(params_ptr);
                    params_vec.push_back(params_ptr);
                }
            }
            return params_vec;
        }
        /**
        * @return A vector of pointers to the unique parameters of the network.
        */
        inline std::vector<Parameters<Scalar>*> get_all_unique_params() {
            std::vector<Parameters<Scalar>*> params_vec;
            std::set<Parameters<Scalar>*> params_set;
            for (auto params_ptr : get_all_params()) {
                if (params_set.find(params_ptr) == params_set.end()) {
                    params_set.insert(params_ptr);
                    params_vec.push_back(params_ptr);
                }
            }
            return params_vec;
        }
        /**
        * Sets all parameters of the network to the specified frozens state.
        *
        * @param frozen Whether the parameters of the network should be frozen i.e. temporarily
        * not optimizable.
        */
        inline virtual void set_frozen(bool frozen) {
            for (auto params_ptr : get_all_unique_params())
                params_ptr->set_frozen(frozen);
        }
        /**
        * Initializes all parameters of the network.
        */
        inline virtual void init() {
            for (auto params_ptr : get_all_unique_params())
                params_ptr->init();
        }
        /**
        * It propagates the input through the neural network and outputs its prediction
        * according to its current parameters.
        *
        * @param input The input to be mapped.
        * @return The inference/prediction of the neural network.
        */
        inline virtual Data infer(Data input) {
            return propagate(std::move(input), false);
        }
        /**
        * It serializes the values of the unique parameters of the network into files in a
        * specified folder.
        *
        * @param dir_path The path to the directory.
        * @param binary Whether the parameters are to be serialized into a binary format.
        * @param file_name_prefix A prefix to the names of the serialized parameter files.
        */
        inline void save_all_unique_params_values(const std::string& dir_path, bool binary = true,
                const std::string& file_name_prefix = PARAM_SERIAL_PREFIX) const {
            std::vector<const Parameters<Scalar>*> params_vec = get_all_unique_params();
            for (std::size_t i = 0; i < params_vec.size(); ++i) {
                const Matrix<Scalar>& values = params_vec[i]->get_values();
                std::string file_path = dir_path + "/" + file_name_prefix + std::to_string(i) + ".prms";
                if (binary)
                    serialize_binary<Scalar>(values, file_path);
                else
                    serialize<Scalar>(values, file_path);
            }
        }
        /**
        * It sets the values of the unique parameters of the network from files containing
        * serialized parameter values.
        *
        * @param dir_path The path to the directory containing the parameter files.
        * @param binary Whether the parameter files binary.
        * @param file_name_prefix The prefix of the names of the parameter files.
        */
        inline void load_all_unique_params_values(const std::string& dir_path, bool binary = true,
                const std::string& file_name_prefix = PARAM_SERIAL_PREFIX) {
            std::vector<Parameters<Scalar>*> params_vec = get_all_unique_params();
            for (std::size_t i = 0; i < params_vec.size(); ++i) {
                std::string file_path = dir_path + "/" + file_name_prefix + std::to_string(i) + ".prms";
                params_vec[i]->set_values(binary ? deserialize_binary<Scalar>(file_path) :
                        deserialize<Scalar>(file_path));
            }
        }
    };

    template<typename Scalar, std::size_t Rank, bool Sequential>
    std::string NeuralNetwork<Scalar,Rank,Sequential>::PARAM_SERIAL_PREFIX =
            "c-attl3_neural_net_params_";


    /**
    * An abstract class template for loss functions. Implementations of this class should be stateless.
    */
    template<typename Scalar, std::size_t Rank, bool Sequential>
    class Loss {
        static_assert(std::is_floating_point<Scalar>::value, "non floating-point scalar type");
        static_assert(Rank > 0 && Rank < 4, "illegal loss rank");
    protected:
        static constexpr std::size_t DATA_RANK = Rank + Sequential + 1;
        typedef Tensor<Scalar,DATA_RANK> Data;
    public:
        virtual ~Loss() = default;
        /**
        * It calculates the error on each sample given the output and the objective tensors.
        *
        * @param out The output tensor.
        * @param obj The objective tensor (of the same dimensionality as the output).
        * @return A column vector containing the loss for each sample.
        */
        virtual ColVector<Scalar> function(Data out, Data obj) const = 0;
        /**
        * It calculates the derivative of the loss function w.r.t. the output.
        *
        * @param out The output tensor.
        * @param obj The objective tensor (of the same dimensionality as the output).
        * @return The derivative of the loss function w.r.t. the output.
        */
        virtual Data d_function(Data out, Data obj) const = 0;
    };


    /**
    * A utility class template containing static methods and variables to help with
    * numerical issues.
    */
    template<typename Scalar>
    class NumericUtils {
        static_assert(std::is_floating_point<Scalar>::value, "non floating-point scalar type");
    private:
    public:
        static constexpr Scalar MIN = std::numeric_limits<Scalar>::lowest();
        static constexpr Scalar MAX = std::numeric_limits<Scalar>::max();
        static constexpr Scalar EPSILON1 = std::numeric_limits<Scalar>::epsilon();
        static constexpr Scalar EPSILON2 = 1e-5;
        NumericUtils() = delete;
        /**
        * Returns whether the two numerals are close enough to be considered equal.
        *
        * @param n1 The first numeral.
        * @param n2 The second numeral.
        * @param abs_epsilon The maximum absolute difference that would still allow
        * them to be considered equal.
        * @param rel_epsilon The maximum relative (to the greater numeral of the two)
        * difference that would still allow them to be considered equal.
        * @return Whether the two numerals can be considered equal.
        */
        inline static bool almost_equal(Scalar n1, Scalar n2, Scalar abs_epsilon = EPSILON1,
                Scalar rel_epsilon = EPSILON1) {
            Scalar diff = std::abs(n1 - n2);
            if (diff <= abs_epsilon)
                return true;
            Scalar max = std::max(std::abs(n1), std::abs(n2));
            return diff <= max * rel_epsilon;
        }
        /**
        * Returns whether a numeral is greater than another one by a margin great enough for
        * them not to be considered almost equal.
        *
        * @param n1 The first numeral.
        * @param n2 The second numeral.
        * @param abs_epsilon The maximum absolute difference that would still allow
        * them to be considered equal.
        * @param rel_epsilon The maximum relative (to the greater numeral of the two)
        * difference that would still allow them to be considered equal.
        * @return Whether the the first numeral is sufficiently greater than the second.
        */
        inline static bool decidedly_greater(Scalar n1, Scalar n2, Scalar abs_epsilon = EPSILON1,
                Scalar rel_epsilon = EPSILON1) {
            return n1 > n2 && !almost_equal(n1, n2, abs_epsilon, rel_epsilon);
        }
        /**
        * Returns whether a numeral is not lesser than another one by a margin great enough for
        * them not to be considered almost equal.
        *
        * @param n1 The first numeral.
        * @param n2 The second numeral.
        * @param abs_epsilon The maximum absolute difference that would still allow
        * them to be considered equal.
        * @param rel_epsilon The maximum relative (to the greater numeral of the two)
        * difference that would still allow them to be considered equal.
        * @return Whether the the first numeral is sufficiently great enough not to be considered
        * decidedly smaller than the second.
        */
        inline static bool greater_or_almost_equal(Scalar n1, Scalar n2, Scalar abs_epsilon = EPSILON1,
                Scalar rel_epsilon = EPSILON1) {
            return n1 > n2 || almost_equal(n1, n2, abs_epsilon, rel_epsilon);
        }
        /**
        * Returns whether a numeral is lesser than another one by a margin great enough for
        * them not to be considered almost equal.
        *
        * @param n1 The first numeral.
        * @param n2 The second numeral.
        * @param abs_epsilon The maximum absolute difference that would still allow
        * them to be considered equal.
        * @param rel_epsilon The maximum relative (to the greater numeral of the two)
        * difference that would still allow them to be considered equal.
        * @return Whether the the first numeral is sufficiently lesser than the second.
        */
        inline static bool decidedly_lesser(Scalar n1, Scalar n2, Scalar abs_epsilon = EPSILON1,
                Scalar rel_epsilon = EPSILON1) {
            return n1 < n2 && !almost_equal(n1, n2, abs_epsilon, rel_epsilon);
        }
        /**
        * Returns whether a numeral is not greater than another one by a margin great enough for
        * them not to be considered almost equal.
        *
        * @param n1 The first numeral.
        * @param n2 The second numeral.
        * @param abs_epsilon The maximum absolute difference that would still allow
        * them to be considered equal.
        * @param rel_epsilon The maximum relative (to the greater numeral of the two)
        * difference that would still allow them to be considered equal.
        * @return Whether the the first numeral is sufficiently small enough not to be considered
        * decidedly greater than the second.
        */
        inline static bool lesser_or_almost_equal(Scalar n1, Scalar n2, Scalar abs_epsilon = EPSILON1,
                Scalar rel_epsilon = EPSILON1) {
            return n1 < n2 || almost_equal(n1, n2, abs_epsilon, rel_epsilon);
        }
    };

    /**
    * An alias for a unique pointer to a loss function of arbitrary rank, scalar type and
    * sequentiality.
    */
    template<typename Scalar, std::size_t Rank, bool Sequential>
    using LossSharedPtr = std::shared_ptr<Loss<Scalar,Rank,Sequential>>;

    /**
    * An abstract class template for neural network optimizer algorithm implementations.
    */
    template<typename Scalar, std::size_t Rank, bool Sequential>
    class Optimizer {
        static_assert(std::is_floating_point<Scalar>::value, "non floating-point scalar type");
        static_assert(Rank > 0 && Rank < 4, "illegal optimizer rank");
    protected:
        typedef NeuralNetwork<Scalar,Rank,Sequential> Net;
        typedef DataProvider<Scalar,Rank,Sequential> Provider;
        typedef Tensor<Scalar,Rank + Sequential + 1> Data;
    public:
        inline Optimizer(LossSharedPtr<Scalar,Rank,Sequential> loss) :
                    loss(loss) {
            assert(loss != nullptr);
        }
        virtual ~Optimizer() = default;
        /**
        * It optimizes the specified neural network using the given data providers according to the
        * optimizers loss function. It also fits the optimizer to the network before the otpimization
        * begins.
        *
        * @param net A reference to the network whose parameters are to be optimized.
        * @param training_prov A reference to the provider of the training data.
        * @param test_prov A reference to the provider of the test data.
        * @param epochs The number of epochs for which the optimization should proceed.
        * @param early_stop An std::size_t integer denoting the number of consecutive loss increases
        * after which the optimization process is to be terminated prematurely. If it is 0, the
        * process is never terminated prematurely.
        * @param target_loss The target test loss value. If the test loss reaches this value or
        * drops below it, the optimization process is terminated.
        * @param verbose Whether the training, test, and regularization losses for each epoch should
        * be printed to the standard out stream.
        * @return The test loss of the last epoch.
        */
        inline Scalar optimize(Net& net, Provider& training_prov, Provider& test_prov, std::size_t epochs,
                std::size_t early_stop = 0, Scalar target_loss = NumericUtils<Scalar>::MIN, bool verbose = true) {
            assert(net.get_input_dims() == training_prov.get_obs_dims());
            assert(net.get_output_dims() == training_prov.get_obj_dims());
            assert(training_prov.get_obs_dims() == test_prov.get_obs_dims());
            assert(training_prov.get_obj_dims() == test_prov.get_obj_dims());
            // Fit the optimizer parameters to the network.
            fit(net);
            Scalar prev_test_loss = NumericUtils<Scalar>::MAX;
            std::size_t cons_loss_inc = 0;
            if (verbose)
                std::cout << "<Optimization>" << std::endl;
            // Start the optimization iterations.
            for (std::size_t i = 0; i <= epochs; ++i) {
                if (verbose)
                    std::cout << "Epoch " << std::setw(3) << i << std::string(28, '-') << std::endl;
                // Train.
                if (i != 0) {
                    training_prov.reset();
                    Scalar train_loss = _train(net, training_prov, i, verbose);
                    if (verbose) {
                        std::cout << std::left << std::setw(20) << "\ttraining loss: " << std::right <<
                                std::to_string(train_loss) << std::endl;
                    }
                }
                // Validate.
                test_prov.reset();
                Scalar test_loss = _test(net, test_prov, i, verbose);
                if (verbose) {
                    std::cout << std::left << std::setw(20) << "\ttest loss: " << std::right <<
                            std::to_string(test_loss);
                }
                if (test_loss >= prev_test_loss) {
                    cons_loss_inc++;
                    if (verbose)
                        std::cout << " *****INCREASED LOSS*****";
                    if (early_stop > 0 && cons_loss_inc >= early_stop)
                        break;
                } else
                    cons_loss_inc = 0;
                if (verbose)
                    std::cout << std::endl << std::endl;
                prev_test_loss = test_loss;
                if (prev_test_loss <= target_loss)
                    break;
            }
            // Reset the providers.
            training_prov.reset();
            test_prov.reset();
            // Empty the network caches.
            net.empty_caches();
            return prev_test_loss;
        }
        /**
        * It trains the specified neural network using the given training data provider according to
        * the optimizers loss function for the specified number of epochs. It does not fit the
        * optimizer to the network, thus the #fit(Net&) method might need to be invoked beforehand.
        *
        * @param net A reference to the network whose parameters are to be optimized.
        * @param prov A reference to the provider of the training data.
        * @param epochs The number of epochs for which the training should proceed.
        * @param early_stop An std::size_t integer denoting the number of consecutive loss increases
        * after which the optimization process is to be terminated prematurely. If it is 0, the
        * process is never terminated prematurely.
        * @param target_loss The target test loss value. If the test loss reaches this value or
        * drops below it, the optimization process is terminated.
        * @param verbose Whether the training losses for of the epochs should be printed to the
        * standard out stream.
        * @return The training loss of the last epoch.
        */
        inline Scalar train(Net& net, Provider& prov, std::size_t epochs, std::size_t early_stop = 0,
                Scalar target_loss = NumericUtils<Scalar>::MIN, bool verbose = false) {
            assert(net.get_input_dims() == prov.get_obs_dims());
            assert(net.get_output_dims() == prov.get_obj_dims());
            Scalar train_loss;
            Scalar prev_train_loss = NumericUtils<Scalar>::MAX;
            std::size_t cons_loss_inc = 0;
            if (verbose)
                std::cout << "<Training>" << std::endl;
            for (std::size_t i = 1; i <= epochs; ++i) {
                if (verbose)
                    std::cout << "Epoch " << std::setw(3) << i << std::string(28, '-') << std::endl;
                prov.reset();
                train_loss = _train(net, prov, i, verbose);
                if (verbose) {
                    std::cout << std::left << std::setw(20) << "\ttraining loss: " << std::right <<
                            std::to_string(train_loss);
                }
                if (train_loss >= prev_train_loss) {
                    cons_loss_inc++;
                    if (verbose)
                        std::cout << " *****INCREASED LOSS*****";
                    if (early_stop > 0 && cons_loss_inc >= early_stop)
                        break;
                } else
                    cons_loss_inc = 0;
                if (verbose)
                    std::cout << std::endl << std::endl;
                prev_train_loss = train_loss;
                if (prev_train_loss <= target_loss)
                    break;
            }
            prov.reset();
            net.empty_caches();
            return prev_train_loss;
        }
        /**
        * It tests the specified neural network using the given test data provider according to the
        * optimizers loss function. It does not fit the optimizer to the network, thus the
        * #fit(Net&) method might need to be invoked beforehand.
        *
        * @param net A reference to the network whose parameters are to be optimized.
        * @param prov A reference to the provider of the training data.
        * @param verbose Whether the training, test, and regularization losses for each epoch should
        * be printed to the standard out stream.
        * @return The test loss of the last epoch.
        */
        inline Scalar test(Net& net, Provider& prov, bool verbose = false) {
            assert(net.get_input_dims() == prov.get_obs_dims());
            assert(net.get_output_dims() == prov.get_obj_dims());
            if (verbose)
                std::cout << "<Testing>" << std::endl;
            prov.reset();
            Scalar test_loss = _test(net, prov, 0, verbose);
            if (verbose) {
                std::cout << std::left << std::setw(20) << "\ttest loss: " << std::right <<
                        std::to_string(test_loss) << std::endl;
            }
            prov.reset();
            net.empty_caches();
            return test_loss;
        }
        /**
        * It fits the optimizer to the neural network. It allows optimizers with individual
        * learning rates for each parameter to set up their necessary internal data structures.
        *
        * @param net A reference to the neural network that is to be optimized.
        */
        virtual void fit(Net& net) = 0;
    protected:
        /**
        * It trains the specified neural network for a single epoch on data provided by the
        * specified data provider.
        *
        * @param net A reference to the neural network to optimize.
        * @param training_prov A reference to the training data provider.
        * @param epoch The index of the current epoch. starting from 1.
        * @param verbose Whether the optimization is performed in verbose mode; i.e. whether
        * information should be printed to the standard out stream.
        * @return The training loss of the epoch.
        */
        virtual Scalar _train(Net& net, Provider& training_prov, std::size_t epoch, bool verbose) = 0;
        /**
        * It tests the specified neural network for a single epoch on the test data
        * provided by the specified data provider.
        *
        * @param net A reference to the neural network to test.
        * @param test_prov A reference to the test data provider.
        * @param epoch The index of the epoch starting from 0.
        * @param verbose Whether the optimization is performed in verbose mode; i.e. whether
        * information should be printed to the standard out stream.
        * @return The test loss of the epoch.
        */
        virtual Scalar _test(Net& net, Provider& test_prov, std::size_t epoch, bool verbose) = 0;
        const LossSharedPtr<Scalar,Rank,Sequential> loss;
    };

    /**
    * An abstract class template for different weight initialization methods for kernel layers.
    */
    template<typename Scalar>
    class ParameterInitialization {
        static_assert(std::is_floating_point<Scalar>::value, "non floating-point scalar type");
    public:
        virtual ~ParameterInitialization() = default;
        /**
        * It initializes the values of the parameters.
        *
        * @param params A reference to the parameter matrix.
        */
        virtual void apply(Matrix<Scalar>& params) const = 0;
    };

    /**
    * An abstract template class for different regularization penalties for neural network
    * layer parameters. Implementations of this class should be stateless.
    */
    template<typename Scalar>
    class ParameterRegularization {
        static_assert(std::is_floating_point<Scalar>::value, "non floating-point scalar type");
    public:
        virtual ~ParameterRegularization() = default;
        /**
        * It computes the regularization penalty for the given parameter values.
        *
        * @param params A constant reference to the parameter matrix.
        * @return The regularization penalty as a single scalar.
        */
        virtual Scalar function(const Matrix<Scalar>& params) const;
        /**
        * It differentiates the regularization function and returns its derivative
        * w.r.t. the parameters.
        *
        * @param params A constant reference to the parameter matrix.
        * @return The gradient matrix.
        */
        virtual Matrix<Scalar> d_function(const Matrix<Scalar>& params) const;
    };


    /**
    * An abstract class template for data preprocessors.
    */
    template<typename Scalar, std::size_t Rank, bool Sequential>
    class Preprocessor {
        static_assert(std::is_floating_point<Scalar>::value, "non floating-point scalar type");
        static_assert(Rank > 0 && Rank < 4, "illegal pre-processor rank");
    public:
        virtual ~Preprocessor() = default;
        /**
        * It fits the preprocessor to the specified data.
        *
        * @param data A constant reference to a data tensor.
        */
        virtual void fit(const Tensor<Scalar,Rank + Sequential + 1>& data) = 0;
        /**
        * It transforms the specified tensor according to the preprocessors current state
        * created by #fit(const Tensor<Scalar,Rank + Sequential + 1>&).
        *
        * @param data A non-constant reference to a data tensor.
        */
        virtual void transform(Tensor<Scalar,Rank + Sequential + 1>& data) const = 0;
    };

    /**
    * An abstract class template for a data provider backed by data on disk in the form of an arbitrary
    * number of files containing both the observations and the objectives. Implementations are responsible
    * for specifying the dimensions of both the observations and the objectives, for reading batches of
    * observation-objective pairs from the file, and for skipping arbitrary number of data instances.
    */
    template<typename Scalar, std::size_t Rank, bool Sequential, bool Binary = false>
    class JointFileDataProvider : public DataProvider<Scalar,Rank,Sequential> {
        typedef DataProvider<Scalar,Rank,Sequential> Base;
    public:
        virtual ~JointFileDataProvider() = default;
        inline const typename Base::Dims& get_obs_dims() const {
            return obs_dims;
        }
        inline const typename Base::Dims& get_obj_dims() const {
            return obj_dims;
        }
        inline bool has_more() {
            if (current_file_stream_has_more())
                return true;
            ++current_file_ind;
            for (; current_file_ind < files.size(); ++current_file_ind) {
                init_current_file_stream();
                if (current_file_stream_has_more())
                    return true;
            }
            return false;
        }
        inline DataPair<Scalar,Rank,Sequential> get_data(std::size_t batch_size) {
            if (!has_more())
                throw std::out_of_range("no more data left to fetch");
            DataPair<Scalar,Rank,Sequential> data_pair = _get_data(files[current_file_ind], current_file_stream,
                    batch_size);
            assert(data_pair.first.dimension(0) == data_pair.second.dimension(0));
            /* If the data contains fewer batches than expected, the end of the file has been reached and the
            * rest of the data should be read from the next file. */
            while (data_pair.first.dimension(0) < batch_size && has_more()) {
                DataPair<Scalar,Rank,Sequential> add_data_pair = _get_data(files[current_file_ind],
                        current_file_stream, batch_size - data_pair.first.dimension(0));
                assert(add_data_pair.first.dimension(0) == add_data_pair.second.dimension(0));
                // It has to be evaluated into a temporary due to the dimension incompatibility.
                typename Base::Data obs_concat = data_pair.first.concatenate(std::move(add_data_pair.first), 0);
                data_pair.first = std::move(obs_concat);
                typename Base::Data obj_concat = data_pair.second.concatenate(std::move(add_data_pair.second), 0);
                data_pair.second = std::move(obj_concat);
            }
            return data_pair;
        }
        inline void reset() {
            current_file_ind = 0;
            init_current_file_stream();
            _set_to_beg(current_file_stream);
        }
        inline void skip(std::size_t instances) {
            if (!has_more())
                return;
            std::size_t skipped = _skip(current_file_stream, instances);
            while (skipped < instances && has_more())
                skipped += _skip(current_file_stream, instances - skipped);
        }
    protected:
        inline JointFileDataProvider(const typename Base::Dims& obs_dims, const typename Base::Dims& obj_dims,
                const std::vector<std::string>& dataset_paths) :
                    obs_dims(obs_dims),
                    obj_dims(obj_dims),
                    files(dataset_paths),
                    current_file_ind(0) {
            assert(!files.empty());
            init_current_file_stream();
        }
        inline JointFileDataProvider(const typename Base::Dims& obs_dims, const typename Base::Dims& obj_dims,
                std::string dataset_path) :
                    JointFileDataProvider(obs_dims, obj_dims, { dataset_path }) { }
        /**
        * It sets the position of the file stream to the beginning of the data set.
        *
        * @param file_stream A reference to the file stream of the data set.
        */
        virtual inline void _set_to_beg(std::ifstream& file_stream) {
            file_stream.seekg(0, std::ios::beg);
        }
        /**
        * It reads at most the specified number of observation-objective pairs from the provided
        * file stream. The file stream can be expected not to have any of its fail flags set
        * initially and to have at least 1 more character left to read.
        *
        * @param file_name The name of the data source file.
        * @param file_stream The input stream of the file.
        * @param batch_size The number of data points to return.
        * @return A pair of tensors containing the data batch.
        */
        virtual DataPair<Scalar,Rank,Sequential> _get_data(const std::string& file_name, std::ifstream& file_stream,
                std::size_t batch_size) = 0;
        /**
        * Skips at most the specified number of instances in the data stream. The file stream can
        * be expected not to have any of its fail flags set initially.
        *
        * @param file_stream A reference to the file stream of the data set.
        * @param instances The number of instances to skip.
        * @return The number of instances actually skipped. It may be less than the specified
        * amount if there are fewer remaining instances in the data stream.
        */
        virtual std::size_t _skip(std::ifstream& file_stream, std::size_t instances) = 0;
        const typename Base::Dims obj_dims, obs_dims;
    private:
        bool current_file_stream_has_more() {
            return current_file_stream && current_file_stream.peek() != EOF;
        }
        void init_current_file_stream() {
            current_file_stream = std::ifstream(files[current_file_ind], Binary ? std::ios::binary : std::ios::in);
            assert(current_file_stream.is_open());
        }
        std::vector<std::string> files;
        std::size_t current_file_ind;
        std::ifstream current_file_stream;
    };

    /**
    * An enum denoting different CIFAR data set types.
    */
    enum CIFARType { CIFAR_10, CIFAR_100 };

    /**
    * A data provider template for the CIFAR-10 and CIFAR-100 data sets.
    *
    * \see https://www.cs.toronto.edu/~kriz/cifar.html
    */
    template<typename Scalar, CIFARType CIFARType = CIFAR_10>
    class CIFARDataProvider : public JointFileDataProvider<Scalar,3,false,true> {
        static_assert(CIFARType == CIFAR_10 || CIFARType == CIFAR_100, "invalid CIFAR type");
        typedef DataProvider<Scalar,3,false> Root;
        typedef JointFileDataProvider<Scalar,3,false,true> Base;
        typedef std::array<std::size_t,4> RankwiseArray;
        static constexpr std::size_t INSTANCE_LENGTH = CIFARType == CIFAR_10 ? 3073 : 3074;
        static constexpr std::size_t NUM_LABELS = CIFARType == CIFAR_10 ? 10 : 100;
    public:
        /**
        * @param file_paths The paths to the data set files.
        */
        inline CIFARDataProvider(std::vector<std::string> file_paths) :
                Base::JointFileDataProvider({ 32u, 32u, 3u }, { NUM_LABELS, 1u, 1u }, file_paths),
                offsets({ 0u, 0u, 0u, 0u }),
                obs_extents(Base::obs_dims.template promote<>()),
                obj_extents(Base::obj_dims.template promote<>()) {
            Base::reset();
        }
        /**
        * @param file_path The path to the data set file.
        */
        inline CIFARDataProvider(std::string file_path) :
                CIFARDataProvider(std::vector<std::string>({ file_path })) { }
    protected:
        inline DataPair<Scalar,3,false> _get_data(const std::string& file_name, std::ifstream& file_stream,
                std::size_t batch_size) {
            typename Root::Data obs(batch_size, Base::obs_dims(0), Base::obs_dims(1), Base::obs_dims(2));
            typename Root::Data obj(batch_size, Base::obj_dims(0), Base::obj_dims(1), Base::obj_dims(2));
            obj.setZero();
            std::size_t i;
            for (i = 0; i < batch_size && file_stream.read(buffer, INSTANCE_LENGTH); ++i) {
                unsigned char* u_buffer = reinterpret_cast<unsigned char*>(buffer);
                std::size_t buffer_ind = 0;
                // Set the label.
                if (CIFARType == CIFAR_100)
                    buffer_ind++;
                obj(i,u_buffer[buffer_ind++],0u,0u) = (Scalar) 1;
                // Set the image.
                for (std::size_t channel = 0; channel < Base::obs_dims(2); ++channel) {
                    for (std::size_t height = 0; height < Base::obs_dims(1); ++height) {
                        for (std::size_t width = 0; width < Base::obj_dims(0); ++width)
                            obs(i,height,width,channel) = (Scalar) u_buffer[buffer_ind++];
                    }
                }
                assert(buffer_ind == INSTANCE_LENGTH);
            }
            if (i == batch_size)
                return std::make_pair(obs, obj);
            obs_extents[0] = i;
            obj_extents[0] = i;
            return std::make_pair(obs.slice(offsets, obs_extents), obj.slice(offsets, obj_extents));
        }
        inline std::size_t _skip(std::ifstream& file_stream, std::size_t instances) {
            std::streampos curr_pos = file_stream.tellg();
            file_stream.seekg(0, std::ios::end);
            std::size_t skip_extent = file_stream.tellg() - curr_pos;
            file_stream.seekg(curr_pos);
            file_stream.ignore(instances * INSTANCE_LENGTH);
            return std::min(instances, skip_extent / INSTANCE_LENGTH);
        }
    private:
        char buffer[INSTANCE_LENGTH];
        RankwiseArray offsets, obs_extents, obj_extents;
    };


    /**
    * An alias for a read-only dictionary mapping words to indices.
    */
    typedef const std::map<std::string,std::size_t> Vocab;

    /**
    * An alias for a shared pointer to a vocabulary.
    */
    typedef std::shared_ptr<Vocab> VocabSharedPtr;

    /**
    * An enumeration for the different objective types to use for the IMDB data set.
    */
    enum IMDBObjType { BINARY, SMOOTH, CATEGORICAL };

    /**
    * A data provider template for the IMDB Large Movie Review Dataset.
    *
    * \see http://ai.stanford.edu/~amaas/data/sentiment/
    */
    template<typename Scalar, IMDBObjType ObjType = BINARY>
    class IMDBDataProvider : public JointFileDataProvider<Scalar,1,true> {
        static_assert(ObjType >= BINARY && ObjType <= CATEGORICAL, "invalid IMDB objective type");
        typedef DataProvider<Scalar,1,true> Root;
        typedef JointFileDataProvider<Scalar,1,true> Base;
        typedef std::array<std::size_t,Root::DATA_RANK> RankwiseArray;
        static constexpr std::size_t PREALLOC_SEQ_LENGTH = 100;
    public:
        static constexpr std::size_t PAD_IND = 0;
        static constexpr std::size_t UNK_IND = 1;
        /**
        * @param pos_reviews_folder_path The path to the folder containing the positive reviews
        * (without a trailing path separator).
        * @param neg_reviews_folder_path The pathe to the folder containing the negative reveies
        * (without a trailing path separator).
        * @param vocab A shared pointer to the vocabulary to use.
        * @param seq_length The sequence length to trim or pad the data to so as to enable batch
        * training. If it is set to 0, no sequence trimming or padding is to be performed (which
        * is likely to make batch training impossible).
        */
        inline IMDBDataProvider(std::string pos_reviews_folder_path, std::string neg_reviews_folder_path,
                VocabSharedPtr vocab, std::size_t seq_length = 0) :
                    Base::JointFileDataProvider({ vocab->size() }, { ObjType == CATEGORICAL ? 10 : 1 },
                            resolve_review_files(pos_reviews_folder_path, neg_reviews_folder_path)),
                    vocab(vocab),
                    seq_length(seq_length) {
            Base::reset();
        }
        /**
        * It populates a dictionary mapping words to indices given the path to a
        * vocabulary file.
        *
        * @param vocab_path The path to the file listing all words appearing in the
        * corpus.
        * @return A map data structure representing the corpus' vocabulary.
        */
        inline static VocabSharedPtr build_vocab(std::string vocab_path) {
            std::ifstream vocab_stream(vocab_path);
            assert(vocab_stream.is_open());
            std::map<std::string,std::size_t> vocab;
            // Reserve the first two indices for padding and unknown words.
            vocab.emplace(std::make_pair("<PAD>", +PAD_IND));
            vocab.emplace(std::make_pair("<UNK>", +UNK_IND));
            std::size_t index = 2;
            std::string word;
            while (vocab_stream >> word)
                vocab.emplace(std::make_pair(word, index++));
            return std::make_shared<const std::map<std::string,std::size_t>>(std::move(vocab));
        }
        /**
        * Converts the embedded text into a string.
        *
        * @param obs A tensor representing the embedded text.
        * @param vocab A shared pointer to the vocabulary.
        * @return The text in the form of a string reconstructed from the provided
        * tensor using the specified vocabulary.
        */
        inline static std::string convert_to_text(const Tensor<Scalar,3>& obs, VocabSharedPtr vocab) {
            std::stringstream text_stream;
            std::string separator("");
            for (std::size_t i = 0; i < obs.dimension(0); ++i) {
                for (std::size_t j = 0; j < obs.dimension(1); ++j) {
                    for (std::size_t k = 0; k < obs.dimension(2); ++k) {
                        if (obs(i,j,k) == 1) {
                            for (auto& entry : *vocab) {
                                if (entry.second == k) {
                                    text_stream << separator << entry.first;
                                    separator = " ";
                                }
                            }
                        }
                    }
                }
            }
            return text_stream.str();
        }
    protected:
        /**
        * @param dir_path The path to the directory.
        * @param file_names A vector to be populated by the paths to all the files
        * contained in the directory.
        */
        inline static void read_files_in_dir(std::string dir_path, std::vector<std::string>& file_names) {
            auto dir_ptr = opendir(dir_path.c_str());
            struct dirent* dir_ent_ptr;
            while ((dir_ent_ptr = readdir(dir_ptr)))
                file_names.push_back(dir_path + "/" + std::string(dir_ent_ptr->d_name));
            closedir(dir_ptr);
        }
        /**
        * @param pos_reviews_folder_path The path to the directory containing
        * the positive movie reviews.
        * @param neg_reviews_folder_path The path to the directory containing
        * the negative movie reviews.
        * @return A randomly shuffled vector of the paths to all files contained
        * in the directory.
        */
        inline static std::vector<std::string> resolve_review_files(std::string pos_reviews_folder_path,
                std::string neg_reviews_folder_path) {
            std::vector<std::string> file_names;
            read_files_in_dir(pos_reviews_folder_path, file_names);
            read_files_in_dir(neg_reviews_folder_path, file_names);
            std::random_shuffle(file_names.begin(), file_names.end());
            return file_names;
        }
        /**
        * It cleans the document by replacing all unsupported characters and words.
        *
        * @param document A string stream to the document to clean
        */
        inline static void clean_document(std::stringstream& document) {
            std::string doc_string = document.str();
            std::transform(doc_string.begin(), doc_string.end(), doc_string.begin(),
                    static_cast<int (*)(int)>(std::tolower));
            // Replace illegal character sequences by white spaces.
            static std::regex illegal_regex("(<br />)+|([^a-zA-Z-'!\?]+)");
            doc_string = std::regex_replace(doc_string, illegal_regex, " ");
            // Add a white space before the supported punctuation marks.
            static std::regex punct_regex("([!\?]{1})");
            doc_string = std::regex_replace(doc_string, punct_regex, " $1");
            document.str(doc_string);
        }
        inline DataPair<Scalar,1,true> _get_data(const std::string& file_name, std::ifstream& file_stream,
                std::size_t batch_size) {
            assert(batch_size > 0);
            const bool fixed_seq_length = seq_length != 0;
            typename Root::Data obs(1, (fixed_seq_length ? seq_length : +PREALLOC_SEQ_LENGTH), Base::obs_dims(0u));
            typename Root::Data obj(1, 1, Base::obj_dims(0u));
            obs.setZero();
            // Parse the rating from the name of the file.
            std::size_t last_under_score = file_name.find_last_of('_');
            std::size_t last_period = file_name.find_last_of('.');
            std::string rating_string = file_name.substr(last_under_score + 1,
                    last_period - last_under_score - 1);
            unsigned rating = (unsigned) std::stoi(rating_string);
            switch (ObjType) {
                case BINARY:
                    obj(0u,0u,0u) = (Scalar) (rating > 5);
                    break;
                case SMOOTH:
                    obj(0u,0u,0u) = ((Scalar) (rating - 1)) / 9;
                    break;
                case CATEGORICAL:
                    obj.setZero();
                    obj(0u,0u,rating) = (Scalar) 1;
                    break;
                default:
                    assert(false);
            }
            // Read the document into a string so it can be pre-processed.
            std::stringstream doc_stream;
            doc_stream << file_stream.rdbuf();
            clean_document(doc_stream);
            // Tokenize the document.
            std::size_t time_step = 0;
            std::string word;
            while (doc_stream >> word && (time_step < seq_length || !fixed_seq_length)) {
                std::size_t ind;
                Vocab::const_iterator val = vocab->find(word);
                ind = (val != vocab->end()) ? val->second : +UNK_IND;
                if (!fixed_seq_length && time_step >= obs.dimension(1)) {
                    typename Root::Data extra_obs(1, +PREALLOC_SEQ_LENGTH, Base::obs_dims(0u));
                    extra_obs.setZero();
                    obs = typename Root::Data(obs.concatenate(std::move(extra_obs), 1));
                }
                obs(0u,time_step++,ind) = (Scalar) 1;
            }
            if (fixed_seq_length) {
                for (; time_step < seq_length; ++time_step)
                    obs(0u,time_step,+PAD_IND) = (Scalar) 1;
            } else {
                if (time_step < obs.dimension(1)) {
                    RankwiseArray offsets({ 0u, 0u, 0u });
                    RankwiseArray extents({ 1u, time_step, Base::obs_dims(0u) });
                    obs = typename Root::Data(obs.slice(offsets, extents));
                }
            }
            return std::make_pair(std::move(obs), std::move(obj));
        }
        inline std::size_t _skip(std::ifstream& file_stream, std::size_t instances) {
            file_stream.seekg(0, std::ios::end);
            return instances - 1;
        }
    private:
        const VocabSharedPtr vocab;
        const std::size_t seq_length;
    };

    /**
    * An abstract class template for a data provider backed by an arbitrary number of file pairs
    * containing the separated observations and the objectives. Implementations are responsible for
    * specifying the dimensions of both the observations and the objectives, for reading batches of
    * observation-objective pairs from the file, and for skipping arbitrary number of data instances.
    */
    template<typename Scalar, std::size_t Rank, bool Sequential, bool ObsBinary = false, bool ObjBinary = false>
    class SplitFileDataProvider : public DataProvider<Scalar,Rank,Sequential> {
        typedef DataProvider<Scalar,Rank,Sequential> Base;
        typedef std::pair<std::string,std::string> FilePair;
        typedef std::pair<std::ifstream,std::ifstream> FileStreamPair;
    public:
        virtual ~SplitFileDataProvider() = default;
        inline const typename Base::Dims& get_obs_dims() const {
            return obs_dims;
        }
        inline const typename Base::Dims& get_obj_dims() const {
            return obj_dims;
        }
        inline bool has_more() {
            if (current_file_stream_pair_has_more())
                return true;
            ++current_file_pair_ind;
            for (; current_file_pair_ind < file_pairs.size(); ++current_file_pair_ind) {
                init_current_file_stream_pair();
                if (current_file_stream_pair_has_more())
                    return true;
            }
            return false;
        }
        inline DataPair<Scalar,Rank,Sequential> get_data(std::size_t batch_size) {
            if (!has_more())
                throw std::out_of_range("no more data left to fetch");
            const FilePair& current_file_pair = file_pairs[current_file_pair_ind];
            DataPair<Scalar,Rank,Sequential> data_pair = _get_data(current_file_pair.first,
                    current_file_stream_pair.first, current_file_pair.second, current_file_stream_pair.second,
                    batch_size);
            assert(data_pair.first.dimension(0) == data_pair.second.dimension(0));
            while (data_pair.first.dimension(0) < batch_size && has_more()) {
                const FilePair& new_current_file_pair = file_pairs[current_file_pair_ind];
                DataPair<Scalar,Rank,Sequential> add_data_pair = _get_data(new_current_file_pair.first,
                        current_file_stream_pair.first, new_current_file_pair.second, current_file_stream_pair.second,
                        batch_size - data_pair.first.dimension(0));
                assert(add_data_pair.first.dimension(0) == add_data_pair.second.dimension(0));
                typename Base::Data obs_concat = data_pair.first.concatenate(std::move(add_data_pair.first), 0);
                data_pair.first = std::move(obs_concat);
                typename Base::Data obj_concat = data_pair.second.concatenate(std::move(add_data_pair.second), 0);
                data_pair.second = std::move(obj_concat);
            }
            return data_pair;
        }
        inline void reset() {
            current_file_pair_ind = 0;
            init_current_file_stream_pair();
            _set_to_beg(current_file_stream_pair.first, current_file_stream_pair.second);
        }
        inline void skip(std::size_t instances) {
            if (!has_more())
                return;
            std::size_t skipped = _skip(current_file_stream_pair.first,current_file_stream_pair.second,
                    instances);
            while (skipped < instances && has_more())
                skipped += _skip(current_file_stream_pair.first, current_file_stream_pair.second,
                        instances - skipped);
        }
    protected:
        inline SplitFileDataProvider(const typename Base::Dims& obs_dims, const typename Base::Dims& obj_dims,
                const std::vector<FilePair>& dataset_path_pairs) :
                    obs_dims(obs_dims),
                    obj_dims(obj_dims),
                    file_pairs(dataset_path_pairs),
                    current_file_pair_ind(0) {
            assert(!dataset_path_pairs.empty());
            init_current_file_stream_pair();
        }
        inline SplitFileDataProvider(const typename Base::Dims& obs_dims, const typename Base::Dims& obj_dims,
                FilePair dataset_path_pair) :
                    SplitFileDataProvider(obs_dims, obj_dims, std::vector<FilePair>({ dataset_path_pair })) { }
        /**
        * It sets the positions of the file streams to the beginning of the observation data set and
        * the objective data set respectively.
        *
        * @param obs_file_stream A reference to the file stream to a file containing observations.
        * @param obj_file_stream A reference to the file stream to a file containing objectives.
        */
        virtual inline void _set_to_beg(std::ifstream& obs_file_stream, std::ifstream& obj_file_stream) {
            obs_file_stream.seekg(0, std::ios::beg);
            obj_file_stream.seekg(0, std::ios::beg);
        }
        /**
        * It reads at most the specified number of observations from the observation-file and at
        * most the specified number of objectives from the objective-file. The file streams can
        * be expected not to have any of their fail flags set initially and to have at least 1
        * more character left to read in each.
        *
        * @param obs_file The name of the observation source file.
        * @param obs_file_stream The input stream of the observation file.
        * @param obj_file The name of the objective source file.
        * @param obj_file_stream The input stream of the objective file.
        * @param batch_size The number of data points to read.
        * @return The paired observations and objectives.
        */
        virtual DataPair<Scalar,Rank,Sequential> _get_data(const std::string& obs_file,
                std::ifstream& obs_file_stream, const std::string& obj_file,
                std::ifstream& obj_file_stream, std::size_t batch_size) = 0;
        /**
        * Skips at most the specified number of instances in the data streams. The file streams can
        * be expected not to have any of their fail flags set initially.
        *
        * @param obs_file_stream A reference to the file stream to a file containing observations.
        * @param obj_file_stream A reference to the file stream to a file containing objectives.
        * @param instances The number of data points to skip.
        * @return The number of actual data points skipped. It may be less than the specified
        * amount if there are fewer remaining instances in the data streams.
        */
        virtual std::size_t _skip(std::ifstream& obs_file_stream, std::ifstream& obj_file_stream,
                std::size_t instances) = 0;
        const typename Base::Dims obj_dims, obs_dims;
    private:
        bool current_file_stream_pair_has_more() {
            return current_file_stream_pair.first && current_file_stream_pair.first.peek() != EOF &&
                    current_file_stream_pair.second && current_file_stream_pair.second.peek() != EOF;
        }
        void init_current_file_stream_pair() {
            const FilePair& file_pair = file_pairs[current_file_pair_ind];
            std::ifstream obs_stream(file_pair.first, ObsBinary ? std::ios::binary : std::ios::in);
            assert(obs_stream.is_open());
            std::ifstream obj_stream(file_pair.second, ObjBinary ? std::ios::binary : std::ios::in);
            assert(obj_stream.is_open());
            current_file_stream_pair = std::make_pair(std::move(obs_stream), std::move(obj_stream));
        }
        std::vector<FilePair> file_pairs;
        std::size_t current_file_pair_ind;
        FileStreamPair current_file_stream_pair;
    };

    /**
    * A data provider template for the MNIST data set.
    *
    * \see http://yann.lecun.com/exdb/mnist/
    */
    template<typename Scalar>
    class MNISTDataProvider : public SplitFileDataProvider<Scalar,3,false,true,true> {
        typedef DataProvider<Scalar,3,false> Root;
        typedef SplitFileDataProvider<Scalar,3,false,true,true> Base;
        typedef std::array<std::size_t,4> RankwiseArray;
        static constexpr std::size_t OBS_OFFSET = 16;
        static constexpr std::size_t LABEL_OFFSET = 8;
        static constexpr std::size_t OBS_INSTANCE_LENGTH = 784;
        static constexpr std::size_t LABEL_INSTANCE_LENGTH = 1;
    public:
        /**
        * @param obs_path The path to the file containing the observations.
        * @param labels_path The path to the file containing the corresponding labels.
        */
        MNISTDataProvider(std::string obs_path, std::string labels_path) :
                Base::SplitFileDataProvider({ 28u, 28u, 1u }, { 10u, 1u, 1u }, std::make_pair(obs_path, labels_path)),
                offsets({ 0u, 0u, 0u, 0u }),
                obs_extents(Base::obs_dims.template promote<>()),
                obj_extents(Base::obs_dims.template promote<>()) {
            Base::reset();
        }
    protected:
        inline void _set_to_beg(std::ifstream& obs_file_stream, std::ifstream& obj_file_stream) {
            Base::_set_to_beg(obs_file_stream, obj_file_stream);
            obs_file_stream.ignore(OBS_OFFSET);
            obj_file_stream.ignore(LABEL_OFFSET);
        }
        inline DataPair<Scalar,3,false> _get_data(const std::string& obs_file,
                std::ifstream& obs_file_stream, const std::string& obj_file,
                std::ifstream& obj_file_stream, std::size_t batch_size) {
            typename Root::Data obs(batch_size, Base::obs_dims(0), Base::obs_dims(1), Base::obs_dims(2));
            typename Root::Data obj(batch_size, Base::obj_dims(0), Base::obj_dims(1), Base::obj_dims(2));
            obj.setZero();
            std::size_t i;
            for (i = 0; i < batch_size && obs_file_stream.read(obs_buffer, OBS_INSTANCE_LENGTH); ++i) {
                // Read and set the label.
                char label;
                obj_file_stream.read(&label, LABEL_INSTANCE_LENGTH);
                obj(i,static_cast<std::size_t>(label),0u,0u) = (Scalar) 1;
                // Set the image.
                unsigned char* u_buffer = reinterpret_cast<unsigned char*>(obs_buffer);
                std::size_t buffer_ind = 0;
                for (std::size_t height = 0; height < Base::obs_dims(0); ++height) {
                    for (std::size_t width = 0; width < Base::obs_dims(1); ++width)
                        obs(i,height,width,0u) = (Scalar) u_buffer[buffer_ind++];
                }
                assert(buffer_ind == OBS_INSTANCE_LENGTH);
            }
            if (i == batch_size)
                return std::make_pair(obs, obj);
            obs_extents[0] = i;
            obj_extents[0] = i;
            return std::make_pair(obs.slice(offsets, obs_extents), obj.slice(offsets, obj_extents));
        }
        inline std::size_t _skip(std::ifstream& obs_file_stream, std::ifstream& obj_file_stream,
                std::size_t instances) {
            // Skip observations.
            std::streampos curr_obs_pos = obs_file_stream.tellg();
            obs_file_stream.seekg(0, std::ios::end);
            std::size_t obs_skip_extent = obs_file_stream.tellg() - curr_obs_pos;
            obs_file_stream.seekg(curr_obs_pos);
            obs_file_stream.ignore(instances * OBS_INSTANCE_LENGTH);
            std::size_t skipped_obs = std::min(instances, obs_skip_extent / OBS_INSTANCE_LENGTH);
            // Skip labels.
            std::streampos curr_label_pos = obj_file_stream.tellg();
            obj_file_stream.seekg(0, std::ios::end);
            std::size_t label_skip_extent = obj_file_stream.tellg() - curr_label_pos;
            obj_file_stream.seekg(curr_label_pos);
            obj_file_stream.ignore(instances * LABEL_INSTANCE_LENGTH);
            std::size_t skipped_labels = std::min(instances, label_skip_extent / LABEL_INSTANCE_LENGTH);
            assert(skipped_obs == skipped_labels);
            return skipped_obs;
        }
    private:
        char obs_buffer[OBS_INSTANCE_LENGTH];
        RankwiseArray offsets, obs_extents, obj_extents;
    };

    /**
    * An alias for a unique pointer to a tensor.
    */
    template<typename Scalar, std::size_t Rank>
    using TensorPtr = std::unique_ptr<Tensor<Scalar,Rank>>;

    /**
    * A data provider that reads from the memory. It requires the entire observation and
    * objective data sets to be loaded into memory, but it fetches pairs faster.
    */
    template<typename Scalar, std::size_t Rank, bool Sequential, bool Shuffle = true>
    class MemoryDataProvider : public DataProvider<Scalar,Rank,Sequential> {
        typedef DataProvider<Scalar,Rank,Sequential> Base;
        typedef TensorPtr<Scalar,Base::DATA_RANK> DataPtr;
        typedef std::array<std::size_t,Base::DATA_RANK> RankwiseArray;
    public:
        /**
        * It constructs a data provider backed by the specified tensors.
        *
        * @param obs A unique pointer to the tensor containing the observations.
        * @param obj A unique pointer to the tensor containing the objectives.
        * shuffled.
        */
        inline MemoryDataProvider(DataPtr obs, DataPtr obj) :
                obs(std::move(obs)),
                obj(std::move(obj)),
                offsets() {
            assert(this->obs != nullptr);
            assert(this->obj != nullptr);
            assert(this->obs->dimension(0) == this->obj->dimension(0) &&
                    "mismatched data and obj tensor row numbers");
            Dimensions<std::size_t,Base::DATA_RANK> obs_dims = this->obs->dimensions();
            Dimensions<std::size_t,Base::DATA_RANK> obj_dims = this->obj->dimensions();
            this->obs_dims = obs_dims.template demote<Sequential + 1>();
            this->obj_dims = obj_dims.template demote<Sequential + 1>();
            instances = (std::size_t) this->obs->dimension(0);
            offsets.fill(0);
            data_extents = obs_dims;
            obj_extents = obj_dims;
            if (Shuffle)
                shuffle_tensor_rows();
        }
        inline const typename Base::Dims& get_obs_dims() const {
            return obs_dims;
        }
        inline const typename Base::Dims& get_obj_dims() const {
            return obj_dims;
        }
        inline bool has_more() {
            return offsets[0] < (int) instances;
        }
        inline DataPair<Scalar,Rank,Sequential> get_data(std::size_t batch_size = std::numeric_limits<std::size_t>::max()) {
            if (!has_more())
                throw std::out_of_range("no more data left to fetch");
            std::size_t max_batch_size = std::min(batch_size, instances - offsets[0]);
            data_extents[0] = max_batch_size;
            obj_extents[0] = max_batch_size;
            typename Base::Data data_batch = obs->slice(offsets, data_extents);
            typename Base::Data obj_batch = obj->slice(offsets, obj_extents);
            offsets[0] = std::min(instances, offsets[0] + max_batch_size);
            return std::make_pair(std::move(data_batch), std::move(obj_batch));
        }
        inline void reset() {
            offsets[0] = 0;
        }
        inline void skip(std::size_t instances) {
            offsets[0] = std::min((int) this->instances, (int) (offsets[0] + instances));
        }
    protected:
        inline void shuffle_tensor_rows() {
            std::size_t rows = obs->dimension(0);
            MatrixMap<Scalar> obs_mat(obs->data(), rows, obs->size() / rows);
            MatrixMap<Scalar> obj_mat(obj->data(), rows, obj->size() / rows);
            // Create an identity matrix.
            PermMatrix perm(rows);
            perm.setIdentity();
            // Shuffle its indices.
            std::random_shuffle(perm.indices().data(), perm.indices().data() + perm.indices().size());
            // And apply the same row permutation transformation to both the observations and the objectives.
            obs_mat = perm * obs_mat;
            obj_mat = perm * obj_mat;
        }
    private:
        DataPtr obs, obj;
        typename Base::Dims obs_dims, obj_dims;
        std::size_t instances;
        RankwiseArray offsets, data_extents, obj_extents;
    };

    /**
    * A wrapper class template for data providers associated with continuous partitions of other data
    * providers. It enables the partitioning of a data provider into training and test data providers
    * by mapping two contiguous blocks of its data to two PartitionedDataProvider instances.
    */
    template<typename Scalar, std::size_t Rank, bool Sequential>
    class PartitionDataProvider : public DataProvider<Scalar,Rank,Sequential> {
        typedef DataProvider<Scalar,Rank,Sequential> Base;
    public:
        inline PartitionDataProvider(Base& orig_provider, std::size_t offset, std::size_t length) :
                orig_provider(orig_provider),
                offset(offset),
                length(length) {
            assert(length > 0);
            reset();
        }
        inline const typename Base::Dims& get_obs_dims() const {
            return orig_provider.get_obs_dims();
        }
        inline const typename Base::Dims& get_obj_dims() const {
            return orig_provider.get_obj_dims();
        }
        inline bool has_more() {
            return instances_read < length && orig_provider.has_more();
        }
        inline DataPair<Scalar,Rank,Sequential> get_data(std::size_t batch_size) {
            std::size_t instances_to_read = std::min(batch_size, length - instances_read);
            instances_read += instances_to_read;
            return orig_provider.get_data(instances_to_read);
        }
        inline void reset() {
            orig_provider.reset();
            orig_provider.skip(offset);
            instances_read = 0;
        }
        inline void skip(std::size_t instances) {
            orig_provider.skip(instances);
        }
    private:
        Base& orig_provider;
        const std::size_t offset, length;
        std::size_t instances_read;
    };

    /**
    * An alias for a shared pointer to a Parameters instance.
    */
    template<typename Scalar>
    using ParamsSharedPtr = std::shared_ptr<Parameters<Scalar>>;

    /**
    * An abstract class template that represents an activation function layer.
    */
    template<typename Scalar, std::size_t Rank>
    class ActivationLayer : public Layer<Scalar,Rank> {
        typedef Layer<Scalar,Rank> Base;
        typedef ActivationLayer<Scalar,Rank> Self;
    public:
        virtual ~ActivationLayer() = default;
        virtual Base* clone() const = 0;
        inline Base* clone_with_shared_params() {
            return clone();
        }
        inline const Base& get_params_owner() const {
            return owner;
        }
        inline const typename Base::Dims& get_input_dims() const {
            return dims;
        }
        inline const typename Base::Dims& get_output_dims() const {
            return dims;
        }
        inline bool is_input_layer() const {
            return input_layer;
        }
        inline void set_input_layer(bool input_layer) {
            this->input_layer = input_layer;
        }
        inline std::vector<const Parameters<Scalar>*> get_params() const {
            return params ? std::vector<const Parameters<Scalar>*>({ params.get() }) :
                    std::vector<const Parameters<Scalar>*>(0);
        }
        inline std::vector<Parameters<Scalar>*> get_params() {
            return params ? std::vector<Parameters<Scalar>*>({ params.get() }) :
                    std::vector<Parameters<Scalar>*>(0);
        }
    protected:
        /**
        * @param dims The input and output dimensions of the layer.
        * @param params The parameters of the layer; it can be null if the layer is not parametric.
        */
        inline ActivationLayer(const typename Base::Dims& dims, ParamsSharedPtr<Scalar> params = nullptr) :
                owner(*this),
                dims(dims),
                params(params),
                input_layer(false) { }
        inline ActivationLayer(const Self& layer, bool share_params = false) :
                owner(share_params && layer.params ? layer.owner : *this),
                dims(layer.dims),
                params(share_params ? layer.params : (!layer.params ? nullptr :
                        ParamsSharedPtr<Scalar>(layer.params->clone()))),
                input_layer(layer.input_layer) { }
        const Self& owner;
        const typename Base::Dims dims;
        ParamsSharedPtr<Scalar> params;
        bool input_layer;
    };    

    /**
    * A class template that represents a binary step activation function that outputs either
    * 1 or 0 based on the signum of its input. This function is not differentiable.
    *
    * \f[
    *   f(x) = \begin{cases}
    *     0 & \text{for } x < 0\\
    *     1 & \text{for } x \geq 0
    *   \end{cases}
    * \f]
    */
    template<typename Scalar, std::size_t Rank>
    class BinaryStepActivationLayer : public ActivationLayer<Scalar,Rank> {
        typedef Layer<Scalar,Rank> Root;
        typedef ActivationLayer<Scalar,Rank> Base;
    public:
        /**
        * @param dims The dimensionality of the input tensor.
        */
        inline BinaryStepActivationLayer(const typename Root::Dims& dims) :
                Base::ActivationLayer(dims) { }
        inline Root* clone() const {
            return new BinaryStepActivationLayer(*this);
        }
        inline void empty_cache() { }
        inline typename Root::Data pass_forward(typename Root::Data in, bool training) {
            assert((Dimensions<std::size_t,Base::DATA_RANK>(in.dimensions()).template demote<>()) == Base::dims);
            assert(in.dimension(0) > 0);
            batch_size = in.dimension(0);
            return in.unaryExpr([](Scalar e) { return (Scalar) (e >= 0); });
        }
        inline typename Root::Data pass_back(typename Root::Data out_grad) {
            assert((Dimensions<std::size_t,Base::DATA_RANK>(out_grad.dimensions()).template demote<>()) == Base::dims);
            assert(out_grad.dimension(0) > 0 && batch_size == out_grad.dimension(0));
            if (Base::is_input_layer())
                return typename Root::Data();
            return out_grad.constant(0);
        }
    private:
        std::size_t batch_size;
    };    

    /**
    * A class template representing an exponential linear unit (ELU) activation function. ELUs
    * apply an exponential (e based) function scaled by alpha to the negative elements of the input.
    * ELU layers are not differentiable.
    *
    * \f[
    *   f(x) = \begin{cases}
    *     \alpha (e^x - 1) & \text{for } x < 0\\
    *     x & \text{for } x \geq 0
    *   \end{cases}
    * \f]
    *
    * \see https://arxiv.org/abs/1511.07289
    */
    template<typename Scalar, std::size_t Rank>
    class ELUActivationLayer : public ActivationLayer<Scalar,Rank> {
        typedef Layer<Scalar,Rank> Root;
        typedef ActivationLayer<Scalar,Rank> Base;
        typedef std::array<std::size_t,Root::DATA_RANK> RankwiseArray;
    public:
        /**
        * @param dims The dimensionality of the input tensor.
        * @param alpha The factor by which negative inputs are to be scaled.
        */
        inline ELUActivationLayer(const typename Root::Dims& dims, Scalar alpha = 1e-1) :
                Base::ActivationLayer(dims),
                alpha(alpha),
                conversion_dims(dims.template promote<>()) { }
        inline Root* clone() const {
            return new ELUActivationLayer(*this);
        }
        inline void empty_cache() {
            in_mat_cache = Matrix<Scalar>();
            out_mat_cache = Matrix<Scalar>();
        }
        inline typename Root::Data pass_forward(typename Root::Data in, bool training) {
            assert((Dimensions<std::size_t,Base::DATA_RANK>(in.dimensions()).template demote<>()) == Base::dims);
            assert(in.dimension(0) > 0);
            conversion_dims[0] = in.dimension(0);
            in_mat_cache = MatrixMap<Scalar>(in.data(), in.dimension(0), Base::dims.get_volume());
            out_mat_cache = in_mat_cache.unaryExpr([this](Scalar e) {
                return (Scalar) (e >= 0 ? e : (alpha * (exp(e) - 1)));
            });
            return TensorMap<Scalar,Root::DATA_RANK>(out_mat_cache.data(), conversion_dims);
        }
        inline typename Root::Data pass_back(typename Root::Data out_grad) {
            assert((Dimensions<std::size_t,Base::DATA_RANK>(out_grad.dimensions()).template demote<>()) == Base::dims);
            assert(out_grad.dimension(0) > 0 && conversion_dims[0] == out_grad.dimension(0));
            if (Base::is_input_layer())
                return typename Root::Data();
            MatrixMap<Scalar> out_grad_mat(out_grad.data(), conversion_dims[0], Base::dims.get_volume());
            Matrix<Scalar> prev_out_grad_mat(in_mat_cache.rows(), in_mat_cache.cols());
            for (int i = 0; i < in_mat_cache.cols(); ++i) {
                for (int j = 0; j < in_mat_cache.rows(); ++j)
                    prev_out_grad_mat(j,i) = (Scalar) ((in_mat_cache(j,i) >= 0 ?
                            1 : (out_mat_cache(j,i) + alpha)) * out_grad_mat(j,i));
            }
            return TensorMap<Scalar,Root::DATA_RANK>(prev_out_grad_mat.data(), conversion_dims);
        }
    private:
        const Scalar alpha;
        RankwiseArray conversion_dims;
        // Staged computation caches.
        Matrix<Scalar> in_mat_cache, out_mat_cache;
    };

    /**
    * A class template representing an identity activation layer that merely outputs
    * its input.
    *
    * \f$f(x) = x\f$
    */
    template<typename Scalar, std::size_t Rank>
    class IdentityActivationLayer : public ActivationLayer<Scalar,Rank> {
        typedef Layer<Scalar,Rank> Root;
        typedef ActivationLayer<Scalar,Rank> Base;
    public:
        /**
        * @param dims The dimensionality of the input tensor.
        */
        inline IdentityActivationLayer(const typename Root::Dims& dims) :
                Base::ActivationLayer(dims) { }
        inline Root* clone() const {
            return new IdentityActivationLayer(*this);
        }
        inline void empty_cache() { }
        inline typename Root::Data pass_forward(typename Root::Data in, bool training) {
            assert((Dimensions<std::size_t,Base::DATA_RANK>(in.dimensions()).template demote<>()) == Base::dims);
            assert(in.dimension(0) > 0);
            batch_size = in.dimension(0);
            return in;
        }
        inline typename Root::Data pass_back(typename Root::Data out_grad) {
            assert((Dimensions<std::size_t,Base::DATA_RANK>(out_grad.dimensions()).template demote<>()) == Base::dims);
            assert(out_grad.dimension(0) > 0 && batch_size == out_grad.dimension(0));
            if (Base::is_input_layer())
                return typename Root::Data();
            return out_grad;
        }
    private:
        std::size_t batch_size;
    };


    /**
    * A class template representing a leaky rectified linear unit activation function. Unlike
    * traditional ReLU layers leaky ReLU layers do not set negative elements of the input to
    * 0 but scale them by a small constant alpha. This function is not differentiable.
    *
    * \f[
    *   f(x) = \begin{cases}
    *     \alpha x & \text{for } x < 0\\
    *     x & \text{for } x \geq 0
    *   \end{cases}
    * \f]
    *
    * \see https://ai.stanford.edu/~amaas/papers/relu_hybrid_icml2013_final.pdf
    */
    template<typename Scalar, std::size_t Rank>
    class LeakyReLUActivationLayer : public ActivationLayer<Scalar,Rank> {
        typedef Layer<Scalar,Rank> Root;
        typedef ActivationLayer<Scalar,Rank> Base;
    public:
        /**
        * @param dims The dimensionality of the input tensor.
        * @param alpha The factor by which negative inputs are to be scaled.
        */
        inline LeakyReLUActivationLayer(const typename Root::Dims& dims, Scalar alpha = 1e-1) :
                Base::ActivationLayer(dims),
                alpha(alpha) { }
        inline Root* clone() const {
            return new LeakyReLUActivationLayer(*this);
        }
        inline void empty_cache() {
            in_cache = typename Root::Data();
        }
        inline typename Root::Data pass_forward(typename Root::Data in, bool training) {
            assert((Dimensions<std::size_t,Base::DATA_RANK>(in.dimensions()).template demote<>()) == Base::dims);
            assert(in.dimension(0) > 0);
            in_cache = std::move(in);
            return in_cache.cwiseMax(in_cache * alpha);
        }
        inline typename Root::Data pass_back(typename Root::Data out_grad) {
            assert((Dimensions<std::size_t,Base::DATA_RANK>(out_grad.dimensions()).template demote<>()) == Base::dims);
            assert(out_grad.dimension(0) > 0 && in_cache.dimension(0) == out_grad.dimension(0));
            if (Base::is_input_layer())
                return typename Root::Data();
            return in_cache.unaryExpr([this](Scalar e) { return (Scalar) (e >= 0 ? 1 : alpha); }) * out_grad;
        }
    private:
        const Scalar alpha;
        typename Root::Data in_cache;
    };

    /**
    * A class template for a parameter initialization that sets all values to a constant.
    */
    template<typename Scalar>
    class ConstantParameterInitialization : public ParameterInitialization<Scalar> {
    public:
        /**
        * @param constant The value to which all elements of the parameter matrix are to be
        * initialized.
        */
        ConstantParameterInitialization(Scalar constant) :
                constant(constant) { }
        inline void apply(Matrix<Scalar>& params) const {
            params.setConstant(params.rows(), params.cols(), constant);
        }
    private:
        Scalar constant;
    };

    /**
    * A class template for a parameter initialization that sets all values to 0.
    */
    template<typename Scalar>
    class ZeroParameterInitialization : public ParameterInitialization<Scalar> {
    public:
        inline void apply(Matrix<Scalar>& params) const {
            params.setZero(params.rows(), params.cols());
        }
    };

    /**
    * A class template for a parameter initialization that sets all values to 1.
    */
    template<typename Scalar>
    class OneParameterInitialization : public ParameterInitialization<Scalar> {
    public:
        inline void apply(Matrix<Scalar>& params) const {
            params.setOnes(params.rows(), params.cols());
        }
    };

    /**
    * An abstract class template representing a random parameter initialization method
    * that samples from a Gaussian distribution.
    */
    template<typename Scalar>
    class GaussianParameterInitialization : public ParameterInitialization<Scalar> {
    public:
        /**
        * @param mean The mean of the distribution.
        * @param sd_scaling_factor The standard deviation scaling factor.
        */
        inline GaussianParameterInitialization(Scalar mean = 0, Scalar sd_scaling_factor = 1) :
                mean(mean),
                sd_scaling_factor(sd_scaling_factor) {
            assert(sd_scaling_factor > 0);
        }
        virtual ~GaussianParameterInitialization() = default;
        inline virtual void apply(Matrix<Scalar>& params) const {
            int rows = params.rows();
            int cols = params.cols();
            std::default_random_engine gen;
            std::normal_distribution<Scalar> dist(mean, sd_scaling_factor * _sd(rows, cols));
            for (int i = 0; i < cols; ++i) {
                for (int j = 0; j < rows; ++j)
                    params(j,i) = dist(gen);
            }
        }
    protected:
        /**
        * It computes the standard deviation of the distribution to sample from.
        *
        * @param fan_ins The input size of the kernel (the number of rows in the weight matrix
        * excluding the bias row).
        * @param fan_outs The output size of the kernel (the number of columns in the weight
        * matrix).
        * @return The standard deviation of the normal distribution from which the values
        * of the initialized weight matrix are to be sampled.
        */
        inline virtual Scalar _sd(unsigned fan_ins, unsigned fan_outs) const {
            return 1;
        }
    private:
        const Scalar mean, sd_scaling_factor;
    };

    /**
    * An abstract class representing the Xavier/Glorot weight initialization method.
    *
    * \f$\sigma = c \sqrt{\frac{2}{fan_{in} + fan_{out}}}\f$
    *
    * \see http://proceedings.mlr.press/v9/glorot10a/glorot10a.pdf
    */
    template<typename Scalar>
    class GlorotParameterInitialization : public GaussianParameterInitialization<Scalar> {
    public:
        /**
        * @param sd_scaling_factor The value by which the randomly initialized weights
        * are to be scaled.
        */
        inline GlorotParameterInitialization(Scalar sd_scaling_factor = 1) :
                GaussianParameterInitialization<Scalar>(0, sd_scaling_factor) { }
    protected:
        inline Scalar _sd(unsigned fan_ins, unsigned fan_outs) const {
            return sqrt(2 / (Scalar) (fan_ins + fan_outs));
        }
    };

    /**
    * An abstract class representing the He weight initialization method.
    *
    * \f$\sigma = c \sqrt{\frac{2}{fan_{in}}}\f$
    *
    * \see https://arxiv.org/abs/1502.01852
    */
    template<typename Scalar>
    class HeParameterInitialization : public GaussianParameterInitialization<Scalar> {
    public:
        /**
        * @param sd_scaling_factor The value by which the randomly initialized weights
        * are to be scaled.
        */
        inline HeParameterInitialization(Scalar sd_scaling_factor = 1) :
                GaussianParameterInitialization<Scalar>(0, sd_scaling_factor) { }
    protected:
        inline Scalar _sd(unsigned fan_ins, unsigned fan_outs) const {
            return sqrt(2 / (Scalar) fan_ins);
        }
    };    

    /**
    * A weight initializer that assigns linearly increasing values to the elements
    * of the parameter matrix. It is meant to be used for testing.
    */
    template<typename Scalar>
    class IncrementalParameterInitialization : public ParameterInitialization<Scalar> {
    public:
        /**
        * @param start The starting value.
        * @param inc The value by which the parameter value is to be incremented.
        */
        inline IncrementalParameterInitialization(Scalar start, Scalar inc) :
                start(start),
                inc(inc) { }
        inline void apply(Matrix<Scalar>& params) const {
            Scalar val = start;
            for (std::size_t i = 0; i < params.cols(); ++i) {
                for (std::size_t j = 0; j < params.rows(); ++j) {
                    params(j,i) = val;
                    val += inc;
                }
            }
        }
    private:
        const Scalar start;
        const Scalar inc;
    };

    /**
    * An abstract class representing the LeCun weight initialization method.
    *
    * \f$\sigma = c \sqrt{\frac{1}{fan_{in}}}\f$
    *
    * \see http://yann.lecun.com/exdb/publis/pdf/lecun-98b.pdf
    */
    template<typename Scalar>
    class LeCunParameterInitialization : public GaussianParameterInitialization<Scalar> {
    public:
        /**
        * @param sd_scaling_factor The value by which the randomly initialized weights
        * are to be scaled.
        */
        inline LeCunParameterInitialization(Scalar sd_scaling_factor = 1) :
                GaussianParameterInitialization<Scalar>(0, sd_scaling_factor) { }
    protected:
        inline Scalar _sd(unsigned fan_ins, unsigned fan_outs) const {
            return sqrt(1 / (Scalar) fan_ins);
        }
    };

    /**
    * A class template representing the orthogonal weight initialization method.
    *
    * \see https://arxiv.org/abs/1312.6120
    */
    template<typename Scalar>
    class OrthogonalParameterInitialization : public GaussianParameterInitialization<Scalar> {
    public:
        /**
        * @param sd The standard deviation of the normal distribution to sample from.
        */
        inline OrthogonalParameterInitialization(Scalar sd = 1) :
                GaussianParameterInitialization<Scalar>(0, sd) { }
        inline void apply(Matrix<Scalar>& params) const {
            GaussianParameterInitialization<Scalar>::apply(params);
            int rows = params.rows();
            int cols = params.cols();
            bool more_rows = rows > cols;
            SVD<Scalar> svd;
            params.block(0, 0, rows, cols) = more_rows ?
                    svd.compute(params, SVDOptions::ComputeFullU).matrixU().block(0, 0, rows, cols) :
                    svd.compute(params, SVDOptions::ComputeFullV).matrixV().block(0, 0, rows, cols);
        }
    };
    /**
    * An alias for a shared pointer to a WeightInitialization implementation instance of
    * an arbitrary scalar type.
    */
    template<typename Scalar>
    using ParamInitSharedPtr = std::shared_ptr<ParameterInitialization<Scalar>>;

    /**
    * An alias for a shared pointer to a regularization penalty of an arbitrary scalar type.
    */
    template<typename Scalar>
    using ParamRegSharedPtr = std::shared_ptr<ParameterRegularization<Scalar>>;

    template<typename Scalar>
    class StandardParameters : public Parameters<Scalar> {
    public:
        /**
        * @param rows The number of rows of the parameter matrix.
        * @param cols The number of columns of the parameter matrix.
        * @param optimizable Whether the parameters are optimizable. Non-optimizable
        * parameters do not maintain gradient information and are thus not regularizable.
        * @param init The parameter value initialization. If it is a null pointer, the
        * values will not be initialized.
        * @param reg The optional regularization to use on the values of the parameter
        * matrix. If it is a null pointer, no regularization is applied.
        * @param value_clip The maximum allowed absolute parameter value. If it is 0
        * or less, no value clipping is performed.
        * @param value_max_l1_norm The maximum allowed L1 parameter value norm. If it
        * is 0 or less, no L1 max norm constraint is enforced.
        * @param value_max_l2_norm The maximum allowed L2 parameter value norm. If it
        * is 0 or less, no L2 max norm constraint is enforced.
        * @param grad_clip The maximum allowed absolute parameter gradient. If it is 0
        * or less, no gradient clipping is performed.
        * @param grad_max_l1_norm The maximum allowed L1 parameter gradient norm. If it
        * is 0 or less, no L1 gradient max norm constraint is enforced.
        * @param grad_max_l2_norm The maximum allowed L2 parameter gradient norm. If it
        * is 0 or less, no L2 gradient max norm constraint is enforced.
        */
        inline StandardParameters(std::size_t rows, std::size_t cols, bool optimizable = true,
                ParamInitSharedPtr<Scalar> init = nullptr, ParamRegSharedPtr<Scalar> reg = nullptr,
                Scalar value_clip = 0, Scalar value_max_l1_norm = 0, Scalar value_max_l2_norm = 0,
                Scalar grad_clip = 0, Scalar grad_max_l1_norm = 0, Scalar grad_max_l2_norm = 0) :
                    rows(rows),
                    cols(cols),
                    optimizable(optimizable),
                    param_init(init),
                    param_reg(reg),
                    value_clip(value_clip),
                    value_max_l1_norm(value_max_l1_norm),
                    value_max_l2_norm(value_max_l2_norm),
                    grad_clip(grad_clip),
                    grad_max_l1_norm(grad_max_l1_norm),
                    grad_max_l2_norm(grad_max_l2_norm),
                    frozen(false) {
            assert(rows > 0 && cols > 0);
        }
        inline Parameters<Scalar>* clone() const {
            return new StandardParameters<Scalar>(*this);
        }
        inline bool are_optimizable() const {
            return optimizable;
        }
        inline std::size_t get_rows() const {
            return rows;
        }
        inline std::size_t get_cols() const {
            return cols;
        }
        inline void init_values() {
            values = Matrix<Scalar>(rows, cols);
            if (param_init)
                param_init->apply(values);
        }
        inline void init_grad() {
            if (optimizable) {
                grad = Matrix<Scalar>(rows, cols);
                reset_grad();
            }
        }
        inline const Matrix<Scalar>& get_values() const {
            return values;
        }
        inline void set_values(Matrix<Scalar> values) {
            assert(values.rows() == rows && values.cols() == cols);
            this->values = std::move(values);
            enforce_clip_constraint(this->values, value_clip);
            enforce_l1_norm_constraint(this->values, value_max_l1_norm);
            enforce_l2_norm_constraint(this->values, value_max_l2_norm);
        }
        inline const Matrix<Scalar>& get_grad() const {
            return grad;
        }
        inline void accumulate_grad(const Matrix<Scalar>& grad) {
            if (!optimizable)
                return;
            assert(grad.rows() == rows && grad.cols() == cols);
            this->grad += grad;
            enforce_clip_constraint(this->grad, grad_clip);
            enforce_l1_norm_constraint(this->grad, grad_max_l1_norm);
            enforce_l2_norm_constraint(this->grad, grad_max_l2_norm);
        }
        inline void reset_grad() {
            if (optimizable)
                grad.setZero();
        }
        inline Scalar get_regularization_penalty() const {
            if (optimizable && param_reg)
                return param_reg->function(values);
            return 0;
        }
        inline void regularize() {
            if (optimizable && param_reg)
                accumulate_grad(param_reg->d_function(values));
        }
        inline bool are_frozen() const {
            return frozen;
        }
        inline void set_frozen(bool frozen) {
            this->frozen = frozen;
        }
    protected:
        /**
        * It clips the values of the matrix falling outside the interval [-clip, clip].
        *
        * @param matrix The matrix on which the clip constraint is to be enforced.
        * @param clip The clipping limit.
        */
        inline static void enforce_clip_constraint(Matrix<Scalar>& matrix, Scalar clip) {
            if (NumericUtils<Scalar>::decidedly_greater(clip, (Scalar) 0)) {
                matrix = matrix.unaryExpr([clip](Scalar m) {
                    return m > clip ? clip : (m < -clip ? - clip : m);
                });
            }
        }
        /**
        * It limits the L1 norm of the matrix by scaling its coefficients.
        *
        * @param matrix The matrix whose L1 norm is to be constrained.
        * @param max_l1_norm The maximum allowed L1 norm.
        */
        inline static void enforce_l1_norm_constraint(Matrix<Scalar>& matrix, Scalar max_l1_norm) {
            if (NumericUtils<Scalar>::decidedly_greater(max_l1_norm, (Scalar) 0)) {
                Scalar l1_norm = matrix.norm();
                if (l1_norm > max_l1_norm)
                    matrix *= (max_l1_norm / l1_norm);
            }
        }
        /**
        * It limits the L2 norm of the matrix by scaling its coefficients.
        *
        * @param matrix The matrix whose L2 norm is to be constrained.
        * @param max_l2_norm The maximum allowed L2 norm.
        */
        inline static void enforce_l2_norm_constraint(Matrix<Scalar>& matrix, Scalar max_l2_norm) {
            if (NumericUtils<Scalar>::decidedly_greater(max_l2_norm, (Scalar) 0)) {
                Scalar l2_norm = matrix.squaredNorm();
                if (l2_norm > max_l2_norm)
                    matrix *= (max_l2_norm / l2_norm);
            }
        }
        const std::size_t rows, cols;
        const bool optimizable;
        const ParamInitSharedPtr<Scalar> param_init;
        const ParamRegSharedPtr<Scalar> param_reg;
        const Scalar value_clip, value_max_l1_norm, value_max_l2_norm;
        const Scalar grad_clip, grad_max_l1_norm, grad_max_l2_norm;
        Matrix<Scalar> values, grad;
        bool frozen;
    };    

    /**
    * A class template representing a parametric rectified linear unit (PReLU) activation function.
    * PReLU layers are Leaky ReLU activation functions with learnable alphas. PReLU activation
    * functions are not differentiable.
    *
    * \f[
    *   f(x) = \begin{cases}
    *     \alpha x & \text{for } x < 0\\
    *     x & \text{for } x \geq 0
    *   \end{cases}
    * \f]
    *
    * \see https://arxiv.org/abs/1502.01852
    */
    template<typename Scalar, std::size_t Rank>
    class PReLUActivationLayer : public ActivationLayer<Scalar,Rank> {
        typedef Layer<Scalar,Rank> Root;
        typedef ActivationLayer<Scalar,Rank> Base;
        typedef std::array<std::size_t,Root::DATA_RANK> RankwiseArray;
    public:
        /**
        * @param dims The dimensionality of the input tensor.
        * @param init_alpha The initial factor by which negative inputs are to be scaled.
        * @param alpha_reg An optional regularization function to apply to the parameters.
        * @param alpha_clip The maximum allowed absolute parameter value. If it is 0 or less, no
        * value clipping is performed.
        * @param alpha_max_l1_norm The maximum allowed L1 alpha value norm. If it is 0 or
        * less, no L1 max norm constraint is enforced.
        * @param alpha_max_l2_norm The maximum allowed L2 alpha value norm. If it is 0 or
        * less, no L2 max norm constraint is enforced.
        * @param alpha_grad_clip The maximum allowed absolute parameter gradient. If it is 0
        * or less, no gradient clipping is performed.
        * @param alpha_grad_max_l1_norm The maximum allowed L1 alpha gradient norm. If it
        * is 0 or less, no L1 gradient max norm constraint is enforced.
        * @param alpha_grad_max_l2_norm The maximum allowed L2 alpha gradient norm. If it
        * is 0 or less, no L2 gradient max norm constraint is enforced.
        */
        inline PReLUActivationLayer(const typename Root::Dims& dims, Scalar init_alpha = 1e-1,
                ParamRegSharedPtr<Scalar> alpha_reg = nullptr, Scalar alpha_clip = 0, Scalar alpha_max_l1_norm = 0,
                Scalar alpha_max_l2_norm = 0, Scalar alpha_grad_clip = 0, Scalar alpha_grad_max_l1_norm = 0,
                Scalar alpha_grad_max_l2_norm = 0) :
                    Base::ActivationLayer(dims, std::make_shared<StandardParameters<Scalar>>(1, dims.get_volume(),
                            true, std::make_shared<ConstantParameterInitialization<Scalar>>(init_alpha),
                            alpha_reg, alpha_clip, alpha_max_l1_norm, alpha_max_l2_norm, alpha_grad_clip,
                            alpha_grad_max_l1_norm, alpha_grad_max_l2_norm)),
                    conversion_dims(dims.template promote<>()) { }
        inline PReLUActivationLayer(const PReLUActivationLayer<Scalar,Rank>& layer, bool share_params = false) :
                Base::ActivationLayer(layer, share_params),
                conversion_dims(layer.conversion_dims),
                in_mat_cache(layer.in_mat_cache) { }
        inline Root* clone() const {
            return new PReLUActivationLayer(*this);
        }
        inline Root* clone_with_shared_params() {
            return new PReLUActivationLayer(*this, true);
        }
        inline void empty_cache() {
            in_mat_cache = Matrix<Scalar>();
        }
        inline typename Root::Data pass_forward(typename Root::Data in, bool training) {
            assert((Dimensions<std::size_t,Base::DATA_RANK>(in.dimensions()).template demote<>()) == Base::dims);
            assert(in.dimension(0) > 0);
            conversion_dims[0] = in.dimension(0);
            in_mat_cache = MatrixMap<Scalar>(in.data(), conversion_dims[0], in.size() / conversion_dims[0]);
            Matrix<Scalar> out_mat = in_mat_cache.cwiseMax(in_mat_cache * Base::params->get_values().asDiagonal());
            return TensorMap<Scalar,Root::DATA_RANK>(out_mat.data(), conversion_dims);
        }
        inline typename Root::Data pass_back(typename Root::Data out_grad) {
            assert((Dimensions<std::size_t,Base::DATA_RANK>(out_grad.dimensions()).template demote<>()) == Base::dims);
            assert(out_grad.dimension(0) > 0 && conversion_dims[0] == out_grad.dimension(0));
            MatrixMap<Scalar> out_grad_mat(out_grad.data(), conversion_dims[0], out_grad.size() / conversion_dims[0]);
            Matrix<Scalar> prev_out_grad_mat = Matrix<Scalar>(in_mat_cache.rows(), in_mat_cache.cols());
            const Matrix<Scalar>& alphas = Base::params->get_values();
            Matrix<Scalar> alphas_grad = Matrix<Scalar>::Zero(1, out_grad_mat.cols());
            for (int i = 0; i < in_mat_cache.cols(); ++i) {
                for (int j = 0; j < in_mat_cache.rows(); ++j) {
                    Scalar in_mat_ji = in_mat_cache(j,i);
                    Scalar out_mat_ji = out_grad_mat(j,i);
                    if (in_mat_ji >= 0)
                        prev_out_grad_mat(j,i) = out_mat_ji;
                    else {
                        prev_out_grad_mat(j,i) = alphas(0,i) * out_mat_ji;
                        alphas_grad(0,i) += in_mat_ji * out_mat_ji;
                    }
                }
            }
            Base::params->accumulate_grad(alphas_grad);
            return TensorMap<Scalar,Root::DATA_RANK>(prev_out_grad_mat.data(), conversion_dims);
        }
    private:
        RankwiseArray conversion_dims;
        // Staged computation cache.
        Matrix<Scalar> in_mat_cache;
    };    

    /**
    * A class template that represents a linearly scaled activation layer.
    *
    * \f$f(x) = c x\f$
    */
    template<typename Scalar, std::size_t Rank>
    class ScaledActivationLayer : public ActivationLayer<Scalar,Rank> {
        typedef Layer<Scalar,Rank> Root;
        typedef ActivationLayer<Scalar,Rank> Base;
    public:
        /**
        * @param dims The dimensionality of the input tensor.
        * @param scale The factor by which the input is to be scaled.
        */
        inline ScaledActivationLayer(const typename Root::Dims& dims, Scalar scale) :
                Base::ActivationLayer(dims),
                scale(scale) { }
        inline Root* clone() const {
            return new ScaledActivationLayer(*this);
        }
        inline void empty_cache() { }
        inline typename Root::Data pass_forward(typename Root::Data in, bool training) {
            assert((Dimensions<std::size_t,Base::DATA_RANK>(in.dimensions()).template demote<>()) == Base::dims);
            assert(in.dimension(0) > 0);
            batch_size = in.dimension(0);
            return in * scale;
        }
        inline typename Root::Data pass_back(typename Root::Data out_grad) {
            assert((Dimensions<std::size_t,Base::DATA_RANK>(out_grad.dimensions()).template demote<>()) == Base::dims);
            assert(out_grad.dimension(0) > 0 && batch_size == out_grad.dimension(0));
            if (Base::is_input_layer())
                return typename Root::Data();
            return out_grad * scale;
        }
    private:
        const Scalar scale;
        std::size_t batch_size;
    };

    /**
    * A class template representing a sigmoid activation function layer.
    *
    * \f$f(x) = \sigma(x) = \frac{1}{1 + e^{-x}}\f$
    */
    template<typename Scalar, std::size_t Rank>
    class SigmoidActivationLayer : public ActivationLayer<Scalar,Rank> {
        typedef Layer<Scalar,Rank> Root;
        typedef ActivationLayer<Scalar,Rank> Base;
    public:
        /**
        * @param dims The dimensionality of the input tensor.
        */
        inline SigmoidActivationLayer(const typename Root::Dims& dims) :
                Base::ActivationLayer(dims) { }
        inline Root* clone() const {
            return new SigmoidActivationLayer(*this);
        }
        inline void empty_cache() {
            out_cache = typename Root::Data();
        }
        inline typename Root::Data pass_forward(typename Root::Data in, bool training) {
            assert((Dimensions<std::size_t,Base::DATA_RANK>(in.dimensions()).template demote<>()) == Base::dims);
            assert(in.dimension(0) > 0);
            out_cache = ((-in).exp() + in.constant(1)).inverse();
            return out_cache;
        }
        inline typename Root::Data pass_back(typename Root::Data out_grad) {
            assert((Dimensions<std::size_t,Base::DATA_RANK>(out_grad.dimensions()).template demote<>()) == Base::dims);
            assert(out_grad.dimension(0) > 0 && out_cache.dimension(0) == out_grad.dimension(0));
            if (Base::is_input_layer())
                return typename Root::Data();
            return (out_cache * (out_cache.constant(1) - out_cache)) * out_grad;
        }
    private:
        // Staged computation cache.
        typename Root::Data out_cache;
    };

    /**
    * A class template for a softmax activation function layer. Unlike most other activation
    * functions, the softmax layer does not represent a simple coefficient-wise function but
    * a multivariate one. The per-sample sums of the elements of the output tensor of the layer
    * are always 1.
    *
    * \f$f(x_i) = \frac{e^{x_i}}{\epsilon + \sum\limits_{j = 1}^J e^{x_j}}\f$
    */
    template<typename Scalar, std::size_t Rank>
    class SoftmaxActivationLayer : public ActivationLayer<Scalar,Rank> {
        typedef Layer<Scalar,Rank> Root;
        typedef ActivationLayer<Scalar,Rank> Base;
        typedef std::array<std::size_t,Root::DATA_RANK> RankwiseArray;
    public:
        /**
        * @param dims The dimensionality of the input tensor.
        * @param epsilon A small constant to maintain numerical stability.
        */
        inline SoftmaxActivationLayer(const typename Root::Dims& dims,
                Scalar epsilon = NumericUtils<Scalar>::EPSILON2) :
                    Base::ActivationLayer(dims),
                    epsilon(epsilon),
                    conversion_dims(dims.template promote<>()) { }
        inline Root* clone() const {
            return new SoftmaxActivationLayer(*this);
        }
        inline void empty_cache() {
            out_mat_cache = Matrix<Scalar>();
        }
        inline typename Root::Data pass_forward(typename Root::Data in, bool training) {
            assert((Dimensions<std::size_t,Base::DATA_RANK>(in.dimensions()).template demote<>()) == Base::dims);
            assert(in.dimension(0) > 0);
            std::size_t rows = in.dimension(0);
            conversion_dims[0] = rows;
            MatrixMap<Scalar> in_mat(in.data(), rows, in.size() / rows);
            /* First subtract the value of the greatest coefficient from each element row-wise
            * to avoid an overflow due to raising e to great powers. */
            Matrix<Scalar> act_in_mat = (in_mat.array().colwise() - in_mat.array().rowwise().maxCoeff()).exp();
            act_in_mat = act_in_mat.array().colwise() / (act_in_mat.array().rowwise().sum() + epsilon);
            out_mat_cache = std::move(act_in_mat);
            return TensorMap<Scalar,Root::DATA_RANK>(out_mat_cache.data(), conversion_dims);
        }
        inline typename Root::Data pass_back(typename Root::Data out_grad) {
            assert((Dimensions<std::size_t,Base::DATA_RANK>(out_grad.dimensions()).template demote<>()) == Base::dims);
            assert(out_grad.dimension(0) > 0 && out_mat_cache.rows() == out_grad.dimension(0));
            if (Base::is_input_layer())
                return typename Root::Data();
            std::size_t rows = out_grad.dimension(0);
            MatrixMap<Scalar> out_grad_mat(out_grad.data(), rows, out_grad.size() / rows);
            Matrix<Scalar> prev_out_grad_mat(rows, out_mat_cache.cols());
            for (int i = 0; i < rows; ++i) {
                RowVector<Scalar> row_i = out_mat_cache.row(i);
                // FIXME Do not evaluate the expressions into a temporary variable.
                Matrix<Scalar> jacobian = row_i.asDiagonal();
                jacobian -= row_i.transpose() * row_i;
                prev_out_grad_mat.row(i) = out_grad_mat.row(i) * jacobian;
            }
            return TensorMap<Scalar,Root::DATA_RANK>(prev_out_grad_mat.data(), conversion_dims);
        }
    private:
        const Scalar epsilon;
        RankwiseArray conversion_dims;
        // Staged computation cache matrix.
        Matrix<Scalar> out_mat_cache;
    };


    /**
    * A class template representing a softplus activation function layer. The softplus activation function
    * is a differentiable function that approximates the rectified linear unit function.
    *
    * \f$f(x) = \ln(1 + e^x)\f$
    */
    template<typename Scalar, std::size_t Rank>
    class SoftplusActivationLayer : public ActivationLayer<Scalar,Rank> {
        typedef Layer<Scalar,Rank> Root;
        typedef ActivationLayer<Scalar,Rank> Base;
    public:
        /**
        * @param dims The dimensionality of the input tensor.
        */
        inline SoftplusActivationLayer(const typename Root::Dims& dims) :
                Base::ActivationLayer(dims) { }
        inline Root* clone() const {
            return new SoftplusActivationLayer(*this);
        }
        inline void empty_cache() {
            in_cache = typename Root::Data();
        }
        inline typename Root::Data pass_forward(typename Root::Data in, bool training) {
            assert((Dimensions<std::size_t,Base::DATA_RANK>(in.dimensions()).template demote<>()) == Base::dims);
            assert(in.dimension(0) > 0);
            in_cache = std::move(in);
            return (in_cache.exp() + in_cache.constant(1)).log();
        }
        inline typename Root::Data pass_back(typename Root::Data out_grad) {
            assert((Dimensions<std::size_t,Base::DATA_RANK>(out_grad.dimensions()).template demote<>()) == Base::dims);
            assert(out_grad.dimension(0) > 0 && in_cache.dimension(0) == out_grad.dimension(0));
            if (Base::is_input_layer())
                return typename Root::Data();
            return ((-in_cache).exp() + in_cache.constant(1)).inverse() * out_grad;
        }
    private:
        typename Root::Data in_cache;
    };

    /**
    * A class template representing a softsign activation function layer, an alternative to the
    * tanh layer.
    *
    * \f$f(x) = \frac{x}{1 + \left|x\right|}\f$
    */
    template<typename Scalar, std::size_t Rank>
    class SoftsignActivationLayer : public ActivationLayer<Scalar,Rank> {
        typedef Layer<Scalar,Rank> Root;
        typedef ActivationLayer<Scalar,Rank> Base;
    public:
        /**
        * @param dims The dimensionality of the input tensor.
        */
        inline SoftsignActivationLayer(const typename Root::Dims& dims) :
                Base::ActivationLayer(dims) { }
        inline Root* clone() const {
            return new SoftsignActivationLayer(*this);
        }
        inline void empty_cache() {
            denom_cache = typename Root::Data();
        }
        inline typename Root::Data pass_forward(typename Root::Data in, bool training) {
            assert((Dimensions<std::size_t,Base::DATA_RANK>(in.dimensions()).template demote<>()) == Base::dims);
            assert(in.dimension(0) > 0);
            denom_cache = in.constant(1) + in.abs();
            return in / denom_cache;
        }
        inline typename Root::Data pass_back(typename Root::Data out_grad) {
            assert((Dimensions<std::size_t,Base::DATA_RANK>(out_grad.dimensions()).template demote<>()) == Base::dims);
            assert(out_grad.dimension(0) > 0 && denom_cache.dimension(0) == out_grad.dimension(0));
            if (Base::is_input_layer())
                return typename Root::Data();
            return denom_cache.square().inverse() * out_grad;
        }
    private:
        // Staged computation cache.
        typename Root::Data denom_cache;
    };

    /**
    * A class template representing the Swish activation function.
    *
    * \f$f(x) = x \sigma(\beta x)\f$
    *
    * \see https://arxiv.org/abs/1710.05941
    */
    template<typename Scalar, std::size_t Rank>
    class SwishActivationLayer : public ActivationLayer<Scalar,Rank> {
        typedef Layer<Scalar,Rank> Root;
        typedef ActivationLayer<Scalar,Rank> Base;
    public:
        /**
        * @param dims The dimensionality of the input tensor.
        * @param beta The factor by which the input of the sigmoid factor of the Swish
        * function is to be scaled.
        */
        inline SwishActivationLayer(const typename Root::Dims& dims, Scalar beta = 1) :
                Base::ActivationLayer(dims),
                beta(beta) { }
        inline Root* clone() const {
            return new SwishActivationLayer(*this);
        }
        inline void empty_cache() {
            in_cache = typename Root::Data();
            sig_out_cache = typename Root::Data();
        }
        inline typename Root::Data pass_forward(typename Root::Data in, bool training) {
            assert((Dimensions<std::size_t,Base::DATA_RANK>(in.dimensions()).template demote<>()) == Base::dims);
            assert(in.dimension(0) > 0);
            sig_out_cache = ((-beta * in).exp() + in.constant(1)).inverse();
            in_cache = std::move(in);
            return in_cache * sig_out_cache;
        }
        inline typename Root::Data pass_back(typename Root::Data out_grad) {
            assert((Dimensions<std::size_t,Base::DATA_RANK>(out_grad.dimensions()).template demote<>()) == Base::dims);
            assert(out_grad.dimension(0) > 0 && sig_out_cache.dimension(0) == out_grad.dimension(0));
            if (Base::is_input_layer())
                return typename Root::Data();
            return sig_out_cache * ((sig_out_cache.constant(1) - sig_out_cache) * beta * in_cache +
                    sig_out_cache.constant(1)) * out_grad;
        }
    private:
        const Scalar beta;
        // Staged computation cache.
        typename Root::Data in_cache, sig_out_cache;
    };    

    /**
    * A class template representing a hyperbolic tangent activation function layer.
    *
    * \f$f(x) = \text{tanh}(x)\f$
    */
    template<typename Scalar, std::size_t Rank>
    class TanhActivationLayer : public ActivationLayer<Scalar,Rank> {
        typedef Layer<Scalar,Rank> Root;
        typedef ActivationLayer<Scalar,Rank> Base;
    public:
        /**
        * @param dims The dimensionality of the input tensor.
        */
        inline TanhActivationLayer(const typename Root::Dims& dims) :
                Base::ActivationLayer(dims) { }
        inline Root* clone() const {
            return new TanhActivationLayer(*this);
        }
        inline void empty_cache() {
            out_cache = typename Root::Data();
        }
        inline typename Root::Data pass_forward(typename Root::Data in, bool training) {
            assert((Dimensions<std::size_t,Base::DATA_RANK>(in.dimensions()).template demote<>()) == Base::dims);
            assert(in.dimension(0) > 0);
            out_cache = in.tanh();
            return out_cache;
        }
        inline typename Root::Data pass_back(typename Root::Data out_grad) {
            assert((Dimensions<std::size_t,Base::DATA_RANK>(out_grad.dimensions()).template demote<>()) == Base::dims);
            assert(out_grad.dimension(0) > 0 && out_cache.dimension(0) == out_grad.dimension(0));
            if (Base::is_input_layer())
                return typename Root::Data();
            return (out_cache.constant(1) - out_cache * out_cache) * out_grad;
        }
    private:
        typename Root::Data out_cache;
    };

    /**
    * An alias for a shared pointer to a Parameters instance.
    */
    template<typename Scalar>
    using ParamsSharedPtr = std::shared_ptr<Parameters<Scalar>>;

    /**
    * An abstract class template that represents a kernel-based linear transformation function layer.
    */
    template<typename Scalar, std::size_t Rank>
    class KernelLayer : public Layer<Scalar,Rank> {
        typedef Layer<Scalar,Rank> Base;
        typedef KernelLayer<Scalar,Rank> Self;
    public:
        virtual ~KernelLayer() = default;
        inline const Base& get_params_owner() const {
            return owner;
        }
        inline const typename Base::Dims& get_input_dims() const {
            return input_dims;
        }
        inline const typename Base::Dims& get_output_dims() const {
            return output_dims;
        }
        inline bool is_input_layer() const {
            return input_layer;
        }
        inline void set_input_layer(bool input_layer) {
            this->input_layer = input_layer;
        }
        inline std::vector<const Parameters<Scalar>*> get_params() const {
            return std::vector<const Parameters<Scalar>*>({ weights.get(), bias.get() });
        }
        inline std::vector<Parameters<Scalar>*> get_params() {
            return std::vector<Parameters<Scalar>*>({ weights.get(), bias.get() });
        }
    protected:
        /**
        * @param input_dims The input dimensions of the layer.
        * @param output_dims The output dimensions of the layer.
        * @param weights The weight parameters; cannot be null.
        * @param bias The bias parameters; cannot be null.
        */
        inline KernelLayer(const typename Base::Dims& input_dims, const typename Base::Dims& output_dims,
                ParamsSharedPtr<Scalar> weights, ParamsSharedPtr<Scalar> bias) :
                    owner(*this),
                    input_dims(input_dims),
                    output_dims(output_dims),
                    weights(weights),
                    bias(bias),
                    input_layer(false) {
            assert(weights && bias);
        }
        inline KernelLayer(const Self& layer, bool share_params = false) :
                owner(share_params ? layer.owner : *this),
                input_dims(layer.input_dims),
                output_dims(layer.output_dims),
                weights(share_params ? layer.weights : ParamsSharedPtr<Scalar>(layer.weights->clone())),
                bias(share_params ? layer.bias : ParamsSharedPtr<Scalar>(layer.bias->clone())),
                input_layer(layer.input_layer) { }
        const Self& owner;
        const typename Base::Dims input_dims, output_dims;
        ParamsSharedPtr<Scalar> weights, bias;
        bool input_layer;
    };

    
    /**
    * An abstract base class template for a 2D convolutional layer.
    */
    template<typename Scalar, std::size_t Rank>
    class ConvKernelLayerBase : public KernelLayer<Scalar,Rank> {
        typedef Layer<Scalar,Rank> Root;
        typedef KernelLayer<Scalar,Rank> Base;
        typedef std::array<std::size_t,4> Array4;
        typedef std::array<std::pair<std::size_t,std::size_t>,4> PaddingsArray4;
    public:
        inline void empty_cache() {
            in_conv_mat_cache = Matrix<Scalar>();
        }
    protected:
        inline ConvKernelLayerBase(const typename Root::Dims& input_dims, std::size_t filters,
                std::size_t receptor_height, std::size_t receptor_width, std::size_t vertical_padding,
                std::size_t horizontal_padding, std::size_t vertical_stride, std::size_t horizontal_stride,
                std::size_t vertical_dilation, std::size_t horizontal_dilation, ParamInitSharedPtr<Scalar> weight_init,
                ParamRegSharedPtr<Scalar> weight_reg, Scalar weight_clip, Scalar weight_max_l1_norm,
                Scalar weight_max_l2_norm, Scalar weight_grad_clip, Scalar weight_grad_max_l1_norm,
                Scalar weight_grad_max_l2_norm, ParamRegSharedPtr<Scalar> bias_reg, Scalar bias_clip,
                Scalar bias_max_l1_norm, Scalar bias_max_l2_norm, Scalar bias_grad_clip,
                Scalar bias_grad_max_l1_norm, Scalar bias_grad_max_l2_norm) :
                    /* For every filter, there is a column in the weight matrix with the same number of
                    * elements as the area of the receptive field (F * F * D). */
                    Base(input_dims, calculate_adjusted_output_dims(input_dims, filters, receptor_height,
                            receptor_width, vertical_padding, horizontal_padding, vertical_stride, horizontal_stride,
                            vertical_dilation, horizontal_dilation),
                            std::make_shared<StandardParameters<Scalar>>(receptor_height * receptor_width *
                                    input_dims.template extend<3 - Rank>()(2), filters, true, weight_init,
                                    weight_reg, weight_clip, weight_max_l1_norm, weight_max_l2_norm,
                                    weight_grad_clip, weight_grad_max_l1_norm, weight_grad_max_l2_norm),
                            std::make_shared<StandardParameters<Scalar>>(1, filters, true,
                                    std::make_shared<ZeroParameterInitialization<Scalar>>(), bias_reg,
                                    bias_clip, bias_max_l1_norm, bias_max_l2_norm, bias_grad_clip,
                                    bias_grad_max_l1_norm, bias_grad_max_l2_norm)),
                    filters(filters),
                    receptor_height(receptor_height),
                    receptor_width(receptor_width),
                    vertical_padding(vertical_padding),
                    horizontal_padding(horizontal_padding),
                    vertical_stride(vertical_stride),
                    horizontal_stride(horizontal_stride),
                    vertical_dilation(vertical_dilation),
                    horizontal_dilation(horizontal_dilation),
                    ext_input_dims(input_dims.template extend<3 - Rank>()),
                    ext_output_dims(calculate_output_dims(ext_input_dims, filters, receptor_height, receptor_width,
                            vertical_padding, horizontal_padding, vertical_stride, horizontal_stride,
                            vertical_dilation, horizontal_dilation)),
                    padded_height(ext_input_dims(0) + 2 * vertical_padding),
                    padded_width(ext_input_dims(1) + 2 * horizontal_padding),
                    dil_receptor_height(receptor_height + (receptor_height - 1) * vertical_dilation),
                    dil_receptor_width(receptor_width + (receptor_width - 1) * horizontal_dilation),
                    patches_per_sample(ext_output_dims(0) * ext_output_dims(1)),
                    out_conversion_dims({ 0u, ext_output_dims(0), ext_output_dims(1), ext_output_dims(2) }),
                    patch_offsets({ 0u, 0u, 0u, 0u }),
                    patch_extents({ 0u, dil_receptor_height, dil_receptor_width, ext_input_dims(2) }),
                    dil_strides({ 1u, vertical_dilation + 1u, horizontal_dilation + 1u, 1u }),
                    no_padding_offsets({ 0u, vertical_padding, horizontal_padding, 0u }),
                    no_padding_extents({ 0u, ext_input_dims(0), ext_input_dims(1), ext_input_dims(2) }),
                    paddings({ std::make_pair(0, 0), std::make_pair(vertical_padding, vertical_padding),
                            std::make_pair(horizontal_padding, horizontal_padding), std::make_pair(0, 0) }) {
            assert(filters > 0);
            assert(receptor_height > 0);
            assert(receptor_width > 0);
            assert(vertical_stride > 0 && horizontal_stride > 0);
            assert(ext_input_dims(0) + 2 * vertical_padding >= dil_receptor_height &&
                    ext_input_dims(1) + 2 * horizontal_padding >= dil_receptor_width);
        }
        inline ConvKernelLayerBase(const ConvKernelLayerBase<Scalar,Rank>& layer, bool share_params = false) :
                Base(layer, share_params),
                filters(layer.filters),
                receptor_height(layer.receptor_height),
                receptor_width(layer.receptor_width),
                vertical_padding(layer.vertical_padding),
                horizontal_padding(layer.horizontal_padding),
                vertical_stride(layer.vertical_stride),
                horizontal_stride(layer.horizontal_stride),
                vertical_dilation(layer.vertical_dilation),
                horizontal_dilation(layer.horizontal_dilation),
                ext_input_dims(layer.ext_input_dims),
                ext_output_dims(layer.ext_output_dims),
                padded_height(layer.padded_height),
                padded_width(layer.padded_width),
                dil_receptor_height(layer.dil_receptor_height),
                dil_receptor_width(layer.dil_receptor_width),
                patches_per_sample(layer.patches_per_sample),
                out_conversion_dims(layer.out_conversion_dims),
                patch_offsets(layer.patch_offsets),
                patch_extents(layer.patch_extents),
                dil_strides(layer.dil_strides),
                no_padding_offsets(layer.no_padding_offsets),
                no_padding_extents(layer.no_padding_extents),
                paddings(layer.paddings),
                in_conv_mat_cache(layer.in_conv_mat_cache) { }
        inline Tensor<Scalar,4> _pass_forward(Tensor<Scalar,4> in, bool training) {
            // Spatial padding.
            if (vertical_padding > 0 || horizontal_padding > 0)
                in = Tensor<Scalar,4>(in.pad(paddings));
            std::size_t rows = in.dimension(0);
            std::size_t total_patches = rows * patches_per_sample;
            std::size_t receptor_vol = Base::weights->get_values().rows();
            /* Flatten the receptor cuboids into row vectors and concatenate them. Each row stands for one
            * stretched out receptor of one sample. The same receptor location along all samples of the
            * batch is represented by a contiguous block of these rows. */
            std::size_t patch_ind = 0;
            patch_extents[0] = rows;
            in_conv_mat_cache = Matrix<Scalar>(total_patches, receptor_vol);
            for (std::size_t i = 0; i <= padded_width - dil_receptor_width; i += horizontal_stride) {
                patch_offsets[2] = i;
                for (std::size_t j = 0; j <= padded_height - dil_receptor_height; j += vertical_stride) {
                    patch_offsets[1] = j;
                    Tensor<Scalar,4> patch;
                    // If the patch is dilated, skip the spatial gaps when flattening it into a matrix.
                    if (vertical_dilation > 0 || horizontal_dilation > 0)
                        patch = in.slice(patch_offsets, patch_extents).stride(dil_strides);
                    else
                        patch = in.slice(patch_offsets, patch_extents);
                    in_conv_mat_cache.block(patch_ind, 0, rows, receptor_vol) = MatrixMap<Scalar>(patch.data(),
                            rows, receptor_vol);
                    patch_ind += rows;
                }
            }
            assert(patch_ind == total_patches);
            Matrix<Scalar> out_mat = (in_conv_mat_cache * Base::weights->get_values()).rowwise() +
                    Base::bias->get_values().row(0);
            out_conversion_dims[0] = rows;
            return TensorMap<Scalar,4>(out_mat.data(), out_conversion_dims);
        }
        inline Tensor<Scalar,4> _pass_back(Tensor<Scalar,4> out_grad) {
            std::size_t rows = out_grad.dimension(0);
            std::size_t total_patches = rows * patches_per_sample;
            std::size_t receptor_vol = Base::weights->get_values().rows();
            MatrixMap<Scalar> out_grad_mat(out_grad.data(), total_patches, filters);
            Base::weights->accumulate_grad(in_conv_mat_cache.transpose() * out_grad_mat);
            Base::bias->accumulate_grad(out_grad_mat.colwise().sum());
            if (Base::is_input_layer())
                return Tensor<Scalar,4>();
            // Compute the gradient of the previous layer's output.
            Matrix<Scalar> prev_out_grad_conv_mat = out_grad_mat * Base::weights->get_values().transpose();
            /* Given the gradient of the stretched out receptor patches, perform a 'backwards' convolution
            * to get the derivative w.r.t. the individual input nodes. */
            Tensor<Scalar,4> prev_out_grad(rows, padded_height, padded_width, ext_input_dims(2));
            prev_out_grad.setZero();
            std::size_t patch_ind = 0;
            patch_extents[0] = rows;
            for (std::size_t i = 0; i <= padded_width - dil_receptor_width; i += horizontal_stride) {
                patch_offsets[2] = i;
                for (std::size_t j = 0; j <= padded_height - dil_receptor_height; j += vertical_stride) {
                    patch_offsets[1] = j;
                    // Accumulate the gradients where the receptor-patch-tensors overlap.
                    Matrix<Scalar> prev_out_grad_conv_mat_block = prev_out_grad_conv_mat.block(patch_ind, 0,
                            rows, receptor_vol);
                    TensorMap<Scalar,4> prev_out_grad_patch(prev_out_grad_conv_mat_block.data(), rows,
                            receptor_height, receptor_width, ext_input_dims(2));
                    if (vertical_dilation > 0 || horizontal_dilation > 0)
                        prev_out_grad.slice(patch_offsets, patch_extents).stride(dil_strides) += prev_out_grad_patch;
                    else
                        prev_out_grad.slice(patch_offsets, patch_extents) += prev_out_grad_patch;
                    patch_ind += rows;
                }
            }
            assert(patch_ind == prev_out_grad_conv_mat.rows());
            if (vertical_padding > 0 || horizontal_padding > 0) {
                // Cut off the padding.
                no_padding_extents[0] = rows;
                return prev_out_grad.slice(no_padding_offsets, no_padding_extents);
            } else
                return prev_out_grad;
        }
        // The defining attributes of the convolutional layer.
        const std::size_t filters, receptor_height, receptor_width, vertical_padding, horizontal_padding,
                vertical_stride, horizontal_stride, vertical_dilation, horizontal_dilation;
    private:
        inline static std::size_t calculate_spatial_output_dim(std::size_t input_dim, std::size_t receptor_size,
                std::size_t padding, std::size_t dilation, std::size_t stride) {
            return (input_dim - receptor_size - (receptor_size - 1) * dilation + 2 * padding) / stride + 1;
        }
        inline static Dimensions<std::size_t,3> calculate_output_dims(const Dimensions<std::size_t,3>& input_dims,
                std::size_t filters, std::size_t receptor_height, std::size_t receptor_width, std::size_t vertical_padding,
                std::size_t horizontal_padding, std::size_t vertical_stride, std::size_t horizontal_stride,
                std::size_t vertical_dilation, std::size_t horizontal_dilation) {
            return { calculate_spatial_output_dim(input_dims(0), receptor_height, vertical_padding, vertical_dilation,
                    vertical_stride), calculate_spatial_output_dim(input_dims(1), receptor_width, horizontal_padding,
                            horizontal_dilation, horizontal_stride), filters };
        }
        inline static Dimensions<std::size_t,Rank> calculate_adjusted_output_dims(
                const Dimensions<std::size_t,Rank>& input_dims, std::size_t filters, std::size_t receptor_height,
                std::size_t receptor_width, std::size_t vertical_padding, std::size_t horizontal_padding,
                std::size_t vertical_stride, std::size_t horizontal_stride, std::size_t vertical_dilation,
                std::size_t horizontal_dilation) {
            auto output_dims = calculate_output_dims(input_dims.template extend<3 - Rank>(), filters, receptor_height,
                    receptor_width, vertical_padding, horizontal_padding, vertical_stride, horizontal_stride,
                    vertical_dilation, horizontal_dilation);
            output_dims(2) /= filters;
            output_dims(Rank - 1) *= filters;
            return output_dims.template contract<3 - Rank>();
        }
        const Dimensions<std::size_t,3> ext_input_dims, ext_output_dims;
        // Pre-computed values to improve propagation-time performance.
        const std::size_t padded_height, padded_width, dil_receptor_height, dil_receptor_width, patches_per_sample;
        Array4 out_conversion_dims, patch_offsets, patch_extents, dil_strides, no_padding_offsets, no_padding_extents;
        PaddingsArray4 paddings;
        // Staged computation caches
        Matrix<Scalar> in_conv_mat_cache;
    };

    /**
    * A class template for a 2D convolutional layer operating on rank-3 data batches (rank-4 tensors).  The results
    * of the convolutions of the filters and the input tensor are concatenated along the highest (4th) rank of the
    * output tensor.
    */
    template<typename Scalar, std::size_t Rank = 3>
    class ConvKernelLayer : public ConvKernelLayerBase<Scalar,Rank> {
        typedef Layer<Scalar,3> Root;
        typedef KernelLayer<Scalar,3> KernelBase;
        typedef ConvKernelLayerBase<Scalar,3> ConvBase;
    public:
        /**
        * @param input_dims The dimensionality of the observations to be processed by the layer.
        * @param filters The number of filters to use.
        * @param weight_init A shared pointer to a weight initialization used to initialize the weights of
        * the layer.
        * @param receptor_height The height of the base of the receptor cuboid.
        * @param receptor_width The width of the base of the receptor cuboid.
        * @param vertical_padding The extent of padding to apply to the input tensor along its height (both
        * at the top and at the bottom).
        * @param horizontal_padding The extent of padding to apply to the input tensor along its width (both
        * at the left and at the right).
        * @param vertical_stride The vertical convolution stride i.e. the number of elements by which the
        * receptor is to be shifted along the height of the input tensor.
        * @param horizontal_stride The horizonzal convolution stride i.e. the number of elements by which the
        * receptor is to be shifted along the width of the input tensor.
        * @param vertical_dilation The extent of vertical dilation to apply to the receptor.
        * @param horizontal_dilation The extent of horizontal dilation to apply to the receptor.
        * @param weight_reg An optional regularization function to apply to the weights.
        * @param weight_clip The maximum allowed absolute weight value. If it is 0 or less, no value
        * clipping is performed.
        * @param weight_max_l1_norm The maximum allowed L1 weight value norm. If it is 0 or less, no L1
        * max norm constraint is enforced.
        * @param weight_max_l2_norm The maximum allowed L2 weight value norm. If it is 0 or less, no L2
        * max norm constraint is enforced.
        * @param weight_grad_clip The maximum allowed absolute weight gradient. If it is 0 or less, no
        * gradient clipping is performed.
        * @param weight_grad_max_l1_norm The maximum allowed L1 weight gradient norm. If it is 0 or less,
        * no L1 gradient max norm constraint is enforced.
        * @param weight_grad_max_l2_norm The maximum allowed L2 weight gradient norm. If it is 0 or less,
        * no L2 gradient max norm constraint is enforced.
        * @param bias_reg An optional regularization function to apply to the bias.
        * @param bias_clip The maximum allowed absolute bias value. If it is 0 or less, no value clipping
        * is performed.
        * @param bias_max_l1_norm The maximum allowed L1 bias value norm. If it is 0 or less, no bias L1
        * max norm constraint is enforced.
        * @param bias_max_l2_norm The maximum allowed L2 bias value norm. If it is 0 or less, no bias L2
        * max norm constraint is enforced.
        * @param bias_grad_clip The maximum allowed absolute bias gradient. If it is 0 or less, no
        * gradient clipping is performed.
        * @param bias_grad_max_l1_norm The maximum allowed L1 bias gradient norm. If it is 0 or less, no
        * bias L1 gradient max norm constraint is enforced.
        * @param bias_grad_max_l2_norm The maximum allowed L2 bias gradient norm. If it is 0 or less, no
        * bias L2 gradient max norm constraint is enforced.
        */
        inline ConvKernelLayer(const typename Root::Dims& input_dims, std::size_t filters,
                ParamInitSharedPtr<Scalar> weight_init, std::size_t receptor_height = 3, std::size_t receptor_width = 3,
                std::size_t vertical_padding = 1, std::size_t horizontal_padding = 1, std::size_t vertical_stride = 1,
                std::size_t horizontal_stride = 1, std::size_t vertical_dilation = 0, std::size_t horizontal_dilation = 0,
                ParamRegSharedPtr<Scalar> weight_reg = nullptr, Scalar weight_clip = 0, Scalar weight_max_l1_norm = 0,
                Scalar weight_max_l2_norm = 0, Scalar weight_grad_clip = 0, Scalar weight_grad_max_l1_norm = 0,
                Scalar weight_grad_max_l2_norm = 0, ParamRegSharedPtr<Scalar> bias_reg = nullptr, Scalar bias_clip = 0,
                Scalar bias_max_l1_norm = 0, Scalar bias_max_l2_norm = 0, Scalar bias_grad_clip = 0,
                Scalar bias_grad_max_l1_norm = 0, Scalar bias_grad_max_l2_norm = 0) :
                    ConvBase::ConvKernelLayerBase(input_dims, filters, receptor_height, receptor_width, vertical_padding,
                            horizontal_padding, vertical_stride, horizontal_stride, vertical_dilation,
                            horizontal_dilation, weight_init, weight_reg, weight_clip, weight_max_l1_norm,
                            weight_max_l2_norm, weight_grad_clip, weight_grad_max_l1_norm, weight_grad_max_l2_norm,
                            bias_reg, bias_clip, bias_max_l1_norm, bias_max_l2_norm, bias_grad_clip,
                            bias_grad_max_l1_norm, bias_grad_max_l2_norm) { }
        inline ConvKernelLayer(const ConvKernelLayer<Scalar,Rank>& layer, bool share_params = false) :
                ConvBase::ConvKernelLayerBase(layer, share_params),
                batch_size(layer.batch_size) { }
        inline Root* clone() const {
            return new ConvKernelLayer(*this);
        }
        inline Root* clone_with_shared_params() {
            return new ConvKernelLayer(*this, true);
        }
        inline typename Root::Data pass_forward(typename Root::Data in, bool training) {
            assert((Dimensions<std::size_t,4>(in.dimensions()).template demote<>()) == KernelBase::input_dims);
            assert(in.dimension(0) > 0);
            batch_size = in.dimension(0);
            return ConvBase::_pass_forward(std::move(in), training);
        }
        inline typename Root::Data pass_back(typename Root::Data out_grad) {
            assert((Dimensions<std::size_t,4>(out_grad.dimensions()).template demote<>()) == KernelBase::output_dims);
            assert(out_grad.dimension(0) > 0 && batch_size == out_grad.dimension(0));
            return ConvBase::_pass_back(std::move(out_grad));
        }
    private:
        std::size_t batch_size;
    };

    /**
    * A class template for a 2D convolutional layer operating on rank-2 data batches (rank-3 tensors).  The results
    * of the convolutions of the filters and the input tensor are concatenated along the highest (3rd) rank of the
    * output tensor.
    */
    template<typename Scalar>
    class ConvKernelLayer<Scalar,2> : public ConvKernelLayerBase<Scalar,2> {
        typedef Layer<Scalar,2> Root;
        typedef KernelLayer<Scalar,2> KernelBase;
        typedef ConvKernelLayerBase<Scalar,2> ConvBase;
    public:
        /**
        * @param input_dims The dimensionality of the observations to be processed by the layer.
        * @param filters The number of filters to use.
        * @param weight_init A shared pointer to a weight initialization used to initialize the weights of
        * the layer.
        * @param receptor_height The height of the base of the receptor cuboid.
        * @param receptor_width The width of the base of the receptor cuboid.
        * @param vertical_padding The extent of padding to apply to the input tensor along its height (both
        * at the top and at the bottom).
        * @param horizontal_padding The extent of padding to apply to the input tensor along its width (both
        * at the left and at the right).
        * @param vertical_stride The vertical convolution stride i.e. the number of elements by which the
        * receptor is to be shifted along the height of the input tensor.
        * @param horizontal_stride The horizonzal convolution stride i.e. the number of elements by which the
        * receptor is to be shifted along the width of the input tensor.
        * @param vertical_dilation The extent of vertical dilation to apply to the receptor.
        * @param horizontal_dilation The extent of horizontal dilation to apply to the receptor.
        * @param weight_reg An optional regularization function to apply to the weights.
        * @param weight_clip The maximum allowed absolute weight value. If it is 0 or less, no value
        * clipping is performed.
        * @param weight_max_l1_norm The maximum allowed L1 weight value norm. If it is 0 or less, no L1
        * max norm constraint is enforced.
        * @param weight_max_l2_norm The maximum allowed L2 weight value norm. If it is 0 or less, no L2
        * max norm constraint is enforced.
        * @param weight_grad_clip The maximum allowed absolute weight gradient. If it is 0 or less, no
        * gradient clipping is performed.
        * @param weight_grad_max_l1_norm The maximum allowed L1 weight gradient norm. If it is 0 or less,
        * no L1 gradient max norm constraint is enforced.
        * @param weight_grad_max_l2_norm The maximum allowed L2 weight gradient norm. If it is 0 or less,
        * no L2 gradient max norm constraint is enforced.
        * @param bias_reg An optional regularization function to apply to the bias.
        * @param bias_clip The maximum allowed absolute bias value. If it is 0 or less, no value clipping
        * is performed.
        * @param bias_max_l1_norm The maximum allowed L1 bias value norm. If it is 0 or less, no bias L1
        * max norm constraint is enforced.
        * @param bias_max_l2_norm The maximum allowed L2 bias value norm. If it is 0 or less, no bias L2
        * max norm constraint is enforced.
        * @param bias_grad_clip The maximum allowed absolute bias gradient. If it is 0 or less, no
        * gradient clipping is performed.
        * @param bias_grad_max_l1_norm The maximum allowed L1 bias gradient norm. If it is 0 or less, no
        * bias L1 gradient max norm constraint is enforced.
        * @param bias_grad_max_l2_norm The maximum allowed L2 bias gradient norm. If it is 0 or less, no
        * bias L2 gradient max norm constraint is enforced.
        */
        inline ConvKernelLayer(const typename Root::Dims& input_dims, std::size_t filters,
                ParamInitSharedPtr<Scalar> weight_init, std::size_t receptor_height = 3, std::size_t receptor_width = 3,
                std::size_t vertical_padding = 1, std::size_t horizontal_padding = 1, std::size_t vertical_stride = 1,
                std::size_t horizontal_stride = 1, std::size_t vertical_dilation = 0, std::size_t horizontal_dilation = 0,
                ParamRegSharedPtr<Scalar> weight_reg = nullptr, Scalar weight_clip = 0, Scalar weight_max_l1_norm = 0,
                Scalar weight_max_l2_norm = 0, Scalar weight_grad_clip = 0, Scalar weight_grad_max_l1_norm = 0,
                Scalar weight_grad_max_l2_norm = 0, ParamRegSharedPtr<Scalar> bias_reg = nullptr, Scalar bias_clip = 0,
                Scalar bias_max_l1_norm = 0, Scalar bias_max_l2_norm = 0, Scalar bias_grad_clip = 0,
                Scalar bias_grad_max_l1_norm = 0, Scalar bias_grad_max_l2_norm = 0) :
                    ConvBase::ConvKernelLayerBase(input_dims, filters, receptor_height, receptor_width, vertical_padding,
                            horizontal_padding, vertical_stride, horizontal_stride, vertical_dilation,
                            horizontal_dilation, weight_init, weight_reg, weight_clip, weight_max_l1_norm,
                            weight_max_l2_norm, weight_grad_clip, weight_grad_max_l1_norm, weight_grad_max_l2_norm,
                            bias_reg, bias_clip, bias_max_l1_norm, bias_max_l2_norm, bias_grad_clip,
                            bias_grad_max_l1_norm, bias_grad_max_l2_norm) { }
        inline ConvKernelLayer(const ConvKernelLayer<Scalar,2>& layer, bool share_params) :
                ConvBase::ConvKernelLayerBase(layer, share_params),
                batch_size(layer.batch_size) { }
        inline Root* clone() const {
            return new ConvKernelLayer(*this);
        }
        inline Root* clone_with_shared_params() {
            return new ConvKernelLayer(*this, true);
        }
        inline typename Root::Data pass_forward(typename Root::Data in, bool training) {
            assert((Dimensions<std::size_t,3>(in.dimensions()).template demote<>()) == KernelBase::input_dims);
            assert(in.dimension(0) > 0);
            batch_size = in.dimension(0);
            return ConvBase::_pass_forward(TensorMap<Scalar,4>(in.data(), { batch_size, in.dimension(1),
                    in.dimension(2), 1u }), training).reshape(std::array<std::size_t,3>({ batch_size,
                            KernelBase::output_dims(0), KernelBase::output_dims(1) }));
        }
        inline typename Root::Data pass_back(typename Root::Data out_grad) {
            assert((Dimensions<std::size_t,3>(out_grad.dimensions()).template demote<>()) == KernelBase::output_dims);
            assert(out_grad.dimension(0) > 0 && batch_size == out_grad.dimension(0));
            Tensor<Scalar,4> prev_out_grad = ConvBase::_pass_back(TensorMap<Scalar,4>(out_grad.data(),
                    { batch_size, KernelBase::output_dims(0), KernelBase::output_dims(1) / ConvBase::filters,
                            ConvBase::filters }));
            if (KernelBase::is_input_layer())
                return typename Root::Data();
            return TensorMap<Scalar,3>(prev_out_grad.data(), { batch_size, KernelBase::input_dims(0),
                    KernelBase::input_dims(1) });
        }
    private:
        std::size_t batch_size;
    };


    /**
    * A class template for a 1D convolutional layer operating on rank-1 data batches (rank-2 tensors).  The results
    * of the convolutions of the filters and the input tensor are concatenated along the highest (2nd) rank of the
    * output tensor.
    */
    template<typename Scalar>
    class ConvKernelLayer<Scalar,1> : public ConvKernelLayerBase<Scalar,1> {
        typedef Layer<Scalar,1> Root;
        typedef KernelLayer<Scalar,1> KernelBase;
        typedef ConvKernelLayerBase<Scalar,1> ConvBase;
    public:
        /**
        * @param input_dims The dimensionality of the observations to be processed by the layer.
        * @param filters The number of filters to use.
        * @param weight_init A shared pointer to a weight initialization used to initialize the weights of
        * the layer.
        * @param receptor_length The length of the receptor.
        * @param padding The extent of padding to apply to the input tensor along its length on both ends.
        * @param stride The convolution stride i.e. the number of elements by which the receptor is to be
        * shifted along the length of the input tensor.
        * @param dilation The extent of dilation to apply to the receptor.
        * @param weight_reg An optional regularization function to apply to the weights.
        * @param weight_clip The maximum allowed absolute weight value. If it is 0 or less, no value
        * clipping is performed.
        * @param weight_max_l1_norm The maximum allowed L1 weight value norm. If it is 0 or less, no L1
        * max norm constraint is enforced.
        * @param weight_max_l2_norm The maximum allowed L2 weight value norm. If it is 0 or less, no L2
        * max norm constraint is enforced.
        * @param weight_grad_clip The maximum allowed absolute weight gradient. If it is 0 or less, no
        * gradient clipping is performed.
        * @param weight_grad_max_l1_norm The maximum allowed L1 weight gradient norm. If it is 0 or less,
        * no L1 gradient max norm constraint is enforced.
        * @param weight_grad_max_l2_norm The maximum allowed L2 weight gradient norm. If it is 0 or less,
        * no L2 gradient max norm constraint is enforced.
        * @param bias_reg An optional regularization function to apply to the bias.
        * @param bias_clip The maximum allowed absolute bias value. If it is 0 or less, no value clipping
        * is performed.
        * @param bias_max_l1_norm The maximum allowed L1 bias value norm. If it is 0 or less, no bias L1
        * max norm constraint is enforced.
        * @param bias_max_l2_norm The maximum allowed L2 bias value norm. If it is 0 or less, no bias L2
        * max norm constraint is enforced.
        * @param bias_grad_clip The maximum allowed absolute bias gradient. If it is 0 or less, no
        * gradient clipping is performed.
        * @param bias_grad_max_l1_norm The maximum allowed L1 bias gradient norm. If it is 0 or less, no
        * bias L1 gradient max norm constraint is enforced.
        * @param bias_grad_max_l2_norm The maximum allowed L2 bias gradient norm. If it is 0 or less, no
        * bias L2 gradient max norm constraint is enforced.
        */
        ConvKernelLayer(const typename Root::Dims& input_dims, std::size_t filters,
                ParamInitSharedPtr<Scalar> weight_init, std::size_t receptor_length = 3, std::size_t padding = 1,
                std::size_t stride = 1, std::size_t dilation = 0, ParamRegSharedPtr<Scalar> weight_reg = nullptr,
                Scalar weight_clip = 0, Scalar weight_max_l1_norm = 0, Scalar weight_max_l2_norm = 0,
                Scalar weight_grad_clip = 0, Scalar weight_grad_max_l1_norm = 0, Scalar weight_grad_max_l2_norm = 0,
                ParamRegSharedPtr<Scalar> bias_reg = nullptr, Scalar bias_clip = 0, Scalar bias_max_l1_norm = 0,
                Scalar bias_max_l2_norm = 0, Scalar bias_grad_clip = 0, Scalar bias_grad_max_l1_norm = 0,
                Scalar bias_grad_max_l2_norm = 0) :
                    ConvBase::ConvKernelLayerBase(input_dims, filters, receptor_length, 1, padding, 0, stride, 1,
                            dilation, 0, weight_init, weight_reg, weight_clip, weight_max_l1_norm, weight_max_l2_norm,
                            weight_grad_clip, weight_grad_max_l1_norm, weight_grad_max_l2_norm, bias_reg, bias_clip,
                            bias_max_l1_norm, bias_max_l2_norm, bias_grad_clip, bias_grad_max_l1_norm,
                            bias_grad_max_l2_norm) { }
        inline ConvKernelLayer(ConvKernelLayer<Scalar,1>& layer, bool share_params) :
                ConvBase::ConvKernelLayerBase(layer, share_params),
                batch_size(layer.batch_size) { }
        inline Root* clone() const {
            return new ConvKernelLayer(*this);
        }
        inline Root* clone_with_shared_params() {
            return new ConvKernelLayer(*this, true);
        }
        inline typename Root::Data pass_forward(typename Root::Data in, bool training) {
            assert((Dimensions<std::size_t,2>(in.dimensions()).template demote<>()) == KernelBase::input_dims);
            assert(in.dimension(0) > 0);
            batch_size = in.dimension(0);
            return ConvBase::_pass_forward(TensorMap<Scalar,4>(in.data(), { batch_size, in.dimension(1), 1u, 1u }),
                    training).reshape(std::array<std::size_t,2>({ batch_size, KernelBase::output_dims(0) }));
        }
        inline typename Root::Data pass_back(typename Root::Data out_grad) {
            assert((Dimensions<std::size_t,2>(out_grad.dimensions()).template demote<>()) == KernelBase::output_dims);
            assert(out_grad.dimension(0) > 0 && batch_size == out_grad.dimension(0));
            Tensor<Scalar,4> prev_out_grad = ConvBase::_pass_back(TensorMap<Scalar,4>(out_grad.data(),
                    { batch_size, KernelBase::output_dims(0) / ConvBase::filters, 1, ConvBase::filters }));
            if (KernelBase::is_input_layer())
                return typename Root::Data();
            return TensorMap<Scalar,2>(prev_out_grad.data(), { batch_size, KernelBase::input_dims(0) });
        }
    private:
        std::size_t batch_size;
    };

    /**
    * A class template representing a fully connected layer.
    */
    template<typename Scalar, std::size_t Rank = 1>
    class DenseKernelLayer : public KernelLayer<Scalar,Rank> {
        typedef Layer<Scalar,Rank> Root;
        typedef KernelLayer<Scalar,Rank> Base;
        typedef std::array<std::size_t,Root::DATA_RANK> RankwiseArray;
    public:
        /**
        * @param input_dims The dimensionality of the observations to be processed by the layer.
        * @param output_size The length of the vector output for each sample.
        * @param weight_init A shared pointer to a weight initialization used to initialize the
        * values of the weights.
        * @param weight_reg An optional regularization function to apply to the weights.
        * @param weight_clip The maximum allowed absolute weight value. If it is 0 or less, no
        * value clipping is performed.
        * @param weight_max_l1_norm The maximum allowed L1 weight value norm. If it is 0 or
        * less, no L1 max norm constraint is enforced.
        * @param weight_max_l2_norm The maximum allowed L2 weight value norm. If it is 0 or
        * less, no L2 max norm constraint is enforced.
        * @param weight_grad_clip The maximum allowed absolute weight gradient. If it is 0
        * or less, no gradient clipping is performed.
        * @param weight_grad_max_l1_norm The maximum allowed L1 weight gradient norm. If it
        * is 0 or less, no L1 gradient max norm constraint is enforced.
        * @param weight_grad_max_l2_norm The maximum allowed L2 weight gradient norm. If it
        * is 0 or less, no L2 gradient max norm constraint is enforced.
        * @param bias_reg An optional regularization function to apply to the bias.
        * @param bias_clip The maximum allowed absolute bias value. If it is 0 or less, no
        * value clipping is performed.
        * @param bias_max_l1_norm The maximum allowed L1 bias value norm. If it is 0 or less,
        * no bias L1 max norm constraint is enforced.
        * @param bias_max_l2_norm The maximum allowed L2 bias value norm. If it is 0 or less,
        * no bias L2 max norm constraint is enforced.
        * @param bias_grad_clip The maximum allowed absolute bias gradient. If it is 0 or
        * less, no gradient clipping is performed.
        * @param bias_grad_max_l1_norm The maximum allowed L1 bias gradient norm. If it is 0
        * or less, no bias L1 gradient max norm constraint is enforced.
        * @param bias_grad_max_l2_norm The maximum allowed L2 bias gradient norm. If it is 0
        * or less, no bias L2 gradient max norm constraint is enforced.
        */
        inline DenseKernelLayer(const typename Root::Dims& input_dims, std::size_t output_size,
                ParamInitSharedPtr<Scalar> weight_init, ParamRegSharedPtr<Scalar> weight_reg = nullptr,
                Scalar weight_clip = 0, Scalar weight_max_l1_norm = 0, Scalar weight_max_l2_norm = 0,
                Scalar weight_grad_clip = 0, Scalar weight_grad_max_l1_norm = 0, Scalar weight_grad_max_l2_norm = 0,
                ParamRegSharedPtr<Scalar> bias_reg = nullptr, Scalar bias_clip = 0, Scalar bias_max_l1_norm = 0,
                Scalar bias_max_l2_norm = 0, Scalar bias_grad_clip = 0, Scalar bias_grad_max_l1_norm = 0,
                Scalar bias_grad_max_l2_norm = 0) :
                    Base(input_dims, { output_size },
                            std::make_shared<StandardParameters<Scalar>>(input_dims.get_volume(), output_size, true,
                                    weight_init, weight_reg, weight_clip, weight_max_l1_norm, weight_max_l2_norm,
                                    weight_grad_clip, weight_grad_max_l1_norm, weight_grad_max_l2_norm),
                            std::make_shared<StandardParameters<Scalar>>(1, output_size, true,
                                    std::make_shared<ZeroParameterInitialization<Scalar>>(), bias_reg, bias_clip,
                                    bias_max_l1_norm, bias_max_l2_norm, bias_grad_clip, bias_grad_max_l1_norm,
                                    bias_grad_max_l2_norm)),
                    out_conversion_dims(Base::output_dims.template promote<>()),
                    prev_out_conversion_dims(Base::input_dims.template promote<>()) { }
        inline DenseKernelLayer(const DenseKernelLayer<Scalar,Rank>& layer, bool share_params = false) :
                    Base(layer, share_params),
                    out_conversion_dims(layer.out_conversion_dims),
                    prev_out_conversion_dims(layer.prev_out_conversion_dims),
                    in_mat_cache(layer.in_mat_cache) { }
        inline Root* clone() const {
            return new DenseKernelLayer(*this);
        }
        inline Root* clone_with_shared_params() {
            return new DenseKernelLayer(*this, true);
        }
        inline void empty_cache() {
            in_mat_cache = Matrix<Scalar>();
        }
        inline typename Root::Data pass_forward(typename Root::Data in, bool training) {
            assert((Dimensions<std::size_t,Root::DATA_RANK>(in.dimensions()).template demote<>()) == Base::input_dims);
            assert(in.dimension(0) > 0);
            std::size_t rows = in.dimension(0);
            out_conversion_dims[0] = rows;
            prev_out_conversion_dims[0] = rows;
            in_mat_cache = MatrixMap<Scalar>(in.data(), in.dimension(0), Base::input_dims.get_volume());
            Matrix<Scalar> out_mat = (in_mat_cache * Base::weights->get_values()).rowwise() +
                    Base::bias->get_values().row(0);
            return TensorMap<Scalar,Root::DATA_RANK>(out_mat.data(), out_conversion_dims);
        }
        inline typename Root::Data pass_back(typename Root::Data out_grad) {
            assert((Dimensions<std::size_t,Root::DATA_RANK>(out_grad.dimensions()).template demote<>()) == Base::output_dims);
            assert(out_grad.dimension(0) > 0 && out_conversion_dims[0] == out_grad.dimension(0));
            // Compute the gradient of the outputs with respect to the weights and the bias.
            MatrixMap<Scalar> out_grad_mat(out_grad.data(), out_grad.dimension(0), Base::output_dims.get_volume());
            Base::weights->accumulate_grad(in_mat_cache.transpose() * out_grad_mat);
            Base::bias->accumulate_grad(out_grad_mat.colwise().sum());
            if (Base::is_input_layer())
                return typename Root::Data();
            // Compute the gradient of the previous layer's output.
            Matrix<Scalar> prev_out_grad_mat = out_grad_mat * Base::weights->get_values().transpose();
            return TensorMap<Scalar,Root::DATA_RANK>(prev_out_grad_mat.data(), prev_out_conversion_dims);
        }
    private:
        RankwiseArray out_conversion_dims, prev_out_conversion_dims;
        // Staged computation caches
        Matrix<Scalar> in_mat_cache;
    };


    /**
    * An abstract base class template for a transposed 2D convolutional layer.
    */
    template<typename Scalar, std::size_t Rank>
    class TransConvKernelLayerBase : public KernelLayer<Scalar,Rank> {
        typedef Layer<Scalar,Rank> Root;
        typedef KernelLayer<Scalar,Rank> Base;
        typedef std::array<std::size_t,4> Array4;
        typedef std::array<std::pair<std::size_t,std::size_t>,4> PaddingsArray4;
    public:
        inline void empty_cache() {
            in_mat_cache = Matrix<Scalar>();
        }
    protected:
        inline TransConvKernelLayerBase(const typename Root::Dims& input_dims, std::size_t filters,
                std::size_t receptor_height, std::size_t receptor_width, std::size_t vertical_padding,
                std::size_t horizontal_padding, std::size_t vertical_stride, std::size_t horizontal_stride,
                std::size_t vertical_dilation, std::size_t horizontal_dilation, ParamInitSharedPtr<Scalar> weight_init,
                ParamRegSharedPtr<Scalar> weight_reg, Scalar weight_clip, Scalar weight_max_l1_norm,
                Scalar weight_max_l2_norm, Scalar weight_grad_clip, Scalar weight_grad_max_l1_norm,
                Scalar weight_grad_max_l2_norm, ParamRegSharedPtr<Scalar> bias_reg, Scalar bias_clip,
                Scalar bias_max_l1_norm, Scalar bias_max_l2_norm, Scalar bias_grad_clip,
                Scalar bias_grad_max_l1_norm, Scalar bias_grad_max_l2_norm) :
                    Base(input_dims, calculate_adjusted_output_dims(input_dims, filters, receptor_height,
                            receptor_width, vertical_padding, horizontal_padding, vertical_stride, horizontal_stride,
                            vertical_dilation, horizontal_dilation),
                            std::make_shared<StandardParameters<Scalar>>(input_dims.template extend<3 - Rank>()(2),
                                    receptor_height * receptor_width * filters, true, weight_init, weight_reg,
                                    weight_clip, weight_max_l1_norm, weight_max_l2_norm, weight_grad_clip,
                                    weight_grad_max_l1_norm, weight_grad_max_l2_norm),
                            std::make_shared<StandardParameters<Scalar>>(1, calculate_adjusted_output_dims(input_dims,
                                    filters, receptor_height, receptor_width, vertical_padding, horizontal_padding,
                                    vertical_stride, horizontal_stride, vertical_dilation, horizontal_dilation)
                                    .get_volume(), true, std::make_shared<ZeroParameterInitialization<Scalar>>(),
                                    bias_reg, bias_clip, bias_max_l1_norm, bias_max_l2_norm, bias_grad_clip,
                                    bias_grad_max_l1_norm, bias_grad_max_l2_norm)),
                    filters(filters),
                    receptor_height(receptor_height),
                    receptor_width(receptor_width),
                    vertical_padding(vertical_padding),
                    horizontal_padding(horizontal_padding),
                    vertical_stride(vertical_stride),
                    horizontal_stride(horizontal_stride),
                    vertical_dilation(vertical_dilation),
                    horizontal_dilation(horizontal_dilation),
                    ext_input_dims(input_dims.template extend<3 - Rank>()),
                    ext_output_dims(calculate_output_dims(ext_input_dims, filters, receptor_height, receptor_width,
                            vertical_padding, horizontal_padding, vertical_stride, horizontal_stride,
                            vertical_dilation, horizontal_dilation)),
                    padded_height(ext_output_dims(0) + 2 * vertical_padding),
                    padded_width(ext_output_dims(1) + 2 * horizontal_padding),
                    dil_receptor_height(receptor_height + (receptor_height - 1) * vertical_dilation),
                    dil_receptor_width(receptor_width + (receptor_width - 1) * horizontal_dilation),
                    patches_per_sample(ext_input_dims(0) * ext_input_dims(1)),
                    prev_out_conversion_dims({ 0u, ext_input_dims(0), ext_input_dims(1), ext_input_dims(2) }),
                    patch_offsets({ 0u, 0u, 0u, 0u }),
                    patch_extents({ 0u, dil_receptor_height, dil_receptor_width, filters }),
                    dil_strides({ 1u, vertical_dilation + 1u, horizontal_dilation + 1u, 1u }),
                    no_padding_offsets({ 0u, vertical_padding, horizontal_padding, 0u }),
                    no_padding_extents({ 0u, ext_output_dims(0), ext_output_dims(1), ext_output_dims(2) }),
                    paddings({ std::make_pair(0, 0), std::make_pair(vertical_padding, vertical_padding),
                            std::make_pair(horizontal_padding, horizontal_padding), std::make_pair(0, 0) }) {
            assert(filters > 0);
            assert(receptor_height > 0);
            assert(receptor_width > 0);
            assert(vertical_stride > 0 && horizontal_stride > 0);
            assert(ext_output_dims(0) + 2 * vertical_padding >= dil_receptor_height &&
                    ext_output_dims(1) + 2 * horizontal_padding >= dil_receptor_width);
        }
        inline TransConvKernelLayerBase(const TransConvKernelLayerBase<Scalar,Rank>& layer, bool share_params = false) :
                Base(layer, share_params),
                filters(layer.filters),
                receptor_height(layer.receptor_height),
                receptor_width(layer.receptor_width),
                vertical_padding(layer.vertical_padding),
                horizontal_padding(layer.horizontal_padding),
                vertical_stride(layer.vertical_stride),
                horizontal_stride(layer.horizontal_stride),
                vertical_dilation(layer.vertical_dilation),
                horizontal_dilation(layer.horizontal_dilation),
                ext_input_dims(layer.ext_input_dims),
                ext_output_dims(layer.ext_output_dims),
                padded_height(layer.padded_height),
                padded_width(layer.padded_width),
                dil_receptor_height(layer.dil_receptor_height),
                dil_receptor_width(layer.dil_receptor_width),
                patches_per_sample(layer.patches_per_sample),
                prev_out_conversion_dims(layer.prev_out_conversion_dims),
                patch_offsets(layer.patch_offsets),
                patch_extents(layer.patch_extents),
                dil_strides(layer.dil_strides),
                no_padding_offsets(layer.no_padding_offsets),
                no_padding_extents(layer.no_padding_extents),
                paddings(layer.paddings),
                in_mat_cache(layer.in_mat_cache) { }
        inline Tensor<Scalar,4> _pass_forward(Tensor<Scalar,4> in, bool training) {
            std::size_t rows = in.dimension(0);
            std::size_t depth = ext_input_dims(2);
            std::size_t receptor_vol = Base::weights->get_values().cols();
            std::size_t total_patches = rows * patches_per_sample;
            in_mat_cache = MatrixMap<Scalar>(in.data(), total_patches, depth);
            Matrix<Scalar> out_conv_mat = in_mat_cache * Base::weights->get_values();
            /* Given the values of the stretched out receptor patches, accumulate them in the output tensor. */
            Tensor<Scalar,4> out(rows, padded_height, padded_width, ext_output_dims(2));
            out.setZero();
            std::size_t patch_ind = 0;
            patch_extents[0] = rows;
            for (std::size_t i = 0; i <= padded_width - dil_receptor_width; i += horizontal_stride) {
                patch_offsets[2] = i;
                for (std::size_t j = 0; j <= padded_height - dil_receptor_height; j += vertical_stride) {
                    patch_offsets[1] = j;
                    // Accumulate the gradients where the receptor-patch-tensors overlap.
                    Matrix<Scalar> out_conv_mat_block = out_conv_mat.block(patch_ind, 0, rows, receptor_vol);
                    TensorMap<Scalar,4> out_patch(out_conv_mat_block.data(), rows, receptor_height,
                            receptor_width, ext_output_dims(2));
                    if (vertical_dilation > 0 || horizontal_dilation > 0)
                        out.slice(patch_offsets, patch_extents).stride(dil_strides) += out_patch;
                    else
                        out.slice(patch_offsets, patch_extents) += out_patch;
                    patch_ind += rows;
                }
            }
            assert(patch_ind == total_patches);
            if (vertical_padding > 0 || horizontal_padding > 0) {
                // Cut off the padding.
                no_padding_extents[0] = rows;
                out = Tensor<Scalar,4>(out.slice(no_padding_offsets, no_padding_extents));
            }
            MatrixMap<Scalar> out_mat(out.data(), rows, ext_output_dims.get_volume());
            out_mat.rowwise() += Base::bias->get_values().row(0);
            return out;
        }
        inline Tensor<Scalar,4> _pass_back(Tensor<Scalar,4> out_grad) {
            std::size_t rows = out_grad.dimension(0);
            std::size_t depth = ext_input_dims(2);
            std::size_t receptor_vol = Base::weights->get_values().cols();
            std::size_t total_patches = rows * patches_per_sample;
            std::size_t patch_ind = 0;
            patch_extents[0] = rows;
            // Compute the gradient of the bias.
            Base::bias->accumulate_grad(MatrixMap<Scalar>(out_grad.data(), rows,
                    ext_output_dims.get_volume()).colwise().sum());
            // Spatial padding.
            if (vertical_padding > 0 || horizontal_padding > 0)
                out_grad = Tensor<Scalar,4>(out_grad.pad(paddings));
            Matrix<Scalar> out_grad_conv_mat(total_patches, receptor_vol);
            for (std::size_t i = 0; i <= padded_width - dil_receptor_width; i += horizontal_stride) {
                patch_offsets[2] = i;
                for (std::size_t j = 0; j <= padded_height - dil_receptor_height; j += vertical_stride) {
                    patch_offsets[1] = j;
                    Tensor<Scalar,4> patch;
                    // If the patch is dilated, skip the spatial gaps when flattening it into a matrix.
                    if (vertical_dilation > 0 || horizontal_dilation > 0)
                        patch = out_grad.slice(patch_offsets, patch_extents).stride(dil_strides);
                    else
                        patch = out_grad.slice(patch_offsets, patch_extents);
                    out_grad_conv_mat.block(patch_ind, 0, rows, receptor_vol) = MatrixMap<Scalar>(patch.data(),
                            rows, receptor_vol);
                    patch_ind += rows;
                }
            }
            assert(patch_ind == total_patches);
            Base::weights->accumulate_grad(in_mat_cache.transpose() * out_grad_conv_mat);
            if (Base::is_input_layer())
                return Tensor<Scalar,4>();
            Matrix<Scalar> prev_out_grad = out_grad_conv_mat * Base::weights->get_values().transpose();
            prev_out_conversion_dims[0] = rows;
            return TensorMap<Scalar,4>(prev_out_grad.data(), prev_out_conversion_dims);
        }
        // The defining attributes of the deconvolutional layer.
        const std::size_t filters;
        const std::size_t receptor_height;
        const std::size_t receptor_width;
        const std::size_t vertical_padding;
        const std::size_t horizontal_padding;
        const std::size_t vertical_stride;
        const std::size_t horizontal_stride;
        const std::size_t vertical_dilation;
        const std::size_t horizontal_dilation;
    private:
        inline static std::size_t calculate_spatial_output_dim(std::size_t input_dim, std::size_t receptor_size,
                std::size_t padding, std::size_t dilation, std::size_t stride) {
            return (input_dim - 1) * stride + receptor_size + (receptor_size - 1) * dilation - 2 * padding;
        }
        inline static Dimensions<std::size_t,3> calculate_output_dims(const Dimensions<std::size_t,3>& input_dims,
                std::size_t filters, std::size_t receptor_height, std::size_t receptor_width,
                std::size_t vertical_padding, std::size_t horizontal_padding, std::size_t vertical_stride,
                std::size_t horizontal_stride, std::size_t vertical_dilation, std::size_t horizontal_dilation) {
            return { calculate_spatial_output_dim(input_dims(0), receptor_height, vertical_padding, vertical_dilation,
                            vertical_stride), calculate_spatial_output_dim(input_dims(1), receptor_width,
                                    horizontal_padding, horizontal_dilation,horizontal_stride), filters };
        }
        inline static Dimensions<std::size_t,Rank> calculate_adjusted_output_dims(const typename Root::Dims& input_dims,
                std::size_t filters, std::size_t receptor_height, std::size_t receptor_width,
                std::size_t vertical_padding, std::size_t horizontal_padding, std::size_t vertical_stride,
                std::size_t horizontal_stride, std::size_t vertical_dilation, std::size_t horizontal_dilation) {
            auto output_dims = calculate_output_dims(input_dims.template extend<3 - Rank>(), filters, receptor_height,
                    receptor_width, vertical_padding, horizontal_padding, vertical_stride, horizontal_stride,
                    vertical_dilation, horizontal_dilation);
            output_dims(2) /= filters;
            output_dims(Rank - 1) *= filters;
            return output_dims.template contract<3 - Rank>();
        }
        const Dimensions<std::size_t,3> ext_input_dims, ext_output_dims;
        // Pre-computed values to improve propagation-time performance.
        const std::size_t padded_height, padded_width, dil_receptor_height, dil_receptor_width, patches_per_sample;
        Array4 prev_out_conversion_dims, patch_offsets, patch_extents, dil_strides,
                no_padding_offsets, no_padding_extents;
        PaddingsArray4 paddings;
        // Staged computation caches
        Matrix<Scalar> in_mat_cache;
    };

    /**
    * A class template for a transposed 2D convolutional layer operating on rank-3 data batches (rank-4 tensors).
    * The results of the convolutions of the filters and the input tensor are concatenated along the highest (4th)
    * rank of the output tensor.
    *
    * \see https://arxiv.org/abs/1603.07285v1
    */
    template<typename Scalar, std::size_t Rank = 3>
    class TransConvKernelLayer : public TransConvKernelLayerBase<Scalar,Rank> {
        typedef Layer<Scalar,3> Root;
        typedef KernelLayer<Scalar,3> KernelBase;
        typedef TransConvKernelLayerBase<Scalar,3> TransConvBase;
    public:
        /**
        * @param input_dims The dimensionality of the observations to be processed by the layer.
        * @param filters The number of filters to use.
        * @param weight_init A shared pointer to a weight initialization used to initialize the weights of
        * the layer.
        * @param receptor_height The height of the base of the receptor cuboid.
        * @param receptor_width The width of the base of the receptor cuboid.
        * @param vertical_padding The extent of padding to apply to the input tensor along its height (both
        * at the top and at the bottom).
        * @param horizontal_padding The extent of padding to apply to the input tensor along its width (both
        * at the left and at the right).
        * @param vertical_stride The vertical transposed convolution stride i.e. the number of elements by
        * which the receptor is to be shifted along the height of the input tensor.
        * @param horizontal_stride The horizonzal transposed convolution stride i.e. the number of elements
        * by which the receptor is to be shifted along the width of the input tensor.
        * @param vertical_dilation The extent of vertical dilation to apply to the receptor.
        * @param horizontal_dilation The extent of horizontal dilation to apply to the receptor.
        * @param weight_reg An optional regularization function to apply to the weights.
        * @param weight_clip The maximum allowed absolute weight value. If it is 0 or less, no value
        * clipping is performed.
        * @param weight_max_l1_norm The maximum allowed L1 weight value norm. If it is 0 or less, no L1
        * max norm constraint is enforced.
        * @param weight_max_l2_norm The maximum allowed L2 weight value norm. If it is 0 or less, no L2
        * max norm constraint is enforced.
        * @param weight_grad_clip The maximum allowed absolute weight gradient. If it is 0 or less, no
        * gradient clipping is performed.
        * @param weight_grad_max_l1_norm The maximum allowed L1 weight gradient norm. If it is 0 or less,
        * no L1 gradient max norm constraint is enforced.
        * @param weight_grad_max_l2_norm The maximum allowed L2 weight gradient norm. If it is 0 or less,
        * no L2 gradient max norm constraint is enforced.
        * @param bias_reg An optional regularization function to apply to the bias.
        * @param bias_clip The maximum allowed absolute bias value. If it is 0 or less, no value clipping
        * is performed.
        * @param bias_max_l1_norm The maximum allowed L1 bias value norm. If it is 0 or less, no bias L1
        * max norm constraint is enforced.
        * @param bias_max_l2_norm The maximum allowed L2 bias value norm. If it is 0 or less, no bias L2
        * max norm constraint is enforced.
        * @param bias_grad_clip The maximum allowed absolute bias gradient. If it is 0 or less, no
        * gradient clipping is performed.
        * @param bias_grad_max_l1_norm The maximum allowed L1 bias gradient norm. If it is 0 or less, no
        * bias L1 gradient max norm constraint is enforced.
        * @param bias_grad_max_l2_norm The maximum allowed L2 bias gradient norm. If it is 0 or less, no
        * bias L2 gradient max norm constraint is enforced.
        */
        inline TransConvKernelLayer(const typename Root::Dims& input_dims, std::size_t filters,
                ParamInitSharedPtr<Scalar> weight_init, std::size_t receptor_height = 3, std::size_t receptor_width = 3,
                std::size_t vertical_padding = 1, std::size_t horizontal_padding = 1, std::size_t vertical_stride = 1,
                std::size_t horizontal_stride = 1, std::size_t vertical_dilation = 0, std::size_t horizontal_dilation = 0,
                ParamRegSharedPtr<Scalar> weight_reg = nullptr, Scalar weight_clip = 0, Scalar weight_max_l1_norm = 0,
                Scalar weight_max_l2_norm = 0, Scalar weight_grad_clip = 0, Scalar weight_grad_max_l1_norm = 0,
                Scalar weight_grad_max_l2_norm = 0, ParamRegSharedPtr<Scalar> bias_reg = nullptr, Scalar bias_clip = 0,
                Scalar bias_max_l1_norm = 0, Scalar bias_max_l2_norm = 0, Scalar bias_grad_clip = 0,
                Scalar bias_grad_max_l1_norm = 0, Scalar bias_grad_max_l2_norm = 0) :
                    TransConvBase::TransConvKernelLayerBase(input_dims, filters, receptor_height, receptor_width,
                            vertical_padding, horizontal_padding, vertical_stride, horizontal_stride, vertical_dilation,
                            horizontal_dilation,weight_init,  weight_reg, weight_clip, weight_max_l1_norm,
                            weight_max_l2_norm, weight_grad_clip, weight_grad_max_l1_norm, weight_grad_max_l2_norm,
                            bias_reg, bias_clip, bias_max_l1_norm, bias_max_l2_norm, bias_grad_clip,
                            bias_grad_max_l1_norm, bias_grad_max_l2_norm) { }
        inline TransConvKernelLayer(const TransConvKernelLayer<Scalar,3>& layer, bool share_params = false) :
                TransConvBase::TransConvKernelLayerBase(layer, share_params),
                batch_size(layer.batch_size) { }
        inline Root* clone() const {
            return new TransConvKernelLayer(*this);
        }
        inline Root* clone_with_shared_params() {
            return new TransConvKernelLayer(*this, true);
        }
        inline typename Root::Data pass_forward(typename Root::Data in, bool training) {
            assert((Dimensions<std::size_t,4>(in.dimensions()).template demote<>()) == KernelBase::input_dims);
            assert(in.dimension(0) > 0);
            batch_size = in.dimension(0);
            return TransConvBase::_pass_forward(std::move(in), training);
        }
        inline typename Root::Data pass_back(typename Root::Data out_grad) {
            assert((Dimensions<std::size_t,4>(out_grad.dimensions()).template demote<>()) == KernelBase::output_dims);
            assert(out_grad.dimension(0) > 0 && batch_size == out_grad.dimension(0));
            return TransConvBase::_pass_back(std::move(out_grad));
        }
    private:
        std::size_t batch_size;
    };

    /**
    * A class template for a transposed 2D convolutional layer operating on rank-2 data batches (rank-3 tensors).
    * The results of the convolutions of the filters and the input tensor are concatenated along the highest (3rd)
    * rank of the output tensor.
    *
    * \see https://arxiv.org/abs/1603.07285v1
    */
    template<typename Scalar>
    class TransConvKernelLayer<Scalar,2> : public TransConvKernelLayerBase<Scalar,2> {
        typedef Layer<Scalar,2> Root;
        typedef KernelLayer<Scalar,2> KernelBase;
        typedef TransConvKernelLayerBase<Scalar,2> TransConvBase;
    public:
        /**
        * @param input_dims The dimensionality of the observations to be processed by the layer.
        * @param filters The number of filters to use.
        * @param weight_init A shared pointer to a weight initialization used to initialize the weights of
        * the layer.
        * @param receptor_height The height of the base of the receptor cuboid.
        * @param receptor_width The width of the base of the receptor cuboid.
        * @param vertical_padding The extent of padding to apply to the input tensor along its height (both
        * at the top and at the bottom).
        * @param horizontal_padding The extent of padding to apply to the input tensor along its width (both
        * at the left and at the right).
        * @param vertical_stride The vertical transposed convolution stride i.e. the number of elements by
        * which the receptor is to be shifted along the height of the input tensor.
        * @param horizontal_stride The horizonzal transposed convolution stride i.e. the number of elements
        * by which the receptor is to be shifted along the width of the input tensor.
        * @param vertical_dilation The extent of vertical dilation to apply to the receptor.
        * @param horizontal_dilation The extent of horizontal dilation to apply to the receptor.
        * @param weight_reg An optional regularization function to apply to the weights.
        * @param weight_clip The maximum allowed absolute weight value. If it is 0 or less, no value
        * clipping is performed.
        * @param weight_max_l1_norm The maximum allowed L1 weight value norm. If it is 0 or less, no L1
        * max norm constraint is enforced.
        * @param weight_max_l2_norm The maximum allowed L2 weight value norm. If it is 0 or less, no L2
        * max norm constraint is enforced.
        * @param weight_grad_clip The maximum allowed absolute weight gradient. If it is 0 or less, no
        * gradient clipping is performed.
        * @param weight_grad_max_l1_norm The maximum allowed L1 weight gradient norm. If it is 0 or less,
        * no L1 gradient max norm constraint is enforced.
        * @param weight_grad_max_l2_norm The maximum allowed L2 weight gradient norm. If it is 0 or less,
        * no L2 gradient max norm constraint is enforced.
        * @param bias_reg An optional regularization function to apply to the bias.
        * @param bias_clip The maximum allowed absolute bias value. If it is 0 or less, no value clipping
        * is performed.
        * @param bias_max_l1_norm The maximum allowed L1 bias value norm. If it is 0 or less, no bias L1
        * max norm constraint is enforced.
        * @param bias_max_l2_norm The maximum allowed L2 bias value norm. If it is 0 or less, no bias L2
        * max norm constraint is enforced.
        * @param bias_grad_clip The maximum allowed absolute bias gradient. If it is 0 or less, no
        * gradient clipping is performed.
        * @param bias_grad_max_l1_norm The maximum allowed L1 bias gradient norm. If it is 0 or less, no
        * bias L1 gradient max norm constraint is enforced.
        * @param bias_grad_max_l2_norm The maximum allowed L2 bias gradient norm. If it is 0 or less, no
        * bias L2 gradient max norm constraint is enforced.
        */
        inline TransConvKernelLayer(const typename Root::Dims& input_dims, std::size_t filters,
                ParamInitSharedPtr<Scalar> weight_init, std::size_t receptor_height = 3, std::size_t receptor_width = 3,
                std::size_t vertical_padding = 1, std::size_t horizontal_padding = 1, std::size_t vertical_stride = 1,
                std::size_t horizontal_stride = 1, std::size_t vertical_dilation = 0, std::size_t horizontal_dilation = 0,
                ParamRegSharedPtr<Scalar> weight_reg = nullptr, Scalar weight_clip = 0, Scalar weight_max_l1_norm = 0,
                Scalar weight_max_l2_norm = 0, Scalar weight_grad_clip = 0, Scalar weight_grad_max_l1_norm = 0,
                Scalar weight_grad_max_l2_norm = 0, ParamRegSharedPtr<Scalar> bias_reg = nullptr, Scalar bias_clip = 0,
                Scalar bias_max_l1_norm = 0, Scalar bias_max_l2_norm = 0, Scalar bias_grad_clip = 0,
                Scalar bias_grad_max_l1_norm = 0, Scalar bias_grad_max_l2_norm = 0) :
                    TransConvBase::TransConvKernelLayerBase(input_dims, filters, receptor_height, receptor_width,
                            vertical_padding, horizontal_padding, vertical_stride, horizontal_stride, vertical_dilation,
                            horizontal_dilation, weight_init, weight_reg, weight_clip, weight_max_l1_norm,
                            weight_max_l2_norm, weight_grad_clip, weight_grad_max_l1_norm, weight_grad_max_l2_norm,
                            bias_reg, bias_clip, bias_max_l1_norm, bias_max_l2_norm, bias_grad_clip,
                            bias_grad_max_l1_norm, bias_grad_max_l2_norm) { }
        inline TransConvKernelLayer(const TransConvKernelLayer<Scalar,2>& layer, bool share_params = false) :
                TransConvBase::TransConvKernelLayerBase(layer, share_params),
                batch_size(layer.batch_size) { }
        inline Root* clone() const {
            return new TransConvKernelLayer(*this);
        }
        inline Root* clone_with_shared_params() {
            return new TransConvKernelLayer(*this, true);
        }
        inline typename Root::Data pass_forward(typename Root::Data in, bool training) {
            assert((Dimensions<std::size_t,3>(in.dimensions()).template demote<>()) == KernelBase::input_dims);
            assert(in.dimension(0) > 0);
            batch_size = in.dimension(0);
            return TransConvBase::_pass_forward(TensorMap<Scalar,4>(in.data(), { batch_size, in.dimension(1),
                    in.dimension(2), 1u }), training).reshape(std::array<std::size_t,3>({ batch_size,
                            KernelBase::output_dims(0), KernelBase::output_dims(1) }));
        }
        inline typename Root::Data pass_back(typename Root::Data out_grad) {
            assert((Dimensions<std::size_t,3>(out_grad.dimensions()).template demote<>()) == KernelBase::output_dims);
            assert(out_grad.dimension(0) > 0 && batch_size == out_grad.dimension(0));
            Tensor<Scalar,4> prev_out_grad = TransConvBase::_pass_back(TensorMap<Scalar,4>(out_grad.data(),
                    { batch_size, KernelBase::output_dims(0), KernelBase::output_dims(1) / TransConvBase::filters,
                            TransConvBase::filters }));
            if (KernelBase::is_input_layer())
                return typename Root::Data();
            return TensorMap<Scalar,3>(prev_out_grad.data(), { batch_size, KernelBase::input_dims(0),
                    KernelBase::input_dims(1) });
        }
    private:
        std::size_t batch_size;
    };

    /**
    * A class template for a transposed 2D convolutional layer operating on rank-1 data batches (rank-2 tensors).
    * The results of the convolutions of the filters and the input tensor are concatenated along the highest (2nd)
    * rank of the output tensor.
    *
    * \see https://arxiv.org/abs/1603.07285v1
    */
    template<typename Scalar>
    class TransConvKernelLayer<Scalar,1> : public TransConvKernelLayerBase<Scalar,1> {
        typedef Layer<Scalar,1> Root;
        typedef KernelLayer<Scalar,1> KernelBase;
        typedef TransConvKernelLayerBase<Scalar,1> TransConvBase;
    public:
        /**
        * @param input_dims The dimensionality of the observations to be processed by the layer.
        * @param filters The number of filters to use.
        * @param weight_init A shared pointer to a weight initialization used to initialize the weights of
        * the layer.
        * @param receptor_length The length of the receptor.
        * @param padding The extent of padding to apply to the input tensor along its length on both ends.
        * @param stride The convolution stride i.e. the number of elements by which the receptor is to be
        * shifted along the length of the input tensor.
        * @param dilation The extent of dilation to apply to the receptor.
        * @param weight_reg An optional regularization function to apply to the weights.
        * @param weight_clip The maximum allowed absolute weight value. If it is 0 or less, no value
        * clipping is performed.
        * @param weight_max_l1_norm The maximum allowed L1 weight value norm. If it is 0 or less, no L1
        * max norm constraint is enforced.
        * @param weight_max_l2_norm The maximum allowed L2 weight value norm. If it is 0 or less, no L2
        * max norm constraint is enforced.
        * @param weight_grad_clip The maximum allowed absolute weight gradient. If it is 0 or less, no
        * gradient clipping is performed.
        * @param weight_grad_max_l1_norm The maximum allowed L1 weight gradient norm. If it is 0 or less,
        * no L1 gradient max norm constraint is enforced.
        * @param weight_grad_max_l2_norm The maximum allowed L2 weight gradient norm. If it is 0 or less,
        * no L2 gradient max norm constraint is enforced.
        * @param bias_reg An optional regularization function to apply to the bias.
        * @param bias_clip The maximum allowed absolute bias value. If it is 0 or less, no value clipping
        * is performed.
        * @param bias_max_l1_norm The maximum allowed L1 bias value norm. If it is 0 or less, no bias L1
        * max norm constraint is enforced.
        * @param bias_max_l2_norm The maximum allowed L2 bias value norm. If it is 0 or less, no bias L2
        * max norm constraint is enforced.
        * @param bias_grad_clip The maximum allowed absolute bias gradient. If it is 0 or less, no
        * gradient clipping is performed.
        * @param bias_grad_max_l1_norm The maximum allowed L1 bias gradient norm. If it is 0 or less, no
        * bias L1 gradient max norm constraint is enforced.
        * @param bias_grad_max_l2_norm The maximum allowed L2 bias gradient norm. If it is 0 or less, no
        * bias L2 gradient max norm constraint is enforced.
        */
        TransConvKernelLayer(const typename Root::Dims& input_dims, std::size_t filters,
                ParamInitSharedPtr<Scalar> weight_init, std::size_t receptor_length = 3, std::size_t padding = 1,
                std::size_t stride = 1, std::size_t dilation = 0, ParamRegSharedPtr<Scalar> weight_reg = nullptr,
                Scalar weight_clip = 0, Scalar weight_max_l1_norm = 0, Scalar weight_max_l2_norm = 0,
                Scalar weight_grad_clip = 0, Scalar weight_grad_max_l1_norm = 0, Scalar weight_grad_max_l2_norm = 0,
                ParamRegSharedPtr<Scalar> bias_reg = nullptr, Scalar bias_clip = 0, Scalar bias_max_l1_norm = 0,
                Scalar bias_max_l2_norm = 0, Scalar bias_grad_clip = 0, Scalar bias_grad_max_l1_norm = 0,
                Scalar bias_grad_max_l2_norm = 0) :
                    TransConvBase::TransConvKernelLayerBase(input_dims, filters, receptor_length, 1, padding, 0,
                            stride, 1, dilation, 0, weight_init, weight_reg, weight_clip, weight_max_l1_norm,
                            weight_max_l2_norm, weight_grad_clip, weight_grad_max_l1_norm, weight_grad_max_l2_norm,
                            bias_reg, bias_clip, bias_max_l1_norm, bias_max_l2_norm, bias_grad_clip,
                            bias_grad_max_l1_norm, bias_grad_max_l2_norm) { }
        inline TransConvKernelLayer(const TransConvKernelLayer<Scalar,1>& layer, bool share_params = false) :
                TransConvBase::TransConvKernelLayerBase(layer, share_params),
                batch_size(layer.batch_size) { }
        inline Root* clone() const {
            return new TransConvKernelLayer(*this);
        }
        inline Root* clone_with_shared_params() {
            return new TransConvKernelLayer(*this, true);
        }
        inline typename Root::Data pass_forward(typename Root::Data in, bool training) {
            assert((Dimensions<std::size_t,2>(in.dimensions()).template demote<>()) == KernelBase::input_dims);
            assert(in.dimension(0) > 0);
            batch_size = in.dimension(0);
            return TransConvBase::_pass_forward(TensorMap<Scalar,4>(in.data(), { batch_size, in.dimension(1), 1u, 1u }),
                    training).reshape(std::array<std::size_t,2>({ batch_size, KernelBase::output_dims(0) }));
        }
        inline typename Root::Data pass_back(typename Root::Data out_grad) {
            assert((Dimensions<std::size_t,2>(out_grad.dimensions()).template demote<>()) == KernelBase::output_dims);
            assert(out_grad.dimension(0) > 0 && batch_size == out_grad.dimension(0));
            Tensor<Scalar,4> prev_out_grad = TransConvBase::_pass_back(TensorMap<Scalar,4>(out_grad.data(),
                    { batch_size, KernelBase::output_dims(0) / TransConvBase::filters, 1, TransConvBase::filters }));
            if (KernelBase::is_input_layer())
                return typename Root::Data();
            return TensorMap<Scalar,2>(prev_out_grad.data(), { batch_size, KernelBase::input_dims(0) });
        }
    private:
        std::size_t batch_size;
    };

    /**
    * An abstract base class template representing a pooling layer.
    */
    template<typename Scalar, std::size_t Rank>
    class PoolLayer : public Layer<Scalar,Rank> {
        typedef Layer<Scalar,Rank> Base;
    public:
        virtual Base* clone() const = 0;
        inline Base* clone_with_shared_params() {
            return clone();
        }
        inline const Base& get_params_owner() const {
            return *this;
        }
        inline const typename Base::Dims& get_input_dims() const {
            return input_dims;
        }
        inline const typename Base::Dims& get_output_dims() const {
            return output_dims;
        }
        inline bool is_input_layer() const {
            return input_layer;
        }
        inline void set_input_layer(bool input_layer) {
            this->input_layer = input_layer;
        }
        inline std::vector<const Parameters<Scalar>*> get_params() const {
            return std::vector<const Parameters<Scalar>*>(0);
        }
        inline std::vector<Parameters<Scalar>*> get_params() {
            return std::vector<Parameters<Scalar>*>(0);
        }
    protected:
        typedef std::array<std::size_t,4> Array4;
        typedef std::array<std::size_t,2> ReductionRanksArray2;
        inline PoolLayer(const typename Base::Dims& input_dims, std::size_t receptor_height, std::size_t receptor_width,
                std::size_t vertical_stride, std::size_t horizontal_stride) :
                    ext_input_dims(input_dims.template extend<3 - Rank>()),
                    ext_output_dims(calculate_output_dims(ext_input_dims, receptor_height, receptor_width,
                            vertical_stride, horizontal_stride)),
                    input_dims(input_dims),
                    output_dims(ext_output_dims.template contract<3 - Rank>()),
                    receptor_height(receptor_height),
                    receptor_width(receptor_width),
                    vertical_stride(vertical_stride),
                    horizontal_stride(horizontal_stride),
                    height_rem(ext_input_dims(0) - receptor_height),
                    width_rem(ext_input_dims(1) - receptor_width),
                    input_layer(false),
                    reduction_ranks({ 1u, 2u }),
                    broadcast({ 1u, receptor_height, receptor_width, 1u }),
                    patch_offsets({ 0u, 0u, 0u, 0u }),
                    patch_extents({ 0u, receptor_height, receptor_width, ext_input_dims(2) }),
                    reduced_patch_offsets({ 0u, 0u, 0u, 0u }),
                    reduced_patch_extents({ 0u, 1u, 1u, ext_input_dims(2) }) {
            assert(receptor_height > 0 && receptor_width > 0);
            assert(vertical_stride > 0 && horizontal_stride > 0);
            assert(ext_input_dims(0) >= receptor_height && ext_input_dims(1) >= receptor_width);
        }
        inline Tensor<Scalar,4> _pass_forward(Tensor<Scalar,4> in, bool training) {
            std::size_t rows = in.dimension(0);
            patch_extents[0] = rows;
            reduced_patch_extents[0] = rows;
            Tensor<Scalar,4> out(rows, ext_output_dims(0), ext_output_dims(1), ext_output_dims(2));
            _init_cache();
            std::size_t patch_ind = 0;
            std::size_t out_i = 0;
            for (std::size_t i = 0; i <= width_rem; i += horizontal_stride, ++out_i) {
                patch_offsets[2] = i;
                reduced_patch_offsets[2] = out_i;
                std::size_t out_j = 0;
                for (std::size_t j = 0; j <= height_rem; j += vertical_stride, ++out_j) {
                    patch_offsets[1] = j;
                    reduced_patch_offsets[1] = out_j;
                    Tensor<Scalar,4> patch = in.slice(patch_offsets, patch_extents);
                    out.slice(reduced_patch_offsets, reduced_patch_extents) = _reduce(patch, patch_ind++);
                }
            }
            return out;
        }
        inline Tensor<Scalar,4> _pass_back(Tensor<Scalar,4> out_grad) {
            Tensor<Scalar,4> prev_out_grad(patch_extents[0], ext_input_dims(0), ext_input_dims(1),  ext_input_dims(2));
            prev_out_grad.setZero();
            std::size_t patch_ind = 0;
            std::size_t out_grad_i = 0;
            for (std::size_t i = 0; i <= width_rem; i += horizontal_stride, ++out_grad_i) {
                patch_offsets[2] = i;
                reduced_patch_offsets[2] = out_grad_i;
                std::size_t out_grad_j = 0;
                for (std::size_t j = 0; j <= height_rem; j += vertical_stride, ++out_grad_j) {
                    patch_offsets[1] = j;
                    reduced_patch_offsets[1] = out_grad_j;
                    Tensor<Scalar,4> reduced_patch_grad = out_grad.slice(reduced_patch_offsets, reduced_patch_extents);
                    // Accumulate the gradients where the patches overlap.
                    prev_out_grad.slice(patch_offsets, patch_extents) += _d_reduce(reduced_patch_grad, patch_ind++);
                }
            }
            return prev_out_grad;
        }
        /**
        * Initializes the cache required for back-propagation.
        */
        virtual void _init_cache() = 0;
        /**
        * Reduces the input tensor patch along the specified ranks.
        *
        * @param patch A tensor representing a spatial patch of the input tensor.
        * @param patch_ind The index of the patch.
        * @return The reduced tensor.
        */
        virtual Tensor<Scalar,4> _reduce(const Tensor<Scalar,4>& patch, std::size_t patch_ind) = 0;
        /**
        * Differentiates the reduction function and returns the derivative of the loss function
        * w.r.t. the non-reduced patch.
        *
        * @param grad The derivative of the loss function w.r.t. the reduced patch.
        * @param patch_ind The index of the patch.
        * @return The derivative of the loss function w.r.t. the non-reduced patch.
        */
        virtual Tensor<Scalar,4> _d_reduce(const Tensor<Scalar,4>& grad, std::size_t patch_ind) = 0;
        const Dimensions<std::size_t,3> ext_input_dims, ext_output_dims;
        const typename Base::Dims input_dims, output_dims;
        const std::size_t receptor_height, receptor_width, vertical_stride, horizontal_stride, height_rem, width_rem;
        // Arrays for tensor manipulation.
        ReductionRanksArray2 reduction_ranks;
        Array4 broadcast, patch_offsets, patch_extents, reduced_patch_offsets, reduced_patch_extents, dil_strides;
    private:
        inline static std::size_t calculate_spatial_output_dim(std::size_t input_dim, std::size_t receptor_size,
                std::size_t stride) {
            return (input_dim - receptor_size) / stride + 1;
        }
        inline static Dimensions<std::size_t,3> calculate_output_dims(const Dimensions<std::size_t,3>& input_dims,
                std::size_t receptor_height, std::size_t receptor_width, std::size_t vertical_stride,
                std::size_t horizontal_stride) {
            return { calculate_spatial_output_dim(input_dims(0), receptor_height, vertical_stride),
                    calculate_spatial_output_dim(input_dims(1), receptor_width, horizontal_stride),
                    input_dims(2) };
        }
        bool input_layer;
    };


    /**
    * An abstract class template representing a pooling layer that reduces patches of the input by taking their
    * maxima.
    */
    template<typename Scalar, std::size_t Rank>
    class MaxPoolLayerBase : public PoolLayer<Scalar,Rank> {
        typedef Layer<Scalar,Rank> Root;
        typedef PoolLayer<Scalar,Rank> Base;
    public:
        inline void empty_cache() {
            max_inds = std::vector<std::vector<unsigned>>(0);
        }
    protected:
        inline MaxPoolLayerBase(const typename Root::Dims& input_dims, std::size_t receptor_height,
                std::size_t receptor_width, std::size_t vertical_stride, std::size_t horizontal_stride) :
                    Base::PoolLayer(input_dims, receptor_height, receptor_width, vertical_stride,
                        horizontal_stride) { }
        inline void _init_cache() {
            max_inds = std::vector<std::vector<unsigned>>(Base::ext_output_dims(0) * Base::ext_output_dims(1));
        }
        inline Tensor<Scalar,4> _reduce(const Tensor<Scalar,4>& patch, std::size_t patch_ind) {
            std::size_t rows = patch.dimension(0);
            std::size_t depth = patch.dimension(3);
            std::vector<unsigned> inds(rows * depth);
            Tensor<Scalar,4> reduced_patch(rows, 1u, 1u, depth);
            for (std::size_t i = 0; i < depth; ++i) {
                for (std::size_t j = 0; j < rows; ++j) {
                    Scalar max = NumericUtils<Scalar>::MIN;
                    unsigned max_height = 0;
                    unsigned max_width = 0;
                    for (std::size_t k = 0; k < Base::receptor_width; ++k) {
                        for (std::size_t l = 0; l < Base::receptor_height; ++l) {
                            Scalar val = patch(j,l,k,i);
                            if (val > max) {
                                max = val;
                                max_height = l;
                                max_width = k;
                            }
                        }
                    }
                    inds[i * rows + j] = max_width * Base::receptor_height + max_height;
                    reduced_patch(j,0u,0u,i) = max;
                }
            }
            max_inds[patch_ind] = inds;
            return reduced_patch;
        }
        inline Tensor<Scalar,4> _d_reduce(const Tensor<Scalar,4>& grad, std::size_t patch_ind) {
            std::size_t rows = grad.dimension(0);
            std::size_t depth = grad.dimension(3);
            Tensor<Scalar,4> patch(rows, Base::receptor_height, Base::receptor_width, depth);
            patch.setZero();
            std::vector<unsigned>& inds = max_inds[patch_ind];
            for (std::size_t i = 0; i < depth; ++i) {
                for (std::size_t j = 0; j < rows; ++j) {
                    unsigned max_ind = inds[i * rows + j];
                    unsigned max_height = max_ind % Base::receptor_height;
                    unsigned max_width = max_ind / Base::receptor_height;
                    patch(j,max_height,max_width,i) = grad(j,0u,0u,i);
                }
            }
            return patch;
        }
    private:
        std::vector<std::vector<unsigned>> max_inds;
    };

    /**
    * A class template representing a 2D max pooling layer operating on rank-3 data.
    *
    * \see http://ais.uni-bonn.de/papers/icann2010_maxpool.pdf
    */
    template<typename Scalar, std::size_t Rank = 3>
    class MaxPoolLayer : public MaxPoolLayerBase<Scalar,Rank> {
        typedef Layer<Scalar,3> Root;
        typedef PoolLayer<Scalar,3> PoolBase;
        typedef MaxPoolLayerBase<Scalar,3> MaxPoolBase;
    public:
        /**
        * @param input_dims The dimensionality of the input tensor.
        * @param receptor_height The height of the pooling receptor.
        * @param receptor_width The width of the pooling receptor.
        * @param vertical_stride The vertical stride at which the input is to be pooled (i.e. the number of
        * elements by which the receptor is to be shifted along the height of the input tensor).
        * @param horizontal_stride The horizontal stride at which the input is to be pooled (i.e. the number
        * of elements by which the receptor is to be shifted along the width of the input tensor).
        */
        inline MaxPoolLayer(const typename Root::Dims& input_dims, std::size_t receptor_height = 2,
                std::size_t receptor_width = 2, std::size_t vertical_stride = 2, std::size_t horizontal_stride = 2) :
                    MaxPoolBase::MaxPoolLayerBase(input_dims, receptor_height, receptor_width, vertical_stride,
                            horizontal_stride) { }
        inline Root* clone() const {
            return new MaxPoolLayer(*this);
        }
        inline typename Root::Data pass_forward(typename Root::Data in, bool training) {
            assert((Dimensions<std::size_t,4>(in.dimensions()).template demote<>()) == PoolBase::input_dims);
            assert(in.dimension(0) > 0);
            batch_size = in.dimension(0);
            return PoolBase::_pass_forward(std::move(in), training);
        }
        inline typename Root::Data pass_back(typename Root::Data out_grad) {
            assert((Dimensions<std::size_t,4>(out_grad.dimensions()).template demote<>()) == PoolBase::output_dims);
            assert(out_grad.dimension(0) > 0 && batch_size == out_grad.dimension(0));
            if (PoolBase::is_input_layer())
                return typename Root::Data();
            return PoolBase::_pass_back(std::move(out_grad));
        }
    private:
        std::size_t batch_size;
    };

    /**
    * A class template representing a 2D max pooling layer operating on rank-2 data.
    *
    * \see http://ais.uni-bonn.de/papers/icann2010_maxpool.pdf
    */
    template<typename Scalar>
    class MaxPoolLayer<Scalar,2> : public MaxPoolLayerBase<Scalar,2> {
        typedef Layer<Scalar,2> Root;
        typedef PoolLayer<Scalar,2> PoolBase;
        typedef MaxPoolLayerBase<Scalar,2> MaxPoolBase;
    public:
        /**
        * @param input_dims The dimensionality of the input tensor.
        * @param receptor_height The height of the pooling receptor.
        * @param receptor_width The width of the pooling receptor.
        * @param vertical_stride The vertical stride at which the input is to be pooled (i.e. the number of
        * elements by which the receptor is to be shifted along the height of the input tensor).
        * @param horizontal_stride The horizontal stride at which the input is to be pooled (i.e. the number
        * of elements by which the receptor is to be shifted along the width of the input tensor).
        */
        inline MaxPoolLayer(const typename Root::Dims& input_dims, std::size_t receptor_height = 2,
                std::size_t receptor_width = 2, std::size_t vertical_stride = 2, std::size_t horizontal_stride = 2) :
                    MaxPoolBase::MaxPoolLayerBase(input_dims, receptor_height, receptor_width, vertical_stride,
                            horizontal_stride) { }
        inline Root* clone() const {
            return new MaxPoolLayer(*this);
        }
        inline typename Root::Data pass_forward(typename Root::Data in, bool training) {
            assert((Dimensions<std::size_t,3>(in.dimensions()).template demote<>()) == PoolBase::input_dims);
            assert(in.dimension(0) > 0);
            batch_size = in.dimension(0);
            return PoolBase::_pass_forward(TensorMap<Scalar,4>(in.data(), { batch_size, in.dimension(1),
                    in.dimension(2), 1u }), training).reshape(std::array<std::size_t,3>({ batch_size,
                            PoolBase::output_dims(0), PoolBase::output_dims(1) }));
        }
        inline typename Root::Data pass_back(typename Root::Data out_grad) {
            assert((Dimensions<std::size_t,3>(out_grad.dimensions()).template demote<>()) == PoolBase::output_dims);
            assert(out_grad.dimension(0) > 0 && batch_size == out_grad.dimension(0));
            if (PoolBase::is_input_layer())
                return typename Root::Data();
            Tensor<Scalar,4> prev_out_grad = PoolBase::_pass_back(TensorMap<Scalar,4>(out_grad.data(),
                    { batch_size, PoolBase::ext_output_dims(0), PoolBase::ext_output_dims(1),
                            PoolBase::ext_output_dims(2) }));
            return TensorMap<Scalar,3>(prev_out_grad.data(), { batch_size, PoolBase::input_dims(0),
                    PoolBase::input_dims(1) });
        }
    private:
        std::size_t batch_size;
    };

    /**
    * A class template representing a 1D max pooling layer.
    *
    * \see http://ais.uni-bonn.de/papers/icann2010_maxpool.pdf
    */
    template<typename Scalar>
    class MaxPoolLayer<Scalar,1> : public MaxPoolLayerBase<Scalar,1> {
        typedef Layer<Scalar,1> Root;
        typedef PoolLayer<Scalar,1> PoolBase;
        typedef MaxPoolLayerBase<Scalar,1> MaxPoolBase;
    public:
        /**
        * @param input_dims The dimensionality of the input tensor.
        * @param receptor_length The length of the pooling receptor.
        * @param stride The stride at which the input is to be pooled (i.e. the number of
        * elements by which the receptor is to be shifted along the length of the input tensor).
        */
        inline MaxPoolLayer(const typename Root::Dims& input_dims, std::size_t receptor_length = 2,
                std::size_t stride = 2) :
                    MaxPoolBase::MaxPoolLayerBase(input_dims, receptor_length, 1, stride, 1) { }
        inline Root* clone() const {
            return new MaxPoolLayer(*this);
        }
        inline typename Root::Data pass_forward(typename Root::Data in, bool training) {
            assert((Dimensions<std::size_t,2>(in.dimensions()).template demote<>()) == PoolBase::input_dims);
            assert(in.dimension(0) > 0);
            batch_size = in.dimension(0);
            return PoolBase::_pass_forward(TensorMap<Scalar,4>(in.data(), { batch_size,
                    in.dimension(1), 1u, 1u }), training).reshape(std::array<std::size_t,2>({ batch_size,
                            PoolBase::output_dims(0) }));
        }
        inline typename Root::Data pass_back(typename Root::Data out_grad) {
            assert((Dimensions<std::size_t,2>(out_grad.dimensions()).template demote<>()) == PoolBase::output_dims);
            assert(out_grad.dimension(0) > 0 && batch_size == out_grad.dimension(0));
            if (PoolBase::is_input_layer())
                return typename Root::Data();
            Tensor<Scalar,4> prev_out_grad = PoolBase::_pass_back(TensorMap<Scalar,4>(out_grad.data(),
                    { batch_size, PoolBase::ext_output_dims(0), PoolBase::ext_output_dims(1),
                            PoolBase::ext_output_dims(2) }));
            return TensorMap<Scalar,2>(prev_out_grad.data(), { batch_size, PoolBase::input_dims(0) });
        }
    private:
        std::size_t batch_size;
    };

    /**
    * An abstract class template representing a pooling layer that reduces patches of the input by taking their
    * means.
    */
    template<typename Scalar, std::size_t Rank>
    class MeanPoolLayerBase : public PoolLayer<Scalar,Rank> {
        typedef Layer<Scalar,Rank> Root;
        typedef PoolLayer<Scalar,Rank> Base;
    public:
        inline void empty_cache() { }
    protected:
        inline MeanPoolLayerBase(const typename Root::Dims& input_dims, std::size_t receptor_height,
                std::size_t receptor_width, std::size_t vertical_stride, std::size_t horizontal_stride) :
                    Base::PoolLayer(input_dims, receptor_height, receptor_width, vertical_stride,
                            horizontal_stride),
                    receptor_area(receptor_height * receptor_width) { }
        inline void _init_cache() { }
        inline Tensor<Scalar,4> _reduce(const Tensor<Scalar,4>& patch, std::size_t patch_ind) {
            Tensor<Scalar,2> reduced_patch = patch.mean(Base::reduction_ranks);
            return TensorMap<Scalar,4>(reduced_patch.data(), Base::reduced_patch_extents);
        }
        inline Tensor<Scalar,4> _d_reduce(const Tensor<Scalar,4>& grad, std::size_t patch_ind) {
            return (grad / (Scalar) receptor_area).broadcast(Base::broadcast);
        }
    private:
        std::size_t receptor_area;
    };

    /**
    * A class template representing a 2D mean pooling layer operating on rank-3 data.
    */
    template<typename Scalar, std::size_t Rank = 3>
    class MeanPoolLayer : public MeanPoolLayerBase<Scalar,Rank> {
        typedef Layer<Scalar,3> Root;
        typedef PoolLayer<Scalar,3> PoolBase;
        typedef MeanPoolLayerBase<Scalar,3> MeanPoolBase;
    public:
        /**
        * @param input_dims The dimensionality of the input tensor.
        * @param receptor_height The height of the pooling receptor.
        * @param receptor_width The width of the pooling receptor.
        * @param vertical_stride The vertical stride at which the input is to be pooled (i.e. the number of
        * elements by which the receptor is to be shifted along the height of the input tensor).
        * @param horizontal_stride The horizontal stride at which the input is to be pooled (i.e. the number
        * of elements by which the receptor is to be shifted along the width of the input tensor).
        */
        inline MeanPoolLayer(const typename Root::Dims& input_dims, std::size_t receptor_height = 2,
                std::size_t receptor_width = 2, std::size_t vertical_stride = 2, std::size_t horizontal_stride = 2) :
                    MeanPoolBase::MeanPoolLayerBase(input_dims, receptor_height, receptor_width, vertical_stride,
                            horizontal_stride) { }
        inline Root* clone() const {
            return new MeanPoolLayer(*this);
        }
        inline typename Root::Data pass_forward(typename Root::Data in, bool training) {
            assert((Dimensions<std::size_t,4>(in.dimensions()).template demote<>()) == PoolBase::input_dims);
            assert(in.dimension(0) > 0);
            batch_size = in.dimension(0);
            return PoolBase::_pass_forward(std::move(in), training);
        }
        inline typename Root::Data pass_back(typename Root::Data out_grad) {
            assert((Dimensions<std::size_t,4>(out_grad.dimensions()).template demote<>()) == PoolBase::output_dims);
            assert(out_grad.dimension(0) > 0 && batch_size == out_grad.dimension(0));
            if (PoolBase::is_input_layer())
                return typename Root::Data();
            return PoolBase::_pass_back(std::move(out_grad));
        }
    private:
        std::size_t batch_size;
    };

    /**
    * A class template representing a 2D mean pooling layer operating on rank-2 data.
    */
    template<typename Scalar>
    class MeanPoolLayer<Scalar,2> : public MeanPoolLayerBase<Scalar,2> {
        typedef Layer<Scalar,2> Root;
        typedef PoolLayer<Scalar,2> PoolBase;
        typedef MeanPoolLayerBase<Scalar,2> MeanPoolBase;
    public:
        /**
        * @param input_dims The dimensionality of the input tensor.
        * @param receptor_height The height of the pooling receptor.
        * @param receptor_width The width of the pooling receptor.
        * @param vertical_stride The vertical stride at which the input is to be pooled (i.e. the number of
        * elements by which the receptor is to be shifted along the height of the input tensor).
        * @param horizontal_stride The horizontal stride at which the input is to be pooled (i.e. the number
        * of elements by which the receptor is to be shifted along the width of the input tensor).
        */
        inline MeanPoolLayer(const typename Root::Dims& input_dims, std::size_t receptor_height = 2,
                std::size_t receptor_width = 2, std::size_t vertical_stride = 2, std::size_t horizontal_stride = 2) :
                    MeanPoolBase::MeanPoolLayerBase(input_dims, receptor_height, receptor_width, vertical_stride,
                            horizontal_stride) { }
        inline Root* clone() const {
            return new MeanPoolLayer(*this);
        }
        inline typename Root::Data pass_forward(typename Root::Data in, bool training) {
            assert((Dimensions<std::size_t,3>(in.dimensions()).template demote<>()) == PoolBase::input_dims);
            assert(in.dimension(0) > 0);
            batch_size = in.dimension(0);
            return PoolBase::_pass_forward(TensorMap<Scalar,4>(in.data(), { batch_size, in.dimension(1),
                    in.dimension(2), 1u }), training).reshape(std::array<std::size_t,3>({ batch_size,
                            PoolBase::output_dims(0), PoolBase::output_dims(1) }));
        }
        inline typename Root::Data pass_back(typename Root::Data out_grad) {
            assert((Dimensions<std::size_t,3>(out_grad.dimensions()).template demote<>()) == PoolBase::output_dims);
            assert(out_grad.dimension(0) > 0 && batch_size == out_grad.dimension(0));
            if (PoolBase::is_input_layer())
                return typename Root::Data();
            Tensor<Scalar,4> prev_out_grad = PoolBase::_pass_back(TensorMap<Scalar,4>(out_grad.data(),
                    { batch_size, PoolBase::ext_output_dims(0), PoolBase::ext_output_dims(1),
                            PoolBase::ext_output_dims(2) }));
            return TensorMap<Scalar,3>(prev_out_grad.data(), { batch_size, PoolBase::input_dims(0),
                    PoolBase::input_dims(1) });
        }
    private:
        std::size_t batch_size;
    };

    /**
    * A class template representing a 1D mean pooling layer.
    */
    template<typename Scalar>
    class MeanPoolLayer<Scalar,1> : public MeanPoolLayerBase<Scalar,1> {
        typedef Layer<Scalar,1> Root;
        typedef PoolLayer<Scalar,1> PoolBase;
        typedef MeanPoolLayerBase<Scalar,1> MeanPoolBase;
    public:
        /**
        * @param input_dims The dimensionality of the input tensor.
        * @param receptor_length The length of the pooling receptor.
        * @param stride The stride at which the input is to be pooled (i.e. the number of
        * elements by which the receptor is to be shifted along the length of the input tensor).
        */
        inline MeanPoolLayer(const typename Root::Dims& input_dims, std::size_t receptor_length = 2,
                std::size_t stride = 2) :
                    MeanPoolBase::MeanPoolLayerBase(input_dims, receptor_length, 1, stride, 1) { }
        inline Root* clone() const {
            return new MeanPoolLayer(*this);
        }
        inline typename Root::Data pass_forward(typename Root::Data in, bool training) {
            assert((Dimensions<std::size_t,2>(in.dimensions()).template demote<>()) == PoolBase::input_dims);
            assert(in.dimension(0) > 0);
            batch_size = in.dimension(0);
            return PoolBase::_pass_forward(TensorMap<Scalar,4>(in.data(), { batch_size, in.dimension(1), 1u, 1u }),
                    training).reshape(std::array<std::size_t,2>({ batch_size, PoolBase::output_dims(0) }));
        }
        inline typename Root::Data pass_back(typename Root::Data out_grad) {
            assert((Dimensions<std::size_t,2>(out_grad.dimensions()).template demote<>()) == PoolBase::output_dims);
            assert(out_grad.dimension(0) > 0 && batch_size == out_grad.dimension(0));
            if (PoolBase::is_input_layer())
                return typename Root::Data();
            Tensor<Scalar,4> prev_out_grad = PoolBase::_pass_back(TensorMap<Scalar,4>(out_grad.data(),
                    { batch_size, PoolBase::ext_output_dims(0), PoolBase::ext_output_dims(1),
                            PoolBase::ext_output_dims(2) }));
            return TensorMap<Scalar,2>(prev_out_grad.data(), { batch_size, PoolBase::input_dims(0) });
        }
    private:
        std::size_t batch_size;
    };

    

    /**
    * A class template for a per-channel batch normalization layer.
    *
    * \see https://arxiv.org/abs/1502.03167
    */
    template<typename Scalar, std::size_t Rank, bool PerLastRank = (Rank == 3)>
    class BatchNormLayer : public Layer<Scalar,Rank> {
        typedef Layer<Scalar,Rank> Base;
        typedef BatchNormLayer<Scalar,Rank,PerLastRank> Self;
        typedef std::array<std::size_t,Base::DATA_RANK> RankwiseArray;
        typedef std::shared_ptr<StandardParameters<Scalar>> StdParamsSharedPtr;
    public:
        /**
        * @param dims The dimensionality of the input tensor.
        * @param norm_avg_decay The decay rate of the maintained means and variances.
        * @param epsilon A small constant used to maintain numerical stability.
        * @param gamma_reg An optional regularization function to apply to the gammas.
        * @param gamma_clip The maximum allowed absolute gamma value. If it is 0 or less, no value
        * clipping is performed.
        * @param gamma_max_l1_norm The maximum allowed L1 gamma value norm. If it is 0 or less, no L1
        * max norm constraint is enforced.
        * @param gamma_max_l2_norm The maximum allowed L2 gamma value norm. If it is 0 or less, no L2
        * max norm constraint is enforced.
        * @param gamma_grad_clip The maximum allowed absolute gamma gradient. If it is 0 or less, no
        * gradient clipping is performed.
        * @param gamma_grad_max_l1_norm The maximum allowed L1 gamma gradient norm. If it is 0 or less,
        * no L1 gradient max norm constraint is enforced.
        * @param gamma_grad_max_l2_norm The maximum allowed L2 gamma gradient norm. If it is 0 or less,
        * no L2 gradient max norm constraint is enforced.
        * @param beta_reg An optional regularization function to apply to the beta.
        * @param beta_clip The maximum allowed absolute beta value. If it is 0 or less, no value clipping
        * is performed.
        * @param beta_max_l1_norm The maximum allowed L1 beta value norm. If it is 0 or less, no beta L1
        * max norm constraint is enforced.
        * @param beta_max_l2_norm The maximum allowed L2 beta value norm. If it is 0 or less, no beta L2
        * max norm constraint is enforced.
        * @param beta_grad_clip The maximum allowed absolute beta gradient. If it is 0 or less, no
        * gradient clipping is performed.
        * @param beta_grad_max_l1_norm The maximum allowed L1 beta gradient norm. If it is 0 or less, no
        * beta L1 gradient max norm constraint is enforced.
        * @param beta_grad_max_l2_norm The maximum allowed L2 beta gradient norm. If it is 0 or less, no
        * beta L2 gradient max norm constraint is enforced.
        */
        inline BatchNormLayer(const typename Base::Dims& dims, Scalar norm_avg_decay = .1,
                Scalar epsilon = NumericUtils<Scalar>::EPSILON2, ParamRegSharedPtr<Scalar> gamma_reg = nullptr,
                Scalar gamma_clip = 0, Scalar gamma_max_l1_norm = 0, Scalar gamma_max_l2_norm = 0,
                Scalar gamma_grad_clip = 0, Scalar gamma_grad_max_l1_norm = 0, Scalar gamma_grad_max_l2_norm = 0,
                ParamRegSharedPtr<Scalar> beta_reg = nullptr, Scalar beta_clip = 0, Scalar beta_max_l1_norm = 0,
                Scalar beta_max_l2_norm = 0, Scalar beta_grad_clip = 0, Scalar beta_grad_max_l1_norm = 0,
                Scalar beta_grad_max_l2_norm = 0) :
                    owner(*this),
                    dims(dims),
                    norm_avg_decay(norm_avg_decay),
                    epsilon(epsilon),
                    channels(dims(Rank - 1)),
                    input_layer(false),
                    offsets(),
                    extents(dims.template promote<>()),
                    memb_vec(channels) {
            assert(norm_avg_decay >= 0 && norm_avg_decay <= 1 &&
                    "norm avg decay must not be less than 0 or greater than 1");
            assert(epsilon > 0 && "epsilon must be greater than 0");
            offsets.fill(0);
            extents[Rank] = 1;
            auto gamma_init = std::make_shared<OneParameterInitialization<Scalar>>();
            auto beta_init = std::make_shared<ZeroParameterInitialization<Scalar>>();
            for (std::size_t i = 0; i < channels; ++i) {
                ChannelSpecificMembers& memb = memb_vec[i];
                memb.avg_means = std::make_shared<StandardParameters<Scalar>>(1, 1, false);
                memb.avg_inv_sds = std::make_shared<StandardParameters<Scalar>>(1, 1, false);
                memb.gammas = std::make_shared<StandardParameters<Scalar>>(1, 1, true, gamma_init, gamma_reg,
                        gamma_clip, gamma_max_l1_norm, gamma_max_l2_norm, gamma_grad_clip, gamma_grad_max_l1_norm,
                        gamma_grad_max_l2_norm);
                memb.betas = std::make_shared<StandardParameters<Scalar>>(1, 1, true, beta_init, beta_reg, beta_clip,
                        beta_max_l1_norm, beta_max_l2_norm, beta_grad_clip, beta_grad_max_l1_norm,
                        beta_grad_max_l2_norm);
                memb.avgs_init = false;
            }
        }
        inline BatchNormLayer(const Self& layer, bool share_params = false) :
                owner(share_params ? layer.owner : *this),
                dims(layer.dims),
                norm_avg_decay(layer.norm_avg_decay),
                epsilon(layer.epsilon),
                channels(layer.channels),
                input_layer(layer.input_layer),
                offsets(layer.offsets),
                extents(layer.extents),
                memb_vec(channels) {
            for (std::size_t i = 0; i < channels; ++i) {
                ChannelSpecificMembers& memb1 = memb_vec[i];
                const ChannelSpecificMembers& memb2 = layer.memb_vec[i];
                if (share_params) {
                    memb1.avg_means = memb2.avg_means;
                    memb1.avg_inv_sds = memb2.avg_inv_sds;
                    memb1.gammas = memb2.gammas;
                    memb1.betas = memb2.betas;
                } else {
                    memb1.avg_means = StdParamsSharedPtr(static_cast<StandardParameters<Scalar>*>(
                            memb2.avg_means->clone()));
                    memb1.avg_inv_sds = StdParamsSharedPtr(static_cast<StandardParameters<Scalar>*>(
                            memb2.avg_inv_sds->clone()));
                    memb1.gammas = StdParamsSharedPtr(static_cast<StandardParameters<Scalar>*>(
                            memb2.gammas->clone()));
                    memb1.betas = StdParamsSharedPtr(static_cast<StandardParameters<Scalar>*>(
                            memb2.betas->clone()));
                }
                memb1.avgs_init = memb2.avgs_init;
                memb1.inv_in_sd_cache = memb2.inv_in_sd_cache;
                memb1.std_in_mat_cache = memb2.std_in_mat_cache;
            }
        }
        inline Base* clone() const {
            return new BatchNormLayer(*this);
        }
        inline Base* clone_with_shared_params() {
            return new BatchNormLayer(*this, true);
        }
        inline const Base& get_params_owner() const {
            return owner;
        }
        inline const typename Base::Dims& get_input_dims() const {
            return dims;
        }
        inline const typename Base::Dims& get_output_dims() const {
            return dims;
        }
        inline bool is_input_layer() const {
            return input_layer;
        }
        inline void set_input_layer(bool input_layer) {
            this->input_layer = input_layer;
        }
        inline std::vector<const Parameters<Scalar>*> get_params() const {
            std::vector<const Parameters<Scalar>*> params_vec;
            populate_params_vector<const Parameters<Scalar>*>(params_vec);
            return params_vec;
        }
        inline std::vector<Parameters<Scalar>*> get_params() {
            std::vector<Parameters<Scalar>*> params_vec;
            populate_params_vector<Parameters<Scalar>*>(params_vec);
            return params_vec;
        }
        inline void empty_cache() {
            for (unsigned i = 0; i < memb_vec.size(); ++i)
                memb_vec[i].std_in_mat_cache = Matrix<Scalar>();
        }
        inline typename Base::Data pass_forward(typename Base::Data in, bool training) {
            assert((Dimensions<std::size_t,Base::DATA_RANK>(in.dimensions()).template demote<>()) == dims);
            assert(in.dimension(0) > 0);
            std::size_t rows = in.dimension(0);
            extents[0] = rows;
            typename Base::Data out;
            if (channels == 1) {
                MatrixMap<Scalar> in_mat(in.data(), rows, in.size() / rows);
                Matrix<Scalar> out_mat = _pass_forward(in_mat, 0, training);
                out = TensorMap<Scalar,Base::DATA_RANK>(out_mat.data(), extents);
            } else {
                out = typename Base::Data(in.dimensions());
                for (std::size_t i = 0; i < channels; ++i) {
                    offsets[Rank] = i;
                    typename Base::Data in_slice = in.slice(offsets, extents);
                    MatrixMap<Scalar> in_slice_mat(in_slice.data(), rows, in_slice.size() / rows);
                    Matrix<Scalar> out_slice_mat = _pass_forward(in_slice_mat, i, training);
                    out.slice(offsets, extents) = TensorMap<Scalar,Base::DATA_RANK>(out_slice_mat.data(), extents);
                }
            }
            return out;
        }
        inline typename Base::Data pass_back(typename Base::Data out_grad) {
            assert((Dimensions<std::size_t,Base::DATA_RANK>(out_grad.dimensions()).template demote<>()) == dims);
            assert(out_grad.dimension(0) > 0 && extents[0] == out_grad.dimension(0));
            std::size_t rows = out_grad.dimension(0);
            typename Base::Data prev_out_grad;
            if (channels == 1) {
                MatrixMap<Scalar> out_grad_mat(out_grad.data(), rows, out_grad.size() / rows);
                if (input_layer) {
                    _pass_back(out_grad_mat, 0);
                    return typename Base::Data();
                } else {
                    Matrix<Scalar> prev_out_grad_mat = _pass_back(out_grad_mat, 0);
                    prev_out_grad = TensorMap<Scalar,Base::DATA_RANK>(prev_out_grad_mat.data(), extents);
                }
            } else {
                prev_out_grad = input_layer ? typename Base::Data() : typename Base::Data(out_grad.dimensions());
                for (std::size_t i = 0; i < channels; ++i) {
                    offsets[Rank] = i;
                    typename Base::Data out_grad_slice = out_grad.slice(offsets, extents);
                    MatrixMap<Scalar> out_grad_slice_mat(out_grad_slice.data(), rows, out_grad_slice.size() / rows);
                    if (input_layer) {
                        _pass_back(out_grad_slice_mat, i);
                        continue;
                    } else {
                        Matrix<Scalar> prev_out_grad_slice_mat = _pass_back(out_grad_slice_mat, i);
                        prev_out_grad.slice(offsets, extents) =
                                TensorMap<Scalar,Base::DATA_RANK>(prev_out_grad_slice_mat.data(), extents);
                    }
                }
            }
            return prev_out_grad;
        }
    private:
        inline Matrix<Scalar> _pass_forward(MatrixMap<Scalar>& in, std::size_t i, bool training) {
            Matrix<Scalar> out_mat;
            ChannelSpecificMembers& memb = memb_vec[i];
            if (training) {
                Scalar mean = in.mean();
                Matrix<Scalar> norm_in_mat = in.array() - mean;
                memb.inv_in_sd_cache = 1 / sqrt(norm_in_mat.array().square().mean() + epsilon);
                memb.std_in_mat_cache = norm_in_mat * memb.inv_in_sd_cache;
                out_mat = memb.std_in_mat_cache;
                if (memb.avgs_init) {
                    memb.avg_means->set_values((1 - norm_avg_decay) * memb.avg_means->get_values().array() +
                            norm_avg_decay * mean);
                    memb.avg_inv_sds->set_values((1 - norm_avg_decay) * memb.avg_inv_sds->get_values().array() +
                            norm_avg_decay * memb.inv_in_sd_cache);
                } else {
                    memb.avg_means->set_values(Matrix<Scalar>::Constant(1, 1, mean));
                    memb.avg_inv_sds->set_values(Matrix<Scalar>::Constant(1, 1, memb.inv_in_sd_cache));
                    memb.avgs_init = true;
                }
            } else {
                assert(memb.avgs_init);
                out_mat = (in.array() - memb.avg_means->get_values()(0,0)) * memb.avg_inv_sds->get_values()(0,0);
            }
            return (out_mat * memb.gammas->get_values()(0,0)).array() + memb.betas->get_values()(0,0);
        }
        inline Matrix<Scalar> _pass_back(MatrixMap<Scalar>& out_grad, std::size_t i) {
            ChannelSpecificMembers& memb = memb_vec[i];
            memb.gammas->accumulate_grad(Matrix<Scalar>::Constant(1, 1,
                    out_grad.cwiseProduct(memb.std_in_mat_cache).sum()));
            memb.betas->accumulate_grad(Matrix<Scalar>::Constant(1, 1, out_grad.sum()));
            if (input_layer)
                return Matrix<Scalar>();
            std::size_t locations = out_grad.size();
            Matrix<Scalar> std_in_grad_mat = out_grad * memb.gammas->get_values()(0,0);
            return (((locations * std_in_grad_mat).array() - std_in_grad_mat.sum()).matrix() -
                    memb.std_in_mat_cache * memb.std_in_mat_cache.cwiseProduct(std_in_grad_mat).sum()) *
                    (((Scalar) 1 / locations) * memb.inv_in_sd_cache);
        }
        template<typename _ParamsPtr>
        inline void populate_params_vector(std::vector<_ParamsPtr>& params_vec) const {
            for (std::size_t i = 0; i < channels; ++i) {
                const ChannelSpecificMembers& memb = memb_vec[i];
                params_vec.push_back(memb.avg_means.get());
                params_vec.push_back(memb.avg_inv_sds.get());
                params_vec.push_back(memb.gammas.get());
                params_vec.push_back(memb.betas.get());
            }
        }
        const Self& owner;
        const typename Base::Dims dims;
        const Scalar norm_avg_decay, epsilon;
        const std::size_t channels;
        bool input_layer;
        RankwiseArray offsets, extents;
        struct ChannelSpecificMembers {
            // Dynamic batch normalization parameters.
            StdParamsSharedPtr avg_means, avg_inv_sds;
            bool avgs_init;
            // The optimizable parameters.
            StdParamsSharedPtr gammas, betas;
            // Staged computation cache.
            Scalar inv_in_sd_cache;
            Matrix<Scalar> std_in_mat_cache;
        };
        std::vector<ChannelSpecificMembers> memb_vec;
    };

    /**
    * A class template for a per-activation batch normalization layer.
    *
    * \see https://arxiv.org/abs/1502.03167
    */
    template<typename Scalar, std::size_t Rank>
    class BatchNormLayer<Scalar,Rank,false> : public Layer<Scalar,Rank> {
        typedef Layer<Scalar,Rank> Base;
        typedef BatchNormLayer<Scalar,Rank,false> Self;
        typedef std::shared_ptr<StandardParameters<Scalar>> StdParamsSharedPtr;
    public:
        /**
        * @param dims The dimensionality of the input tensor.
        * @param norm_avg_decay The decay rate of the maintained means and variances.
        * @param epsilon A small constant used to maintain numerical stability.
        * @param gamma_reg An optional regularization function to apply to the gammas.
        * @param gamma_clip The maximum allowed absolute gamma value. If it is 0 or less, no value
        * clipping is performed.
        * @param gamma_max_l1_norm The maximum allowed L1 gamma value norm. If it is 0 or less, no L1
        * max norm constraint is enforced
        * @param gamma_max_l2_norm The maximum allowed L2 gamma value norm. If it is 0 or less, no L2
        * max norm constraint is enforced.
        * @param gamma_grad_clip The maximum allowed absolute gamma gradient. If it is 0 or less, no
        * gradient clipping is performed.
        * @param gamma_grad_max_l1_norm The maximum allowed L1 gamma gradient norm. If it is 0 or less,
        * no L1 gradient max norm constraint is enforced.
        * @param gamma_grad_max_l2_norm The maximum allowed L2 gamma gradient norm. If it is 0 or less,
        * no L2 gradient max norm constraint is enforced.
        * @param beta_reg An optional regularization function to apply to the beta.
        * @param beta_clip The maximum allowed absolute beta value. If it is 0 or less, no value clipping
        * is performed.
        * @param beta_max_l1_norm The maximum allowed L1 beta value norm. If it is 0 or less, no beta L1
        * max norm constraint is enforced.
        * @param beta_max_l2_norm The maximum allowed L2 beta value norm. If it is 0 or less, no beta L2
        * max norm constraint is enforced.
        * @param beta_grad_clip The maximum allowed absolute beta gradient. If it is 0 or less, no
        * gradient clipping is performed.
        * @param beta_grad_max_l1_norm The maximum allowed L1 beta gradient norm. If it is 0 or less, no
        * beta L1 gradient max norm constraint is enforced.
        * @param beta_grad_max_l2_norm The maximum allowed L2 beta gradient norm. If it is 0 or less, no
        * beta L2 gradient max norm constraint is enforced.
        */
        inline BatchNormLayer(const typename Base::Dims& dims, Scalar norm_avg_decay = .1,
                Scalar epsilon = NumericUtils<Scalar>::EPSILON2, ParamRegSharedPtr<Scalar> gamma_reg = nullptr,
                Scalar gamma_clip = 0, Scalar gamma_max_l1_norm = 0, Scalar gamma_max_l2_norm = 0,
                Scalar gamma_grad_clip = 0, Scalar gamma_grad_max_l1_norm = 0, Scalar gamma_grad_max_l2_norm = 0,
                ParamRegSharedPtr<Scalar> beta_reg = nullptr, Scalar beta_clip = 0, Scalar beta_max_l1_norm = 0,
                Scalar beta_max_l2_norm = 0, Scalar beta_grad_clip = 0, Scalar beta_grad_max_l1_norm = 0,
                Scalar beta_grad_max_l2_norm = 0) :
                    owner(*this),
                    dims(dims),
                    norm_avg_decay(norm_avg_decay),
                    epsilon(epsilon),
                    avg_means(std::make_shared<StandardParameters<Scalar>>(1, dims.get_volume(), false)),
                    avg_inv_sds(std::make_shared<StandardParameters<Scalar>>(1, dims.get_volume(), false)),
                    avgs_init(false),
                    gammas(std::make_shared<StandardParameters<Scalar>>(1, dims.get_volume(), true,
                            std::make_shared<OneParameterInitialization<Scalar>>(), gamma_reg, gamma_clip,
                            gamma_max_l1_norm, gamma_max_l2_norm, gamma_grad_clip, gamma_grad_max_l1_norm,
                            gamma_grad_max_l2_norm)),
                    betas(std::make_shared<StandardParameters<Scalar>>(1, dims.get_volume(), true,
                            std::make_shared<ZeroParameterInitialization<Scalar>>(), beta_reg, beta_clip,
                            beta_max_l1_norm, beta_max_l2_norm, beta_grad_clip, beta_grad_max_l1_norm,
                            beta_grad_max_l2_norm)),
                    input_layer(false) {
            assert(norm_avg_decay >= 0 && norm_avg_decay <= 1 &&
                    "norm avg decay must not be less than 0 or greater than 1");
            assert(epsilon > 0 && "epsilon must be greater than 0");
        }
        inline BatchNormLayer(const Self& layer, bool share_params = false) :
                owner(share_params ? layer.owner : *this),
                dims(layer.dims),
                norm_avg_decay(layer.norm_avg_decay),
                epsilon(layer.epsilon),
                input_layer(layer.input_layer),
                avg_means(share_params ? layer.avg_means :
                        StdParamsSharedPtr(static_cast<StandardParameters<Scalar>*>(layer.avg_means->clone()))),
                avg_inv_sds(share_params ? layer.avg_inv_sds :
                        StdParamsSharedPtr(static_cast<StandardParameters<Scalar>*>(layer.avg_inv_sds->clone()))),
                avgs_init(layer.avgs_init),
                gammas(share_params ? layer.gammas :
                        StdParamsSharedPtr(static_cast<StandardParameters<Scalar>*>(layer.gammas->clone()))),
                betas(share_params ? layer.betas :
                        StdParamsSharedPtr(static_cast<StandardParameters<Scalar>*>(layer.betas->clone()))),
                inv_in_sd_cache(layer.inv_in_sd_cache),
                std_in_mat_cache(layer.std_in_mat_cache) { }
        inline Base* clone() const {
            return new BatchNormLayer(*this);
        }
        inline Base* clone_with_shared_params() {
            return new BatchNormLayer(*this, true);
        }
        inline const Base& get_params_owner() const {
            return owner;
        }
        inline const typename Base::Dims& get_input_dims() const {
            return dims;
        }
        inline const typename Base::Dims& get_output_dims() const {
            return dims;
        }
        inline bool is_input_layer() const {
            return input_layer;
        }
        inline void set_input_layer(bool input_layer) {
            this->input_layer = input_layer;
        }
        inline std::vector<Parameters<Scalar>*> get_params() {
            return std::vector<Parameters<Scalar>*>({ avg_means.get(), avg_inv_sds.get(),
                    gammas.get(), betas.get() });
        }
        inline std::vector<const Parameters<Scalar>*> get_params() const  {
            return std::vector<const Parameters<Scalar>*>({ avg_means.get(), avg_inv_sds.get(),
                    gammas.get(), betas.get() });
        }
        inline void empty_cache() {
            inv_in_sd_cache = RowVector<Scalar>();
            std_in_mat_cache = Matrix<Scalar>();
        }
        inline typename Base::Data pass_forward(typename Base::Data in, bool training) {
            assert((Dimensions<std::size_t,Base::DATA_RANK>(in.dimensions()).template demote<>()) == dims);
            assert(in.dimension(0) > 0);
            std::size_t rows = in.dimension(0);
            MatrixMap<Scalar> in_mat(in.data(), rows, in.size() / rows);
            if (training) {
                RowVector<Scalar> mean_vec = in_mat.colwise().mean();
                Matrix<Scalar> norm_in_mat = in_mat.rowwise() - mean_vec;
                inv_in_sd_cache = (norm_in_mat.array().square().colwise().mean() + epsilon).sqrt().inverse();
                std_in_mat_cache = norm_in_mat * inv_in_sd_cache.asDiagonal();
                in_mat = std_in_mat_cache;
                // Maintain a moving average of means and variances for testing.
                if (avgs_init) {
                    avg_means->set_values((1 - norm_avg_decay) * avg_means->get_values() + norm_avg_decay * mean_vec);
                    avg_inv_sds->set_values((1 - norm_avg_decay) * avg_inv_sds->get_values() +
                            norm_avg_decay * inv_in_sd_cache);
                } else {
                    avg_means->set_values(mean_vec);
                    avg_inv_sds->set_values(inv_in_sd_cache);
                    avgs_init = true;
                }
            } else {
                // For testing, use the moving averages.
                assert(avgs_init);
                in_mat = (in_mat.rowwise() - avg_means->get_values().row(0)) * avg_inv_sds->get_values().asDiagonal();
            }
            Matrix<Scalar> out_mat = (in_mat * gammas->get_values().asDiagonal()).rowwise() +
                    betas->get_values().row(0);
            return TensorMap<Scalar,Base::DATA_RANK>(out_mat.data(), in.dimensions());
        }
        inline typename Base::Data pass_back(typename Base::Data out_grad) {
            assert((Dimensions<std::size_t,Base::DATA_RANK>(out_grad.dimensions()).template demote<>()) == dims);
            assert(out_grad.dimension(0) > 0 && std_in_mat_cache.rows() == out_grad.dimension(0));
            std::size_t rows = out_grad.dimension(0);
            MatrixMap<Scalar> out_grad_mat(out_grad.data(), rows, out_grad.size() / rows);
            gammas->accumulate_grad(out_grad_mat.cwiseProduct(std_in_mat_cache).colwise().sum());
            betas->accumulate_grad(out_grad_mat.colwise().sum());
            if (input_layer)
                return typename Base::Data();
            Matrix<Scalar> std_in_grad_mat = out_grad_mat * gammas->get_values().asDiagonal();
            Matrix<Scalar> prev_out_grad_mat = (((rows * std_in_grad_mat).rowwise() -
                    std_in_grad_mat.colwise().sum()) - std_in_mat_cache *
                    (std_in_mat_cache.cwiseProduct(std_in_grad_mat).colwise().sum().asDiagonal())) *
                    (((Scalar) 1 / rows) * inv_in_sd_cache).asDiagonal();
            return TensorMap<Scalar,Base::DATA_RANK>(prev_out_grad_mat.data(), out_grad.dimensions());
        }
    private:
        const Self& owner;
        const typename Base::Dims dims;
        const Scalar norm_avg_decay, epsilon;
        bool input_layer;
        // Dynamic batch normalization parameters.
        StdParamsSharedPtr avg_means, avg_inv_sds;
        bool avgs_init;
        // Betas and gammas
        StdParamsSharedPtr gammas, betas;
        // Staged computation caches.
        RowVector<Scalar> inv_in_sd_cache;
        Matrix<Scalar> std_in_mat_cache;
    };


    /**
    * A class template representing a broadcasting layer that repeats the contents of its input tensors
    * along its ranks.
    */
    template<typename Scalar, std::size_t Rank>
    class BroadcastLayer : public Layer<Scalar,Rank> {
        typedef Layer<Scalar,Rank> Base;
        typedef std::array<std::size_t,Base::DATA_RANK> RankwiseArray;
    public:
        /**
        * @param input_dims The nominal input dimensions of the layer.
        * @param broadcast The number of times the input tensor's contents are
        * repeated along each rank. All elements should be greater than 0.
        */
        inline BroadcastLayer(const typename Base::Dims& input_dims, const typename Base::Dims& broadcast) :
                    input_dims(input_dims),
                    output_dims(input_dims * broadcast),
                    input_layer(false),
                    broadcast(broadcast.template promote<>()) {
            slice_offsets.fill(0);
            for (std::size_t i = 0; i < Rank; ++i)
                assert(broadcast(i) > 0);
        }
        inline Base* clone() const {
            return new BroadcastLayer(*this);
        }
        inline Base* clone_with_shared_params() {
            return clone();
        }
        inline const Base& get_params_owner() const {
            return *this;
        }
        inline const typename Base::Dims& get_input_dims() const {
            return input_dims;
        }
        inline const typename Base::Dims& get_output_dims() const {
            return output_dims;
        }
        inline bool is_input_layer() const {
            return input_layer;
        }
        inline void set_input_layer(bool input_layer) {
            this->input_layer = input_layer;
        }
        inline std::vector<const Parameters<Scalar>*> get_params() const {
            return std::vector<const Parameters<Scalar>*>(0);
        }
        inline std::vector<Parameters<Scalar>*> get_params() {
            return std::vector<Parameters<Scalar>*>(0);
        }
        inline void empty_cache() { }
        inline typename Base::Data pass_forward(typename Base::Data in, bool training) {
            assert((Dimensions<std::size_t,Base::DATA_RANK>(in.dimensions()).template demote<>()) == input_dims);
            assert(in.dimension(0) > 0);
            rows = in.dimension(0);
            return in.broadcast(broadcast);
        }
        inline typename Base::Data pass_back(typename Base::Data out_grad) {
            assert((Dimensions<std::size_t,Base::DATA_RANK>(out_grad.dimensions()).template demote<>()) == output_dims);
            assert(out_grad.dimension(0) > 0 && rows == out_grad.dimension(0));
            if (input_layer)
                return typename Base::Data();
            typename Base::Data prev_out_grad = std::move(out_grad);
            slice_offsets.fill(0);
            slice_extents = output_dims.template promote<>();
            slice_extents[0] = rows;
            for (std::size_t i = 0; i < Rank; ++i) {
                if (broadcast[i + 1] <= 1)
                    continue;
                slice_extents[i + 1] = input_dims(i);
                typename Base::Data work_tensor(slice_extents);
                work_tensor.setZero();
                for (std::size_t j = 0; j < broadcast[i + 1]; ++j) {
                    work_tensor += prev_out_grad.slice(slice_offsets, slice_extents);
                    slice_offsets[i + 1] += input_dims(i);
                }
                slice_offsets[i + 1] = 0;
                prev_out_grad = std::move(work_tensor);
            }
            return prev_out_grad;
        }
    private:
        const typename Base::Dims input_dims, output_dims;
        RankwiseArray broadcast, slice_offsets, slice_extents;
        std::size_t rows;
        bool input_layer;
    };

    /**
    * A class template representing a drop-out layer.
    *
    * \see https://arxiv.org/abs/1207.0580
    * \see http://www.cs.toronto.edu/~rsalakhu/papers/srivastava14a.pdf
    */
    template<typename Scalar, std::size_t Rank>
    class DropoutLayer : public Layer<Scalar,Rank> {
        typedef Layer<Scalar,Rank> Base;
    public:
        /**
        * @param dims The dimensionality of the input tensor.
        * @param dropout_prob The probability of an element of the input tensor being set to 0.
        * @param epsilon A small constant used to maintain numerical stability.
        */
        inline DropoutLayer(const typename Base::Dims& dims, Scalar dropout_prob,
                Scalar epsilon = NumericUtils<Scalar>::EPSILON2) :
                    dims(dims),
                    dropout_prob(dropout_prob),
                    epsilon(epsilon),
                    input_layer(false) {
            assert(dropout_prob > 0 && dropout_prob <= 1 &&
                    "dropout probability must be greater than 0 and no greater than 1");
            assert(epsilon > 0 && "epsilon must be greater than 0");
        }
        inline Base* clone() const {
            return new DropoutLayer(*this);
        }
        inline Base* clone_with_shared_params() {
            return clone();
        }
        inline const Base& get_params_owner() const {
            return *this;
        }
        inline const typename Base::Dims& get_input_dims() const {
            return dims;
        }
        inline const typename Base::Dims& get_output_dims() const {
            return dims;
        }
        inline bool is_input_layer() const {
            return input_layer;
        }
        inline void set_input_layer(bool input_layer) {
            this->input_layer = input_layer;
        }
        inline std::vector<const Parameters<Scalar>*> get_params() const {
            return std::vector<const Parameters<Scalar>*>(0);
        }
        inline std::vector<Parameters<Scalar>*> get_params() {
            return std::vector<Parameters<Scalar>*>(0);
        }
        inline void empty_cache() {
            dropout_mask = typename Base::Data();
        }
        inline typename Base::Data pass_forward(typename Base::Data in, bool training) {
            assert((Dimensions<std::size_t,Base::DATA_RANK>(in.dimensions()).template demote<>()) == dims);
            assert(in.dimension(0) > 0);
            if (training) {
                // Inverted dropout.
                Scalar scaling_factor = (Scalar) 1 / (1 - dropout_prob + epsilon);
                dropout_mask = in.random().unaryExpr([this,scaling_factor](Scalar e) {
                    return (Scalar) (e <= dropout_prob ? 0 : scaling_factor);
                });
                return in * dropout_mask;
            }
            return std::move(in);
        }
        inline typename Base::Data pass_back(typename Base::Data out_grad) {
            assert((Dimensions<std::size_t,Base::DATA_RANK>(out_grad.dimensions()).template demote<>()) == dims);
            assert(out_grad.dimension(0) > 0 && dropout_mask.dimension(0) == out_grad.dimension(0));
            if (input_layer)
                return typename Base::Data();
            // The derivative of the dropout function.
            return out_grad * dropout_mask;
        }
    private:
        const typename Base::Dims dims;
        const Scalar dropout_prob, epsilon;
        bool input_layer;
        // Staged computation cache.
        typename Base::Data dropout_mask;
    };

    /**
    * A class template representing a reshaping layer that outputs a reshaped copy of the input
    * tensor with the same volume. The data backing the tensor is not shifted in any way.
    */
    template<typename Scalar, std::size_t Rank>
    class ReshapeLayer : public Layer<Scalar,Rank> {
        typedef Layer<Scalar,Rank> Base;
        typedef std::array<std::size_t,Base::DATA_RANK> RankwiseArray;
    public:
        /**
        * @param input_dims The nominal input dimensions of the layer.
        * @param output_dims The dimensions of the reshaped tensor. The output tensor must have
        * the same volume as the input tensor.
        */
        inline ReshapeLayer(const typename Base::Dims& input_dims, const typename Base::Dims& output_dims) :
                    input_dims(input_dims),
                    output_dims(output_dims),
                    input_layer(false),
                    input_conversion_dims(output_dims.template promote<>()),
                    output_conversion_dims(input_dims.template promote<>()) {
            assert(input_dims.get_volume() == output_dims.get_volume());
        }
        inline Base* clone() const {
            return new ReshapeLayer(*this);
        }
        inline Base* clone_with_shared_params() {
            return clone();
        }
        inline const Base& get_params_owner() const {
            return *this;
        }
        inline const typename Base::Dims& get_input_dims() const {
            return input_dims;
        }
        inline const typename Base::Dims& get_output_dims() const {
            return output_dims;
        }
        inline bool is_input_layer() const {
            return input_layer;
        }
        inline void set_input_layer(bool input_layer) {
            this->input_layer = input_layer;
        }
        inline std::vector<const Parameters<Scalar>*> get_params() const {
            return std::vector<const Parameters<Scalar>*>(0);
        }
        inline std::vector<Parameters<Scalar>*> get_params() {
            return std::vector<Parameters<Scalar>*>(0);
        }
        inline void empty_cache() { }
        inline typename Base::Data pass_forward(typename Base::Data in, bool training) {
            assert((Dimensions<std::size_t,Base::DATA_RANK>(in.dimensions()).template demote<>()) == input_dims);
            assert(in.dimension(0) > 0);
            input_conversion_dims[0] = in.dimension(0);
            return in.reshape(input_conversion_dims);
        }
        inline typename Base::Data pass_back(typename Base::Data out_grad) {
            assert((Dimensions<std::size_t,Base::DATA_RANK>(out_grad.dimensions()).template demote<>()) == output_dims);
            assert(out_grad.dimension(0) > 0 && input_conversion_dims[0] == out_grad.dimension(0));
            if (input_layer)
                return typename Base::Data();
            output_conversion_dims[0] = input_conversion_dims[0];
            return out_grad.reshape(output_conversion_dims);
        }
    private:
        const typename Base::Dims input_dims, output_dims;
        RankwiseArray input_conversion_dims, output_conversion_dims;
        bool input_layer;
    };

    /**
    * An abstract class template for loss functions for both sequential and non-sequential data.
    * Implementations of this class should be stateless.
    */
    template<typename Scalar, std::size_t Rank, bool Sequential>
    class UniversalLoss : public Loss<Scalar,Rank,Sequential> {
        typedef Loss<Scalar,Rank,Sequential> Base;
    protected:
        typedef std::array<std::size_t,Base::DATA_RANK> RankwiseArray;
        /**
        * It computes the loss of a batch of non-sequential data.
        *
        * @param out The output tensor.
        * @param obj The objective tensor.
        * @return A column vector representing the losses of the samples in the batch.
        */
        virtual ColVector<Scalar> _function(typename Base::Data out, typename Base::Data obj) const = 0;
        /**
        * It computes the gradient of the output batch.
        *
        * @param out The output tensor.
        * @param obj The objective tensor.
        * @param grad_dims The dimensions of the gradient tensor.
        * @return The gradient tensor of the output batch.
        */
        virtual typename Base::Data _d_function(typename Base::Data out, typename Base::Data obj,
                const RankwiseArray& grad_dims) const = 0;
    public:
        virtual ~UniversalLoss() = default;
        inline ColVector<Scalar> function(typename Base::Data out, typename Base::Data obj) const {
            assert(out.dimensions() == obj.dimensions());
            return _function(std::move(out), std::move(obj));
        }
        inline typename Base::Data d_function(typename Base::Data out, typename Base::Data obj) const {
            assert(out.dimensions() == obj.dimensions());
            RankwiseArray dims = out.dimensions();
            return _d_function(std::move(out), std::move(obj), dims);
        }
    };

    /**
    * Partial template specialization for sequential data. Implementations
    * of this class should be stateless.
    */
    template<typename Scalar, std::size_t Rank>
    class UniversalLoss<Scalar,Rank,true> : public Loss<Scalar,Rank,true> {
        typedef Loss<Scalar,Rank,true> Base;
    protected:
        typedef std::array<std::size_t,Base::DATA_RANK> RankwiseArray;
        /**
        * It computes the loss of a single time step in a batch. The total loss of the batch is the sum of the losses
        * of all its time steps.
        *
        * @param out The output tensor.
        * @param obj The objective tensor.
        * @return A column vector representing the losses of the samples in the batch for the given time step.
        */
        virtual ColVector<Scalar> _function(typename Base::Data out, typename Base::Data obj) const = 0;
        /**
        * It computes the gradient of a single time step of the output sequence batch.
        *
        * @param out The output tensor.
        * @param obj The objective tensor.
        * @param grad_dims The dimensions of the gradient tensor.
        * @return The gradient tensor of the provided time step of the output batch.
        */
        virtual typename Base::Data _d_function(typename Base::Data out, typename Base::Data obj,
                const RankwiseArray& grad_dims) const = 0;
    public:
        virtual ~UniversalLoss() = default;
        inline ColVector<Scalar> function(typename Base::Data out, typename Base::Data obj) const {
            assert(out.dimensions() == obj.dimensions());
            int time_steps = out.dimension(1);
            if (time_steps == 1)
                return _function(std::move(out), std::move(obj));
            RankwiseArray offsets;
            RankwiseArray extents = out.dimensions();
            offsets.fill(0);
            extents[1] = 1;
            ColVector<Scalar> loss = ColVector<Scalar>::Zero(out.dimension(0), 1);
            for (int i = 0; i < time_steps; ++i) {
                offsets[1] = i;
                typename Base::Data out_i = out.slice(offsets, extents);
                typename Base::Data obj_i = obj.slice(offsets, extents);
                loss += _function(std::move(out_i), std::move(obj_i));
            }
            return loss;
        }
        inline typename Base::Data d_function(const typename Base::Data out, const typename Base::Data obj) const {
            assert(out.dimensions() == obj.dimensions());
            int time_steps = out.dimension(1);
            if (time_steps == 1)
                return _d_function(std::move(out), std::move(obj), out.dimensions());
            RankwiseArray offsets;
            RankwiseArray extents = out.dimensions();
            offsets.fill(0);
            typename Base::Data grads(extents);
            extents[1] = 1;
            grads.setZero();
            for (int i = 0; i < time_steps; ++i) {
                offsets[1] = i;
                typename Base::Data out_i = out.slice(offsets, extents);
                typename Base::Data obj_i = obj.slice(offsets, extents);
                grads.slice(offsets, extents) += _d_function(std::move(out_i), std::move(obj_i), extents);
            }
            return grads;
        }
    };

    /**
    * A template class representing the absolute error (L1) loss function.
    *
    * \f$L_i = \left|\hat{y_i} - y_i\right|\f$
    */
    template<typename Scalar, std::size_t Rank, bool Sequential>
    class AbsoluteLoss : public UniversalLoss<Scalar,Rank,Sequential> {
        typedef Loss<Scalar,Rank,Sequential> Root;
        typedef UniversalLoss<Scalar,Rank,Sequential> Base;
    protected:
        inline ColVector<Scalar> _function(typename Root::Data out, typename Root::Data obj) const {
            std::size_t rows = out.dimension(0);
            std::size_t cols = out.size() / rows;
            return (MatrixMap<Scalar>(out.data(), rows, cols) - MatrixMap<Scalar>(obj.data(), rows, cols))
                    .array().abs().rowwise().sum();
        }
        inline typename Root::Data _d_function(typename Root::Data out, typename Root::Data obj,
                const typename Base::RankwiseArray& grad_dims) const {
            typename Root::Data diff = out - obj;
            return diff.unaryExpr([this](Scalar i) { return (Scalar) (i >= 0 ? 1 : -1); });
        }
    };


    /**
    * A template class representing the binary cross entropy loss function. The objective
    * is expected to be a size-1 tensor with values in the range [0, 1].
    *
    * \f$L_i = -(y_i \ln(\hat{y_i} + \epsilon) + (1 - y_i) \ln(1 + \epsilon - \hat{y_i}))\f$
    */
    template<typename Scalar, std::size_t Rank, bool Sequential>
    class BinaryCrossEntropyLoss : public UniversalLoss<Scalar,Rank,Sequential> {
        typedef Loss<Scalar,Rank,Sequential> Root;
        typedef UniversalLoss<Scalar,Rank,Sequential> Base;
    public:
        /**
        * @param epsilon A small constant used to maintain numerical stability.
        */
        BinaryCrossEntropyLoss(Scalar epsilon = NumericUtils<Scalar>::EPSILON2) :
                epsilon(epsilon) { };
    protected:
        inline ColVector<Scalar> _function(typename Root::Data out, typename Root::Data obj) const {
            assert(out.size() == out.dimension(0));
            typename Root::Data loss = -(obj * (out + out.constant(epsilon)).log() +
                    (obj.constant(1) - obj) * (out.constant(1 + epsilon) - out).log());
            return MatrixMap<Scalar>(loss.data(), out.dimension(0), 1);
        }
        inline typename Root::Data _d_function(typename Root::Data out, typename Root::Data obj,
                const typename Base::RankwiseArray& grad_dims) const {
            assert(out.size() == out.dimension(0));
            return -(obj / (out + out.constant(epsilon)) -
                    (obj.constant(1) - obj) / (out.constant(1 + epsilon) - out));
        }
    private:
        Scalar epsilon;
    };

    /**
    * A template class representing the cross entropy loss function. This class assumes the objective
    * values for each sample (and time step) to be in the range [0, 1].
    *
    * \f$L_i = -\ln(\hat{y_i} + \epsilon) y_i\f$
    */
    template<typename Scalar, std::size_t Rank, bool Sequential>
    class CrossEntropyLoss : public UniversalLoss<Scalar,Rank,Sequential> {
        typedef Loss<Scalar,Rank,Sequential> Root;
        typedef UniversalLoss<Scalar,Rank,Sequential> Base;
    public:
        /**
        * @param epsilon A small constant used to maintain numerical stability.
        */
        CrossEntropyLoss(Scalar epsilon = NumericUtils<Scalar>::EPSILON2) :
                epsilon(epsilon) { };
    protected:
        inline ColVector<Scalar> _function(typename Root::Data out, typename Root::Data obj) const {
            std::size_t rows = out.dimension(0);
            std::size_t cols = out.size() / rows;
            return -((MatrixMap<Scalar>(out.data(), rows, cols).array() + epsilon).log() *
                    MatrixMap<Scalar>(obj.data(), rows, cols).array()).matrix().rowwise().sum();
        }
        inline typename Root::Data _d_function(typename Root::Data out, typename Root::Data obj,
                const typename Base::RankwiseArray& grad_dims) const {
            return -obj / (out + epsilon);
        }
    private:
        Scalar epsilon;
    };


    /**
    * A template class representing the hinge loss function. This class assumes the objectives for
    * each sample (and time step) to be a one-hot vector (tensor rank).
    *
    * \f$L_i = \sum_{j \neq y_i} \max(0, \hat{y_i}_j - \hat{y_i}_{y_i} + 1)\f$ or
    * \f$L_i = \sum_{j \neq y_i} \max(0, \hat{y_i}_j - \hat{y_i}_{y_i} + 1)^2\f$
    */
    template<typename Scalar, std::size_t Rank, bool Sequential, bool Squared = false>
    class HingeLoss : public UniversalLoss<Scalar,Rank,Sequential> {
        typedef Loss<Scalar,Rank,Sequential> Root;
        typedef UniversalLoss<Scalar,Rank,Sequential> Base;
    protected:
        inline ColVector<Scalar> _function(typename Root::Data out, typename Root::Data obj) const {
            std::size_t rows = out.dimension(0);
            std::size_t cols = out.size() / rows;
            MatrixMap<Scalar> out_mat(out.data(), rows, cols);
            MatrixMap<Scalar> obj_mat(obj.data(), rows, cols);
            ColVector<Scalar> loss(rows);
            for (int i = 0; i < rows; ++i) {
                unsigned ones = 0;
                int correct_class_ind = -1;
                for (int j = 0; j < cols; ++j) {
                    Scalar obj_ij = obj_mat(i,j);
                    assert((NumericUtils<Scalar>::almost_equal(obj_ij, (Scalar) 0) ||
                            NumericUtils<Scalar>::almost_equal(obj_ij, (Scalar) 1)));
                    if (NumericUtils<Scalar>::almost_equal(obj_ij, (Scalar) 1)) {
                        ones++;
                        correct_class_ind = j;
                    }
                }
                assert(ones == 1);
                Scalar loss_i = 0;
                Scalar correct_class_score = out_mat(i,correct_class_ind);
                for (int j = 0; j < cols; ++j) {
                    if (j == correct_class_ind)
                        continue;
                    Scalar loss_ij = std::max((Scalar) 0, (Scalar) (out_mat(i,j) - correct_class_score + 1));
                    loss_i += Squared ? loss_ij * loss_ij : loss_ij;
                }
                loss(i) = loss_i;
            }
            return loss;
        }
        inline typename Root::Data _d_function(typename Root::Data out, typename Root::Data obj,
                const typename Base::RankwiseArray& grad_dims) const {
            std::size_t rows = out.dimension(0);
            std::size_t cols = out.size() / rows;
            MatrixMap<Scalar> out_mat(out.data(), rows, cols);
            MatrixMap<Scalar> obj_mat(obj.data(), rows, cols);
            Matrix<Scalar> out_grad(rows, cols);
            for (int i = 0; i < rows; ++i) {
                unsigned ones = 0;
                int correct_class_ind = -1;
                for (int j = 0; j < cols; ++j) {
                    Scalar obj_ij = obj_mat(i,j);
                    assert((NumericUtils<Scalar>::almost_equal(obj_ij, (Scalar) 0) ||
                            NumericUtils<Scalar>::almost_equal(obj_ij, (Scalar) 1)));
                    if (NumericUtils<Scalar>::almost_equal(obj_ij, (Scalar) 1)) {
                        ones++;
                        correct_class_ind = j;
                    }
                }
                assert(ones == 1);
                Scalar total_out_grad = 0;
                Scalar correct_class_score = out_mat(i,correct_class_ind);
                for (int j = 0; j < cols; ++j) {
                    if (j == correct_class_ind)
                        continue;
                    Scalar out_ij = out_mat(i,j);
                    Scalar margin = out_ij - correct_class_score + 1;
                    if (NumericUtils<Scalar>::decidedly_greater(margin, (Scalar) 0)) {
                        Scalar out_grad_ij = Squared ? 2 * margin : 1;
                        total_out_grad += out_grad_ij;
                        out_grad(i,j) = out_grad_ij;
                    } else
                        out_grad(i,j) = 0;
                }
                out_grad(i,correct_class_ind) = -total_out_grad;
            }
            return TensorMap<Scalar,Root::DATA_RANK>(out_grad.data(), grad_dims);
        }
    };

    /**
    * A template class representing the cross Kullback-Leibler divergence loss function.
    *
    * \f$L_i = -\ln(\frac{-\hat{y_i}}{y_i + \epsilon} + \epsilon) y_i\f$
    */
    template<typename Scalar, std::size_t Rank, bool Sequential>
    class KullbackLeiblerLoss : public UniversalLoss<Scalar,Rank,Sequential> {
        typedef Loss<Scalar,Rank,Sequential> Root;
        typedef UniversalLoss<Scalar,Rank,Sequential> Base;
    public:
        /**
        * @param epsilon A small constant used to maintain numerical stability.
        */
        KullbackLeiblerLoss(Scalar epsilon = NumericUtils<Scalar>::EPSILON2) :
                epsilon(epsilon) { };
    protected:
        inline ColVector<Scalar> _function(typename Root::Data out, typename Root::Data obj) const {
            std::size_t rows = out.dimension(0);
            std::size_t cols = out.size() / rows;
            MatrixMap<Scalar> obj_mat(obj.data(), rows, cols);
            return -((MatrixMap<Scalar>(out.data(), rows, cols).array() /
                    (obj_mat.array() + epsilon) + epsilon).log() *
                    obj_mat.array()).matrix().rowwise().sum();
        }
        inline typename Root::Data _d_function(typename Root::Data out, typename Root::Data obj,
                const typename Base::RankwiseArray& grad_dims) const {
            return -obj / (out + epsilon);
        }
    private:
        Scalar epsilon;
    };

    /**
    * A class template representing the hinge loss function for multi-label objectives. True labels
    * are expected to have the value 1, while false labels are expected to correspond to the value -1.
    *
    * \f$L_i = \sum_j \max(0, 1 - {y_i}_j \hat{y_i}_j)\f$ or
    * \f$L_i = \sum_j \max(0, 1 - {y_i}_j \hat{y_i}_j)^2\f$
    */
    template<typename Scalar, std::size_t Rank, bool Sequential, bool Squared = false>
    class MultiLabelHingeLoss : public UniversalLoss<Scalar,Rank,Sequential> {
        typedef Loss<Scalar,Rank,Sequential> Root;
        typedef UniversalLoss<Scalar,Rank,Sequential> Base;
    protected:
        inline ColVector<Scalar> _function(typename Root::Data out, typename Root::Data obj) const {
            std::size_t rows = out.dimension(0);
            std::size_t cols = out.size() / rows;
            MatrixMap<Scalar> out_mat(out.data(), rows, cols);
            MatrixMap<Scalar> obj_mat(obj.data(), rows, cols);
            ColVector<Scalar> loss(rows);
            for (int i = 0; i < rows; ++i) {
                Scalar loss_i = 0;
                for (int j = 0; j < cols; ++j) {
                    Scalar obj_ij = obj_mat(i,j);
                    assert((NumericUtils<Scalar>::almost_equal(obj_ij, (Scalar) -1) ||
                            NumericUtils<Scalar>::almost_equal(obj_ij, (Scalar) 1)));
                    Scalar loss_ij = std::max((Scalar) 0, (Scalar) 1 - obj_ij * out_mat(i,j));
                    loss_i += Squared ? loss_ij * loss_ij : loss_ij;
                }
                loss(i) = loss_i;
            }
            return loss;
        }
        inline typename Root::Data _d_function(typename Root::Data out, typename Root::Data obj,
                const typename Base::RankwiseArray& grad_dims) const {
            std::size_t rows = out.dimension(0);
            std::size_t cols = out.size() / rows;
            MatrixMap<Scalar> out_mat(out.data(), rows, cols);
            MatrixMap<Scalar> obj_mat(obj.data(), rows, cols);
            Matrix<Scalar> out_grad(rows, cols);
            for (int i = 0; i < cols; ++i) {
                for (int j = 0; j < rows; ++j) {
                    Scalar obj_ji = obj_mat(j,i);
                    assert((NumericUtils<Scalar>::almost_equal(obj_ji, (Scalar) -1) ||
                            NumericUtils<Scalar>::almost_equal(obj_ji, (Scalar) 1)));
                    Scalar out_ji = out_mat(j,i);
                    Scalar margin = 1 - obj_ji * out_ji;
                    if (NumericUtils<Scalar>::decidedly_greater(margin, (Scalar) 0))
                        out_grad(j,i) = Squared ? 2 * out_ji - 2 * obj_ji : -obj_ji;
                    else
                        out_grad(j,i) = 0;
                }
            }
            return TensorMap<Scalar,Root::DATA_RANK>(out_grad.data(), grad_dims);
        }
    };


    /**
    * A class template representing the logarithmic loss function for multi-label objectives. True
    * labels are expected to have the value 1, while false labels are expected to correspond to the
    * value 0.
    *
    * \f$L_i = \sum_j {y_i}_j \ln(\hat{y_i}_j + \epsilon) + (1 - {y_i}_j) \ln(1 + \epsilon - \hat{y_i}_j)\f$
    */
    template<typename Scalar, std::size_t Rank, bool Sequential>
    class MultiLabelLogLoss : public UniversalLoss<Scalar,Rank,Sequential> {
        typedef Loss<Scalar,Rank,Sequential> Root;
        typedef UniversalLoss<Scalar,Rank,Sequential> Base;
    public:
        /**
        * @param epsilon A small constant used to maintain numerical stability.
        */
        MultiLabelLogLoss(Scalar epsilon = NumericUtils<Scalar>::EPSILON2) :
            epsilon(epsilon) { };
    protected:
        inline ColVector<Scalar> _function(typename Root::Data out, typename Root::Data obj) const {
            std::size_t rows = out.dimension(0);
            std::size_t cols = out.size() / rows;
            MatrixMap<Scalar> out_mat(out.data(), rows, cols);
            MatrixMap<Scalar> obj_mat(obj.data(), rows, cols);
            ColVector<Scalar> loss(rows);
            for (int i = 0; i < rows; ++i) {
                Scalar loss_i = 0;
                for (int j = 0; j < cols; ++j) {
                    Scalar obj_ij = obj_mat(i,j);
                    assert(NumericUtils<Scalar>::almost_equal(obj_ij, (Scalar) 0) ||
                            NumericUtils<Scalar>::almost_equal(obj_ij, (Scalar) 1));
                    Scalar out_ij = out_mat(i,j);
                    loss_i += obj_ij * log(out_ij + epsilon) + (1 - obj_ij) * log(1 + epsilon - out_ij);
                }
                loss(i) = loss_i;
            }
            return loss;
        }
        inline typename Root::Data _d_function(typename Root::Data out, typename Root::Data obj,
                const typename Base::RankwiseArray& grad_dims) const {
            std::size_t rows = out.dimension(0);
            std::size_t cols = out.size() / rows;
            MatrixMap<Scalar> out_mat(out.data(), rows, cols);
            MatrixMap<Scalar> obj_mat(obj.data(), rows, cols);
            Matrix<Scalar> out_grad(rows, cols);
            for (int i = 0; i < cols; ++i) {
                for (int j = 0; j < rows; ++j) {
                    Scalar obj_ji = obj_mat(j,i);
                    assert(NumericUtils<Scalar>::almost_equal(obj_ji, (Scalar) 0) ||
                            NumericUtils<Scalar>::almost_equal(obj_ji, (Scalar) 1));
                    Scalar out_ji = out_mat(j,i);
                    Scalar denominator = out_ji * (1 - out_ji);
                    if (out_ji == 0)
                        out_ji += epsilon;
                    out_grad(j,i) = (obj_ji - out_ji) / denominator;
                }
            }
            return TensorMap<Scalar,Root::DATA_RANK>(out_grad.data(), grad_dims);
        }
    private:
        Scalar epsilon;
    };    

    /**
    * An alias for a unique pointer to a loss function of arbitrary rank, scalar type and
    * sequentiality.
    */
    template<typename Scalar, std::size_t Rank, bool Sequential>
    using LossSharedPtr = std::shared_ptr<Loss<Scalar,Rank,Sequential>>;

    /**
    * A wrapper class template for negating losses and thus allowing for their maximization
    * via the standard optimization methods.
    */
    template<typename Scalar, std::size_t Rank, bool Sequential>
    class NegatedLoss : public Loss<Scalar,Rank,Sequential> {
        typedef Loss<Scalar,Rank,Sequential> Base;
    public:
        /**
        * @param loss A shared pointer to the loss instance to negate.
        */
        NegatedLoss(LossSharedPtr<Scalar,Rank,Sequential> loss) :
                loss(loss) {
            assert(loss);
        }
        inline ColVector<Scalar> function(typename Base::Data out, typename Base::Data obj) const {
            return -(loss->function(std::move(out), std::move(obj)));
        }
        inline typename Base::Data d_function(typename Base::Data out, typename Base::Data obj) const {
            return -(loss->d_function(std::move(out), std::move(obj)));
        }
    private:
        LossSharedPtr<Scalar,Rank,Sequential> loss;
    };    

    /**
    * A loss function template that applies the softmax function to its input before calculating the cross
    * entropy loss. This allows for increased numerical stability and faster computation.
    *
    * \f$L_i = -\ln(\text{softmax}(\hat{y_i}) + \epsilon) y_i\f$
    */
    template<typename Scalar, std::size_t Rank, bool Sequential>
    class SoftmaxCrossEntropyLoss : public UniversalLoss<Scalar,Rank,Sequential> {
        typedef Loss<Scalar,Rank,Sequential> Root;
        typedef UniversalLoss<Scalar,Rank,Sequential> Base;
    public:
        /**
        * @param epsilon A small constant used to maintain numerical stability.
        */
        SoftmaxCrossEntropyLoss(Scalar epsilon = NumericUtils<Scalar>::EPSILON2) :
                epsilon(epsilon) { };
    protected:
        inline ColVector<Scalar> _function(typename Root::Data out, typename Root::Data obj) const {
            std::size_t rows = out.dimension(0);
            std::size_t cols = out.size() / rows;
            MatrixMap<Scalar> out_mat(out.data(), rows, cols);
            Matrix<Scalar> out_exp = (out_mat.array().colwise() - out_mat.array().rowwise().maxCoeff()).exp();
            return -((out_exp.array().colwise() / (out_exp.array().rowwise().sum() + epsilon)).log() *
                    MatrixMap<Scalar>(obj.data(), rows, cols).array()).matrix().rowwise().sum();
        }
        inline typename Root::Data _d_function(typename Root::Data out, typename Root::Data obj,
                const typename Base::RankwiseArray& grad_dims) const {
            std::size_t rows = out.dimension(0);
            std::size_t cols = out.size() / rows;
            MatrixMap<Scalar> out_mat(out.data(), rows, cols);
            Matrix<Scalar> out_exp = (out_mat.array().colwise() - out_mat.array().rowwise().maxCoeff()).exp();
            Matrix<Scalar> grads = (out_exp.array().colwise() / (out_exp.array().rowwise().sum() + epsilon)) -
                    MatrixMap<Scalar>(obj.data(), rows, cols).array();
            return TensorMap<Scalar,Root::DATA_RANK>(grads.data(), grad_dims);
        }
    private:
        Scalar epsilon;
    };

    /**
    * A template class representing the squared error (L2) loss function.
    *
    * \f$L_i = (\hat{y_i} - y_i)^2\f$
    */
    template<typename Scalar, std::size_t Rank, bool Sequential>
    class SquaredLoss : public UniversalLoss<Scalar,Rank,Sequential> {
        typedef Loss<Scalar,Rank,Sequential> Root;
        typedef UniversalLoss<Scalar,Rank,Sequential> Base;
    protected:
        inline ColVector<Scalar> _function(typename Root::Data out, typename Root::Data obj) const {
            std::size_t rows = out.dimension(0);
            std::size_t cols = out.size() / rows;
            return (MatrixMap<Scalar>(out.data(), rows, cols) - MatrixMap<Scalar>(obj.data(), rows, cols))
                    .array().square().rowwise().sum();
        }
        inline typename Root::Data _d_function(typename Root::Data out, typename Root::Data obj,
                const typename Base::RankwiseArray& grad_dims) const {
            return 2 * (out - obj);
        }
    };

    /**
    * An alias for a unique pointer to a neural network of arbitrary scalar type, rank,
    * and sequentiality.
    */
    template<typename Scalar, std::size_t Rank, bool Sequential>
    using NeuralNetPtr = std::unique_ptr<NeuralNetwork<Scalar,Rank,Sequential>>;

    /**
    * A class template for composite neural networks consisting of one or more neural
    * network modules.
    */
    template<typename Scalar, std::size_t Rank, bool Sequential, typename Module>
    class CompositeNeuralNetwork : public NeuralNetwork<Scalar,Rank,Sequential> {
    public:
        /**
        * @return A vector of pointers pointing to the sub-modules of the composite
        * network instance. The ownership of the modules is not transferred to the
        * caller of the method.
        */
        virtual std::vector<Module*> get_modules() = 0;
    };

    /**
    * An abstract class template for unidirectional recurrent neural networks.
    */
    template<typename Scalar, std::size_t Rank>
    class UnidirectionalNeuralNetwork : public NeuralNetwork<Scalar,Rank,true> {
        typedef NeuralNetwork<Scalar,Rank,true> Base;
    public:
        virtual ~UnidirectionalNeuralNetwork() = default;
        /**
        * @return Whether the direction along the time-step rank in which the network processes
        * its inputs is reversed.
        */
        virtual bool is_reversed() const=0;
        /**
        * Flips the direction along the time-step rank in which the network processes its inputs
        * is reversed.
        */
        virtual void reverse()=0;
        /**
        * Reverses a tensor along its time axis.
        *
        * @param tensor The tensor to reverse.
        */
        inline static void reverse_along_time_axis(typename Base::Data& tensor) {
            std::array<bool,Base::DATA_RANK> reverse;
            reverse.fill(false);
            reverse[1] = true;
            tensor = tensor.reverse(reverse);
        }
    };


    /**
 * An enumeration type for the different ways the outputs of sub-modules of neural networks
 * may be merged.
 */
enum BidirectionalOutputMergeType { BIDIRECTIONAL_CONCAT_LO_RANK, BIDIRECTIONAL_CONCAT_HI_RANK,
		BIDIRECTIONAL_SUM, BIDIRECTIONAL_MUL };

/**
 * An alias for unidirectional recurrent neural network of arbitrary scalar type and rank.
 */
template<typename Scalar, std::size_t Rank>
using UnidirNeuralNetPtr = std::unique_ptr<UnidirectionalNeuralNetwork<Scalar,Rank>>;

    /**
    * A class template for a bidirectional neural network that takes a unidirectional recurrent
    * network, clones it, reverses the clone's processing direction, and uses the two networks
    * as its parallel sub-modules. The outputs of the two sub-networks can be merged by summation
    * or concatenation either along the lowest (the 3rd after the sample and time-step ranks) or
    * highest rank.
    *
    * \see https://pdfs.semanticscholar.org/4b80/89bc9b49f84de43acc2eb8900035f7d492b2.pdf
    */
    template<typename Scalar, std::size_t Rank, BidirectionalOutputMergeType MergeType = BIDIRECTIONAL_CONCAT_LO_RANK>
    class BidirectionalNeuralNetwork :
            public CompositeNeuralNetwork<Scalar,Rank,true,UnidirectionalNeuralNetwork<Scalar,Rank>> {
        typedef NeuralNetwork<Scalar,Rank,true> Base;
        typedef BidirectionalNeuralNetwork<Scalar,Rank,MergeType> Self;
        typedef UnidirNeuralNetPtr<Scalar,Rank> UnidirNet;
        typedef std::array<std::size_t,Base::DATA_RANK> RankwiseArray;
        static_assert(MergeType >= BIDIRECTIONAL_CONCAT_LO_RANK && MergeType <= BIDIRECTIONAL_MUL,
                "illegal merge type value");
        static constexpr std::size_t CONCAT_RANK = MergeType == BIDIRECTIONAL_CONCAT_HI_RANK ? Rank - 1 : 0;
        static constexpr std::size_t CONCAT_BATCH_RANK = CONCAT_RANK + 2;
    public:
        /**
        * @param network A unique pointer to a unidirectional recurrent neural network that,
        * along with its reversed clone, will constitute the bidirectional network.
        * @param foremost Whether the network is to function as a foremost network.
        */
        inline BidirectionalNeuralNetwork(UnidirNet&& network, bool foremost = true) :
                net(std::move(network)),
                foremost(foremost) {
            assert(this->net);
            net_rev = UnidirNet(static_cast<UnidirectionalNeuralNetwork<Scalar,Rank>*>(this->net->clone()));
            net_rev->reverse();
            input_dims = this->net->get_input_dims();
            output_dims = this->net->get_output_dims();
            if (MergeType == BIDIRECTIONAL_CONCAT_LO_RANK || MergeType == BIDIRECTIONAL_CONCAT_HI_RANK)
                output_dims(+CONCAT_RANK) *= 2;
        }
        inline BidirectionalNeuralNetwork(const Self& network) :
                net(UnidirNet(static_cast<UnidirectionalNeuralNetwork<Scalar,Rank>*>(network.net->clone()))),
                net_rev(UnidirNet(static_cast<UnidirectionalNeuralNetwork<Scalar,Rank>*>(network.net_rev->clone()))),
                foremost(network.foremost),
                input_dims(network.input_dims),
                output_dims(network.output_dims),
                output(network.output),
                output_rev(network.output_rev) { }
        inline BidirectionalNeuralNetwork(Self&& network) {
            swap(*this, network);
        }
        ~BidirectionalNeuralNetwork() = default;
        inline Self& operator=(Self network) {
            swap(*this, network);
            return *this;
        }
        inline Base* clone() const {
            return new BidirectionalNeuralNetwork(*this);
        }
        inline const typename Base::Dims& get_input_dims() const {
            return input_dims;
        }
        inline const typename Base::Dims& get_output_dims() const {
            return output_dims;
        }
        inline std::vector<const Layer<Scalar,Rank>*> get_layers() const {
            std::vector<const Layer<Scalar,Rank>*> layer_ptrs;
            populate_layer_vector<const Layer<Scalar,Rank>*>(layer_ptrs);
            return layer_ptrs;
        }
        inline std::vector<Layer<Scalar,Rank>*> get_layers() {
            std::vector<Layer<Scalar,Rank>*> layer_ptrs;
            populate_layer_vector<Layer<Scalar,Rank>*>(layer_ptrs);
            return layer_ptrs;
        }
        inline std::vector<UnidirectionalNeuralNetwork<Scalar,Rank>*> get_modules() {
            std::vector<UnidirectionalNeuralNetwork<Scalar,Rank>*> modules;
            modules.push_back(net.get());
            modules.push_back(net_rev.get());
            return modules;
        }
        inline bool is_foremost() const {
            return foremost;
        }
        inline void set_foremost(bool foremost) {
            net->set_foremost(foremost);
            net_rev->set_foremost(foremost);
            this->foremost = foremost;
        }
        inline void empty_caches() {
            net->empty_caches();
            net_rev->empty_caches();
            output = typename Base::Data();
            output_rev = typename Base::Data();
        }
        inline typename Base::Data propagate(typename Base::Data input, bool training) {
            assert(input_dims == (Dimensions<std::size_t,Base::DATA_RANK>(input.dimensions()).template demote<2>()));
            pthread_attr_t attr;
            pthread_t helper_thread;
            int pthread_state;
            pthread_state = pthread_attr_init(&attr);
            assert(!pthread_state);
            pthread_state = pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
            assert(!pthread_state);
            PropArgs args;
            args.obj = this;
            args.training = training;
            args.in = &input;
            pthread_state = pthread_create(&helper_thread, &attr, propagate, &args);
            assert(!pthread_state);
            typename Base::Data forward_out = net->propagate(input, training);
            pthread_state = pthread_join(helper_thread, nullptr);
            assert(!pthread_state);
            pthread_state = pthread_attr_destroy(&attr);
            assert(!pthread_state);
            assert(forward_out.dimension(1) == args.out.dimension(1));
            if (MergeType == BIDIRECTIONAL_SUM)
                return forward_out + args.out;
            else if (MergeType == BIDIRECTIONAL_MUL) {
                output = std::move(forward_out);
                output_rev = std::move(args.out);
                return output * output_rev;
            } else
                return forward_out.concatenate(args.out, +CONCAT_BATCH_RANK);
        }
        inline typename Base::Data backpropagate(typename Base::Data out_grad) {
            Dimensions<std::size_t,Base::DATA_RANK> dims(out_grad.dimensions());
            assert(output_dims == dims.template demote<2>());
            pthread_attr_t attr;
            pthread_t helper_thread;
            int pthread_state;
            pthread_state = pthread_attr_init(&attr);
            assert(!pthread_state);
            pthread_state = pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
            assert(!pthread_state);
            BackpropArgs args;
            args.obj = this;
            typename Base::Data forward_prev_out_grad;
            if (MergeType == BIDIRECTIONAL_SUM) {
                args.out_grad = &out_grad;
                pthread_state = pthread_create(&helper_thread, &attr, backpropagate, &args);
                assert(!pthread_state);
                forward_prev_out_grad = net->backpropagate(out_grad);
                pthread_state = pthread_join(helper_thread, nullptr);
                assert(!pthread_state);
                out_grad = typename Base::Data();
            } else if (MergeType == BIDIRECTIONAL_MUL) {
                typename Base::Data out_grad_rev = output * out_grad;
                args.out_grad = &out_grad_rev;
                pthread_state = pthread_create(&helper_thread, &attr, backpropagate, &args);
                assert(!pthread_state);
                out_grad *= output_rev;
                forward_prev_out_grad = net->backpropagate(std::move(out_grad));
                pthread_state = pthread_join(helper_thread, nullptr);
                assert(!pthread_state);
            } else {
                RankwiseArray offsets;
                RankwiseArray extents = dims;
                offsets.fill(0);
                extents[+CONCAT_BATCH_RANK] /= 2;
                offsets[+CONCAT_BATCH_RANK] += extents[+CONCAT_BATCH_RANK];
                typename Base::Data backward_slice = out_grad.slice(offsets, extents);
                args.out_grad = &backward_slice;
                pthread_state = pthread_create(&helper_thread, &attr, backpropagate, &args);
                assert(!pthread_state);
                offsets[+CONCAT_BATCH_RANK] -= extents[+CONCAT_BATCH_RANK];
                typename Base::Data forward_slice = out_grad.slice(offsets, extents);
                out_grad = typename Base::Data();
                forward_prev_out_grad = net->backpropagate(std::move(forward_slice));
                // Make sure that backward_slice does not go out of scope before the thread terminates.
                pthread_state = pthread_join(helper_thread, nullptr);
                assert(!pthread_state);
            }
            pthread_state = pthread_attr_destroy(&attr);
            assert(!pthread_state);
            return forward_prev_out_grad + args.prev_out_grad;
        }
        inline friend void swap(Self& network1, Self& network2) {
            using std::swap;
            swap(network1.net, network2.net);
            swap(network1.net_rev, network2.net_rev);
            swap(network1.foremost, network2.foremost);
            swap(network1.input_dims, network2.input_dims);
            swap(network1.output_dims, network2.output_dims);
            swap(network1.output, network2.output);
            swap(network1.output_rev, network2.output_rev);
        }
    private:
        /**
        * The propagation function executed in a different thread for each lane of a
        * parallel network.
        *
        * @param args_ptr The propagation argument struct containing all necessary
        * information.
        */
        inline static void* propagate(void* args_ptr) {
            PropArgs& args = *((PropArgs*) args_ptr);
            args.out = args.obj->net_rev->propagate(*args.in, args.training);
            return nullptr;
        }
        /**
        * The back-propagation function executed in a different thread for each lane of a
        * parallel network.
        *
        * @param args_ptr The back-propagation argument struct containing all necessary
        * information.
        */
        inline static void* backpropagate(void* args_ptr) {
            BackpropArgs& args = *((BackpropArgs*) args_ptr);
            args.prev_out_grad = args.obj->net_rev->backpropagate(*args.out_grad);
            return nullptr;
        }
        template<typename _LayerPtr>
        inline void populate_layer_vector(std::vector<_LayerPtr>& layer_ptrs) const {
            std::vector<Layer<Scalar,Rank>*> net_layer_ptrs = net->get_layers();
            for (std::size_t i = 0; i < net_layer_ptrs.size(); ++i)
                layer_ptrs.push_back(net_layer_ptrs[i]);
            std::vector<Layer<Scalar,Rank>*> net_rev_layer_ptrs = net_rev->get_layers();
            for (std::size_t i = 0; i < net_rev_layer_ptrs.size(); ++i)
                layer_ptrs.push_back(net_rev_layer_ptrs[i]);
        }
        UnidirNet net, net_rev;
        bool foremost;
        typename Base::Dims input_dims, output_dims;
        typename Base::Data output, output_rev;
        /**
        * A struct containing the data required for propagation.
        */
        struct PropArgs {
            Self* obj;
            bool training;
            typename Base::Data* in;
            typename Base::Data out;
        };
        /**
        * A struct containing the data require for back-propagation.
        */
        struct BackpropArgs {
            Self* obj;
            typename Base::Data* out_grad;
            typename Base::Data prev_out_grad;
        };
    };

    /**
    * An enumeration type for the different ways the input of a layer in a dense network may be concatenated
    * to its output.
    */
    enum DenseConcatType { DENSE_LOWEST_RANK, DENSE_HIGHEST_RANK };

    /**
    * A class template for DenseNet architectures. These networks consist of sub-modules that are all
    * 'connected' to each other as in the input of each module is concatenated to its output and then
    * propagated to the next module as its input. The input is concatenated to the output either along
    * its lowest or highest rank.
    *
    * \see https://arxiv.org/abs/1608.06993
    */
    template<typename Scalar, std::size_t Rank, DenseConcatType ConcatType = DENSE_HIGHEST_RANK>
    class DenseNeuralNetwork :
            public CompositeNeuralNetwork<Scalar,Rank,false,NeuralNetwork<Scalar,Rank,false>> {
        typedef NeuralNetwork<Scalar,Rank,false> Base;
        typedef NeuralNetPtr<Scalar,Rank,false> Module;
        typedef DenseNeuralNetwork<Scalar,Rank,ConcatType> Self;
        typedef std::array<std::size_t,Base::DATA_RANK> RankwiseArray;
        static_assert(ConcatType >= DENSE_LOWEST_RANK && ConcatType <= DENSE_HIGHEST_RANK, "illegal merge type value");
        static constexpr std::size_t CONCAT_RANK = ConcatType == DENSE_HIGHEST_RANK ? Rank - 1 : 0;
        static constexpr std::size_t CONCAT_BATCH_RANK = CONCAT_RANK + 1;
    public:
        /**
        * @param modules A vector of dense modules.
        * @param foremost Whether the network is to function as a foremost network.
        */
        inline DenseNeuralNetwork(std::vector<Module>&& modules, bool foremost = true) :
                modules(std::move(modules)),
                foremost(foremost) {
            assert(this->modules.size() > 0 && "modules must contain at least 1 element");
            Base& first_module = *this->modules[0];
            input_dims = first_module.get_input_dims();
            typename Base::Dims output_dims = first_module.get_output_dims();
            output_dims(+CONCAT_RANK) += input_dims(+CONCAT_RANK);
            if (ConcatType == DENSE_LOWEST_RANK) {
                for (std::size_t i = Rank - 1; i > +CONCAT_RANK; --i)
                    assert(input_dims(i) == output_dims(i));
            } else {
                for (std::size_t i = 0; i < +CONCAT_RANK; ++i)
                    assert(input_dims(i) == output_dims(i));
            }
            for (std::size_t i = 1; i < this->modules.size(); ++i) {
                Base& module = *this->modules[i];
                const typename Base::Dims& module_input_dims = module.get_input_dims();
                assert(module_input_dims == output_dims && "incompatible module dimensions");
                output_dims(+CONCAT_RANK) += module.get_output_dims()(+CONCAT_RANK);
                module.set_foremost(false);
            }
            this->output_dims = output_dims;
            first_module.set_foremost(foremost);
        }
        /**
        * @param module A single dense module.
        * @param foremost Whether the network is to function as a foremost network.
        */
        inline DenseNeuralNetwork(Module&& module, bool foremost = true) :
                DenseNeuralNetwork(create_vector(std::move(module), foremost)) { }
        inline DenseNeuralNetwork(const Self& network) :
                modules(network.modules.size()),
                foremost(network.foremost),
                input_dims(network.input_dims),
                output_dims(network.output_dims) {
            for (std::size_t i = 0; i < modules.size(); ++i)
                modules[i] = Module(network.modules[i]->clone());
        }
        inline DenseNeuralNetwork(Self&& network) {
            swap(*this, network);
        }
        ~DenseNeuralNetwork() = default;
        inline Self& operator=(Self network) {
            swap(*this, network);
            return *this;
        }
        inline Base* clone() const {
            return new DenseNeuralNetwork(*this);
        }
        inline const typename Base::Dims& get_input_dims() const {
            return input_dims;
        }
        inline const typename Base::Dims& get_output_dims() const {
            return output_dims;
        }
        inline std::vector<const Layer<Scalar,Rank>*> get_layers() const {
            std::vector<const Layer<Scalar,Rank>*> layer_ptrs;
            populate_layer_vector<const Layer<Scalar,Rank>*>(layer_ptrs);
            return layer_ptrs;
        }
        inline std::vector<Layer<Scalar,Rank>*> get_layers() {
            std::vector<Layer<Scalar,Rank>*> layer_ptrs;
            populate_layer_vector<Layer<Scalar,Rank>*>(layer_ptrs);
            return layer_ptrs;
        }
        inline std::vector<Base*> get_modules() {
            std::vector<Base*> modules;
            for (std::size_t i = 0; i < this->modules.size(); ++i)
                modules.push_back(this->modules[i].get());
            return modules;
        }
        inline bool is_foremost() const {
            return foremost;
        }
        inline void set_foremost(bool foremost) {
            modules[0]->set_foremost(foremost);
            this->foremost = foremost;
        }
        inline void empty_caches() {
            for (std::size_t i = 0; i < modules.size(); ++i)
                modules[i]->empty_caches();
        }
        inline typename Base::Data propagate(typename Base::Data input, bool training) {
            assert(input_dims == (Dimensions<std::size_t,Base::DATA_RANK>(input.dimensions()).template demote<>()));
            for (std::size_t i = 0; i < modules.size(); ++i) {
                typename Base::Data concat = input.concatenate(modules[i]->propagate(input, training), +CONCAT_BATCH_RANK);
                input = std::move(concat);
            }
            return input;
        }
        inline typename Base::Data backpropagate(typename Base::Data out_grad) {
            assert(output_dims == (Dimensions<std::size_t,Base::DATA_RANK>(out_grad.dimensions()).template demote<>()));
            RankwiseArray offsets;
            RankwiseArray extents = input_dims.template promote<>();
            offsets.fill(0);
            extents[0] = out_grad.dimension(0);
            for (int i = modules.size() - 1; i >= 0; --i) {
                Base& module = *modules[i];
                int layer_input_concat_rank_dim = module.get_input_dims()(+CONCAT_RANK);
                int layer_output_concat_rank_dim = module.get_output_dims()(+CONCAT_RANK);
                offsets[+CONCAT_BATCH_RANK] = layer_input_concat_rank_dim;
                extents[+CONCAT_BATCH_RANK] = layer_output_concat_rank_dim;
                typename Base::Data out_grad_i = out_grad.slice(offsets, extents);
                offsets[+CONCAT_BATCH_RANK] = 0;
                extents[+CONCAT_BATCH_RANK] = layer_input_concat_rank_dim;
                if (foremost && i == 0)
                    module.backpropagate(std::move(out_grad_i));
                else
                    out_grad = typename Base::Data(out_grad.slice(offsets, extents) +
                            module.backpropagate(std::move(out_grad_i)));
            }
            return out_grad;
        }
        inline friend void swap(Self& network1, Self& network2) {
            using std::swap;
            swap(network1.modules, network2.modules);
            swap(network1.foremost, network2.foremost);
            swap(network1.input_dims, network2.input_dims);
            swap(network1.output_dims, network2.output_dims);
        }
    private:
        inline static std::vector<Module> create_vector(Module&& module) {
            std::vector<Module> vec(1);
            vec[0] = std::move(module);
            return vec;
        }
        template<typename _LayerPtr>
        inline void populate_layer_vector(std::vector<_LayerPtr>& layer_ptrs) const {
            for (std::size_t i = 0; i < modules.size(); ++i) {
                std::vector<Layer<Scalar,Rank>*> internal_layer_ptrs = modules[i]->get_layers();
                for (std::size_t j = 0; j < internal_layer_ptrs.size(); ++j)
                    layer_ptrs.push_back(internal_layer_ptrs[j]);
            }
        }
        std::vector<Module> modules;
        bool foremost;
        typename Base::Dims input_dims, output_dims;
    };

    /**
    * An alias for a unique pointer to a layer of arbitrary rank and scalar type.
    */
    template<typename Scalar, std::size_t Rank>
    using LayerPtr = std::unique_ptr<Layer<Scalar,Rank>>;

    /**
    * A class template representing a simple feed-forward neural network.
    */
    template<typename Scalar, std::size_t Rank>
    class FeedforwardNeuralNetwork : public NeuralNetwork<Scalar,Rank,false> {
        typedef NeuralNetwork<Scalar,Rank,false> Base;
        typedef FeedforwardNeuralNetwork<Scalar,Rank> Self;
    public:
        /**
        * @param layers A vector of unique smart pointers to the layers that constitute the neural network.
        * @param foremost Whether the network directly receives its input. If it is set to false, back-propagation
        * returns an empty tensor.
        */
        inline FeedforwardNeuralNetwork(std::vector<LayerPtr<Scalar,Rank>>&& layers, bool foremost = true) :
                layers(std::move(layers)),
                foremost(foremost) {
            assert(this->layers.size() > 0 && "layers must contain at least 1 element");
            assert(this->layers[0] != nullptr);
            Layer<Scalar,Rank>& first_layer = *this->layers[0];
            input_dims = first_layer.get_input_dims();
            output_dims = this->layers[this->layers.size() - 1]->get_output_dims();
            typename Base::Dims prev_dims = first_layer.get_output_dims();
            for (std::size_t i = 1; i < this->layers.size(); ++i) {
                assert(this->layers[i] != nullptr && "layers contains null pointers");
                assert(prev_dims == this->layers[i]->get_input_dims() && "incompatible layer dimensions");
                prev_dims = this->layers[i]->get_output_dims();
            }
            first_layer.set_input_layer(foremost);
        }
        /**
        * @param layer A unique pointer to the single layer of the network.
        * @param foremost Whether the network is to function as a foremost network.
        */
        inline FeedforwardNeuralNetwork(LayerPtr<Scalar,Rank>&& layer, bool foremost = true) :
                FeedforwardNeuralNetwork(create_vector(std::move(layer)), foremost) { }
        // Copy constructor.
        inline FeedforwardNeuralNetwork(const Self& network) :
                layers(network.layers.size()),
                foremost(network.foremost),
                input_dims(network.input_dims),
                output_dims(network.output_dims) {
            for (std::size_t i = 0; i < layers.size(); ++i)
                layers[i] = LayerPtr<Scalar,Rank>(network.layers[i]->clone());
        }
        // Move constructor.
        inline FeedforwardNeuralNetwork(Self&& network) {
            swap(*this, network);
        }
        // The smart pointers take care of deleting the layers.
        ~FeedforwardNeuralNetwork() = default;
        /* The assignment uses the move or copy constructor to pass the parameter
        * based on whether it is an rvalue or an lvalue. */
        inline Self& operator=(Self network) {
            swap(*this, network);
            return *this;
        }
        inline Base* clone() const {
            return new FeedforwardNeuralNetwork(*this);
        }
        inline const typename Base::Dims& get_input_dims() const {
            return input_dims;
        }
        inline const typename Base::Dims& get_output_dims() const {
            return output_dims;
        }
        inline std::vector<const Layer<Scalar,Rank>*> get_layers() const {
            std::vector<const Layer<Scalar,Rank>*> layer_ptrs(layers.size());
            populate_layer_vector<const Layer<Scalar,Rank>*>(layer_ptrs);
            return layer_ptrs;
        }
        inline std::vector<Layer<Scalar,Rank>*> get_layers() {
            std::vector<Layer<Scalar,Rank>*> layer_ptrs(layers.size());
            populate_layer_vector<Layer<Scalar,Rank>*>(layer_ptrs);
            return layer_ptrs;
        }
        inline bool is_foremost() const {
            return foremost;
        }
        inline void set_foremost(bool foremost) {
            layers[0]->set_input_layer(foremost);
            this->foremost = foremost;
        }
        inline virtual void empty_caches() {
            std::vector<Layer<Scalar,Rank>*> layers = get_layers();
            for (std::size_t i = 0; i < layers.size(); ++i)
                layers[i]->empty_cache();
        }
        inline typename Base::Data propagate(typename Base::Data input, bool training) {
            assert(input_dims == (Dimensions<std::size_t,Base::DATA_RANK>(input.dimensions()).template demote<>()));
            for (std::size_t i = 0; i < layers.size(); ++i) {
                Layer<Scalar,Rank>& layer = *layers[i];
                input = layer.pass_forward(std::move(input), training);
            }
            return input;
        }
        inline typename Base::Data backpropagate(typename Base::Data out_grad) {
            assert(output_dims == (Dimensions<std::size_t,Base::DATA_RANK>(out_grad.dimensions()).template demote<>()));
            for (int i = layers.size() - 1; i >= 0; --i) {
                Layer<Scalar,Rank>& layer = *layers[i];
                out_grad = layer.pass_back(std::move(out_grad));
            }
            return out_grad;
        }
        // For the copy-and-swap idiom.
        inline friend void swap(Self& network1, Self& network2) {
            using std::swap;
            swap(network1.layers, network2.layers);
            swap(network1.foremost, network2.foremost);
            swap(network1.input_dims, network2.input_dims);
            swap(network1.output_dims, network2.output_dims);
        }
    private:
        inline static std::vector<LayerPtr<Scalar,Rank>> create_vector(LayerPtr<Scalar,Rank>&& layer) {
            std::vector<LayerPtr<Scalar,Rank>> vec(1);
            vec[0] = std::move(layer);
            return vec;
        }
        template<typename _LayerPtr>
        inline void populate_layer_vector(std::vector<_LayerPtr>& layer_ptrs) const {
            for (std::size_t i = 0; i < layers.size(); ++i)
                layer_ptrs[i] = layers[i].get();
        }
        std::vector<LayerPtr<Scalar,Rank>> layers;
        bool foremost;
        typename Base::Dims input_dims, output_dims;
    };

    /**
    * An alias for a unique pointer to a kernel layer of arbitrary rank and scalar type.
    */
    template<typename Scalar, std::size_t Rank>
    using KernelPtr = std::unique_ptr<KernelLayer<Scalar,Rank>>;

    /**
    * An alias for a unique pointer to an activation layer of arbitrary rank and scalar type.
    */
    template<typename Scalar, std::size_t Rank>
    using ActivationPtr = std::unique_ptr<ActivationLayer<Scalar,Rank>>;

    /**
    * A class template representing a long-short term memory (LSTM) recurrent neural network. The network
    * can use multiplicative integration to combine its linearly transformed inputs and its linearly
    * transformed hidden outputs. A stateful network retains its hidden state across sequences as long as
    * the batch size is constant.
    *
    * \see http://www.bioinf.jku.at/publications/older/2604.pdf
    */
    template<typename Scalar, std::size_t Rank, bool MulInt = false, bool Stateful = false>
    class LSTMNeuralNetwork : public UnidirectionalNeuralNetwork<Scalar,Rank> {
        typedef NeuralNetwork<Scalar,Rank,true> Root;
        typedef UnidirectionalNeuralNetwork<Scalar,Rank> Base;
        typedef LSTMNeuralNetwork<Scalar,Rank,MulInt,Stateful> Self;
        typedef std::array<std::size_t,Root::DATA_RANK> RankwiseArray;
        typedef std::function<std::pair<std::size_t,std::size_t>(std::size_t)> OutputSeqSizeFunc;
        typedef Tensor<Scalar,Rank + 1> TimeStepData;
    public:
        /**
        * @param input_forget_kernel The forget kernel to apply to the input of the network.
        * @param output_forget_kernel The forget kernel to apply to the hidden output of the network
        * at the previous time step.
        * @param input_write_kernel The write kernel to apply to the input of the network.
        * @param output_write_kernel The write kernel to apply to the hidden output of the network
        * at the previous time step.
        * @param input_candidate_kernel The candidate kernel to apply to the input of the network.
        * @param output_candidate_kernel The candidate kernel to apply to the hidden output of the
        * network at the previous time step.
        * @param input_read_kernel The read kernel to apply to the input of the network.
        * @param output_read_kernel The read kernel to apply to the hidden output of the network
        * at the previous time step.
        * @param forget_act The activation layer of the forget gate. Usually a sigmoid activation
        * function.
        * @param write_act The activation layer of the filter of the write gate. Usually a sigmoid
        * activation function.
        * @param candidate_act The activation layer of the candidates of the write gate. Usually
        * a hyperbolic tangent activation function.
        * @param state_act The activation layer of the state at the read gate. Usually a hyperbolic
        * tangent activation function.
        * @param read_act The activation layer of the read filter. Usually a sigmoid activation
        * function.
        * @param output_seq_size_func A function parameterized by the input sequence length that
        * determines the output sequence delay and length. The output of the function is a pair of unsigned
        * integers where the first element is the sequence length and the second element is the sequence
        * delay.
        * @param reversed Whether the network is to reverse its inputs along the time-step rank.
        * @param foremost Whether the network is to function as a foremost network.
        */
        inline LSTMNeuralNetwork(KernelPtr<Scalar,Rank>&& input_forget_kernel,
                KernelPtr<Scalar,Rank>&& output_forget_kernel, KernelPtr<Scalar,Rank>&& input_write_kernel,
                KernelPtr<Scalar,Rank>&& output_write_kernel, KernelPtr<Scalar,Rank>&& input_candidate_kernel,
                KernelPtr<Scalar,Rank>&& output_candidate_kernel, KernelPtr<Scalar,Rank>&& input_read_kernel,
                KernelPtr<Scalar,Rank>&& output_read_kernel, ActivationPtr<Scalar,Rank>&& forget_act,
                ActivationPtr<Scalar,Rank>&& write_act, ActivationPtr<Scalar,Rank>&& candidate_act,
                ActivationPtr<Scalar,Rank>&& state_act, ActivationPtr<Scalar,Rank>&& read_act,
                OutputSeqSizeFunc output_seq_size_func, bool reversed = false, bool foremost = true) :
                    main_cell(),
                    output_seq_size_func(output_seq_size_func),
                    reversed(reversed),
                    foremost(foremost),
                    cells(0),
                    batch_size(-1),
                    input_seq_length(-1),
                    output_seq_length(-1),
                    output_seq_delay(-1) {
            assert(output_forget_kernel && input_forget_kernel && output_write_kernel && input_write_kernel &&
                    output_candidate_kernel && input_candidate_kernel && output_read_kernel && input_read_kernel &&
                    forget_act && write_act && candidate_act && state_act && read_act);
            typename Root::Dims in_forget_kernel_input_dims = input_forget_kernel->get_input_dims();
            typename Root::Dims out_forget_kernel_input_dims = output_forget_kernel->get_input_dims();
            assert(out_forget_kernel_input_dims == input_forget_kernel->get_output_dims() &&
                    in_forget_kernel_input_dims == input_write_kernel->get_input_dims() &&
                    in_forget_kernel_input_dims == input_candidate_kernel->get_input_dims() &&
                    in_forget_kernel_input_dims == input_write_kernel->get_input_dims() &&
                    out_forget_kernel_input_dims == output_write_kernel->get_input_dims() &&
                    out_forget_kernel_input_dims == output_candidate_kernel->get_input_dims() &&
                    out_forget_kernel_input_dims == output_write_kernel->get_input_dims() &&
                    out_forget_kernel_input_dims == output_forget_kernel->get_output_dims() &&
                    out_forget_kernel_input_dims == input_write_kernel->get_output_dims() &&
                    out_forget_kernel_input_dims == output_write_kernel->get_output_dims() &&
                    out_forget_kernel_input_dims == input_candidate_kernel->get_output_dims() &&
                    out_forget_kernel_input_dims == output_candidate_kernel->get_output_dims() &&
                    out_forget_kernel_input_dims == input_read_kernel->get_output_dims() &&
                    out_forget_kernel_input_dims == output_read_kernel->get_output_dims() &&
                    out_forget_kernel_input_dims == forget_act->get_input_dims() &&
                    out_forget_kernel_input_dims == write_act->get_input_dims() &&
                    out_forget_kernel_input_dims == candidate_act->get_input_dims() &&
                    out_forget_kernel_input_dims == state_act->get_input_dims() &&
                    out_forget_kernel_input_dims == read_act->get_input_dims());
            main_cell.input_forget_kernel = std::move(input_forget_kernel);
            main_cell.output_forget_kernel = std::move(output_forget_kernel);
            main_cell.input_write_kernel = std::move(input_write_kernel);
            main_cell.output_write_kernel = std::move(output_write_kernel);
            main_cell.input_candidate_kernel = std::move(input_candidate_kernel);
            main_cell.output_candidate_kernel = std::move(output_candidate_kernel);
            main_cell.input_read_kernel = std::move(input_read_kernel);
            main_cell.output_read_kernel = std::move(output_read_kernel);
            main_cell.forget_act = std::move(forget_act);
            main_cell.write_act = std::move(write_act);
            main_cell.candidate_act = std::move(candidate_act);
            main_cell.state_act = std::move(state_act);
            main_cell.read_act = std::move(read_act);
            input_dims = std::move(in_forget_kernel_input_dims);
            output_dims = std::move(out_forget_kernel_input_dims);
            set_foremost(foremost);
        }
        inline LSTMNeuralNetwork(const Self& network) :
                main_cell(network.main_cell),
                output_seq_size_func(network.output_seq_size_func),
                reversed(network.reversed),
                foremost(network.foremost),
                input_dims(network.input_dims),
                output_dims(network.output_dims),
                cells(0),
                state(network.state),
                batch_size(network.batch_size),
                input_seq_length(network.input_seq_length),
                output_seq_length(network.output_seq_length),
                output_seq_delay(network.output_seq_delay) {
            for (std::size_t i = 0; i < network.cells.size(); i++)
                cells.push_back(Cell(network.cells[i]));
        }
        inline LSTMNeuralNetwork(Self&& network) {
            swap(*this, network);
        }
        ~LSTMNeuralNetwork() = default;
        inline Self& operator=(Self network) {
            swap(*this, network);
            return *this;
        }
        inline Root* clone() const {
            return new LSTMNeuralNetwork(*this);
        }
        inline bool is_reversed() const {
            return reversed;
        }
        inline void reverse() {
            reversed = !reversed;
        }
        inline const typename Root::Dims& get_input_dims() const {
            return input_dims;
        }
        inline const typename Root::Dims& get_output_dims() const {
            return output_dims;
        }
        inline std::vector<const Layer<Scalar,Rank>*> get_layers() const {
            std::vector<const Layer<Scalar,Rank>*> layer_ptrs(13);
            populate_layer_vector<const Layer<Scalar,Rank>*>(layer_ptrs);
            return layer_ptrs;
        }
        inline std::vector<Layer<Scalar,Rank>*> get_layers() {
            std::vector<Layer<Scalar,Rank>*> layer_ptrs(13);
            populate_layer_vector<Layer<Scalar,Rank>*>(layer_ptrs);
            return layer_ptrs;
        }
        inline bool is_foremost() const {
            return foremost;
        }
        inline void set_foremost(bool foremost) {
            main_cell.input_forget_kernel->set_input_layer(foremost);
            main_cell.input_write_kernel->set_input_layer(foremost);
            main_cell.input_candidate_kernel->set_input_layer(foremost);
            main_cell.input_read_kernel->set_input_layer(foremost);
            this->foremost = foremost;
        }
        inline void empty_caches() {
            main_cell.input_read_kernel->empty_cache();
            main_cell.input_candidate_kernel->empty_cache();
            main_cell.input_write_kernel->empty_cache();
            main_cell.input_forget_kernel->empty_cache();
            main_cell.output_read_kernel->empty_cache();
            main_cell.output_candidate_kernel->empty_cache();
            main_cell.output_write_kernel->empty_cache();
            main_cell.output_forget_kernel->empty_cache();
            main_cell.write_act->empty_cache();
            main_cell.forget_act->empty_cache();
            main_cell.candidate_act->empty_cache();
            main_cell.state_act->empty_cache();
            main_cell.read_act->empty_cache();
            main_cell.forget_filter_cache = TimeStepData();
            main_cell.prev_state_cache = TimeStepData();
            main_cell.write_filter_cache = TimeStepData();
            main_cell.candidate_cache = TimeStepData();
            main_cell.read_filter_cache = TimeStepData();
            main_cell.activated_state_cache = TimeStepData();
            main_cell.weighted_input_forget_cache = TimeStepData();
            main_cell.weighted_output_forget_cache = TimeStepData();
            main_cell.weighted_input_write_cache = TimeStepData();
            main_cell.weighted_output_write_cache = TimeStepData();
            main_cell.weighted_input_candidate_cache = TimeStepData();
            main_cell.weighted_output_candidate_cache = TimeStepData();
            main_cell.weighted_input_read_cache = TimeStepData();
            main_cell.weighted_output_read_cache = TimeStepData();
            // Clear the state as well.
            batch_size = -1;
            state = TimeStepData();
            input_seq_length = -1;
            output_seq_length = -1;
            output_seq_delay = -1;
            cells = std::vector<Cell>(0);
        }
        inline typename Root::Data propagate(typename Root::Data input, bool training) {
            Dimensions<std::size_t,Root::DATA_RANK> data_dims = input.dimensions();
            assert(input_dims == data_dims.template demote<2>());
            int samples = data_dims(0);
            int input_seq_length = data_dims(1);
            std::pair<std::size_t,std::size_t> output_seq_info = output_seq_size_func((std::size_t) input_seq_length);
            int output_seq_length = (int) output_seq_info.first;
            int output_seq_delay = (int) output_seq_info.second;
            assert(output_seq_length > 0);
            if (reversed)
                Base::reverse_along_time_axis(input);
            int output_end = output_seq_length + output_seq_delay;
            int time_steps = std::max(input_seq_length, output_end);
            // Only unroll the network in training mode and if the sequence alignment has changed.
            if (training && (input_seq_length != this->input_seq_length ||
                    output_seq_length != this->output_seq_length || output_seq_delay != this->output_seq_delay))
                unroll_network(time_steps, input_seq_length);
            setup_hidden_state(samples);
            RankwiseArray input_offsets, output_offsets;
            RankwiseArray input_extents = data_dims;
            RankwiseArray output_extents = output_dims.template promote<2>();
            input_offsets.fill(0);
            output_offsets.fill(0);
            input_extents[1] = 1;
            output_extents[0] = samples;
            typename Root::Data out;
            if (output_seq_length > 1) {
                output_extents[1] = output_seq_length;
                out = typename Root::Data(output_extents);
            }
            output_extents[1] = 1;
            Dimensions<std::size_t,Rank + 1> input_time_step_dims = input_dims.template promote<>();
            input_time_step_dims(0) = samples;
            TimeStepData hidden_out;
            for (int i = 0; i < time_steps; ++i) {
                Cell& cell = !training || i == 0 ? main_cell : cells[i - 1];
                TimeStepData input_res;
                // State update.
                if (i < input_seq_length) {
                    if (input_seq_length > 1) {
                        typename Root::Data input_slice = input.slice(input_offsets, input_extents);
                        input_offsets[1] += 1;
                        input_res = TensorMap<Scalar,Rank + 1>(input_slice.data(), input_time_step_dims);
                    } else
                        input_res = TensorMap<Scalar,Rank + 1>(input.data(), input_time_step_dims);
                    if (i == 0) {
                        // There must be an input at this time step and there cannot be an output from the previous one.
                        TimeStepData weighted_input_forget = cell.input_forget_kernel->pass_forward(input_res, training);
                        // Cache the factors of the multiplication for the backward pass.
                        cell.forget_filter_cache = cell.forget_act->pass_forward(std::move(weighted_input_forget),
                                training);
                        cell.prev_state_cache = std::move(state);
                        // Selective remembrance.
                        state = cell.forget_filter_cache * cell.prev_state_cache;
                        TimeStepData weighted_input_write = cell.input_write_kernel->pass_forward(input_res, training);
                        cell.write_filter_cache = cell.write_act->pass_forward(std::move(weighted_input_write), training);
                        TimeStepData weighted_input_candidates = cell.input_candidate_kernel->pass_forward(input_res,
                                training);
                        cell.candidate_cache = cell.candidate_act->pass_forward(std::move(weighted_input_candidates),
                                training);
                        state += cell.write_filter_cache * cell.candidate_cache;
                    } else {
                        // There is both an input and an output from the previous time step.
                        TimeStepData weighted_input_forget = cell.input_forget_kernel->pass_forward(input_res, training);
                        TimeStepData weighted_output_forget = cell.output_forget_kernel->pass_forward(hidden_out,
                                training);
                        TimeStepData weighted_forget;
                        if (MulInt) {
                            if (training) {
                                cell.weighted_input_forget_cache = std::move(weighted_input_forget);
                                cell.weighted_output_forget_cache = std::move(weighted_output_forget);
                                weighted_forget = cell.weighted_input_forget_cache * cell.weighted_output_forget_cache;
                            } else
                                weighted_forget = weighted_input_forget * weighted_output_forget;
                        } else
                            weighted_forget = weighted_input_forget + weighted_output_forget;
                        cell.forget_filter_cache = cell.forget_act->pass_forward(std::move(weighted_forget), training);
                        cell.prev_state_cache = std::move(state);
                        state = cell.forget_filter_cache * cell.prev_state_cache;
                        TimeStepData weighted_input_write = cell.input_write_kernel->pass_forward(input_res, training);
                        TimeStepData weighted_output_write = cell.output_write_kernel->pass_forward(hidden_out, training);
                        TimeStepData weighted_write;
                        if (MulInt) {
                            if (training) {
                                cell.weighted_input_write_cache = std::move(weighted_input_write);
                                cell.weighted_output_write_cache = std::move(weighted_output_write);
                                weighted_write = cell.weighted_input_write_cache *
                                        cell.weighted_output_write_cache;
                            } else
                                weighted_write = weighted_input_write * weighted_output_write;
                        } else
                            weighted_write = weighted_input_write + weighted_output_write;
                        cell.write_filter_cache = cell.write_act->pass_forward(std::move(weighted_write), training);
                        TimeStepData weighted_input_candidates = cell.input_candidate_kernel->
                                pass_forward(input_res, training);
                        TimeStepData weighted_output_candidates = cell.output_candidate_kernel->
                                pass_forward(hidden_out, training);
                        TimeStepData weighted_candidates;
                        if (MulInt) {
                            if (training) {
                                cell.weighted_input_candidate_cache = std::move(weighted_input_candidates);
                                cell.weighted_output_candidate_cache = std::move(weighted_output_candidates);
                                weighted_candidates = cell.weighted_input_candidate_cache *
                                        cell.weighted_output_candidate_cache;
                            } else
                                weighted_candidates = weighted_input_candidates * weighted_output_candidates;
                        } else
                            weighted_candidates = weighted_input_candidates + weighted_output_candidates;
                        cell.candidate_cache = cell.candidate_act->pass_forward(std::move(weighted_candidates), training);
                        state += cell.write_filter_cache * cell.candidate_cache;
                    }
                } else {
                    // There is only the output from the previous time step and no new input (i must be greater than 0).
                    TimeStepData weighted_output_forget = cell.output_forget_kernel->pass_forward(hidden_out, training);
                    cell.forget_filter_cache = cell.forget_act->pass_forward(std::move(weighted_output_forget), training);
                    cell.prev_state_cache = std::move(state);
                    state = cell.forget_filter_cache * cell.prev_state_cache;
                    TimeStepData weighted_output_write = cell.output_write_kernel->pass_forward(hidden_out, training);
                    cell.write_filter_cache = cell.write_act->pass_forward(std::move(weighted_output_write), training);
                    TimeStepData weighted_output_candidates = cell.output_candidate_kernel->pass_forward(hidden_out,
                            training);
                    cell.candidate_cache = cell.candidate_act->pass_forward(std::move(weighted_output_candidates),
                            training);
                    state += cell.write_filter_cache * cell.candidate_cache;
                }
                // Output computation.
                TimeStepData weighted_read;
                if (i < input_seq_length) {
                    if (i == 0)
                        weighted_read = cell.input_read_kernel->pass_forward(input_res, training);
                    else {
                        TimeStepData weighted_input_read = cell.input_read_kernel->pass_forward(input_res, training);
                        TimeStepData weighted_output_read = cell.output_read_kernel->pass_forward(hidden_out, training);
                        if (MulInt) {
                            if (training) {
                                cell.weighted_input_read_cache = std::move(weighted_input_read);
                                cell.weighted_output_read_cache = std::move(weighted_output_read);
                                weighted_read = cell.weighted_input_read_cache *
                                        cell.weighted_output_read_cache;
                            } else
                                weighted_read = weighted_input_read * weighted_output_read;
                        } else
                            weighted_read = weighted_input_read + weighted_output_read;
                    }
                } else
                    weighted_read = cell.output_read_kernel->pass_forward(hidden_out, training);
                cell.read_filter_cache = cell.read_act->pass_forward(std::move(weighted_read), training);
                cell.activated_state_cache = cell.state_act->pass_forward(state, training);
                hidden_out = cell.read_filter_cache * cell.activated_state_cache;
                // If there is a non-hidden output at this time step...
                if (i >= output_seq_delay && i < output_end) {
                    TensorMap<Scalar,Root::DATA_RANK> out_i_seq(hidden_out.data(), output_extents);
                    if (output_seq_length > 1) {
                        out.slice(output_offsets, output_extents) = out_i_seq;
                        output_offsets[1] += 1;
                    } else
                        out = out_i_seq;
                }
            }
            batch_size = samples;
            this->input_seq_length = input_seq_length;
            this->output_seq_length = output_seq_length;
            this->output_seq_delay = output_seq_delay;
            return out;
        }
        inline typename Root::Data backpropagate(typename Root::Data out_grad) {
            Dimensions<std::size_t,Root::DATA_RANK> data_dims = out_grad.dimensions();
            assert(output_dims == data_dims.template demote<2>() && batch_size == data_dims(0) &&
                    output_seq_length == data_dims(1));
            RankwiseArray output_offsets;
            RankwiseArray output_extents = data_dims;
            RankwiseArray input_offsets;
            RankwiseArray input_extents = input_dims.template promote<2>();
            output_offsets.fill(0);
            output_offsets[1] = output_seq_length - 1;
            input_offsets.fill(0);
            input_offsets[1] = input_seq_length - 1;
            output_extents[1] = 1;
            input_extents[0] = batch_size;
            typename Root::Data prev_out_grad;
            if (input_seq_length > 1) {
                input_extents[1] = input_seq_length;
                prev_out_grad = typename Root::Data(input_extents);
            }
            input_extents[1] = 1;
            TimeStepData state_grad(state.dimensions());
            TimeStepData hidden_out_grad(state.dimensions());
            state_grad.setZero();
            hidden_out_grad.setZero();
            Dimensions<std::size_t,Rank + 1> out_time_step_dims = output_dims.template promote<>();
            out_time_step_dims(0) = batch_size;
            int output_end = output_seq_length + output_seq_delay;
            int time_steps = std::max(input_seq_length, output_end);
            for (int i = time_steps - 1; i >= 0; --i) {
                Cell& cell = i == 0 ? main_cell : cells[i - 1];
                // If there was a non-hidden output at the time step, let the gradients flow into the hidden output gradients.
                if (i >= output_seq_delay && i < output_end) {
                    if (output_seq_length == 1) {
                        hidden_out_grad += TensorMap<Scalar,Rank + 1>(out_grad.data(), out_time_step_dims);
                    } else {
                        typename Root::Data out_grad_seq = out_grad.slice(output_offsets, output_extents);
                        output_offsets[1] -= 1;
                        hidden_out_grad += TensorMap<Scalar,Rank + 1>(out_grad_seq.data(), out_time_step_dims);
                    }
                }
                state_grad += cell.state_act->pass_back(cell.read_filter_cache * hidden_out_grad);
                TimeStepData weighted_read_grad = cell.read_act->pass_back(cell.activated_state_cache * hidden_out_grad);
                TimeStepData candidate_grad = cell.candidate_act->pass_back(cell.write_filter_cache * state_grad);
                TimeStepData weighted_write_grad = cell.write_act->pass_back(cell.candidate_cache * state_grad);
                TimeStepData weighted_forget_grad = cell.forget_act->pass_back(cell.prev_state_cache * state_grad);
                state_grad *= cell.forget_filter_cache;
                if (i < input_seq_length) {
                    TimeStepData prev_out_grad_i;
                    if (MulInt) {
                        if (i != 0) {
                            // Calculate the previous hidden output gradients.
                            hidden_out_grad = cell.output_read_kernel->pass_back(
                                            cell.weighted_input_read_cache * weighted_read_grad) +
                                    cell.output_candidate_kernel->pass_back(
                                            cell.weighted_input_candidate_cache * candidate_grad) +
                                    cell.output_write_kernel->pass_back(
                                            cell.weighted_input_write_cache * weighted_write_grad) +
                                    cell.output_forget_kernel->pass_back(
                                            cell.weighted_input_forget_cache * weighted_forget_grad);
                            // Calculate the input gradients.
                            prev_out_grad_i = cell.input_read_kernel->pass_back(
                                            cell.weighted_output_read_cache * weighted_read_grad) +
                                    cell.input_candidate_kernel->pass_back(
                                            cell.weighted_output_candidate_cache * candidate_grad) +
                                    cell.input_write_kernel->pass_back(
                                            cell.weighted_output_write_cache * weighted_write_grad) +
                                    cell.input_forget_kernel->pass_back(
                                            cell.weighted_output_forget_cache * weighted_forget_grad);
                        } else {
                            prev_out_grad_i = cell.input_read_kernel->pass_back(std::move(weighted_read_grad)) +
                                    cell.input_candidate_kernel->pass_back(std::move(candidate_grad)) +
                                    cell.input_write_kernel->pass_back(std::move(weighted_write_grad)) +
                                    cell.input_forget_kernel->pass_back(std::move(weighted_forget_grad));
                        }
                    } else {
                        if (i != 0) {
                            hidden_out_grad = cell.output_read_kernel->pass_back(weighted_read_grad) +
                                    cell.output_candidate_kernel->pass_back(candidate_grad) +
                                    cell.output_write_kernel->pass_back(weighted_write_grad) +
                                    cell.output_forget_kernel->pass_back(weighted_forget_grad);
                        }
                        prev_out_grad_i = cell.input_read_kernel->pass_back(std::move(weighted_read_grad)) +
                                cell.input_candidate_kernel->pass_back(std::move(candidate_grad)) +
                                cell.input_write_kernel->pass_back(std::move(weighted_write_grad)) +
                                cell.input_forget_kernel->pass_back(std::move(weighted_forget_grad));
                    }
                    if (!foremost) {
                        TensorMap<Scalar,Root::DATA_RANK> prev_out_grad_i_seq(prev_out_grad_i.data(), input_extents);
                        if (input_seq_length > 1) {
                            prev_out_grad.slice(input_offsets, input_extents) = prev_out_grad_i_seq;
                            input_offsets[1] -= 1;
                        } else
                            prev_out_grad = prev_out_grad_i_seq;
                    }
                } else {
                    hidden_out_grad = cell.output_read_kernel->pass_back(std::move(weighted_read_grad)) +
                            cell.output_candidate_kernel->pass_back(std::move(candidate_grad)) +
                            cell.output_write_kernel->pass_back(std::move(weighted_write_grad)) +
                            cell.output_forget_kernel->pass_back(std::move(weighted_forget_grad));
                }
            }
            return prev_out_grad;
        }
        inline friend void swap(Self& network1, Self& network2) {
            using std::swap;
            swap(network1.main_cell, network2.main_cell);
            swap(network1.output_seq_size_func, network2.output_seq_size_func);
            swap(network1.reversed, network2.reversed);
            swap(network1.foremost, network2.foremost);
            swap(network1.input_dims, network2.input_dims);
            swap(network1.output_dims, network2.output_dims);
            swap(network1.cells, network2.cells);
            swap(network1.state, network2.state);
            swap(network1.batch_size, network2.batch_size);
            swap(network1.input_seq_length, network2.input_seq_length);
            swap(network1.output_seq_length, network2.output_seq_length);
            swap(network1.output_seq_delay, network2.output_seq_delay);
        }
    private:
        template<typename _LayerPtr>
        inline void populate_layer_vector(std::vector<_LayerPtr>& layer_ptrs) const {
            layer_ptrs[0] = main_cell.input_forget_kernel.get();
            layer_ptrs[1] = main_cell.output_forget_kernel.get();
            layer_ptrs[2] = main_cell.input_write_kernel.get();
            layer_ptrs[3] = main_cell.output_write_kernel.get();
            layer_ptrs[4] = main_cell.input_candidate_kernel.get();
            layer_ptrs[5] = main_cell.output_candidate_kernel.get();
            layer_ptrs[6] = main_cell.input_read_kernel.get();
            layer_ptrs[7] = main_cell.output_read_kernel.get();
            layer_ptrs[8] = main_cell.forget_act.get();
            layer_ptrs[9] = main_cell.write_act.get();
            layer_ptrs[10] = main_cell.candidate_act.get();
            layer_ptrs[11] = main_cell.read_act.get();
            layer_ptrs[12] = main_cell.state_act.get();
        }
        inline void unroll_network(std::size_t time_steps, std::size_t input_seq_length) {
            if (time_steps > 1) {
                empty_caches();
                cells = std::vector<Cell>(time_steps - 1);
                for (int j = 1; j < time_steps; ++j) {
                    Cell& cell = cells[j - 1];
                    cell.output_forget_kernel = KernelPtr<Scalar,Rank>(static_cast<KernelLayer<Scalar,Rank>*>(
                            main_cell.output_forget_kernel->clone_with_shared_params()));
                    cell.output_write_kernel = KernelPtr<Scalar,Rank>(static_cast<KernelLayer<Scalar,Rank>*>(
                            main_cell.output_write_kernel->clone_with_shared_params()));
                    cell.output_candidate_kernel = KernelPtr<Scalar,Rank>(static_cast<KernelLayer<Scalar,Rank>*>(
                            main_cell.output_candidate_kernel->clone_with_shared_params()));
                    cell.output_read_kernel = KernelPtr<Scalar,Rank>(static_cast<KernelLayer<Scalar,Rank>*>(
                            main_cell.output_read_kernel->clone_with_shared_params()));
                    cell.write_act = ActivationPtr<Scalar,Rank>(static_cast<ActivationLayer<Scalar,Rank>*>(
                            main_cell.write_act->clone_with_shared_params()));
                    cell.forget_act = ActivationPtr<Scalar,Rank>(static_cast<ActivationLayer<Scalar,Rank>*>(
                            main_cell.forget_act->clone_with_shared_params()));
                    cell.candidate_act = ActivationPtr<Scalar,Rank>(static_cast<ActivationLayer<Scalar,Rank>*>(
                            main_cell.candidate_act->clone_with_shared_params()));
                    cell.state_act = ActivationPtr<Scalar,Rank>(static_cast<ActivationLayer<Scalar,Rank>*>(
                            main_cell.state_act->clone_with_shared_params()));
                    cell.read_act = ActivationPtr<Scalar,Rank>(static_cast<ActivationLayer<Scalar,Rank>*>(
                            main_cell.read_act->clone_with_shared_params()));
                    if (j < input_seq_length) {
                        cell.input_forget_kernel = KernelPtr<Scalar,Rank>(static_cast<KernelLayer<Scalar,Rank>*>(
                                main_cell.input_forget_kernel->clone_with_shared_params()));
                        cell.input_write_kernel = KernelPtr<Scalar,Rank>(static_cast<KernelLayer<Scalar,Rank>*>(
                                main_cell.input_write_kernel->clone_with_shared_params()));
                        cell.input_candidate_kernel = KernelPtr<Scalar,Rank>(static_cast<KernelLayer<Scalar,Rank>*>(
                                main_cell.input_candidate_kernel->clone_with_shared_params()));
                        cell.input_read_kernel = KernelPtr<Scalar,Rank>(static_cast<KernelLayer<Scalar,Rank>*>(
                                main_cell.input_read_kernel->clone_with_shared_params()));
                    }
                }
            } else
                cells = std::vector<Cell>(0);
        }
        inline void setup_hidden_state(std::size_t samples) {
            if (!Stateful || batch_size == -1) {
                Dimensions<std::size_t,Rank + 1> dims = main_cell.forget_act->get_output_dims().template promote<>();
                dims(0) = samples;
                state = TimeStepData(dims);
                state.setZero();
            } else if (samples != batch_size) {
                std::array<std::size_t,Rank + 1> offsets;
                std::array<std::size_t,Rank + 1> extents = main_cell.forget_act->get_output_dims().template promote<>();
                offsets.fill(0);
                extents[0] = samples;
                TimeStepData new_state;
                if (samples > batch_size) {
                    new_state = TimeStepData(extents);
                    new_state.setZero();
                    extents[0] = batch_size;
                    new_state.slice(offsets, extents) = state;
                } else
                    new_state = state.slice(offsets, extents);
                state = std::move(new_state);
            }
        }
        /**
        * A struct representing a cell in the unrolled LSTM.
        */
        struct Cell {
            inline Cell() { }
            inline Cell(const Cell& cell) :
                    input_forget_kernel(KernelPtr<Scalar,Rank>(static_cast<KernelLayer<Scalar,Rank>*>(
                            cell.input_forget_kernel->clone()))),
                    output_forget_kernel(KernelPtr<Scalar,Rank>(static_cast<KernelLayer<Scalar,Rank>*>(
                            cell.output_forget_kernel->clone()))),
                    input_write_kernel(KernelPtr<Scalar,Rank>(static_cast<KernelLayer<Scalar,Rank>*>(
                            cell.input_write_kernel->clone()))),
                    output_write_kernel(KernelPtr<Scalar,Rank>(static_cast<KernelLayer<Scalar,Rank>*>(
                            cell.output_write_kernel->clone()))),
                    input_candidate_kernel(KernelPtr<Scalar,Rank>(static_cast<KernelLayer<Scalar,Rank>*>(
                            cell.input_candidate_kernel->clone()))),
                    output_candidate_kernel(KernelPtr<Scalar,Rank>(static_cast<KernelLayer<Scalar,Rank>*>(
                            cell.output_candidate_kernel->clone()))),
                    input_read_kernel(KernelPtr<Scalar,Rank>(static_cast<KernelLayer<Scalar,Rank>*>(
                            cell.input_read_kernel->clone()))),
                    output_read_kernel(KernelPtr<Scalar,Rank>(static_cast<KernelLayer<Scalar,Rank>*>(
                            cell.output_read_kernel->clone()))),
                    write_act(ActivationPtr<Scalar,Rank>(static_cast<ActivationLayer<Scalar,Rank>*>(
                            cell.write_act->clone()))),
                    forget_act(ActivationPtr<Scalar,Rank>(static_cast<ActivationLayer<Scalar,Rank>*>(
                            cell.forget_act->clone()))),
                    candidate_act(ActivationPtr<Scalar,Rank>(static_cast<ActivationLayer<Scalar,Rank>*>(
                            cell.candidate_act->clone()))),
                    state_act(ActivationPtr<Scalar,Rank>(static_cast<ActivationLayer<Scalar,Rank>*>(
                            cell.state_act->clone()))),
                    read_act(ActivationPtr<Scalar,Rank>(static_cast<ActivationLayer<Scalar,Rank>*>(
                            cell.read_act->clone()))),
                    forget_filter_cache(cell.forget_filter_cache),
                    prev_state_cache(cell.prev_state_cache),
                    write_filter_cache(cell.write_filter_cache),
                    candidate_cache(cell.candidate_cache),
                    read_filter_cache(cell.read_filter_cache),
                    activated_state_cache(cell.activated_state_cache),
                    weighted_input_forget_cache(cell.weighted_input_forget_cache),
                    weighted_output_forget_cache(cell.weighted_output_forget_cache),
                    weighted_input_write_cache(cell.weighted_input_write_cache),
                    weighted_output_write_cache(cell.weighted_output_write_cache),
                    weighted_input_candidate_cache(cell.weighted_input_candidate_cache),
                    weighted_output_candidate_cache(cell.weighted_output_candidate_cache),
                    weighted_input_read_cache(cell.weighted_input_read_cache),
                    weighted_output_read_cache(cell.weighted_output_read_cache) { }
            KernelPtr<Scalar,Rank> input_forget_kernel, output_forget_kernel, input_write_kernel, output_write_kernel,
                    input_candidate_kernel, output_candidate_kernel, input_read_kernel, output_read_kernel;
            ActivationPtr<Scalar,Rank> forget_act, write_act, candidate_act, state_act, read_act;
            // Caches for the derivation of multiplicative filtering operations.
            TimeStepData forget_filter_cache, prev_state_cache, write_filter_cache, candidate_cache,
                    read_filter_cache, activated_state_cache;
            // Caches for the derivation of multiplicative integration operations.
            TimeStepData weighted_input_forget_cache, weighted_output_forget_cache, weighted_input_write_cache,
                    weighted_output_write_cache, weighted_input_candidate_cache, weighted_output_candidate_cache,
                    weighted_input_read_cache, weighted_output_read_cache;
        };
        Cell main_cell;
        OutputSeqSizeFunc output_seq_size_func;
        bool reversed, foremost;
        typename Root::Dims input_dims, output_dims;
        // State.
        std::vector<Cell> cells;
        TimeStepData state;
        int batch_size, input_seq_length, output_seq_length, output_seq_delay;
    };


    /**
    * An enumeration type for the different ways the outputs of sub-modules of neural
    * networks may be merged.
    */
    enum ParallelOutputMergeType { PARALLEL_CONCAT_LO_RANK, PARALLEL_CONCAT_HI_RANK, PARALLEL_SUM, PARALLEL_MUL };

    /**
    * A class template for a parallel neural network that consists of one or more
    * lanes of non-sequential neural networks with the same input dimensions. Inputs
    * and gradients are propagated through the lanes simultaneously using multithreading.
    * The outputs of the lanes are merged by concatenation (either along the lowest
    * or hightest rank), summation, or multiplication.
    *
    * \see https://arxiv.org/abs/1409.4842
    */
    template<typename Scalar, std::size_t Rank, ParallelOutputMergeType MergeType = PARALLEL_CONCAT_HI_RANK>
    class ParallelNeuralNetwork :
            public CompositeNeuralNetwork<Scalar,Rank,false,NeuralNetwork<Scalar,Rank,false>> {
        typedef NeuralNetwork<Scalar,Rank,false> Base;
        typedef ParallelNeuralNetwork<Scalar,Rank,MergeType> Self;
        typedef NeuralNetPtr<Scalar,Rank,false> Lane;
        typedef std::array<std::size_t,Base::DATA_RANK> RankwiseArray;
        static_assert(MergeType >= PARALLEL_CONCAT_LO_RANK && MergeType <= PARALLEL_MUL, "illegal merge type value");
        static constexpr std::size_t CONCAT_RANK = MergeType == PARALLEL_CONCAT_HI_RANK ? Rank - 1 : 0;
        static constexpr std::size_t CONCAT_BATCH_RANK = CONCAT_RANK + 1;
    public:
        /**
        * @param lanes A vector of unique pointers to non-sequential neural networks.
        * @param foremost Whether the network is to function as a foremost network.
        */
        inline ParallelNeuralNetwork(std::vector<Lane>&& lanes, bool foremost = true) :
                lanes(std::move(lanes)),
                foremost(foremost),
                outputs(this->lanes.size()) {
            assert(this->lanes.size() > 0 && "lanes must contain at least 1 element");
            assert(this->lanes[0] != nullptr && "lanes contains null pointers");
            Base& first_lane = *this->lanes[0];
            const typename Base::Dims& input_dims = first_lane.get_input_dims();
            typename Base::Dims output_dims = first_lane.get_output_dims();
            for (std::size_t i = 1; i < this->lanes.size(); ++i) {
                assert(this->lanes[i] != nullptr && "lanes contains null pointers");
                Base& lane = *this->lanes[i];
                assert(input_dims == lane.get_input_dims());
                const typename Base::Dims& lane_output_dims = lane.get_output_dims();
                if (MergeType == PARALLEL_CONCAT_HI_RANK || MergeType == PARALLEL_CONCAT_LO_RANK) {
                    if (MergeType == PARALLEL_CONCAT_HI_RANK) {
                        for (std::size_t i = 0; i < +CONCAT_RANK; ++i)
                            assert(output_dims(i) == lane_output_dims(i));
                    } else {
                        for (std::size_t i = Rank - 1; i > +CONCAT_RANK; --i)
                            assert(output_dims(i) == lane_output_dims(i));
                    }
                    output_dims(+CONCAT_RANK) += lane_output_dims(+CONCAT_RANK);
                } else
                    assert(output_dims == lane_output_dims);
            }
            set_foremost(foremost);
            this->input_dims = first_lane.get_input_dims();
            this->output_dims = output_dims;
        }
        /**
        * @param lane A unique pointer to a non-sequential neural network.
        * @param foremost Whether the network is to function as a foremost network.
        */
        inline ParallelNeuralNetwork(Base&& lane, bool foremost = true) :
                ParallelNeuralNetwork(create_vector(std::move(lane)), foremost) { }
        inline ParallelNeuralNetwork(const Self& network) :
                lanes(network.lanes.size()),
                foremost(network.foremost),
                input_dims(network.input_dims),
                output_dims(network.output_dims),
                outputs(network.outputs) {
            for (std::size_t i = 0; i < lanes.size(); ++i)
                lanes[i] = Lane(network.lanes[i]->clone());
        }
        inline ParallelNeuralNetwork(Self&& network) {
            swap(*this, network);
        }
        ~ParallelNeuralNetwork() = default;
        inline Self& operator=(Self network) {
            swap(*this, network);
            return *this;
        }
        inline Base* clone() const {
            return new ParallelNeuralNetwork(*this);
        }
        inline const typename Base::Dims& get_input_dims() const {
            return input_dims;
        }
        inline const typename Base::Dims& get_output_dims() const {
            return output_dims;
        }
        inline std::vector<const Layer<Scalar,Rank>*> get_layers() const {
            std::vector<const Layer<Scalar,Rank>*> layer_ptrs;
            populate_layer_vector<const Layer<Scalar,Rank>*>(layer_ptrs);
            return layer_ptrs;
        }
        inline std::vector<Layer<Scalar,Rank>*> get_layers() {
            std::vector<Layer<Scalar,Rank>*> layer_ptrs;
            populate_layer_vector<Layer<Scalar,Rank>*>(layer_ptrs);
            return layer_ptrs;
        }
        inline std::vector<Base*> get_modules() {
            std::vector<Base*> modules;
            for (std::size_t i = 0; i < lanes.size(); ++i)
                modules.push_back(lanes[i].get());
            return modules;
        }
        inline bool is_foremost() const {
            return foremost;
        }
        inline void set_foremost(bool foremost) {
            for (std::size_t i = 0; i < lanes.size(); ++i)
                lanes[i]->set_foremost(foremost);
            this->foremost = foremost;
        }
        inline void empty_caches() {
            for (std::size_t i = 0; i < lanes.size(); ++i) {
                lanes[i]->empty_caches();
                outputs[i] = typename Base::Data();
            }
        }
        inline typename Base::Data propagate(typename Base::Data input, bool training) {
            assert(input_dims == (Dimensions<std::size_t,Base::DATA_RANK>(input.dimensions()).template demote<>()));
            std::size_t rows = input.dimension(0);
            int pthread_state;
            typename Base::Data out;
            std::size_t lane_num = lanes.size();
            std::size_t helper_thread_num = lane_num - 1;
            pthread_t threads[helper_thread_num];
            pthread_attr_t attr;
            if (helper_thread_num > 0) {
                pthread_state = pthread_attr_init(&attr);
                assert(!pthread_state);
                pthread_state = pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
                assert(!pthread_state);
            }
            PropArgs args_arr[lane_num];
            for (int i = helper_thread_num; i >= 0; --i) {
                PropArgs args;
                args.obj = this;
                args.lane_id = i;
                args.training = training;
                args.in = &input;
                args_arr[i] = args;
                // Leave the first lane to the main thread.
                if (i == 0)
                    propagate(&args_arr[i]);
                else {
                    pthread_state = pthread_create(&threads[i - 1], &attr, propagate, &args_arr[i]);
                    assert(!pthread_state);
                }
            }
            for (std::size_t i = 0; i < lane_num; ++i) {
                if (i == 0) {
                    out = std::move(args_arr[i].out);
                    if (MergeType == PARALLEL_MUL)
                        outputs[i] = out;
                } else {
                    pthread_state = pthread_join(threads[i - 1], nullptr);
                    assert(!pthread_state);
                    if (MergeType == PARALLEL_SUM)
                        out += args_arr[i].out;
                    else if (MergeType == PARALLEL_MUL) {
                        outputs[i] = std::move(args_arr[i].out);
                        out *= outputs[i];
                    } else {
                        // Must be evaluated first due to the dimension difference.
                        typename Base::Data concat = out.concatenate(std::move(args_arr[i].out), +CONCAT_BATCH_RANK);
                        out = std::move(concat);
                    }
                }
            }
            if (helper_thread_num > 0) {
                pthread_state = pthread_attr_destroy(&attr);
                assert(!pthread_state);
            }
            return out;
        }
        inline typename Base::Data backpropagate(typename Base::Data out_grad) {
            assert(output_dims == (Dimensions<std::size_t,Base::DATA_RANK>(out_grad.dimensions()).template demote<>()));
            typename Base::Data prev_out_grad;
            if (foremost)
                prev_out_grad = typename Base::Data();
            else {
                RankwiseArray dims = input_dims.template promote<>();
                dims[0] = out_grad.dimension(0);
                prev_out_grad = typename Base::Data(dims);
                prev_out_grad.setZero();
            }
            int pthread_state;
            std::size_t lane_num = lanes.size();
            std::size_t helper_thread_num = lane_num - 1;
            pthread_t threads[helper_thread_num];
            pthread_attr_t attr;
            if (helper_thread_num > 0) {
                pthread_state = pthread_attr_init(&attr);
                assert(!pthread_state);
                pthread_state = pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
                assert(!pthread_state);
            }
            BackpropArgs args_arr[lane_num];
            int concat_rank_offset = out_grad.dimension(+CONCAT_BATCH_RANK);
            for (int i = helper_thread_num; i >= 0; --i) {
                concat_rank_offset -= lanes[i]->get_output_dims()(+CONCAT_RANK);
                BackpropArgs args;
                args.obj = this;
                args.lane_id = i;
                args.concat_rank_offset = concat_rank_offset;
                args.out_grad = &out_grad;
                args_arr[i] = args;
                // Leave the first lane to the main thread.
                if (i == 0)
                    backpropagate(&args_arr[i]);
                else {
                    pthread_state = pthread_create(&threads[i - 1], &attr, backpropagate, &args_arr[i]);
                    assert(!pthread_state);
                }
            }
            for (std::size_t i = 0; i < lanes.size(); ++i) {
                if (i != 0) {
                    pthread_state = pthread_join(threads[i - 1], nullptr);
                    assert(!pthread_state);
                }
                if (!foremost)
                    prev_out_grad += args_arr[i].prev_out_grad;
            }
            if (helper_thread_num > 0) {
                pthread_state = pthread_attr_destroy(&attr);
                assert(!pthread_state);
            }
            return prev_out_grad;
        }
        inline friend void swap(Self& network1, Self& network2) {
            using std::swap;
            swap(network1.lanes, network2.lanes);
            swap(network1.foremost, network2.foremost);
            swap(network1.input_dims, network2.input_dims);
            swap(network1.output_dims, network2.output_dims);
            swap(network1.outputs, network2.outputs);
        }
    private:
        inline static std::vector<Lane> create_vector(Lane&& net) {
            std::vector<Lane> vec(1);
            vec[0] = std::move(net);
            return vec;
        }
        /**
        * The propagation function executed in a different thread for each lane of a
        * parallel network.
        *
        * @param args_ptr The propagation argument struct containing all necessary
        * information.
        */
        inline static void* propagate(void* args_ptr) {
            PropArgs& args = *((PropArgs*) args_ptr);
            args.out = args.obj->lanes[args.lane_id]->propagate(*args.in, args.training);
            return nullptr;
        }
        /**
        * The back-propagation function executed in a different thread for each lane of a
        * parallel network.
        *
        * @param args_ptr The back-propagation argument struct containing all necessary
        * information.
        */
        inline static void* backpropagate(void* args_ptr) {
            BackpropArgs& args = *((BackpropArgs*) args_ptr);
            Base& lane = *args.obj->lanes[args.lane_id];
            if (MergeType == PARALLEL_SUM)
                args.prev_out_grad = lane.backpropagate(*args.out_grad);
            else if (MergeType == PARALLEL_MUL) {
                typename Base::Data out_grad = *args.out_grad;
                for (std::size_t i = 0; i < args.obj->lanes.size(); ++i) {
                    if (i != (std::size_t) args.lane_id)
                        out_grad *= args.obj->outputs[i];
                }
                args.prev_out_grad = lane.backpropagate(std::move(out_grad));
            } else {
                RankwiseArray offsets;
                RankwiseArray extents = lane.get_output_dims().template promote<>();
                offsets.fill(0);
                offsets[+CONCAT_BATCH_RANK] = args.concat_rank_offset;
                extents[0] = args.out_grad->dimension(0);
                typename Base::Data out_grad_slice = args.out_grad->slice(offsets, extents);
                args.prev_out_grad = lane.backpropagate(std::move(out_grad_slice));
            }
            return nullptr;
        }
        template<typename _LayerPtr>
        inline void populate_layer_vector(std::vector<_LayerPtr>& layer_ptrs) const {
            for (std::size_t i = 0; i < lanes.size(); ++i) {
                std::vector<Layer<Scalar,Rank>*> internal_layer_ptrs = lanes[i]->get_layers();
                for (std::size_t j = 0; j < internal_layer_ptrs.size(); ++j)
                    layer_ptrs.push_back(internal_layer_ptrs[j]);
            }
        }
        std::vector<Lane> lanes;
        bool foremost;
        typename Base::Dims input_dims, output_dims;
        std::vector<typename Base::Data> outputs;
        /**
        * A struct containing the data required for propagation.
        */
        struct PropArgs {
            Self* obj;
            int lane_id;
            bool training;
            typename Base::Data* in;
            typename Base::Data out;
        };
        /**
        * A struct containing the data require for back-propagation.
        */
        struct BackpropArgs {
            Self* obj;
            int lane_id;
            int concat_rank_offset;
            typename Base::Data* out_grad;
            typename Base::Data prev_out_grad;
        };
    };

    /**
    * An alias for a unique pointer to a kernel layer of arbitrary rank and scalar type.
    */
    template<typename Scalar, std::size_t Rank>
    using KernelPtr = std::unique_ptr<KernelLayer<Scalar,Rank>>;

    /**
    * An alias for a unique pointer to an activation layer of arbitrary rank and scalar type.
    */
    template<typename Scalar, std::size_t Rank>
    using ActivationPtr = std::unique_ptr<ActivationLayer<Scalar,Rank>>;

    /**
    * A class template for a simple recurrent neural network (RNN). The network can use multiplicative
    * integration to combine its linearly transformed input and its linearly transformed previous hidden
    * state. A stateful network retains its hidden state across sequences as long as the batch size is
    * constant.
    */
    template<typename Scalar, std::size_t Rank, bool MulInt = false, bool Stateful = false>
    class RecurrentNeuralNetwork : public UnidirectionalNeuralNetwork<Scalar,Rank> {
        typedef NeuralNetwork<Scalar,Rank,true> Root;
        typedef UnidirectionalNeuralNetwork<Scalar,Rank> Base;
        typedef RecurrentNeuralNetwork<Scalar,Rank,MulInt,Stateful> Self;
        typedef std::array<std::size_t,Root::DATA_RANK> RankwiseArray;
        typedef std::function<std::pair<std::size_t,std::size_t>(std::size_t)> OutputSeqSizeFunc;
        typedef Tensor<Scalar,Rank + 1> TimeStepData;
    public:
        /**
        * @param input_kernel The linear layer applied to the input of the network at each time step
        * with an input.
        * @param state_kernel The linear layer applied to the previous hidden state of the network at
        * each time step.
        * @param output_kernel The linear layer applied to the hidden state of the network at each time
        * step with an output.
        * @param state_act The activation function applied to the hidden state at each time step.
        * @param output_act The activation function applied to the linearly transformed hidden state
        * of the network at each time step with an output.
        * @param output_seq_size_func A function parameterized by the input sequence length that
        * determines the output sequence delay and length. The output of the function is a pair of unsigned
        * integers where the first element is the sequence length and the second element is the sequence
        * delay.
        * @param reversed Whether the network is to reverse its inputs along the time-step rank.
        * @param foremost Whether the network is to function as a foremost network.
        */
        inline RecurrentNeuralNetwork(KernelPtr<Scalar,Rank>&& input_kernel, KernelPtr<Scalar,Rank>&& state_kernel,
                KernelPtr<Scalar,Rank>&& output_kernel, ActivationPtr<Scalar,Rank>&& state_act,
                ActivationPtr<Scalar,Rank>&& output_act, OutputSeqSizeFunc output_seq_size_func, bool reversed = false,
                bool foremost = true) :
                    main_cell(),
                    output_seq_size_func(output_seq_size_func),
                    reversed(reversed),
                    foremost(foremost),
                    cells(0),
                    batch_size(-1),
                    input_seq_length(-1),
                    output_seq_length(-1),
                    output_seq_delay(-1) {
            assert(input_kernel && state_kernel && output_kernel && state_act && output_act);
            typename Root::Dims input_layer_input_dims = input_kernel->get_input_dims();
            typename Root::Dims input_layer_output_dims = input_kernel->get_output_dims();
            typename Root::Dims output_layer_output_dims = output_kernel->get_output_dims();
            assert(input_layer_output_dims == state_kernel->get_output_dims() &&
                    input_layer_output_dims == output_kernel->get_input_dims() &&
                    input_layer_output_dims == state_act->get_input_dims() &&
                    output_layer_output_dims == output_act->get_input_dims() &&
                    state_kernel->get_input_dims() == state_kernel->get_output_dims());
            main_cell.input_kernel = std::move(input_kernel);
            main_cell.state_kernel = std::move(state_kernel);
            main_cell.output_kernel = std::move(output_kernel);
            main_cell.state_act = std::move(state_act);
            main_cell.output_act = std::move(output_act);
            input_dims = std::move(input_layer_input_dims);
            output_dims = std::move(output_layer_output_dims);
            set_foremost(foremost);
        }
        // Copy constructor.
        inline RecurrentNeuralNetwork(const Self& network) :
                main_cell(network.main_cell),
                output_seq_size_func(network.output_seq_size_func),
                reversed(network.reversed),
                foremost(network.foremost),
                input_dims(network.input_dims),
                output_dims(network.output_dims),
                cells(0),
                state(network.state),
                batch_size(network.batch_size),
                input_seq_length(network.input_seq_length),
                output_seq_length(network.output_seq_length),
                output_seq_delay(network.output_seq_delay) {
            for (std::size_t i = 0; i < network.cells.size(); i++)
                cells.push_back(Cell(network.cells[i]));
        }
        inline RecurrentNeuralNetwork(Self&& network) {
            swap(*this, network);
        }
        ~RecurrentNeuralNetwork() = default;
        inline Self& operator=(Self network) {
            swap(*this, network);
            return *this;
        }
        inline Root* clone() const {
            return new RecurrentNeuralNetwork(*this);
        }
        inline bool is_reversed() const {
            return reversed;
        }
        inline void reverse() {
            reversed = !reversed;
        }
        inline const typename Root::Dims& get_input_dims() const {
            return input_dims;
        }
        inline const typename Root::Dims& get_output_dims() const {
            return output_dims;
        }
        inline std::vector<const Layer<Scalar,Rank>*> get_layers() const {
            std::vector<const Layer<Scalar,Rank>*> layer_ptrs(5);
            populate_layer_vector<const Layer<Scalar,Rank>*>(layer_ptrs);
            return layer_ptrs;
        }
        inline std::vector<Layer<Scalar,Rank>*> get_layers() {
            std::vector<Layer<Scalar,Rank>*> layer_ptrs(5);
            populate_layer_vector<Layer<Scalar,Rank>*>(layer_ptrs);
            return layer_ptrs;
        }
        inline bool is_foremost() const {
            return foremost;
        }
        inline void set_foremost(bool foremost) {
            main_cell.input_kernel->set_input_layer(foremost);
            this->foremost = foremost;
        }
        inline void empty_caches() {
            main_cell.input_kernel->empty_cache();
            main_cell.state_kernel->empty_cache();
            main_cell.output_kernel->empty_cache();
            main_cell.state_act->empty_cache();
            main_cell.output_act->empty_cache();
            main_cell.state_kernel_cache = TimeStepData();
            main_cell.input_kernel_cache = TimeStepData();
            // Clear the state as well.
            batch_size = -1;
            state = TimeStepData();
            input_seq_length = -1;
            output_seq_length = -1;
            output_seq_delay = -1;
            cells = std::vector<Cell>(0);
        }
        inline typename Root::Data propagate(typename Root::Data input, bool training) {
            Dimensions<std::size_t,Root::DATA_RANK> data_dims = input.dimensions();
            assert(input_dims == data_dims.template demote<2>());
            int samples = data_dims(0);
            int input_seq_length = data_dims(1);
            // Calculate the output sequence length and delay based on the input sequence length.
            std::pair<std::size_t,std::size_t> output_seq_info = output_seq_size_func((std::size_t) input_seq_length);
            int output_seq_length = (int) output_seq_info.first;
            int output_seq_delay = (int) output_seq_info.second;
            assert(output_seq_length > 0);
            if (reversed)
                Base::reverse_along_time_axis(input);
            int output_end = output_seq_length + output_seq_delay;
            int time_steps = std::max(input_seq_length, output_end);
            // If in training mode, unroll the network (unless it has already been unrolled for the same alignment).
            if (training && (input_seq_length != this->input_seq_length ||
                    output_seq_length != this->output_seq_length || output_seq_delay != this->output_seq_delay))
                unroll_network(time_steps, input_seq_length, output_seq_delay, output_end);
            setup_hidden_state(samples);
            RankwiseArray input_offsets, output_offsets;
            RankwiseArray input_extents = data_dims;
            RankwiseArray output_extents = output_dims.template promote<2>();
            input_offsets.fill(0);
            output_offsets.fill(0);
            input_extents[1] = 1;
            output_extents[0] = samples;
            typename Root::Data out;
            // If the output is a single time step prediction, there is no need to create an output tensor.
            if (output_seq_length > 1) {
                output_extents[1] = output_seq_length;
                out = typename Root::Data(output_extents);
            }
            output_extents[1] = 1;
            Dimensions<std::size_t,Rank + 1> input_time_step_dims = input_dims.template promote<>();
            input_time_step_dims(0) = samples;
            for (int i = 0; i < time_steps; ++i) {
                // In inference mode, do not unroll the network.
                Cell& cell = !training || i == 0 ? main_cell : cells[i - 1];
                // Always apply the state kernel.
                state = cell.state_kernel->pass_forward(std::move(state), training);
                // If there is an input for the time step...
                if (i < input_seq_length) {
                    typename Root::Data in_i_seq;
                    if (input_seq_length == 1)
                        in_i_seq = std::move(input);
                    else {
                        in_i_seq = input.slice(input_offsets, input_extents);
                        input_offsets[1] += 1;
                    }
                    TensorMap<Scalar,Rank + 1> in_i(in_i_seq.data(), input_time_step_dims);
                    if (MulInt) {
                        if (training) {
                            /* If multiplicative integration is enabled, cache the factors of the multiplication so that
                            * the function can be differentiated in the backward pass. */
                            cell.state_kernel_cache = state;
                            cell.input_kernel_cache = cell.input_kernel->pass_forward(in_i, training);
                            state *= cell.input_kernel_cache;
                        } else
                            state *= cell.input_kernel->pass_forward(in_i, training);
                    } else
                        state += cell.input_kernel->pass_forward(in_i, training);
                }
                state = cell.state_act->pass_forward(std::move(state), training);
                // If there is an output for the time step...
                if (i >= output_seq_delay && i < output_end) {
                    // If the output is a single time step prediction, just return it.
                    TimeStepData act_out_i = cell.output_act->pass_forward(
                            cell.output_kernel->pass_forward(state, training), training);
                    TensorMap<Scalar,Root::DATA_RANK> out_i_seq(act_out_i.data(), output_extents);
                    if (output_seq_length == 1)
                        out = out_i_seq;
                    else {
                        out.slice(output_offsets, output_extents) = out_i_seq;
                        output_offsets[1] += 1;
                    }
                }
            }
            batch_size = samples;
            this->input_seq_length = input_seq_length;
            this->output_seq_length = output_seq_length;
            this->output_seq_delay = output_seq_delay;
            return out;
        }
        inline typename Root::Data backpropagate(typename Root::Data out_grad) {
            Dimensions<std::size_t,Root::DATA_RANK> data_dims = out_grad.dimensions();
            assert(output_dims == data_dims.template demote<2>() && batch_size == data_dims(0) &&
                    output_seq_length == data_dims(1));
            RankwiseArray output_offsets, input_offsets;
            RankwiseArray output_extents = data_dims;
            RankwiseArray input_extents = input_dims.template promote<2>();
            output_offsets.fill(0);
            input_offsets.fill(0);
            output_offsets[1] = output_seq_length - 1;
            input_offsets[1] = input_seq_length - 1;
            output_extents[1] = 1;
            input_extents[0] = batch_size;
            typename Root::Data prev_out_grad;
            if (input_seq_length > 1) {
                input_extents[1] = input_seq_length;
                prev_out_grad = typename Root::Data(input_extents);
            }
            input_extents[1] = 1;
            TimeStepData state_grad(state.dimensions());
            state_grad.setZero();
            Dimensions<std::size_t,Rank + 1> out_time_step_dims = output_dims.template promote<>();
            out_time_step_dims(0) = batch_size;
            int output_end = output_seq_length + output_seq_delay;
            int time_steps = std::max(input_seq_length, output_end);
            for (int i = time_steps - 1; i >= 0; --i) {
                Cell& cell = i == 0 ? main_cell : cells[i - 1];
                // If there was an output at the time step...
                if (i >= output_seq_delay && i < output_end) {
                    typename Root::Data out_grad_seq_i;
                    if (output_seq_length == 1)
                        out_grad_seq_i = std::move(out_grad);
                    else {
                        out_grad_seq_i = out_grad.slice(output_offsets, output_extents);
                        output_offsets[1] -= 1;
                    }
                    state_grad += cell.output_kernel->pass_back(cell.output_act->pass_back(
                            TensorMap<Scalar,Rank + 1>(out_grad_seq_i.data(), out_time_step_dims)));
                }
                // Always back-propagate the state gradient.
                state_grad = cell.state_act->pass_back(std::move(state_grad));
                // If there was an input at the time step...
                if (i < input_seq_length) {
                    // If it is the foremost layer, the gradients do not need to be propagated further back.
                    if (foremost) {
                        if (MulInt) // Multiplicative integration.
                            cell.input_kernel->pass_back(cell.state_kernel_cache * state_grad);
                        else // Additive integration.
                            cell.input_kernel->pass_back(state_grad);
                    } else if (input_seq_length == 1) {
                        TimeStepData input_i;
                        if (MulInt)
                            input_i = cell.input_kernel->pass_back(cell.state_kernel_cache * state_grad);
                        else
                            input_i = cell.input_kernel->pass_back(state_grad);
                        prev_out_grad = TensorMap<Scalar,Root::DATA_RANK>(input_i.data(), input_extents);
                    } else {
                        TimeStepData input_i;
                        if (MulInt)
                            input_i = cell.input_kernel->pass_back(cell.state_kernel_cache * state_grad);
                        else
                            input_i = cell.input_kernel->pass_back(state_grad);
                        prev_out_grad.slice(input_offsets, input_extents) =
                                TensorMap<Scalar,Root::DATA_RANK>(input_i.data(), input_extents);
                        input_offsets[1] -= 1;
                    }
                    // Compute the the state kernel's gradient.
                    if (MulInt)
                        state_grad = cell.state_kernel->pass_back(cell.input_kernel_cache * state_grad);
                    else
                        state_grad = cell.state_kernel->pass_back(std::move(state_grad));
                } else
                    state_grad = cell.state_kernel->pass_back(std::move(state_grad));
            }
            return prev_out_grad;
        }
        // For the copy-and-swap idiom.
        inline friend void swap(Self& network1, Self& network2) {
            using std::swap;
            swap(network1.main_cell, network2.main_cell);
            swap(network1.output_seq_size_func, network2.output_seq_size_func);
            swap(network1.reversed, network2.reversed);
            swap(network1.foremost, network2.foremost);
            swap(network1.input_dims, network2.input_dims);
            swap(network1.output_dims, network2.output_dims);
            swap(network1.cells, network2.cells);
            swap(network1.state, network2.state);
            swap(network1.batch_size, network2.batch_size);
            swap(network1.input_seq_length, network2.input_seq_length);
            swap(network1.output_seq_length, network2.output_seq_length);
            swap(network1.output_seq_delay, network2.output_seq_delay);
        }
    private:
        template<typename _LayerPtr>
        inline void populate_layer_vector(std::vector<_LayerPtr>& layer_ptrs) const {
            layer_ptrs[0] = main_cell.input_kernel.get();
            layer_ptrs[1] = main_cell.state_kernel.get();
            layer_ptrs[2] = main_cell.output_kernel.get();
            layer_ptrs[3] = main_cell.state_act.get();
            layer_ptrs[4] = main_cell.output_act.get();
        }
        inline void unroll_network(std::size_t time_steps, std::size_t input_seq_length,
                std::size_t output_seq_delay, std::size_t output_end) {
            if (time_steps > 1) {
                // Empty the caches of the main cell to reduce the amount of data to copy.
                empty_caches();
                // Emptying the caches also clears the cell vector, thus it has to be recreated afterwards.
                cells = std::vector<Cell>(time_steps - 1);
                // Unroll the network by creating n -1 copies of the main cell;
                for (int j = 1; j < time_steps; ++j) {
                    Cell& cell = cells[j - 1];
                    cell.state_kernel = KernelPtr<Scalar,Rank>(static_cast<KernelLayer<Scalar,Rank>*>(
                            main_cell.state_kernel->clone_with_shared_params()));
                    cell.state_act = ActivationPtr<Scalar,Rank>(static_cast<ActivationLayer<Scalar,Rank>*>(
                            main_cell.state_act->clone_with_shared_params()));
                    // Only copy the kernels and activations that will actually be used.
                    if (j < input_seq_length) {
                        cell.input_kernel = KernelPtr<Scalar,Rank>(static_cast<KernelLayer<Scalar,Rank>*>(
                                main_cell.input_kernel->clone_with_shared_params()));
                    }
                    if (j >= output_seq_delay && j < output_end) {
                        cell.output_kernel = KernelPtr<Scalar,Rank>(static_cast<KernelLayer<Scalar,Rank>*>(
                                main_cell.output_kernel->clone_with_shared_params()));
                        cell.output_act = ActivationPtr<Scalar,Rank>(static_cast<ActivationLayer<Scalar,Rank>*>(
                                main_cell.output_act->clone_with_shared_params()));
                    }
                }
            } else
                cells = std::vector<Cell>(0);
        }
        inline void setup_hidden_state(std::size_t samples) {
            if (!Stateful || batch_size == -1) {
                Dimensions<std::size_t,Rank + 1> dims = main_cell.input_kernel->get_output_dims().template promote<>();
                dims(0) = samples;
                state = TimeStepData(dims);
                state.setZero();
            } else if (samples != batch_size) {
                // If the network is stateful, retain the state.
                std::array<std::size_t,Rank + 1> offsets;
                std::array<std::size_t,Rank + 1> extents = main_cell.input_kernel->get_output_dims().template promote<>();
                offsets.fill(0);
                extents[0] = samples;
                TimeStepData new_state;
                if (samples > batch_size) {
                    new_state = TimeStepData(extents);
                    new_state.setZero();
                    extents[0] = batch_size;
                    new_state.slice(offsets, extents) = state;
                } else
                    new_state = state.slice(offsets, extents);
                state = std::move(new_state);
            }
        }
        /**
        * A struct representing a cell in the unrolled RNN.
        */
        struct Cell {
            inline Cell() { }
            inline Cell(const Cell& cell) :
                    input_kernel(KernelPtr<Scalar,Rank>(static_cast<KernelLayer<Scalar,Rank>*>(
                            cell.input_kernel->clone()))),
                    state_kernel(KernelPtr<Scalar,Rank>(static_cast<KernelLayer<Scalar,Rank>*>(
                            cell.state_kernel->clone()))),
                    output_kernel(KernelPtr<Scalar,Rank>(static_cast<KernelLayer<Scalar,Rank>*>(
                            cell.output_kernel->clone()))),
                    state_act(ActivationPtr<Scalar,Rank>(static_cast<ActivationLayer<Scalar,Rank>*>(
                            cell.state_act->clone()))),
                    output_act(ActivationPtr<Scalar,Rank>(static_cast<ActivationLayer<Scalar,Rank>*>(
                            cell.output_act->clone()))),
                    state_kernel_cache(cell.state_kernel_cache),
                    input_kernel_cache(cell.input_kernel_cache) { }
            KernelPtr<Scalar,Rank> input_kernel, state_kernel, output_kernel;
            ActivationPtr<Scalar,Rank> state_act, output_act;
            // State and input caches for multiplicative integration.
            TimeStepData state_kernel_cache, input_kernel_cache;
        };
        Cell main_cell;
        OutputSeqSizeFunc output_seq_size_func;
        bool reversed, foremost;
        typename Root::Dims input_dims, output_dims;
        // State.
        std::vector<Cell> cells;
        TimeStepData state;
        int batch_size, input_seq_length, output_seq_length, output_seq_delay;
    };    

    /**
    * A class template for ResNets. These networks take vectors of neural networks as their
    * sub-modules.
    *
    * \see https://arxiv.org/abs/1512.03385
    */
    template<typename Scalar, std::size_t Rank>
    class ResidualNeuralNetwork :
            public CompositeNeuralNetwork<Scalar,Rank,false,NeuralNetwork<Scalar,Rank,false>> {
        typedef NeuralNetwork<Scalar,Rank,false> Base;
        typedef NeuralNetPtr<Scalar,Rank,false> Module;
        typedef ResidualNeuralNetwork<Scalar,Rank> Self;
    public:
        /**
        * @param modules A vector of residual modules.
        * @param foremost Whether the network is to function as a foremost network.
        */
        inline ResidualNeuralNetwork(std::vector<Module>&& modules, bool foremost = true) :
                modules(std::move(modules)),
                foremost(foremost) {
            assert(this->modules.size() > 0 && "modules must contain at least 1 element");
            Base& first_module = *this->modules[0];
            input_dims = first_module.get_input_dims();
            output_dims = this->modules[this->modules.size() - 1]->get_output_dims();
            first_module.set_foremost(foremost);
            typename Base::Dims prev_dims = input_dims;
            for (std::size_t i = 0; i < this->modules.size(); ++i) {
                Base& module = *this->modules[i];
                if (i != 0)
                    module.set_foremost(false);
                assert(module.get_input_dims() == module.get_output_dims() &&
                        "residual module input-output dimension discrepancy");
                assert(prev_dims == module.get_input_dims() && "incompatible module dimensions");
                prev_dims = module.get_output_dims();
            }
        }
        /**
        * @param module A single residual module.
        * @param foremost Whether the network is to function as a foremost network.
        */
        inline ResidualNeuralNetwork(Module&& module, bool foremost = true) :
                ResidualNeuralNetwork(create_vector(std::move(module), foremost)) { }
        inline ResidualNeuralNetwork(const Self& network) :
                modules(network.modules.size()),
                foremost(network.foremost),
                input_dims(network.input_dims),
                output_dims(network.output_dims) {
            for (std::size_t i = 0; i < modules.size(); ++i)
                modules[i] = Module(network.modules[i]->clone());
        }
        inline ResidualNeuralNetwork(Self&& network) {
            swap(*this, network);
        }
        ~ResidualNeuralNetwork() = default;
        inline Self& operator=(Self network) {
            swap(*this, network);
            return *this;
        }
        inline Base* clone() const {
            return new ResidualNeuralNetwork(*this);
        }
        inline const typename Base::Dims& get_input_dims() const {
            return input_dims;
        }
        inline const typename Base::Dims& get_output_dims() const {
            return output_dims;
        }
        inline std::vector<const Layer<Scalar,Rank>*> get_layers() const {
            std::vector<const Layer<Scalar,Rank>*> layer_ptrs;
            populate_layer_vector<const Layer<Scalar,Rank>*>(layer_ptrs);
            return layer_ptrs;
        }
        inline std::vector<Layer<Scalar,Rank>*> get_layers() {
            std::vector<Layer<Scalar,Rank>*> layer_ptrs;
            populate_layer_vector<Layer<Scalar,Rank>*>(layer_ptrs);
            return layer_ptrs;
        }
        inline std::vector<Base*> get_modules() {
            std::vector<Base*> modules;
            for (std::size_t i = 0; i < this->modules.size(); ++i)
                modules.push_back(this->modules[i].get());
            return modules;
        }
        inline bool is_foremost() const {
            return foremost;
        }
        inline void set_foremost(bool foremost) {
            modules[0]->set_foremost(foremost);
            this->foremost = foremost;
        }
        inline void empty_caches() {
            for (std::size_t i = 0; i < modules.size(); ++i)
                modules[i]->empty_caches();
        }
        inline typename Base::Data propagate(typename Base::Data input, bool training) {
            assert(input_dims == (Dimensions<std::size_t,Base::DATA_RANK>(input.dimensions()).template demote<>()));
            for (std::size_t i = 0; i < modules.size(); ++i)
                input += modules[i]->propagate(input, training);
            return input;
        }
        inline typename Base::Data backpropagate(typename Base::Data out_grad) {
            assert(output_dims == (Dimensions<std::size_t,Base::DATA_RANK>(out_grad.dimensions()).template demote<>()));
            for (int i = modules.size() - 1; i >= 0; --i) {
                if (foremost && i == 0)
                    return modules[i]->backpropagate(std::move(out_grad));
                else
                    out_grad += modules[i]->backpropagate(out_grad);
            }
            return out_grad;
        }
        inline friend void swap(Self& network1, Self& network2) {
            using std::swap;
            swap(network1.modules, network2.modules);
            swap(network1.foremost, network2.foremost);
            swap(network1.input_dims, network2.input_dims);
            swap(network1.output_dims, network2.output_dims);
        }
    private:
        inline static std::vector<Module> create_vector(Module&& module) {
            std::vector<Module> vec(1);
            vec[0] = std::move(module);
            return vec;
        }
        template<typename _LayerPtr>
        inline void populate_layer_vector(std::vector<_LayerPtr>& layer_ptrs) const {
            for (std::size_t i = 0; i < modules.size(); ++i) {
                std::vector<Layer<Scalar,Rank>*> internal_layer_ptrs = modules[i]->get_layers();
                for (std::size_t j = 0; j < internal_layer_ptrs.size(); ++j)
                    layer_ptrs.push_back(internal_layer_ptrs[j]);
            }
        }
        std::vector<Module> modules;
        bool foremost;
        typename Base::Dims input_dims, output_dims;
    };    

    /**
    * A class template for a wrapper neural network that enables the use of non-sequential networks on
    * sequential data by joining the 'samples' and 'time steps' ranks of the tensors and splitting them
    * again once the internal, non-sequential network is done processing them.
    */
    template<typename Scalar, std::size_t Rank>
    class SequentialNeuralNetwork :
            public CompositeNeuralNetwork<Scalar,Rank,true,NeuralNetwork<Scalar,Rank,false>> {
        typedef NeuralNetwork<Scalar,Rank,true> Base;
        typedef SequentialNeuralNetwork<Scalar,Rank> Self;
        typedef NeuralNetPtr<Scalar,Rank,false> Net;
        typedef std::array<std::size_t,Base::DATA_RANK> RankwiseArray;
    public:
        /**
        * @param network A unique pointer to a non-sequential neural network to wrap.
        * @param foremost Whether the network is to function as a foremost network.
        */
        inline SequentialNeuralNetwork(Net&& network, bool foremost = true) :
                net(std::move(network)),
                foremost(foremost) {
            assert(net);
            input_dims = net->get_input_dims();
            output_dims = net->get_output_dims();
            joint_input_dims = input_dims.template promote<>();
            joint_output_dims = output_dims.template promote<>();
            split_input_dims = input_dims.template promote<2>();
            split_output_dims = output_dims.template promote<2>();
            set_foremost(foremost);
        }
        inline SequentialNeuralNetwork(const Self& network) :
                net(Net(network.net->clone())),
                foremost(network.foremost),
                input_dims(network.input_dims),
                output_dims(network.output_dims),
                joint_input_dims(network.joint_input_dims),
                joint_output_dims(network.joint_output_dims),
                split_input_dims(network.split_input_dims),
                split_output_dims(network.split_output_dims) { }
        inline SequentialNeuralNetwork(Self&& network) {
            swap(*this, network);
        }
        ~SequentialNeuralNetwork() = default;
        inline Self& operator=(Self network) {
            swap(*this, network);
            return *this;
        }
        inline Base* clone() const {
            return new SequentialNeuralNetwork(*this);
        }
        inline const typename Base::Dims& get_input_dims() const {
            return input_dims;
        }
        inline const typename Base::Dims& get_output_dims() const {
            return output_dims;
        }
        inline std::vector<const Layer<Scalar,Rank>*> get_layers() const {
            return ((const NeuralNetwork<Scalar,Rank,false>&) *net).get_layers();
        }
        inline std::vector<Layer<Scalar,Rank>*> get_layers() {
            return net->get_layers();
        }
        inline std::vector<NeuralNetwork<Scalar,Rank,false>*> get_modules() {
            std::vector<NeuralNetwork<Scalar,Rank,false>*> modules;
            modules.push_back(net.get());
            return modules;
        }
        inline bool is_foremost() const {
            return foremost;
        }
        inline void set_foremost(bool foremost) {
            net->set_foremost(foremost);
            this->foremost = foremost;
        }
        inline void empty_caches() {
            net->empty_caches();
        }
        inline typename Base::Data propagate(typename Base::Data input, bool training) {
            assert(input_dims == (Dimensions<std::size_t,Base::DATA_RANK>(input.dimensions()).template demote<2>()));
            std::size_t batch_size = input.dimension(0);
            std::size_t seq_length = input.dimension(1);
            joint_input_dims[0] = batch_size * seq_length;
            split_output_dims[0] = batch_size;
            split_output_dims[1] = seq_length;
            TensorMap<Scalar,Rank + 1> joint_input(input.data(), joint_input_dims);
            Tensor<Scalar,Rank + 1> out = net->propagate(joint_input, training);
            return TensorMap<Scalar,Rank + 2>(out.data(), split_output_dims);
        }
        inline typename Base::Data backpropagate(typename Base::Data out_grad) {
            assert(output_dims == (Dimensions<std::size_t,Base::DATA_RANK>(out_grad.dimensions()).template demote<2>()));
            assert(split_output_dims[0] == out_grad.dimension(0));
            std::size_t batch_size = out_grad.dimension(0);
            std::size_t seq_length = out_grad.dimension(1);
            joint_output_dims[0] = batch_size * seq_length;
            TensorMap<Scalar,Rank + 1> joint_out_grad(out_grad.data(), joint_output_dims);
            if (foremost) {
                net->backpropagate(joint_out_grad);
                return typename Base::Data();
            } else {
                Tensor<Scalar,Rank + 1> prev_out_grad = net->backpropagate(joint_out_grad);
                split_input_dims[0] = batch_size;
                split_input_dims[1] = seq_length;
                return TensorMap<Scalar,Rank + 2>(prev_out_grad.data(), split_input_dims);
            }
        }
        inline friend void swap(Self& network1, Self& network2) {
            using std::swap;
            swap(network1.net, network2.net);
            swap(network1.foremost, network2.foremost);
            swap(network1.input_dims, network2.input_dims);
            swap(network1.output_dims, network2.output_dims);
            swap(network1.joint_input_dims, network2.joint_input_dims);
            swap(network1.joint_output_dims, network2.joint_output_dims);
            swap(network1.split_input_dims, network2.split_input_dims);
            swap(network1.split_output_dims, network2.split_output_dims);
        }
    private:
        Net net;
        bool foremost;
        typename Base::Dims input_dims, output_dims;
        std::array<std::size_t,Rank + 1> joint_input_dims, joint_output_dims;
        std::array<std::size_t,Rank + 2> split_input_dims, split_output_dims;
    };

    /**
    * A class template for a composite neural network that consists of a set of
    * serially stacked neural network sub-modules.
    */
    template<typename Scalar, std::size_t Rank, bool Sequential>
    class StackedNeuralNetwork :
            public CompositeNeuralNetwork<Scalar,Rank,Sequential,NeuralNetwork<Scalar,Rank,Sequential>> {
        typedef NeuralNetwork<Scalar,Rank,Sequential> Base;
        typedef StackedNeuralNetwork<Scalar,Rank,Sequential> Self;
        typedef NeuralNetPtr<Scalar,Rank,Sequential> Block;
    public:
        /**
        * @param blocks A vector of unique pointers to neural networks.
        * @param foremost Whether the network is to function as a foremost network.
        */
        inline StackedNeuralNetwork(std::vector<Block>&& blocks, bool foremost = true) :
                blocks(std::move(blocks)),
                foremost(foremost) {
            assert(this->blocks.size() > 0 && "blocks must contain at least 1 element");
            assert(this->blocks[0] != nullptr && "blocks contains null pointers");
            Base& first_block = *this->blocks[0];
            input_dims = first_block.get_input_dims();
            output_dims = this->blocks[this->blocks.size() - 1]->get_output_dims();
            typename Base::Dims prev_dims = first_block.get_output_dims();
            for (std::size_t i = 1; i < this->blocks.size(); ++i) {
                assert(this->blocks[i] != nullptr && "blocks contains null pointers");
                Base& block = *this->blocks[i];
                assert(prev_dims == block.get_input_dims() && "incompatible network dimensions");
                block.set_foremost(false);
                prev_dims = block.get_output_dims();
            }
            first_block.set_foremost(foremost);
        }
        /**
        * @param block A unique pointer to a neural network.
        * @param foremost Whether the network is to function as a foremost network.
        */
        inline StackedNeuralNetwork(Block&& block, bool foremost = true) :
                StackedNeuralNetwork(create_vector(std::move(block)), foremost) { }
        inline StackedNeuralNetwork(const Self& network) :
                blocks(network.blocks.size()),
                foremost(network.foremost),
                input_dims(network.input_dims),
                output_dims(network.output_dims) {
            for (std::size_t i = 0; i < blocks.size(); ++i)
                blocks[i] = Block(network.blocks[i]->clone());
        }
        inline StackedNeuralNetwork(Self&& network) {
            swap(*this, network);
        }
        ~StackedNeuralNetwork() = default;
        inline Self& operator=(Self network) {
            swap(*this, network);
            return *this;
        }
        inline Base* clone() const {
            return new StackedNeuralNetwork(*this);
        }
        inline const typename Base::Dims& get_input_dims() const {
            return input_dims;
        }
        inline const typename Base::Dims& get_output_dims() const {
            return output_dims;
        }
        inline std::vector<const Layer<Scalar,Rank>*> get_layers() const {
            std::vector<const Layer<Scalar,Rank>*> layer_ptrs;
            populate_layer_vector<const Layer<Scalar,Rank>*>(layer_ptrs);
            return layer_ptrs;
        }
        inline std::vector<Layer<Scalar,Rank>*> get_layers() {
            std::vector<Layer<Scalar,Rank>*> layer_ptrs;
            populate_layer_vector<Layer<Scalar,Rank>*>(layer_ptrs);
            return layer_ptrs;
        }
        inline std::vector<Base*> get_modules() {
            std::vector<Base*> modules;
            for (std::size_t i = 0; i < blocks.size(); ++i)
                modules.push_back(blocks[i].get());
            return modules;
        }
        inline bool is_foremost() const {
            return foremost;
        }
        inline void set_foremost(bool foremost) {
            blocks[0]->set_foremost(foremost);
            this->foremost = foremost;
        }
        inline void empty_caches() {
            for (std::size_t i = 0; i < blocks.size(); ++i)
                blocks[i]->empty_caches();
        }
        inline typename Base::Data propagate(typename Base::Data input, bool training) {
            assert(input_dims == (Dimensions<std::size_t,Base::DATA_RANK>(input.dimensions()).template demote<>()));
            for (std::size_t i = 0; i < blocks.size(); ++i)
                input = blocks[i]->propagate(std::move(input), training);
            return input;
        }
        inline typename Base::Data backpropagate(typename Base::Data out_grad) {
            assert(output_dims == (Dimensions<std::size_t,Base::DATA_RANK>(out_grad.dimensions()).template demote<>()));
            for (int i = blocks.size() - 1; i >= 0; --i)
                out_grad = blocks[i]->backpropagate(std::move(out_grad));
            return out_grad;
        }
        inline friend void swap(Self& network1, Self& network2) {
            using std::swap;
            swap(network1.blocks, network2.blocks);
            swap(network1.foremost, network2.foremost);
            swap(network1.input_dims, network2.input_dims);
            swap(network1.output_dims, network2.output_dims);
        }
    private:
        inline static std::vector<Block> create_vector(Block&& net) {
            std::vector<Block> vec(1);
            vec[0] = std::move(net);
            return vec;
        }
        template<typename _LayerPtr>
        inline void populate_layer_vector(std::vector<_LayerPtr>& layer_ptrs) const {
            for (std::size_t i = 0; i < blocks.size(); ++i) {
                std::vector<Layer<Scalar,Rank>*> internal_layer_ptrs = blocks[i]->get_layers();
                for (std::size_t j = 0; j < internal_layer_ptrs.size(); ++j)
                    layer_ptrs.push_back(internal_layer_ptrs[j]);
            }
        }
        std::vector<Block> blocks;
        bool foremost;
        typename Base::Dims input_dims, output_dims;
    };

    /**
    * An abstract class template for stochastic gradient descent (SGD) optimizers.
    */
    template<typename Scalar, std::size_t Rank, bool Sequential>
    class SGDOptimizer : public Optimizer<Scalar,Rank,Sequential> {
        typedef Optimizer<Scalar,Rank,Sequential> Base;
    public:
        inline SGDOptimizer(LossSharedPtr<Scalar,Rank,Sequential> loss, std::size_t batch_size) :
                Base::Optimizer(loss),
                    batch_size(batch_size) {
            assert(batch_size > 0);
        }
        virtual ~SGDOptimizer() = default;
        inline void fit(typename Base::Net& net) {
            timestep = 0;
            target_net_ptr = &net;
            _fit(net.get_all_unique_params());
        }
    protected:
        inline Scalar _train(typename Base::Net& net, typename Base::Provider& training_prov, std::size_t epoch,
                bool verbose) {
            assert(target_net_ptr == &net);
            Scalar obj_loss = 0;
            Scalar reg_loss = 0;
            std::size_t instances = 0;
            std::size_t updates = 0;
            // Get all the optimizable parameters.
            std::vector<Parameters<Scalar>*> params_vec = net.get_all_unique_params();
            // Perform an entire training epoch.
            while (training_prov.has_more()) {
                DataPair<Scalar,Rank,Sequential> data_pair = training_prov.get_data(batch_size);
                instances += data_pair.first.dimension(0);
                typename Base::Data out = net.propagate(std::move(data_pair.first), true);
                obj_loss += Base::loss->function(out, data_pair.second).sum();
                /* Divide the gradient by the batch size to decouple the learning rate and the batch
                * size hyper-parameters. Use the nominal batch size as the denominator even if the
                * actual batch size is different (in case the number of samples in the data set is
                * not divisible by the batch size and the last batch of the epoch contains fewer
                * instances than the others) to make sure that the magnitude of the gradient is
                * proportional to the batch size (just like its 'accuracy' is). */
                net.backpropagate(Base::loss->d_function(std::move(out),
                        std::move(data_pair.second)) / (Scalar) batch_size);
                // Update the values of the parameters.
                std::size_t i = 0;
                for (auto params_ptr : params_vec) {
                    if (!params_ptr->are_optimizable() || params_ptr->are_frozen())
                        continue;
                    reg_loss += params_ptr->get_regularization_penalty();
                    params_ptr->regularize();
                }
                _update_params(params_vec, epoch - 1, timestep);
                ++updates;
                ++timestep;
                // Reset the gradients of the optimizable parameters.
                for (auto params_ptr : params_vec)
                    params_ptr->reset_grad();
            }
            Scalar mean_obj_loss = obj_loss / instances;
            Scalar mean_reg_loss = reg_loss / updates;
            if (verbose) {
                std::cout << std::left << std::setw(20) << "\ttraining obj loss: " << std::right <<
                        std::to_string(mean_obj_loss) << std::endl;
                std::cout << std::left << std::setw(20) << "\ttraining reg loss: " << std::right <<
                        std::to_string(mean_reg_loss) << std::endl;
            }
            return mean_obj_loss + mean_reg_loss;
        }
        inline Scalar _test(typename Base::Net& net, typename Base::Provider& test_prov, std::size_t epoch,
                bool verbose) {
            assert(target_net_ptr == &net);
            Scalar obj_loss = 0;
            Scalar instances = 0;
            std::vector<Parameters<Scalar>*> params_vec = net.get_all_unique_params();
            // Perform an entire test epoch.
            while (test_prov.has_more()) {
                DataPair<Scalar,Rank,Sequential> data_pair = test_prov.get_data(batch_size);
                instances += data_pair.first.dimension(0);
                obj_loss += Base::loss->function(net.infer(std::move(data_pair.first)),
                        std::move(data_pair.second)).sum();
            }
            Scalar mean_obj_loss = obj_loss / instances;
            Scalar reg_loss = 0;
            for (auto params_ptr : params_vec) {
                if (params_ptr->are_optimizable() && !params_ptr->are_frozen())
                    reg_loss += params_ptr->get_regularization_penalty();
            }
            if (verbose) {
                std::cout << std::left << std::setw(20) << "\ttest obj loss: " << std::right <<
                        std::to_string(mean_obj_loss) << std::endl;
                std::cout << std::left << std::setw(20) << "\ttest reg loss: " << std::right <<
                        std::to_string(reg_loss) << std::endl;
            }
            return mean_obj_loss + reg_loss;
        }
        /**
        * It fits the optimizer to the provided parameters.
        *
        * @param params_vec All the unique parameters of the network.
        */
        virtual void _fit(const std::vector<Parameters<Scalar>*>& params_vec) = 0;
        /**
        * It updates the parameters based on their gradients after back-propagation.
        *
        * @param params_vec All the unique parameters of the network.
        * @param epoch The index of the epoch.
        * @param timestep The index of the update (time step).
        */
        virtual void _update_params(const std::vector<Parameters<Scalar>*>& params_vec, std::size_t epoch, std::size_t timestep) = 0;
        const std::size_t batch_size;
    private:
        std::size_t timestep;
        const typename Base::Net* target_net_ptr;
    };
    /**
    * A class template for the Adam optimization algorithm.
    *
    * \see https://arxiv.org/abs/1412.6980
    */
    template<typename Scalar, std::size_t Rank, bool Sequential>
    class AdamOptimizer : public SGDOptimizer<Scalar,Rank,Sequential> {
        typedef Optimizer<Scalar,Rank,Sequential> Root;
    public:
        /**
        * @param loss A shared pointer to the loss function to use.
        * @param batch_size The batch size to use for training and testing. It is expected to
        * be greater than 0.
        * @param learning_rate The learning rate (a.k.a. step size) to use. It is expected to
        * be greater than 0.
        * @param l1_decay The decay rate of the accumulated parameter gradients. It is expected
        * to be in the range [0,1]. The greater it is, the faster the accumulated gradients
        * decay.
        * @param l2_decay The decay rate of the accumulated squared parameter gradients. It is
        * expected to be in the range [0,1]. The greater it is, the faster the accumulated
        * squared gradients decay.
        * @param epsilon A small constant used to maintain numerical stability.
        */
        inline AdamOptimizer(LossSharedPtr<Scalar,Rank,Sequential> loss, std::size_t batch_size = 1,
                Scalar learning_rate = 1e-3, Scalar l1_decay = 1e-1, Scalar l2_decay = 1e-3,
                Scalar epsilon = NumericUtils<Scalar>::EPSILON2) :
                    SGDOptimizer<Scalar,Rank,Sequential>::SGDOptimizer(loss, batch_size),
                    learning_rate(learning_rate),
                    l1_decay(l1_decay),
                    l2_decay(l2_decay),
                    epsilon(epsilon) {
            assert(learning_rate > 0);
            assert(l1_decay >= 0 && l1_decay <= 1);
            assert(l2_decay >= 0 && l2_decay <= 1);
            assert(epsilon > 0);
        }
        virtual ~AdamOptimizer() = default;
    protected:
        inline void _fit(const std::vector<Parameters<Scalar>*>& params_vec) {
            pgn_vec = std::vector<ParamsGradNorms>();
            for (auto params_ptr : params_vec) {
                if (!params_ptr->are_optimizable() || params_ptr->are_frozen())
                    continue;
                ParamsGradNorms pgn;
                pgn.params_grad_l1 = Matrix<Scalar>::Zero(params_ptr->get_rows(), params_ptr->get_cols());
                pgn.params_grad_l2 = Matrix<Scalar>::Zero(params_ptr->get_rows(), params_ptr->get_cols());
                pgn_vec.push_back(std::move(pgn));
            }
        }
        inline void _update_params(const std::vector<Parameters<Scalar>*>& params_vec, std::size_t epoch, std::size_t timestep) {
            Scalar l1_corr = (Scalar) 1 / (1 - pow(1 - l1_decay, timestep + 1) + epsilon);
            Scalar l2_corr = (Scalar) 1 / (1 - pow(1 - l2_decay, timestep + 1) + epsilon);
            std::size_t i = 0;
            for (auto params_ptr : params_vec) {
                if (!params_ptr->are_optimizable() || params_ptr->are_frozen())
                    continue;
                ParamsGradNorms& grad_norms = pgn_vec[i++];
                const Matrix<Scalar>& params_grad = params_ptr->get_grad();
                grad_norms.params_grad_l1 = grad_norms.params_grad_l1 * (1 - l1_decay) + params_grad * l1_decay;
                grad_norms.params_grad_l2 = grad_norms.params_grad_l2 * (1 - l2_decay) +
                        params_grad.cwiseProduct(params_grad) * l2_decay;
                params_ptr->set_values(params_ptr->get_values() -
                        ((grad_norms.params_grad_l1 * (learning_rate * l1_corr)).array() /
                        ((grad_norms.params_grad_l2 * l2_corr).array() + epsilon).sqrt()).matrix());
            }
        }
        const Scalar learning_rate, l1_decay, l2_decay, epsilon;
        /**
        * A struct containing the moving averages of the first and second norms of the parameter gradients
        * of a layer.
        */
        struct ParamsGradNorms {
            Matrix<Scalar> params_grad_l1, params_grad_l2;
        };
        std::vector<ParamsGradNorms> pgn_vec;
    };

    /**
    * A class template for the AMSGrad optimization algorithm.
    *
    * \see https://openreview.net/pdf?id=ryQu7f-RZ
    */
    template<typename Scalar, std::size_t Rank, bool Sequential>
    class AMSGradOptimizer : public AdamOptimizer<Scalar,Rank,Sequential> {
        typedef Optimizer<Scalar,Rank,Sequential> Root;
        typedef AdamOptimizer<Scalar,Rank,Sequential> Base;
    public:
        /**
        * @param loss A shared pointer to the loss function to use.
        * @param batch_size The batch size to use for training and testing. It is expected to
        * be greater than 0.
        * @param learning_rate The learning rate (a.k.a. step size) to use. It is expected to
        * be greater than 0.
        * @param l1_decay The decay rate of the accumulated parameter gradients. It is expected
        * to be in the range [0,1]. The greater it is, the faster the accumulated gradients
        * decay.
        * @param l2_decay The decay rate of the accumulated squared parameter gradients. It is
        * expected to be in the range [0,1]. The greater it is, the faster the accumulated
        * squared gradients decay.
        * @param epsilon A small constant used to maintain numerical stability.
        */
        inline AMSGradOptimizer(LossSharedPtr<Scalar,Rank,Sequential> loss, std::size_t batch_size = 1,
                Scalar learning_rate = 1e-3, Scalar l1_decay = 1e-1, Scalar l2_decay = 1e-3,
                Scalar epsilon = NumericUtils<Scalar>::EPSILON2) :
                    Base::AdamOptimizer(loss, batch_size, learning_rate, l1_decay, l2_decay, epsilon) { }
    protected:
        inline void _fit(const std::vector<Parameters<Scalar>*>& params_vec) {
            Base::_fit(params_vec);
            params_grad_l2_max = std::vector<Matrix<Scalar>>(Base::pgn_vec.size());
            for (std::size_t i = 0; i < params_grad_l2_max.size(); ++i)
                params_grad_l2_max[i] = Base::pgn_vec[i].params_grad_l2;
        }
        inline void _update_params(const std::vector<Parameters<Scalar>*>& params_vec, std::size_t epoch, std::size_t timestep) {
            std::size_t i = 0;
            for (auto params_ptr : params_vec) {
                if (!params_ptr->are_optimizable() || params_ptr->are_frozen())
                    continue;
                typename Base::ParamsGradNorms& grad_norms = Base::pgn_vec[i];
                const Matrix<Scalar>& params_grad = params_ptr->get_grad();
                grad_norms.params_grad_l1 = grad_norms.params_grad_l1 * (1 - Base::l1_decay) +
                        params_grad * Base::l1_decay;
                grad_norms.params_grad_l2 = grad_norms.params_grad_l2 * (1 - Base::l2_decay) +
                                params_grad.cwiseProduct(params_grad) * Base::l2_decay;
                Matrix<Scalar>& params_grad_l2_max = this->params_grad_l2_max[i++];
                params_grad_l2_max = grad_norms.params_grad_l2.cwiseMax(params_grad_l2_max);
                params_ptr->set_values(params_ptr->get_values() -
                        (grad_norms.params_grad_l1.array() * Base::learning_rate /
                        (params_grad_l2_max.array() + Base::epsilon).sqrt()).matrix());
            }
        }
    private:
        std::vector<Matrix<Scalar>> params_grad_l2_max;
    };

    /**
    * A class template for the ADADELTA optimization algorithm.
    *
    * \see https://arxiv.org/abs/1212.5701
    */
    template<typename Scalar, std::size_t Rank, bool Sequential>
    class AdaDeltaOptimizer : public SGDOptimizer<Scalar,Rank,Sequential> {
        typedef Optimizer<Scalar,Rank,Sequential> Root;
    public:
        /**
        * @param loss A shared pointer to the loss function to use.
        * @param batch_size The batch size to use for training and testing. It is expected to
        * be greater than 0.
        * @param decay The decay rate of the accelerated accumulated parameter gradients.
        * It is expected to be in the range [0,1]. The greater it is, the faster the accumulated
        * gradients decay.
        * @param epsilon A small constant used to maintain numerical stability.
        */
        inline AdaDeltaOptimizer(LossSharedPtr<Scalar,Rank,Sequential> loss, std::size_t batch_size = 1,
                Scalar decay = 5e-2, Scalar epsilon = NumericUtils<Scalar>::EPSILON2) :
                    SGDOptimizer<Scalar,Rank,Sequential>::SGDOptimizer(loss, batch_size),
                    decay(decay),
                    epsilon(epsilon) {
            assert(decay >= 0 && decay <= 1);
            assert(epsilon > 0);
        }
    protected:
        inline void _fit(const std::vector<Parameters<Scalar>*>& params_vec) {
            pgus_vec = std::vector<ParamsGradAndUpdateSqrs>();
            for (auto params_ptr : params_vec) {
                if (!params_ptr->are_optimizable() || params_ptr->are_frozen())
                    continue;
                ParamsGradAndUpdateSqrs pgus;
                pgus.params_grad = Matrix<Scalar>::Zero(params_ptr->get_rows(), params_ptr->get_cols());
                pgus.params_update = Matrix<Scalar>::Zero(params_ptr->get_rows(), params_ptr->get_cols());
                pgus_vec.push_back(std::move(pgus));
            }
        }
        inline void _update_params(const std::vector<Parameters<Scalar>*>& params_vec, std::size_t epoch, std::size_t timestep) {
            std::size_t i = 0;
            for (auto params_ptr : params_vec) {
                if (!params_ptr->are_optimizable() || params_ptr->are_frozen())
                    continue;
                ParamsGradAndUpdateSqrs& pgus = pgus_vec[i++];
                const Matrix<Scalar>& params_grad = params_ptr->get_grad();
                pgus.params_grad = pgus.params_grad * (1 - decay) + params_grad.cwiseProduct(params_grad) * decay;
                Matrix<Scalar> weight_updates = -params_grad.array() * (pgus.params_update.array() + epsilon).sqrt() /
                        (pgus.params_grad.array() + epsilon).sqrt();
                params_ptr->set_values(params_ptr->get_values() + weight_updates);
                pgus.params_update = pgus.params_update * (1 - decay) +
                        weight_updates.cwiseProduct(weight_updates) * decay;
            }
        }
        const Scalar decay, epsilon;
        /**
        * A struct containing the moving averages of the squared gradients and squared updates of a layer.
        */
        struct ParamsGradAndUpdateSqrs {
            Matrix<Scalar> params_grad, params_update;
        };
        std::vector<ParamsGradAndUpdateSqrs> pgus_vec;
    };


    /**
    * A class template for the AdaGrad optimization algorithm.
    *
    * \see http://www.jmlr.org/papers/volume12/duchi11a/duchi11a.pdf
    */
    template<typename Scalar, std::size_t Rank, bool Sequential>
    class AdaGradOptimizer : public SGDOptimizer<Scalar,Rank,Sequential> {
        typedef Optimizer<Scalar,Rank,Sequential> Root;
    public:
        /**
        * @param loss A shared pointer to the loss function to use.
        * @param batch_size The batch size to use for training and testing. It is expected to
        * be greater than 0.
        * @param learning_rate The learning rate (a.k.a. step size) to use. It is expected to
        * be greater than 0.
        * @param epsilon A small constant used to maintain numerical stability.
        */
        inline AdaGradOptimizer(LossSharedPtr<Scalar,Rank,Sequential> loss, std::size_t batch_size = 1,
                Scalar learning_rate = 1e-2, Scalar epsilon = NumericUtils<Scalar>::EPSILON2) :
                    SGDOptimizer<Scalar,Rank,Sequential>::SGDOptimizer(loss, batch_size),
                    learning_rate(learning_rate),
                    epsilon(epsilon) {
            assert(learning_rate > 0);
            assert(epsilon > 0);
        }
        virtual ~AdaGradOptimizer() = default;
    protected:
        inline void _fit(const std::vector<Parameters<Scalar>*>& params_vec) {
            params_grad_sqrs_vec = std::vector<Matrix<Scalar>>();
            for (auto params_ptr : params_vec) {
                if (params_ptr->are_optimizable() && !params_ptr->are_frozen())
                    params_grad_sqrs_vec.push_back(Matrix<Scalar>::Zero(params_ptr->get_rows(), params_ptr->get_cols()));
            }
        }
        inline void _update_params(const std::vector<Parameters<Scalar>*>& params_vec, std::size_t epoch, std::size_t timestep) {
            std::size_t i = 0;
            for (auto params_ptr : params_vec) {
                if (!params_ptr->are_optimizable() || params_ptr->are_frozen())
                    continue;
                const Matrix<Scalar>& params_grad = params_ptr->get_grad();
                Matrix<Scalar>& params_grad_sqrs = params_grad_sqrs_vec[i++];
                _update_acc_params_grad_sqrs(params_grad_sqrs, params_grad);
                params_ptr->set_values(params_ptr->get_values() - (params_grad.array() * learning_rate /
                        (params_grad_sqrs.array().sqrt() + epsilon)).matrix());
            }
        }
        /**
        * It updates the accumulated squared parameter gradients.
        *
        * @param acc_params_grad_sqrs The accumulated squared parameter gradients.
        * @param params_grad The new parameter gradients.
        */
        inline virtual void _update_acc_params_grad_sqrs(Matrix<Scalar>& acc_params_grad_sqrs,
                const Matrix<Scalar>& params_grad) {
            // Accumulate the squares of the gradients.
            acc_params_grad_sqrs += params_grad.cwiseProduct(params_grad);
        }
        const Scalar learning_rate, epsilon;
        std::vector<Matrix<Scalar>> params_grad_sqrs_vec;
    };

    /**
    * A class template for the AdaMax optimization algorithm.
    *
    * \see https://arxiv.org/abs/1412.6980
    */
    template<typename Scalar, std::size_t Rank, bool Sequential>
    class AdaMaxOptimizer : public AdamOptimizer<Scalar,Rank,Sequential> {
        typedef Optimizer<Scalar,Rank,Sequential> Root;
        typedef AdamOptimizer<Scalar,Rank,Sequential> Base;
    public:
        /**
        * @param loss A shared pointer to the loss function to use.
        * @param batch_size The batch size to use for training and testing. It is expected to
        * be greater than 0.
        * @param learning_rate The learning rate (a.k.a. step size) to use. It is expected to
        * be greater than 0.
        * @param l1_decay The decay rate of the accumulated parameter gradients. It is expected
        * to be in the range [0,1]. The greater it is, the faster the accumulated gradients
        * decay.
        * @param l2_decay The decay rate of the accumulated squared parameter gradients. It is
        * expected to be in the range [0,1]. The greater it is, the faster the accumulated
        * squared gradients decay.
        * @param epsilon A small constant used to maintain numerical stability.
        */
        inline AdaMaxOptimizer(LossSharedPtr<Scalar,Rank,Sequential> loss, std::size_t batch_size = 1,
                Scalar learning_rate = 1e-3, Scalar l1_decay = 1e-1, Scalar l2_decay = 1e-3,
                Scalar epsilon = NumericUtils<Scalar>::EPSILON2) :
                    Base::AdamOptimizer(loss, batch_size, learning_rate, l1_decay, l2_decay, epsilon) { }
    protected:
        inline void _update_params(const std::vector<Parameters<Scalar>*>& params_vec, std::size_t epoch, std::size_t timestep) {
            Scalar l1_corr = (Scalar) 1 / (1 - pow(1 - Base::l1_decay, timestep + 1) + Base::epsilon);
            std::size_t i = 0;
            for (auto params_ptr : params_vec) {
                if (!params_ptr->are_optimizable() || params_ptr->are_frozen())
                    continue;
                typename Base::ParamsGradNorms& grad_norms = Base::pgn_vec[i++];
                const Matrix<Scalar>& params_grad = params_ptr->get_grad();
                grad_norms.params_grad_l1 = grad_norms.params_grad_l1 * (1 - Base::l1_decay) +
                        params_grad * Base::l1_decay;
                grad_norms.params_grad_l2 = (grad_norms.params_grad_l2 * (1 - Base::l2_decay))
                        .cwiseMax(params_grad.cwiseAbs());
                params_ptr->set_values(params_ptr->get_values() -
                        ((grad_norms.params_grad_l1 * (Base::learning_rate * l1_corr)).array() /
                        (grad_norms.params_grad_l2.array() + Base::epsilon)).matrix());
            }
        }
    };

    /**
    * A class template for a momentum accelerated SGD optimizer.
    */
    template<typename Scalar, std::size_t Rank, bool Sequential>
    class MomentumSGDOptimizer : public SGDOptimizer<Scalar,Rank,Sequential> {
        typedef Optimizer<Scalar,Rank,Sequential> Root;
        typedef SGDOptimizer<Scalar,Rank,Sequential> Base;
    public:
        /**
        * @param loss A shared pointer to the loss function to use.
        * @param batch_size The batch size to use for training and testing. It is expected to
        * be greater than 0.
        * @param init_learning_rate The initial learning rate (a.k.a. step size) to use. It
        * is expected to be greater than 0.
        * @param annealing_rate The rate at which the learning rate is to be annealed. It is
        * expected to be greater than or equal to 0. The greater it is, the faster the learning
        * rate decreases.
        * @param momentum The momentum rate to use. The greater the momentum, the lesser the
        * effect of newer gradients. It is expected to be greater than 0 and less than 1.
        */
        inline MomentumSGDOptimizer(LossSharedPtr<Scalar,Rank,Sequential> loss, std::size_t batch_size = 1,
                Scalar init_learning_rate = 1e-3, Scalar annealing_rate = 1e-3, Scalar momentum = .9) :
                    SGDOptimizer<Scalar,Rank,Sequential>::SGDOptimizer(loss, batch_size),
                    init_learning_rate(init_learning_rate),
                    annealing_rate(annealing_rate),
                    momentum(momentum) {
            assert(init_learning_rate > 0);
            assert(annealing_rate >= 0);
            assert(momentum > 0 && momentum < 1);
        }
        virtual ~MomentumSGDOptimizer() = default;
    protected:
        inline void _fit(const std::vector<Parameters<Scalar>*>& params_vec) {
            params_grad_vec = std::vector<Matrix<Scalar>>();
            for (auto params_ptr : params_vec) {
                if (params_ptr->are_optimizable() && !params_ptr->are_frozen())
                    params_grad_vec.push_back(Matrix<Scalar>::Zero(params_ptr->get_rows(), params_ptr->get_cols()));
            }
        }
        inline void _update_params(const std::vector<Parameters<Scalar>*>& params_vec, std::size_t epoch, std::size_t timestep) {
            Scalar learning_rate = calculate_learning_rate(epoch);
            std::size_t i = 0;
            for (auto params_ptr : params_vec) {
                if (!params_ptr->are_optimizable() || params_ptr->are_frozen())
                    continue;
                params_grad_vec[i] = params_grad_vec[i] * momentum + params_ptr->get_grad() * learning_rate;
                params_ptr->set_values(params_ptr->get_values() - params_grad_vec[i++]);
            }
        }
        /**
        * It calculates the annealed learning rate as a function of the epoch index.
        *
        * @param epoch The epoch index.
        * @return The learning rate to use.
        */
        Scalar calculate_learning_rate(std::size_t epoch) {
            return init_learning_rate / (1 + annealing_rate * epoch);
        }
        const Scalar init_learning_rate, annealing_rate, momentum;
        std::vector<Matrix<Scalar>> params_grad_vec;
    };

    /**
    * A class template for the Nesterov accelerated Adam (Nadam) optimization algorithm.
    *
    * \see http://cs229.stanford.edu/proj2015/054_report.pdf
    */
    template<typename Scalar, std::size_t Rank, bool Sequential>
    class NadamOptimizer : public AdamOptimizer<Scalar,Rank,Sequential> {
        typedef Optimizer<Scalar,Rank,Sequential> Root;
        typedef AdamOptimizer<Scalar,Rank,Sequential> Base;
    public:
        /**
        * @param loss A shared pointer to the loss function to use.
        * @param batch_size The batch size to use for training and testing. It is expected to
        * be greater than 0.
        * @param learning_rate The learning rate (a.k.a. step size) to use. It is expected to
        * be greater than 0.
        * @param l1_decay The decay rate of the accumulated parameter gradients. It is expected
        * to be in the range [0,1]. The greater it is, the faster the accumulated gradients
        * decay.
        * @param l2_decay The decay rate of the accumulated squared parameter gradients. It is
        * expected to be in the range [0,1]. The greater it is, the faster the accumulated
        * squared gradients decay.
        * @param epsilon A small constant used to maintain numerical stability.
        */
        inline NadamOptimizer(LossSharedPtr<Scalar,Rank,Sequential> loss, std::size_t batch_size = 1,
                Scalar learning_rate = 1e-3, Scalar l1_decay = 1e-1, Scalar l2_decay = 1e-3,
                Scalar epsilon = NumericUtils<Scalar>::EPSILON2) :
                    Base::AdamOptimizer(loss, batch_size, learning_rate, l1_decay, l2_decay, epsilon) { }
    protected:
        inline void _update_params(const std::vector<Parameters<Scalar>*>& params_vec, std::size_t epoch, std::size_t timestep) {
            Scalar l1_corr = (Scalar) 1 / (1 - pow(1 - Base::l1_decay, timestep + 1) + Base::epsilon);
            Scalar l1_next_corr = (Scalar) 1 / (1 - pow(1 - Base::l1_decay, timestep + 2) + Base::epsilon);
            Scalar l2_corr = (Scalar) 1 / (1 - pow(1 - Base::l2_decay, timestep + 1) + Base::epsilon);
            std::size_t i = 0;
            for (auto params_ptr : params_vec) {
                if (!params_ptr->are_optimizable() || params_ptr->are_frozen())
                    continue;
                typename Base::ParamsGradNorms& grad_norms = Base::pgn_vec[i++];
                const Matrix<Scalar>& params_grad = params_ptr->get_grad();
                grad_norms.params_grad_l1 = grad_norms.params_grad_l1 * (1 - Base::l1_decay) +
                        params_grad * Base::l1_decay;
                grad_norms.params_grad_l2 = grad_norms.params_grad_l2 * (1 - Base::l2_decay) +
                        params_grad.cwiseProduct(params_grad) * Base::l2_decay;
                params_ptr->set_values(params_ptr->get_values() -
                        ((params_grad * (Base::l1_decay * l1_corr) + grad_norms.params_grad_l1 *
                        ((1 - Base::l1_decay) * l1_next_corr)).array() * Base::learning_rate /
                        ((grad_norms.params_grad_l2 * l2_corr).array() + Base::epsilon).sqrt()).matrix());
            }
        }
    };    

    /**
    * A class template for Nesterov momentum accelerated SGD optimizers.
    *
    * \see https://arxiv.org/abs/1212.0901
    */
    template<typename Scalar, std::size_t Rank, bool Sequential>
    class NesterovMomentumSGDOptimizer : public MomentumSGDOptimizer<Scalar,Rank,Sequential> {
        typedef Optimizer<Scalar,Rank,Sequential> Root;
        typedef MomentumSGDOptimizer<Scalar,Rank,Sequential> Base;
    public:
        /**
        * @param loss A shared pointer to the loss function to use.
        * @param batch_size The batch size to use for training and testing. It is expected to
        * be greater than 0.
        * @param init_learning_rate The initial learning rate (a.k.a. step size) to use. It is
        * expected to be greater than 0.
        * @param annealing_rate The rate at which the learning rate is to be annealed. It is
        * expected to be greater than or equal to 0. The greater it is, the faster the learning
        * rate decreases.
        * @param momentum The momentum rate to use. The greater the momentum, the lesser the
        * effect of newer gradients. It is expected to be greater than 0 and less than 1.
        */
        inline NesterovMomentumSGDOptimizer(LossSharedPtr<Scalar,Rank,Sequential> loss,
                std::size_t batch_size = 1, Scalar init_learning_rate = 1e-3, Scalar annealing_rate = 1e-3,
                Scalar momentum = .9) :
                    Base::MomentumSGDOptimizer(loss, batch_size, init_learning_rate, annealing_rate,
                            momentum) { };
    protected:
        inline void _update_params(const std::vector<Parameters<Scalar>*>& params_vec, std::size_t epoch, std::size_t timestep) {
            Scalar learning_rate = Base::calculate_learning_rate(epoch);
            std::size_t i = 0;
            for (auto params_ptr : params_vec) {
                if (!params_ptr->are_optimizable() || params_ptr->are_frozen())
                    continue;
                Matrix<Scalar> old_acc_params_grad = Base::params_grad_vec[i];
                Base::params_grad_vec[i] = old_acc_params_grad * Base::momentum -
                        params_ptr->get_grad() * learning_rate;
                params_ptr->set_values(params_ptr->get_values() + old_acc_params_grad * -Base::momentum +
                        Base::params_grad_vec[i++] * (1 + Base::momentum));
            }
        }
    };    

    /**
    * A class template for the RMSProp optimizer.
    *
    * \see https://www.cs.toronto.edu/~tijmen/csc321/slides/lecture_slides_lec6.pdf
    */
    template<typename Scalar, std::size_t Rank, bool Sequential>
    class RMSPropOptimizer : public AdaGradOptimizer<Scalar,Rank,Sequential> {
    public:
        /**
        * @param loss A shared pointer to the loss function to use.
        * @param batch_size The batch size to use for training and testing. It is expected to
        * be greater than 0.
        * @param learning_rate The learning rate (a.k.a. step size) to use. It is expected to
        * be greater than 0.
        * @param l2_decay The decay rate of the accumulated squared parameter gradients.
        * It is expected to be in the range [0,1]. The greater it is, the faster the accumulated
        * gradients decay.
        * @param epsilon A small constant used to maintain numerical stability.
        */
        inline RMSPropOptimizer(LossSharedPtr<Scalar,Rank,Sequential> loss, std::size_t batch_size = 1,
                Scalar learning_rate = 1e-3, Scalar l2_decay = 1e-1, Scalar epsilon = NumericUtils<Scalar>::EPSILON2) :
                    AdaGradOptimizer<Scalar,Rank,Sequential>::AdaGradOptimizer(loss, batch_size,
                            learning_rate, epsilon),
                    l2_decay(l2_decay) {
            assert(l2_decay >= 0 && l2_decay <= 1);
        }
    protected:
        inline void _update_acc_params_grad_sqrs(Matrix<Scalar>& acc_params_grad_sqrs,
                const Matrix<Scalar>& params_grad) {
            acc_params_grad_sqrs = acc_params_grad_sqrs * (1 - l2_decay) +
                    params_grad.cwiseProduct(params_grad) * l2_decay;
        }
        const Scalar l2_decay;
    };

    /**
    * A class template for a vanilla SGD optimizer.
    */
    template<typename Scalar, std::size_t Rank, bool Sequential>
    class VanillaSGDOptimizer : public SGDOptimizer<Scalar,Rank,Sequential> {
        typedef Optimizer<Scalar,Rank,Sequential> Root;
        typedef SGDOptimizer<Scalar,Rank,Sequential> Base;
    public:
        /**
        * @param loss A shared pointer to the loss function to use.
        * @param batch_size The batch size to use for training and testing. It is expected to
        * be greater than 0.
        * @param learning_rate The learning rate (a.k.a. step size) to use. It is expected to
        * be greater than 0.
        */
        inline VanillaSGDOptimizer(LossSharedPtr<Scalar,Rank,Sequential> loss, std::size_t batch_size = 1,
                Scalar learning_rate = 1e-3) :
                    SGDOptimizer<Scalar,Rank,Sequential>::SGDOptimizer(loss, batch_size),
                    learning_rate(learning_rate) {
            assert(learning_rate > 0);
        }
    protected:
        inline void _fit(const std::vector<Parameters<Scalar>*>& params_vec) { }
        inline void _update_params(const std::vector<Parameters<Scalar>*>& params_vec, std::size_t epoch, std::size_t timestep) {
            for (auto params_ptr : params_vec) {
                if (params_ptr->are_optimizable() && !params_ptr->are_frozen())
                    params_ptr->set_values(params_ptr->get_values() - params_ptr->get_grad() * learning_rate);
            }
        }
        const Scalar learning_rate;
    };

    /**
    * A class template for the elastic net regularization penalty which is a combination of
    * the L1 and L2 regularization penalties.
    *
    * \f$P = \lambda_1 \sum\limits_{i = 1}^n \left|w_i\right| + \frac{\lambda_2}{2} \sum\limits_{i = 1}^n w_i^2\f$
    *
    * \see http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.124.4696
    */
    template<typename Scalar>
    class ElasticNetParameterRegularization : public ParameterRegularization<Scalar> {
    public:
        /**
        * @param l1_lambda The constant by which the L1 penalty is to be scaled.
        * @param l2_lambda The constant by which the L2 penalty is to be scaled.
        */
        inline ElasticNetParameterRegularization(Scalar l1_lambda = 1e-2, Scalar l2_lambda = 1e-2) :
                l1_lambda(l1_lambda),
                l2_lambda(l2_lambda) { }
        inline Scalar function(const Matrix<Scalar>& params) const {
            return params.array().abs().sum() * l1_lambda + params.squaredNorm() * ((Scalar) .5 * l2_lambda);
        }
        inline Matrix<Scalar> d_function(const Matrix<Scalar>& params) const {
            return params.unaryExpr([this](Scalar e) { return e >= 0 ? l1_lambda : -l1_lambda; }) + params * l2_lambda;
        }
    private:
        const Scalar l1_lambda;
        const Scalar l2_lambda;
    };

    /**
    * A class template for an L1 (first-norm) regularization penalty.
    *
    * \f$P = \lambda \sum\limits_{i = 1}^n \left|w_i\right|\f$
    */
    template<typename Scalar>
    class L1ParameterRegularization : public ParameterRegularization<Scalar> {
    public:
        /**
        * @param lambda The constant by which the penalty is to be scaled.
        */
        inline L1ParameterRegularization(Scalar lambda = 1e-2) :
                lambda(lambda) { }
        inline Scalar function(const Matrix<Scalar>& params) const {
            return params.cwiseAbs().sum() * lambda;
        }
        inline Matrix<Scalar> d_function(const Matrix<Scalar>& params) const {
            return params.unaryExpr([this](Scalar e) { return e >= 0 ? lambda : -lambda; });
        }
    private:
        const Scalar lambda;
    };

    /**
    * A class template for an L2 (second-norm) regularization penalty.
    *
    * \f$P = \frac{\lambda_2}{2} \sum\limits_{i = 1}^n w_i^2\f$
    */
    template<typename Scalar>
    class L2ParameterRegularization : public ParameterRegularization<Scalar> {
    public:
        /**
        * @param lambda The constant by which the penalty is to be scaled.
        */
        inline L2ParameterRegularization(Scalar lambda = 1e-2) :
                lambda(lambda) { }
        inline Scalar function(const Matrix<Scalar>& params) const {
            return params.squaredNorm() * ((Scalar) .5 * lambda);
        }
        inline Matrix<Scalar> d_function(const Matrix<Scalar>& params) const {
            return params * lambda;
        }
    private:
        const Scalar lambda;
    };

    /**
    * A class template for a normalization (mean-subtraction) preprocessor that optionally also
    * standardizes the data (divides it by the standard deviation).
    */
    template<typename Scalar, std::size_t Rank, bool Standardize = false, bool PerLastRank = (Rank == 3)>
    class NormalizationPreprocessor : public Preprocessor<Scalar,Rank,false> {
    public:
        /**
        * @param epsilon A small constant used to maintain numerical stability.
        */
        inline NormalizationPreprocessor(Scalar epsilon = NumericUtils<Scalar>::EPSILON2) :
                epsilon(epsilon) {
            assert(epsilon > 0 && "epsilon must be greater than 0");
        }
        virtual ~NormalizationPreprocessor() = default;
        inline virtual void fit(const Tensor<Scalar,Rank + 1>& data) {
            auto rows = data.dimension(0);
            assert(rows > 0);
            dims = (Dimensions<std::size_t,Rank + 1>(data.dimensions()).template demote<>());
            channels = dims(Rank - 1);
            means = Matrix<Scalar>(channels, dims.get_volume() / channels);
            sd = Matrix<Scalar>(means.rows(), means.cols());
            if (channels == 1) {
                Tensor<Scalar,Rank + 1> data_copy(data);
                MatrixMap<Scalar> data_mat(data_copy.data(), rows, data.size() / rows);
                means.row(0) = data_mat.colwise().mean();
                if (Standardize)
                    sd.row(0) = (data_mat.rowwise() - means.row(0)).array().square().colwise().mean().sqrt();
            } else {
                std::array<std::size_t,Rank + 1> offsets;
                offsets.fill(0);
                auto extents = data.dimensions();
                extents[Rank] = 1;
                for (std::size_t i = 0; i < channels; ++i) {
                    offsets[Rank] = i;
                    Tensor<Scalar,Rank + 1> data_slice = data.slice(offsets, extents);
                    MatrixMap<Scalar> data_slice_mat(data_slice.data(), rows, data_slice.size() / rows);
                    means.row(i) = data_slice_mat.colwise().mean();
                    if (Standardize)
                        sd.row(i) = (data_slice_mat.rowwise() - means.row(i)).array().square().colwise().mean().sqrt();
                }
            }
        }
        inline virtual void transform(Tensor<Scalar,Rank + 1>& data) const {
            auto rows = data.dimension(0);
            assert(rows > 0);
            assert((Dimensions<std::size_t,Rank + 1>(data.dimensions()).template demote<>()) == dims &&
                    "mismatched fit and transform input tensor dimensions");
            if (channels == 1) {
                MatrixMap<Scalar> data_mat(data.data(), rows, data.size() / rows);
                data_mat = data_mat.rowwise() - means.row(0);
                if (Standardize)
                    data_mat *= (sd.row(0).array() + epsilon).inverse().matrix().asDiagonal();
            } else {
                std::array<std::size_t,Rank + 1> offsets;
                offsets.fill(0);
                auto extents = data.dimensions();
                extents[Rank] = 1;
                for (std::size_t i = 0; i < channels; ++i) {
                    offsets[Rank] = i;
                    Tensor<Scalar,Rank + 1> data_slice = data.slice(offsets, extents);
                    MatrixMap<Scalar> data_slice_mat(data_slice.data(), rows, data_slice.size() / rows);
                    data_slice_mat = data_slice_mat.rowwise() - means.row(i);
                    if (Standardize)
                        data_slice_mat *= (sd.row(i).array() + epsilon).inverse().matrix().asDiagonal();
                    data.slice(offsets, extents) = std::move(data_slice);
                }
            }
        }
    protected:
        const Scalar epsilon;
        Dimensions<std::size_t,Rank> dims;
        Matrix<Scalar> means;
        Matrix<Scalar> sd;
        std::size_t channels;
    };

    /**
    * Partial template specialization for pre-activation normalization.
    */
    template<typename Scalar, std::size_t Rank, bool Standardize>
    class NormalizationPreprocessor<Scalar,Rank,Standardize,false> : public Preprocessor<Scalar,Rank,false> {
    public:
        /**
        * @param epsilon A small constant used to maintain numerical stability.
        */
        inline NormalizationPreprocessor(Scalar epsilon = NumericUtils<Scalar>::EPSILON2) :
                epsilon(epsilon) {
            assert(epsilon > 0 && "epsilon must be greater than 0");
        }
        virtual ~NormalizationPreprocessor() = default;
        inline virtual void fit(const Tensor<Scalar,Rank + 1>& data) {
            std::size_t rows = data.dimension(0);
            assert(rows > 0);
            dims = (Dimensions<std::size_t,Rank + 1>(data.dimensions()).template demote<>());
            Tensor<Scalar,Rank + 1> data_copy(data);
            MatrixMap<Scalar> data_mat(data_copy.data(), rows, data.size() / rows);
            means = data_mat.colwise().mean();
            if (Standardize)
                sd = (data_mat.rowwise() - means).array().square().colwise().mean().sqrt();
        }
        inline virtual void transform(Tensor<Scalar,Rank + 1>& data) const {
            std::size_t rows = data.dimension(0);
            assert(rows > 0);
            assert((Dimensions<std::size_t,Rank + 1>(data.dimensions()).template demote<>()) == dims &&
                    "mismatched fit and transform input tensor dimensions");
            MatrixMap<Scalar> data_mat(data.data(), rows, data.size() / rows);
            data_mat = data_mat.rowwise() - means;
            if (Standardize)
                data_mat *= (sd.array() + epsilon).inverse().matrix().asDiagonal();
        }
    protected:
        const Scalar epsilon;
        Dimensions<std::size_t,Rank> dims;
        RowVector<Scalar> means;
        RowVector<Scalar> sd;
    };

    /**
    * An abstract base class template for a principal component analysis (PCA) preprocessor that can also
    * standardize and whiten the data.
    */
    template<typename Scalar, std::size_t Rank, bool Standardize, bool Whiten, bool PerLastRank>
    class PCAPreprocessorBase : public NormalizationPreprocessor<Scalar,Rank,Standardize,PerLastRank> {
    protected:
        typedef NormalizationPreprocessor<Scalar,Rank,Standardize,PerLastRank> Base;
        inline PCAPreprocessorBase(Scalar min_rel_var_to_retain, Scalar epsilon) :
                    Base::NormalizationPreprocessor(epsilon),
                    min_rel_var_to_retain(min_rel_var_to_retain),
                    reduce_dims(NumericUtils<Scalar>::decidedly_lesser(min_rel_var_to_retain, (Scalar) 1)) {
            assert(min_rel_var_to_retain > 0 && min_rel_var_to_retain <= 1 &&
                    "the minimum relative variance to be retained must be greater "
                    "then 0 and less than or equal to 1");
        }
        inline void _fit(Tensor<Scalar,Rank + 1> data, int i) {
            std::size_t rows = data.dimension(0);
            MatrixMap<Scalar> data_mat(data.data(), rows, data.size() / rows);
            Matrix<Scalar> normalized_data = data_mat.rowwise() - Base::means.row(i);
            if (Standardize)
                normalized_data *= (Base::sd.row(i).array() + Base::epsilon).inverse().matrix().asDiagonal();
            // Compute the covariance matrix.
            Matrix<Scalar> cov = normalized_data.transpose() * normalized_data / normalized_data.rows();
            // Eigen decomposition.
            EigenSolver<Scalar> eigen_solver(cov);
            // Determine the number of components to retain.
            const ColVector<Scalar>& eigen_values = eigen_solver.eigenvalues();
            int dims_to_retain = 0;
            if (reduce_dims) {
                const Scalar min_var_to_retain = eigen_values.sum() * min_rel_var_to_retain;
                Scalar var = 0;
                for (; dims_to_retain < eigen_values.rows(); ++dims_to_retain) {
                    // The eigen values are in ascending order.
                    var += eigen_values(eigen_values.rows() - (1 + dims_to_retain));
                    if (NumericUtils<Scalar>::decidedly_greater(var, min_var_to_retain))
                        break;
                }
            } else
                dims_to_retain = eigen_values.rows();
            // The eigen vectors are sorted by the magnitude of their corresponding eigen values.
            ed_vec[i].eigen_basis = eigen_solver.eigenvectors().rightCols(dims_to_retain);
            if (Whiten) // The eigen values are only needed if whitening is enabled.
                ed_vec[i].eigen_values = eigen_values.bottomRows(dims_to_retain).transpose();
        }
        inline Tensor<Scalar,Rank + 1> _transform(Tensor<Scalar,Rank + 1> data, int i) const {
            std::size_t rows = data.dimension(0);
            MatrixMap<Scalar> data_mat(data.data(), rows, data.size() / rows);
            if (reduce_dims) {
                Dimensions<std::size_t,Rank + 1> output_dims;
                output_dims(0) = rows;
                Matrix<Scalar> transformed_data_mat = data_mat * ed_vec[i].eigen_basis;
                if (Whiten)
                    transformed_data_mat *= (ed_vec[i].eigen_values.array() + Base::epsilon).sqrt().inverse().matrix().asDiagonal();
                output_dims(1) = transformed_data_mat.cols();
                return TensorMap<Scalar,Rank + 1>(transformed_data_mat.data(), output_dims);
            } else {
                Dimensions<std::size_t,Rank + 1> output_dims = data.dimensions();
                data_mat *= ed_vec[i].eigen_basis;
                if (Whiten)
                    data_mat *= (ed_vec[i].eigen_values.array() + Base::epsilon).sqrt().inverse().matrix().asDiagonal();
                return TensorMap<Scalar,Rank + 1>(data_mat.data(), output_dims);
            }
        }
        Scalar min_rel_var_to_retain;
        bool reduce_dims;
        struct EigenDecomposition {
            Matrix<Scalar> eigen_basis;
            RowVector<Scalar> eigen_values;
        };
        std::vector<EigenDecomposition> ed_vec;
    };

    /**
    * A class template for a PCA preprocessor that can also be used to standardize and whiten data across
    * multiple channels.
    */
    template<typename Scalar, std::size_t Rank, bool Standardize = false, bool Whiten = false,
            bool PerLastRank = (Rank == 3)>
    class PCAPreprocessor : public PCAPreprocessorBase<Scalar,Rank,Standardize,Whiten,PerLastRank> {
        typedef PCAPreprocessorBase<Scalar,Rank,Standardize,Whiten,PerLastRank> Base;
        typedef typename Base::Base Root;
    public:
        /**
        * @param min_rel_var_to_retain The minimum relative variance in the data
        * to retain. It is expected to be within the range (0,1]. If it is 1,
        * the dimensionality of the preprocessed data is guaranteed not to be
        * reduced. If it is less than 1, the data cannot be a multi-channel
        * tensor.
        * @param epsilon A small consant used to maintain numerical stability.
        */
        inline PCAPreprocessor(Scalar min_rel_var_to_retain = 1, Scalar epsilon = NumericUtils<Scalar>::EPSILON2) :
                    Base::PCAPreprocessorBase(min_rel_var_to_retain, epsilon) { }
        inline void fit(const Tensor<Scalar,Rank + 1>& data) {
            assert((!Base::reduce_dims || data.dimension(Rank) == 1) && "cannot reduce the dimensionality of multi-channel data");
            Root::fit(data);
            channels = data.dimension(Rank);
            Base::ed_vec = std::vector<typename Base::EigenDecomposition>(channels);
            if (channels == 1)
                Base::_fit(data, 0);
            else {
                std::array<std::size_t,Rank + 1> offsets;
                offsets.fill(0);
                auto extents = data.dimensions();
                extents[Rank] = 1;
                for (std::size_t i = 0; i < channels; ++i) {
                    offsets[Rank] = i;
                    Tensor<Scalar,Rank + 1> data_slice_i = data.slice(offsets, extents);
                    Base::_fit(std::move(data_slice_i), i);
                }
            }
        }
        inline void transform(Tensor<Scalar,Rank + 1>& data) const {
            Root::transform(data);
            if (Base::reduce_dims || channels == 1)
                data = Base::_transform(std::move(data), 0);
            else {
                std::array<std::size_t,Rank + 1> offsets;
                offsets.fill(0);
                auto extents = data.dimensions();
                extents[Rank] = 1;
                for (std::size_t i = 0; i < channels; ++i) {
                    offsets[Rank] = i;
                    Tensor<Scalar,Rank + 1> data_slice_i = data.slice(offsets, extents);
                    data.slice(offsets, extents) = Base::_transform(std::move(data_slice_i), i);
                }
            }
        }
    private:
        std::size_t channels;
    };

    /**
    * Partial template specialization of the PCA preprocessor for single channel data.
    */
    template<typename Scalar, std::size_t Rank, bool Standardize, bool Whiten>
    class PCAPreprocessor<Scalar,Rank,Standardize,Whiten,false> :
            public PCAPreprocessorBase<Scalar,Rank,Standardize,Whiten,false> {
        typedef PCAPreprocessorBase<Scalar,Rank,Standardize,Whiten,false> Base;
        typedef typename Base::Base Root;
    public:
        /**
        * @param min_rel_var_to_retain The minimum relative variance in the data
        * to retain. It is expected to be within the range (0,1]. If it is 1,
        * the dimensionality of the preprocessed data is guaranteed not to be
        * reduced.
        * @param epsilon A small consant used to maintain numerical stability.
        */
        inline PCAPreprocessor(Scalar min_rel_var_to_retain = 1, Scalar epsilon = NumericUtils<Scalar>::EPSILON2) :
                    Base::PCAPreprocessorBase(min_rel_var_to_retain, epsilon) { }
        inline void fit(const Tensor<Scalar,Rank + 1>& data) {
            Root::fit(data);
            Base::ed_vec = std::vector<typename Base::EigenDecomposition>(1);
            Base::_fit(data, 0);
        }
        inline void transform(Tensor<Scalar,Rank + 1>& data) const {
            Root::transform(data);
            data = Base::_transform(std::move(data), 0);
        }
    };

        /**
    * A utility class for performing gradient checks on neural network and loss function
    * implementations.
    */
    template<typename Scalar, std::size_t Rank, bool Sequential>
    class GradientCheck {
        static_assert(std::is_floating_point<Scalar>::value, "non floating-point scalar type");
        static_assert(Rank > 0 && Rank < 4, "illegal optimizer rank");
        typedef DataProvider<Scalar,Rank,Sequential> Provider;
        typedef NeuralNetwork<Scalar,Rank,Sequential> Net;
        typedef Loss<Scalar,Rank,Sequential> LossType;
    public:
        /**
        * It performs a gradient check to verify the correctness of the neural network and layer implementations.
        * It is recommended to use double precision floating points.
        *
        * @param provider A reference to the data provider to use for the gradient check.
        * @param net A reference to the network on which the gradient check is to be performed.
        * @param loss The loss function to use for the gradient check.
        * @param verbose Whether the analytic and numerical derivatives of the variables should be printed to the
        * standard out stream.
        * @param step_size The step size for numerical differentiation.
        * @param abs_epsilon The maximum acceptable absolute difference between the numerical and analytic
        * gradients.
        * @param rel_epsilon The maximum acceptable relative (to the greater out of the two) difference between
        * the numerical and analytic gradients.
        * @return Whether the gradient check has been passed or failed.
        */
        inline static bool verify_gradients(Provider& provider, Net& net, const LossType& loss, bool verbose = true,
                Scalar step_size = NumericUtils<Scalar>::EPSILON2 / 2,
                Scalar abs_epsilon = NumericUtils<Scalar>::EPSILON2,
                Scalar rel_epsilon = NumericUtils<Scalar>::EPSILON2) {
            assert(net.get_input_dims() == provider.get_obs_dims());
            assert(net.get_output_dims() == provider.get_obj_dims());
            assert(step_size > 0);
            assert(abs_epsilon >= 0 && rel_epsilon > 0);
            bool failure = false;
            DataPair<Scalar,Rank,Sequential> data_pair = provider.get_data(std::numeric_limits<std::size_t>::max());
            std::size_t instances = data_pair.first.dimension(0);
            provider.reset();
            /* As the loss to minimize is the mean of the losses for all the training observations, the gradient to
            * back-propagate is to be divided by the number of observations in the batch. */
            net.backpropagate(loss.d_function(net.propagate(data_pair.first, true),
                    data_pair.second) / (Scalar) instances);
            std::size_t i = 0;
            for (auto params_ptr : net.get_all_unique_params()) {
                if (!params_ptr->are_optimizable())
                    continue;
                Parameters<Scalar>& params = *params_ptr;
                if (verbose) {
                    std::cout << "Parameter Set " << std::setw(3) << std::to_string(i) <<
                            std::string(28, '-') << std::endl;
                }
                /* Add the derivative of the regularization function w.r.t. to the parameters of the layer to the
                * parameters' gradient. */
                params.regularize();
                Matrix<Scalar> params_values = params.get_values();
                const Matrix<Scalar>& params_grad = params.get_grad();
                for (int j = 0; j < params_values.rows(); ++j) {
                    for (int k = 0; k < params_values.cols(); ++k) {
                        if (verbose)
                            std::cout << "\tparam[" << i << "," << j << "," << k << "]:" << std::endl;
                        Scalar ana_grad = params_grad(j,k);
                        if (verbose)
                            std::cout << "\t\tanalytic gradient = " << ana_grad << std::endl;
                        Scalar param = params_values(j,k);
                        params_values(j,k) = param + step_size;
                        params.set_values(params_values);
                        /* Compute the numerical gradients in training mode to ensure that the means and standard
                        * deviations used for batch normalization are the same as those used during the analytic
                        * gradient computation. */
                        Scalar loss_inc = loss.function(net.propagate(data_pair.first, true),
                                data_pair.second).mean();
                        /* Calculate the new regularization penalty as its derivative w.r.t. the layer's
                        * parameters is included in the gradient. */
                        Scalar reg_pen_inc = params.get_regularization_penalty();
                        params_values(j,k) = param - step_size;
                        params.set_values(params_values);
                        Scalar loss_dec = loss.function(net.propagate(data_pair.first, true),
                                data_pair.second).mean();
                        Scalar reg_pen_dec = params.get_regularization_penalty();
                        params_values(j,k) = param;
                        params.set_values(params_values);
                        // Include the regularization penalty as well.
                        Scalar num_grad = (loss_inc + reg_pen_inc - (loss_dec + reg_pen_dec)) / (2 * step_size);
                        if (verbose)
                            std::cout << "\t\tnumerical gradient = " << num_grad;
                        if (!NumericUtils<Scalar>::almost_equal(ana_grad, num_grad, abs_epsilon, rel_epsilon)) {
                            if (verbose)
                                std::cout << " *****FAIL*****";
                            failure = true;
                        }
                        if (verbose)
                            std::cout << std::endl;
                    }
                }
                params.reset_grad();
                ++i;
            }
            // Empty the network caches.
            net.empty_caches();
            return !failure;
        }
    };

}
