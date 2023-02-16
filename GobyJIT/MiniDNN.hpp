#include <Eigen/Core>
#include <vector>
#include <map>
#include <stdexcept>
#include <map>       // std::map
#include <string>    // std::string
#include <sstream>   // std::ostringstream
#include <fstream>   // std::ofstream, std::ifstream
#include <iterator>  // std::ostream_iterator, std::istreambuf_iterator, std::back_inserter
#include <cstdlib>   // atoi
#include <iostream>

#ifdef _WIN32
    #include <direct.h>     // _mkdir
#else
    #include <sys/stat.h> // mkdir
#endif





namespace MiniDNN
{

// Floating-point number type
#ifndef MDNN_SCALAR
typedef double Scalar;
#else
typedef MDNN_SCALAR Scalar;
#endif


///
/// The interface and default implementation of the random number generator (%RNG).
/// The %RNG is used to shuffle data and to initialize parameters of hidden layers.
/// This default implementation is based on the public domain code by Ray Gardner
/// <http://stjarnhimlen.se/snippets/rg_rand.c>.
///
class RNG
{
    private:
        const unsigned int m_a;     // multiplier
        const unsigned long m_max;  // 2^31 - 1
        long m_rand;

        inline long next_long_rand(long seed)
        {
            unsigned long lo, hi;
            lo = m_a * (long)(seed & 0xFFFF);
            hi = m_a * (long)((unsigned long)seed >> 16);
            lo += (hi & 0x7FFF) << 16;

            if (lo > m_max)
            {
                lo &= m_max;
                ++lo;
            }

            lo += hi >> 15;

            if (lo > m_max)
            {
                lo &= m_max;
                ++lo;
            }

            return (long)lo;
        }
    public:
        RNG(unsigned long init_seed) :
            m_a(16807),
            m_max(2147483647L),
            m_rand(init_seed ? (init_seed & m_max) : 1)
        {}

        virtual ~RNG() {}

        virtual void seed(unsigned long seed)
        {
            m_rand = (seed ? (seed & m_max) : 1);
        }

        virtual Scalar rand()
        {
            m_rand = next_long_rand(m_rand);
            return Scalar(m_rand) / Scalar(m_max);
        }
};

namespace internal
{
    ///
    /// Convert a number to a string in C++98
    ///
    /// \tparam NumberType     Type of the number
    /// \param  num            The number to be converted
    /// \return                An std::string containing the number
    ///
    template <class NumberType>
    inline std::string to_string(const NumberType& num)
    {
        std::ostringstream convert;
        convert << num;
        return convert.str();
    }

    ///
    /// Create a directory
    ///
    /// \param dir     Name of the directory to be created
    /// \return        \c true if the directory is successfully created
    ///
    inline bool create_directory(const std::string& dir)
    {
    #ifdef _WIN32
        return 0 == _mkdir(dir.c_str());
    #else
        return 0 == mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
    #endif
    }

    ///
    /// Write an std::vector<Scalar> vector to file
    ///
    /// \param vec          The vector to be written to file
    /// \param filename     The filename of the output
    ///
    inline void write_vector_to_file(
        const std::vector<Scalar>& vec, const std::string& filename
    )
    {
        std::ofstream ofs(filename.c_str(), std::ios::out | std::ios::binary);
        if (ofs.fail())
            throw std::runtime_error("Error while opening file");

        std::ostream_iterator<char> osi(ofs);
        const char* begin_byte = reinterpret_cast<const char*>(&vec[0]);
        const char* end_byte = begin_byte + vec.size() * sizeof(Scalar);
        std::copy(begin_byte, end_byte, osi);
    }

    ///
    /// Write the parameters of an NN model to file
    ///
    /// \param folder       The folder where the parameter files are stored
    /// \param filename     The filename prefix of the parameter files
    /// \param params       The parameters of the NN model
    ///
    inline void write_parameters(
        const std::string& folder, const std::string& filename,
        const std::vector< std::vector< Scalar> >& params
    )
    {
        const int nfiles = params.size();
        for (int i = 0; i < nfiles; i++)
        {
            write_vector_to_file(params[i], folder + "/" + filename + to_string(i));
        }
    }

    ///
    /// Read in an std::vector<Scalar> vector from file
    ///
    /// \param filename     The filename of the input
    /// \return             The vector that has been read
    ///
    inline std::vector<Scalar> read_vector_from_file(const std::string& filename)
    {

        std::ifstream ifs(filename.c_str(), std::ios::in | std::ifstream::binary);
        if (ifs.fail())
            throw std::runtime_error("Error while opening file");

        std::vector<char> buffer;
        std::istreambuf_iterator<char> iter(ifs);
        std::istreambuf_iterator<char> end;
        std::copy(iter, end, std::back_inserter(buffer));
        std::vector<Scalar> vec(buffer.size() / sizeof(Scalar));
        std::copy(&buffer[0], &buffer[0] + buffer.size(), reinterpret_cast<char*>(&vec[0]));
        return vec;
    }

    ///
    /// Read in parameters of an NN model from file
    ///
    /// \param folder       The folder where the parameter files are stored
    /// \param filename     The filename prefix of the parameter files
    /// \param nlayer       Number of layers in the NN model
    /// \return             A vector of vectors that contains the NN parameters
    ///
    inline std::vector< std::vector< Scalar> > read_parameters(
        const std::string& folder, const std::string& filename, int nlayer
    )
    {
        std::vector< std::vector< Scalar> > params;
        params.reserve(nlayer);

        for (int i = 0; i < nlayer; i++)
        {
            params.push_back(read_vector_from_file(folder + "/" + filename + to_string(i)));
        }

        return params;
    }

    ///
    /// Write a map object to file
    ///
    /// \param filename     The filename of the output
    /// \param map          The map object to be exported
    ///
    inline void write_map(const std::string& filename, const std::map<std::string, int>& map)
    {
        if (map.empty())
            return;

        std::ofstream ofs(filename.c_str(), std::ios::out);
        if (ofs.fail())
            throw std::runtime_error("Error while opening file");

        for (std::map<std::string, int>::const_iterator it = map.begin(); it != map.end(); it++)
        {
            ofs << it->first << "=" << it->second << std::endl;
        }
    }

    ///
    /// Read in a map object from file
    ///
    /// \param filename     The filename of the input
    /// \param map          The output map object
    ///
    inline void read_map(const std::string& filename, std::map<std::string, int>& map)
    {
        std::ifstream ifs(filename.c_str(), std::ios::in);
        if (ifs.fail())
            throw std::runtime_error("Error while opening file");

        map.clear();
        std::string buf;
        while (std::getline(ifs, buf))
        {
            std::size_t sep = buf.find('=');
            if (sep == std::string::npos)
                throw std::invalid_argument("File format error");

            std::string key = buf.substr(0, sep);
            std::string value = buf.substr(sep + 1, buf.length() - sep - 1);
            map[key] = atoi(value.c_str());
        }
    }

    // Shuffle the integer array
    inline void shuffle(int* arr, const int n, RNG& rng)
    {
        for (int i = n - 1; i > 0; i--)
        {
            // A random non-negative integer <= i
            const int j = int(rng.rand() * (i + 1));
            // Swap arr[i] and arr[j]
            const int tmp = arr[i];
            arr[i] = arr[j];
            arr[j] = tmp;
        }
    }

    template <typename DerivedX, typename DerivedY, typename XType, typename YType>
    inline int create_shuffled_batches(
        const Eigen::MatrixBase<DerivedX>& x, const Eigen::MatrixBase<DerivedY>& y,
        int batch_size, RNG& rng,
        std::vector<XType>& x_batches, std::vector<YType>& y_batches
    )
    {
        const int nobs = x.cols();
        const int dimx = x.rows();
        const int dimy = y.rows();

        if (y.cols() != nobs)
        {
            throw std::invalid_argument("Input X and Y have different number of observations");
        }

        // Randomly shuffle the IDs
        Eigen::VectorXi id = Eigen::VectorXi::LinSpaced(nobs, 0, nobs - 1);
        shuffle(id.data(), id.size(), rng);

        // Compute batch size
        if (batch_size > nobs)
        {
            batch_size = nobs;
        }

        const int nbatch = (nobs - 1) / batch_size + 1;
        const int last_batch_size = nobs - (nbatch - 1) * batch_size;
        // Create shuffled data
        x_batches.clear();
        y_batches.clear();
        x_batches.reserve(nbatch);
        y_batches.reserve(nbatch);

        for (int i = 0; i < nbatch; i++)
        {
            const int bsize = (i == nbatch - 1) ? last_batch_size : batch_size;
            x_batches.push_back(XType(dimx, bsize));
            y_batches.push_back(YType(dimy, bsize));
            // Copy data
            const int offset = i * batch_size;

            for (int j = 0; j < bsize; j++)
            {
                x_batches[i].col(j).noalias() = x.col(id[offset + j]);
                y_batches[i].col(j).noalias() = y.col(id[offset + j]);
            }
        }

        return nbatch;
    }

    // Fill array with N(mu, sigma^2) random numbers
    inline void set_normal_random(Scalar* arr, const int n, RNG& rng,
                                const Scalar& mu = Scalar(0),
                                const Scalar& sigma = Scalar(1))
    {
        // For simplicity we use Box-Muller transform to generate normal random variates
        const double two_pi = 6.283185307179586476925286766559;

        for (int i = 0; i < n - 1; i += 2)
        {
            const double t1 = sigma * std::sqrt(-2 * std::log(rng.rand()));
            const double t2 = two_pi * rng.rand();
            arr[i]     = t1 * std::cos(t2) + mu;
            arr[i + 1] = t1 * std::sin(t2) + mu;
        }

        if (n % 2 == 1)
        {
            const double t1 = sigma * std::sqrt(-2 * std::log(rng.rand()));
            const double t2 = two_pi * rng.rand();
            arr[n - 1] = t1 * std::cos(t2) + mu;
        }
    }
    
    // We assume the following memory layout:
    // There are 'n_obs' images, each with 'in_channels' channels
    // Each channel has 'channel_rows' rows and 'channel_cols' columns
    // The data starts from 'src'
    // If 'image_outer_loop == true', the data first iterates on channels,
    // and then images:
    /*
    * ###############################################################
    * #           #           #           #           #
    * # channel 1 # channel 2 # channel 3 # channel 1 # ...
    * #           #           #           #           #
    * ###############################################################
    * |<------------ image 1 ------------>|<------------ image 2 ----
    */
    // If 'image_outer_loop == false', the layout looks like below:
    /*
    * ###############################################################
    * #           #           #           #           #
    * #  image 1  #  image 2  #  image 3  #  image 1  # ...
    * # channel 1 # channel 1 # channel 1 # channel 2 #
    * #           #           #           #           #
    * ###############################################################
    * |<----------- channel 1 ----------->|<----------- channel 2 ----
    */
    //
    // Then we assume there are 'out_channels' output channels, so in total
    // we have 'in_channels * out_channels' filters for each image
    // Each filter has 'filter_rows' rows and 'filter_cols' columns
    // Filters start from 'filter_data'. The layout looks like below, with each
    // block consisting of 'filter_rows * filter_cols' elements
    /*
    * #########################################################################
    * #               #               #               #               #
    * # out channel 1 # out channel 2 # out channel 3 # out channel 1 # ...
    * #               #               #               #               #
    * #########################################################################
    * |<---------------- in channel 1 --------------->|<---------------- in channel 2 ----
    */
    //
    // Convolution results from different input channels are summed up to produce the
    // result for each output channel
    // Convolution results for different output channels are concatenated to preoduce
    // the result for each image
    //
    // The final result is written to the memory pointed by 'dest', with a similar
    // layout to 'src' in the 'image_outer_loop == true' case
    //
    // Memory efficient convolution (MEC)
    // Algorithm is based on https://arxiv.org/abs/1706.06873
    //
    // First define a simple structure to store the various dimensions of convolution
    struct ConvDims
    {
        // Input parameters
        const int in_channels;
        const int out_channels;
        const int channel_rows;
        const int channel_cols;
        const int filter_rows;
        const int filter_cols;
        // Image dimension -- one observation with all channels
        const int img_rows;
        const int img_cols;
        // Dimension of the convolution result for each output channel
        const int conv_rows;
        const int conv_cols;

        ConvDims(
            const int in_channels_, const int out_channels_,
            const int channel_rows_, const int channel_cols_,
            const int filter_rows_, const int filter_cols_
        ) :
            in_channels(in_channels_), out_channels(out_channels_),
            channel_rows(channel_rows_), channel_cols(channel_cols_),
            filter_rows(filter_rows_), filter_cols(filter_cols_),
            img_rows(channel_rows_), img_cols(in_channels_ * channel_cols_),
            conv_rows(channel_rows_ - filter_rows_ + 1),
            conv_cols(channel_cols_ - filter_cols_ + 1)
        {}
    };
    // Transform original matrix to "lower" form as described in the MEC paper
    // I feel that it is better called the "flat" form
    //
    // Helper function to "flatten" source images
    // 'flat_mat' will be overwritten
    // We focus on one channel, and let 'stride' be the distance between two images
    inline void flatten_mat(
        const ConvDims& dim, const Scalar* src, const int stride, const int n_obs,
        Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& flat_mat
    )
    {
        // Number of bytes in the segment that will be copied at one time
        const int& segment_size = dim.filter_rows;
        const std::size_t copy_bytes = sizeof(Scalar) * segment_size;
        Scalar* writer = flat_mat.data();
        const int channel_size = dim.channel_rows * dim.channel_cols;

        for (int i = 0; i < n_obs; i++, src += stride)
        {
            const Scalar* reader_row = src;
            const Scalar* const reader_row_end = src + dim.conv_rows;

            for (; reader_row < reader_row_end; reader_row++)
            {
                const Scalar* reader = reader_row;
                const Scalar* const reader_end = reader + channel_size;

                for (; reader < reader_end; reader += dim.channel_rows, writer += segment_size)
                {
                    std::memcpy(writer, reader, copy_bytes);
                }
            }
        }
    }
    // A special matrix product. We select a window from 'mat1' and calculates its product with 'mat2',
    // and progressively move the window to the right
    inline void moving_product(
        const int step,
        const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>&
        mat1,
        Eigen::Map< const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> >& mat2,
        Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>& res
    )
    {
        const int row1 = mat1.rows();
        const int col1 = mat1.cols();
        const int row2 = mat2.rows();
        const int col2 = mat2.cols();
        const int col_end = col1 - row2;
        int res_start_col = 0;

        for (int left_end = 0; left_end <= col_end;
                left_end += step, res_start_col += col2)
        {
            res.block(0, res_start_col, row1, col2).noalias() += mat1.block(0, left_end,
                    row1, row2) * mat2;
        }
    }
    // The main convolution function using the "valid" rule
    inline void convolve_valid(
        const ConvDims& dim,
        const Scalar* src, const bool image_outer_loop, const int n_obs,
        const Scalar* filter_data,
        Scalar* dest)
    {
        typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Matrix;
        typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
        RMatrix;
        typedef Eigen::Map<const Matrix> ConstMapMat;
        // Flat matrix
        const int flat_rows = dim.conv_rows * n_obs;
        const int flat_cols = dim.filter_rows * dim.channel_cols;
        const int channel_size = dim.channel_rows * dim.channel_cols;
        // Distance between two images
        const int img_stride = image_outer_loop ? (dim.img_rows * dim.img_cols) :
                            channel_size;
        // Distance between two channels
        const int channel_stride = image_outer_loop ? channel_size :
                                (channel_size * n_obs);
        RMatrix flat_mat(flat_rows, flat_cols);
        // Convolution results
        const int& res_rows = flat_rows;
        const int res_cols = dim.conv_cols * dim.out_channels;
        Matrix res = Matrix::Zero(res_rows, res_cols);
        const int& step = dim.filter_rows;
        const int filter_size = dim.filter_rows * dim.filter_cols;
        const int filter_stride = filter_size * dim.out_channels;

        for (int i = 0; i < dim.in_channels;
                i++, src += channel_stride, filter_data += filter_stride)
        {
            // Flatten source image
            flatten_mat(dim, src, img_stride, n_obs, flat_mat);
            // Compute the convolution result
            ConstMapMat filter(filter_data, filter_size, dim.out_channels);
            moving_product(step, flat_mat, filter, res);
        }

        // The layout of 'res' is very complicated
        /*
        * obs0_out0[0, 0] obs0_out1[0, 0] obs0_out2[0, 0] obs0_out0[0, 1] obs0_out1[0, 1] obs0_out2[0, 1] ...
        * obs0_out0[1, 0] obs0_out1[1, 0] obs0_out2[1, 0] obs0_out0[1, 1] obs0_out1[1, 1] obs0_out2[1, 1] ...
        * obs0_out0[2, 0] obs0_out1[2, 0] obs0_out2[2, 0] obs0_out0[2, 1] obs0_out1[2, 1] obs0_out2[2, 1] ...
        * obs1_out0[0, 0] obs1_out1[0, 0] obs1_out2[0, 0] obs1_out0[0, 1] obs1_out1[0, 1] obs1_out2[0, 1] ...
        * obs1_out0[1, 0] obs1_out1[1, 0] obs1_out2[1, 0] obs1_out0[1, 1] obs1_out1[1, 1] obs1_out2[1, 1] ...
        * obs1_out0[2, 0] obs1_out1[2, 0] obs1_out2[2, 0] obs1_out0[2, 1] obs1_out1[2, 1] obs1_out2[2, 1] ...
        * ...
        *
        */
        // obs<k>_out<l> means the convolution result of the k-th image on the l-th output channel
        // [i, j] gives the matrix indices
        // The destination has the layout
        /*
        * obs0_out0[0, 0] obs0_out0[0, 1] obs0_out0[0, 2] obs0_out1[0, 0] obs0_out1[0, 1] obs0_out1[0, 2] ...
        * obs0_out0[1, 0] obs0_out0[1, 1] obs0_out0[1, 2] obs0_out1[1, 0] obs0_out1[1, 1] obs0_out1[1, 2] ...
        * obs0_out0[2, 0] obs0_out0[2, 1] obs0_out0[2, 2] obs0_out1[2, 0] obs0_out1[2, 1] obs0_out1[2, 2] ...
        *
        */
        // which in a larger scale looks like
        // [obs0_out0 obs0_out1 obs0_out2 obs1_out0 obs1_out1 obs1_out2 obs2_out0 ...]
        // Copy data to destination
        // dest[a, b] corresponds to obs<k>_out<l>[i, j]
        // where k = b / (conv_cols * out_channels),
        //       l = (b % (conv_cols * out_channels)) / conv_cols
        //       i = a,
        //       j = b % conv_cols
        // and then obs<k>_out<l>[i, j] corresponds to res[c, d]
        // where c = k * conv_rows + i,
        //       d = j * out_channels + l
        const int dest_rows = dim.conv_rows;
        const int dest_cols = res_cols * n_obs;
        const Scalar* res_data = res.data();
        const std::size_t copy_bytes = sizeof(Scalar) * dest_rows;

        for (int b = 0; b < dest_cols; b++, dest += dest_rows)
        {
            const int k = b / res_cols;
            const int l = (b % res_cols) / dim.conv_cols;
            const int j = b % dim.conv_cols;
            const int d = j * dim.out_channels + l;
            const int res_col_head = d * res_rows;
            std::memcpy(dest, res_data + res_col_head + k * dim.conv_rows, copy_bytes);
        }
    }



    // The moving_product() function for the "full" rule
    inline void moving_product(
        const int padding, const int step,
        const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>&
        mat1,
        const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>& mat2,
        Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>& res
    )
    {
        const int row1 = mat1.rows();
        const int col1 = mat1.cols();
        const int row2 = mat2.rows();
        const int col2 = mat2.cols();
        int res_start_col = 0;
        // Left padding
        int left_end = -padding;
        int right_end = step;

        for (; left_end < 0
                && right_end <= col1;
                left_end += step, right_end += step, res_start_col += col2)
        {
            res.block(0, res_start_col, row1, col2).noalias() += mat1.leftCols(right_end) *
                    mat2.bottomRows(right_end);
        }

        // Main part
        for (; right_end <= col1;
                left_end += step, right_end += step, res_start_col += col2)
        {
            res.block(0, res_start_col, row1, col2).noalias() += mat1.block(0, left_end,
                    row1, row2) * mat2;
        }

        // Right padding
        for (; left_end < col1; left_end += step, res_start_col += col2)
        {
            if (left_end <= 0)
            {
                res.block(0, res_start_col, row1, col2).noalias() += mat1 * mat2.block(0,
                        -left_end, col1, row2);
            }
            else
            {
                const int overlap = col1 - left_end;
                res.block(0, res_start_col, row1, col2).noalias() += mat1.rightCols(overlap) *
                        mat2.topRows(overlap);
            }
        }
    }
    // The main convolution function for the "full" rule
    inline void convolve_full(
        const ConvDims& dim,
        const Scalar* src, const int n_obs, const Scalar* filter_data,
        Scalar* dest)
    {
        typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Matrix;
        typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
        RMatrix;
        typedef Eigen::Map<const Matrix> ConstMapMat;
        // Padding sizes
        const int padding_top = dim.filter_rows - 1;
        const int padding_left = dim.filter_cols - 1;
        // Dimension of convolution result using "full" rule
        const int conv_rows = dim.channel_rows + padding_top;
        const int conv_cols = dim.channel_cols + padding_left;
        // Add (top and bottom) padding to source images
        const int pad_rows = dim.img_rows + padding_top * 2;
        const int pad_cols = dim.img_cols * n_obs;
        Matrix pad_mat(pad_rows, pad_cols);
        ConstMapMat src_mat(src, dim.img_rows, pad_cols);
        pad_mat.topRows(padding_top).setZero();
        pad_mat.bottomRows(padding_top).setZero();
        pad_mat.block(padding_top, 0, dim.img_rows, pad_cols).noalias() = src_mat;
        src = pad_mat.data();
        ConvDims pad_dim(dim.in_channels, dim.out_channels, pad_rows, dim.channel_cols,
                        dim.filter_rows, dim.filter_cols);
        // Flat matrix
        const int flat_rows = conv_rows * n_obs;
        const int flat_cols = dim.filter_rows * dim.channel_cols;
        const int img_stride = pad_rows * dim.img_cols;
        const int channel_stride = pad_rows * dim.channel_cols;
        RMatrix flat_mat(flat_rows, flat_cols);
        // The processing of filters are different from the "valid" rule in two ways:
        // 1. The layout of input channels and output channels are switched
        // 2. The filters need to be rotated, which is equivalent to reversing the vector of each filter
        // We also separate filters that belong to different input channels
        std::vector<Matrix> filters_in(dim.in_channels);
        const int filter_size = dim.filter_rows * dim.filter_cols;
        const int nfilter = dim.in_channels * dim.out_channels;

        for (int i = 0; i < dim.in_channels; i++)
        {
            filters_in[i].resize(filter_size, dim.out_channels);
        }

        const Scalar* reader = filter_data;

        for (int i = 0; i < nfilter; i++, reader += filter_size)
        {
            Scalar* writer = filters_in[i % dim.in_channels].data() +
                            (i / dim.in_channels) * filter_size;
            std::reverse_copy(reader, reader + filter_size, writer);
        }

        // Convolution results
        const int& res_rows = flat_rows;
        const int res_cols = conv_cols * dim.out_channels;
        Matrix res = Matrix::Zero(res_rows, res_cols);
        const int& step = dim.filter_rows;
        const int filter_padding = padding_left * dim.filter_rows;

        for (int i = 0; i < dim.in_channels; i++, src += channel_stride)
        {
            // Flatten source image
            flatten_mat(pad_dim, src, img_stride, n_obs, flat_mat);
            // Compute the convolution result
            moving_product(filter_padding, step, flat_mat, filters_in[i], res);
        }

        // Copy results to destination
        const int& dest_rows = conv_rows;
        const int  dest_cols = res_cols * n_obs;
        const Scalar* res_data = res.data();
        const std::size_t copy_bytes = sizeof(Scalar) * dest_rows;

        for (int b = 0; b < dest_cols; b++, dest += dest_rows)
        {
            const int k = b / res_cols;
            const int l = (b % res_cols) / conv_cols;
            const int j = b % conv_cols;
            const int d = j * dim.out_channels + l;
            const int res_col_head = d * res_rows;
            std::memcpy(dest, res_data + res_col_head + k * conv_rows, copy_bytes);
        }
    }

    // Enumerations for hidden layers
    enum LAYER_ENUM
    {
        FULLY_CONNECTED = 0,
        CONVOLUTIONAL,
        MAX_POOLING
    };

    // Convert a hidden layer type string to an integer
    inline int layer_id(const std::string& type)
    {
        if (type == "FullyConnected")
            return FULLY_CONNECTED;
        if (type == "Convolutional")
            return CONVOLUTIONAL;
        if (type == "MaxPooling")
            return MAX_POOLING;

        throw std::invalid_argument("[function layer_id]: Layer is not of a known type");
        return -1;
    }

    // Enumerations for activation functions
    enum ACTIVATION_ENUM
    {
        IDENTITY = 0,
        RELU,
        SIGMOID,
        SOFTMAX,
        TANH,
        MISH
    };

    // Convert an activation type string to an integer
    inline int activation_id(const std::string& type)
    {
        if (type == "Identity")
            return IDENTITY;
        if (type == "ReLU")
            return RELU;
        if (type == "Sigmoid")
            return SIGMOID;
        if (type == "Softmax")
            return SOFTMAX;
        if (type == "Tanh")
            return TANH;
        if (type == "Mish")
            return MISH;

        throw std::invalid_argument("[function activation_id]: Activation is not of a known type");
        return -1;
    }

    // Enumerations for output layers
    enum OUTPUT_ENUM
    {
        REGRESSION_MSE = 0,
        BINARY_CLASS_ENTROPY,
        MULTI_CLASS_ENTROPY
    };

    // Convert an output layer type string to an integer
    inline int output_id(const std::string& type)
    {
        if (type == "RegressionMSE")
            return REGRESSION_MSE;
        if (type == "MultiClassEntropy")
            return BINARY_CLASS_ENTROPY;
        if (type == "BinaryClassEntropy")
            return MULTI_CLASS_ENTROPY;

        throw std::invalid_argument("[function output_id]: Output is not of a known type");
        return -1;
    }

    // Find the location of the maximum element in x[0], x[1], ..., x[n-1]
    // Special cases for small n using recursive template
    // N is assumed to be >= 2
    template <int N>
    inline int find_max(const Scalar* x)
    {
        const int loc = find_max < N - 1 > (x);
        return (x[N - 1] > x[loc]) ? (N - 1) : loc;
    }

    template <>
    inline int find_max<2>(const Scalar* x)
    {
        return int(x[1] > x[0]);
    }

    // n is assumed be >= 2
    inline int find_max(const Scalar* x, const int n)
    {
        switch (n)
        {
            case 2:
                return find_max<2>(x);

            case 3:
                return find_max<3>(x);

            case 4:
                return find_max<4>(x);

            case 5:
                return find_max<5>(x);
        }

        int loc = find_max<6>(x);

        for (int i = 6; i < n; i++)
        {
            loc = (x[i] > x[loc]) ? i : loc;
        }

        return loc;
    }

    // Find the maximum element in the block x[0:(nrow-1), 0:(ncol-1)]
    // col_stride is the distance between x[0, 0] and x[0, 1]
    // Special cases for small n
    inline Scalar find_block_max(const Scalar* x, const int nrow, const int ncol,
                                const int col_stride, int& loc)
    {
        // Max element in the first column
        loc = find_max(x, nrow);
        Scalar val = x[loc];
        // 2nd column
        x += col_stride;
        int loc_next = find_max(x, nrow);
        Scalar val_next = x[loc_next];

        if (val_next > val)
        {
            loc = col_stride + loc_next;
            val = val_next;
        }

        if (ncol == 2)
        {
            return val;
        }

        // 3rd column
        x += col_stride;
        loc_next = find_max(x, nrow);
        val_next = x[loc_next];

        if (val_next > val)
        {
            loc = 2 * col_stride + loc_next;
            val = val_next;
        }

        if (ncol == 3)
        {
            return val;
        }

        // 4th column
        x += col_stride;
        loc_next = find_max(x, nrow);
        val_next = x[loc_next];

        if (val_next > val)
        {
            loc = 3 * col_stride + loc_next;
            val = val_next;
        }

        if (ncol == 4)
        {
            return val;
        }

        // 5th column
        x += col_stride;
        loc_next = find_max(x, nrow);
        val_next = x[loc_next];

        if (val_next > val)
        {
            loc = 4 * col_stride + loc_next;
            val = val_next;
        }

        if (ncol == 5)
        {
            return val;
        }

        // Other columns
        for (int i = 5; i < ncol; i++)
        {
            x += col_stride;
            loc_next = find_max(x, nrow);
            val_next = x[loc_next];

            if (val_next > val)
            {
                loc = i * col_stride + loc_next;
                val = val_next;
            }
        }

        return val;
    }

} // namespace internal


    ///
    /// \defgroup Optimizers Optimization Algorithms
    ///

    ///
    /// \ingroup Optimizers
    ///
    /// The interface of optimization algorithms
    ///
    class Optimizer
    {
        protected:
            typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
            typedef Vector::ConstAlignedMapType ConstAlignedMapVec;
            typedef Vector::AlignedMapType AlignedMapVec;

        public:
            virtual ~Optimizer() {}

            ///
            /// Reset the optimizer to clear all historical information
            ///
            virtual void reset() {};

            ///
            /// Update the parameter vector using its gradient
            ///
            /// It is assumed that the memory addresses of `dvec` and `vec` do not
            /// change during the training process. This is used to implement optimization
            /// algorithms that have "memories". See the AdaGrad algorithm for an example.
            ///
            /// \param dvec The gradient of the parameter. Read-only
            /// \param vec  On entering, the current parameter vector. On exit, the
            ///             updated parameters.
            ///
            virtual void update(ConstAlignedMapVec& dvec, AlignedMapVec& vec) = 0;
    };


    ///
    /// \defgroup Layers Hidden Layers
    ///

    ///
    /// \ingroup Layers
    ///
    /// The interface of hidden layers in a neural network. It defines some common
    /// operations of hidden layers such as initialization, forward and backward
    /// propogation, and also functions to get/set parameters of the layer.
    ///
    class Layer
    {
        protected:
            typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Matrix;
            typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
            typedef std::map<std::string, int> MetaInfo;

            const int m_in_size;  // Size of input units
            const int m_out_size; // Size of output units

        public:
            ///
            /// Constructor.
            ///
            /// \param in_size  Number of input units of this hidden Layer. It must be
            ///                 equal to the number of output units of the previous layer.
            /// \param out_size Number of output units of this hidden layer. It must be
            ///                 equal to the number of input units of the next layer.
            ///
            Layer(const int in_size, const int out_size) :
                m_in_size(in_size), m_out_size(out_size)
            {}

            ///
            /// Virtual destructor.
            ///
            virtual ~Layer() {}

            ///
            /// Get the number of input units of this hidden layer.
            ///
            int in_size() const
            {
                return m_in_size;
            }
            ///
            /// Get the number of output units of this hidden layer.
            ///
            int out_size() const
            {
                return m_out_size;
            }

            ///
            /// Initialize layer parameters using \f$N(\mu, \sigma^2)\f$ distribution.
            ///
            /// \param mu    Mean of the normal distribution.
            /// \param sigma Standard deviation of the normal distribution.
            /// \param rng   The random number generator of type RNG.
            virtual void init(const Scalar& mu, const Scalar& sigma, RNG& rng) = 0;

            ///
            /// Initialize layer parameters without arguments. It is used when the layer is
            /// read from file. This function will typically set the sizes of member
            /// matrices and vectors.
            ///
            virtual void init() = 0;

            ///
            /// Compute the output of this layer.
            ///
            /// The purpose of this function is to let the hidden layer compute information
            /// that will be passed to the next layer as the input. The concrete behavior
            /// of this function is subject to the implementation, with the only
            /// requirement that after calling this function, the Layer::output() member
            /// function will return a reference to the output values.
            ///
            /// \param prev_layer_data The output of previous layer, which is also the
            ///                        input of this layer. `prev_layer_data` should have
            ///                        `in_size` rows as in the constructor, and each
            ///                        column of `prev_layer_data` is an observation.
            ///
            virtual void forward(const Matrix& prev_layer_data) = 0;

            ///
            /// Obtain the output values of this layer
            ///
            /// This function is assumed to be called after Layer::forward() in each iteration.
            /// The output are the values of output hidden units after applying activation function.
            /// The main usage of this function is to provide the `prev_layer_data` parameter
            /// in Layer::forward() of the next layer.
            ///
            /// \return A reference to the matrix that contains the output values. The
            ///         matrix should have `out_size` rows as in the constructor,
            ///         and have number of columns equal to that of `prev_layer_data` in the
            ///         Layer::forward() function. Each column represents an observation.
            ///
            virtual const Matrix& output() const = 0;

            ///
            /// Compute the gradients of parameters and input units using back-propagation
            ///
            /// The purpose of this function is to compute the gradient of input units,
            /// which can be retrieved by Layer::backprop_data(), and the gradient of
            /// layer parameters, which could later be used by the Layer::update() function.
            ///
            /// \param prev_layer_data The output of previous layer, which is also the
            ///                        input of this layer. `prev_layer_data` should have
            ///                        `in_size` rows as in the constructor, and each
            ///                        column of `prev_layer_data` is an observation.
            /// \param next_layer_data The gradients of the input units of the next layer,
            ///                        which is also the gradients of the output units of
            ///                        this layer. `next_layer_data` should have
            ///                        `out_size` rows as in the constructor, and the same
            ///                        number of columns as `prev_layer_data`.
            ///
            virtual void backprop(const Matrix& prev_layer_data,
                                const Matrix& next_layer_data) = 0;

            ///
            /// Obtain the gradient of input units of this layer
            ///
            /// This function provides the `next_layer_data` parameter in Layer::backprop()
            /// of the previous layer, since the derivative of the input of this layer is also the derivative
            /// of the output of previous layer.
            ///
            virtual const Matrix& backprop_data() const = 0;

            ///
            /// Update parameters after back-propagation
            ///
            /// \param opt The optimization algorithm to be used. See the Optimizer class.
            ///
            virtual void update(Optimizer& opt) = 0;

            ///
            /// Get serialized values of parameters
            ///
            virtual std::vector<Scalar> get_parameters() const = 0;
            ///
            /// Set the values of layer parameters from serialized data
            ///
            virtual void set_parameters(const std::vector<Scalar>& param) {};

            ///
            /// Get serialized values of the gradient of parameters
            ///
            virtual std::vector<Scalar> get_derivatives() const = 0;

            ///
            /// Return the layer type. It is used to export the NN model.
            ///
            virtual std::string layer_type() const = 0;

            ///
            /// Return the activation layer type. It is used to export the NN model.
            ///
            virtual std::string activation_type() const = 0;

            ///
            /// Fill in the meta information of this layer, such as layer type, input
            /// and output sizes, etc. It is used to export layer to file.
            ///
            /// \param map   A key-value map that contains the meta information of the NN model.
            /// \param index The index of this layer in the NN model. It is used to generate
            ///              the key. For example, the layer may insert {"Layer1": 2},
            ///              where 1 is the index, "Layer1" is the key, and 2 is the value.
            ///
            virtual void fill_meta_info(MetaInfo& map, int index) const = 0;
    };

    ///
    /// \defgroup Outputs Output Layers
    ///

    ///
    /// \ingroup Outputs
    ///
    /// The interface of the output layer of a neural network model. The output
    /// layer is a special layer that associates the last hidden layer with the
    /// target response variable.
    ///
    class Output
    {
        protected:
            typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Matrix;
            typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
            typedef Eigen::RowVectorXi IntegerVector;

        public:
            virtual ~Output() {}

            // Check the format of target data, e.g. in classification problems the
            // target data should be binary (either 0 or 1)
            virtual void check_target_data(const Matrix& target) {}

            // Another type of target data where each element is a class label
            // This version may not be sensible for regression tasks, so by default
            // we raise an exception
            virtual void check_target_data(const IntegerVector& target)
            {
                throw std::invalid_argument("[class Output]: This output type cannot take class labels as target data");
            }

            // A combination of the forward stage and the back-propagation stage for the output layer
            // The computed derivative of the input should be stored in this layer, and can be retrieved by
            // the backprop_data() function
            virtual void evaluate(const Matrix& prev_layer_data, const Matrix& target) = 0;

            // Another type of target data where each element is a class label
            // This version may not be sensible for regression tasks, so by default
            // we raise an exception
            virtual void evaluate(const Matrix& prev_layer_data,
                                const IntegerVector& target)
            {
                throw std::invalid_argument("[class Output]: This output type cannot take class labels as target data");
            }

            // The derivative of the input of this layer, which is also the derivative
            // of the output of previous layer
            virtual const Matrix& backprop_data() const = 0;

            // Return the loss function value after the evaluation
            // This function can be assumed to be called after evaluate(), so that it can make use of the
            // intermediate result to save some computation
            virtual Scalar loss() const = 0;

            // Return the output layer type. It is used to export the NN model.
            virtual std::string output_type() const = 0;
    };


    class Network;

    ///
    /// \defgroup Callbacks Callback Functions
    ///

    ///
    /// \ingroup Callbacks
    ///
    /// The interface and default implementation of the callback function during
    /// model fitting. The purpose of this class is to allow users printing some
    /// messages in each epoch or mini-batch training, for example the time spent,
    /// the loss function values, etc.
    ///
    /// This default implementation is a silent version of the callback function
    /// that basically does nothing. See the VerboseCallback class for a verbose
    /// version that prints the loss function value in each mini-batch.
    ///
    class Callback
    {
        protected:
            typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Matrix;
            typedef Eigen::RowVectorXi IntegerVector;

        public:
            // Public members that will be set by the network during the training process
            int m_nbatch;   // Number of total batches
            int m_batch_id; // The index for the current mini-batch (0, 1, ..., m_nbatch-1)
            int m_nepoch;   // Total number of epochs (one run on the whole data set) in the training process
            int m_epoch_id; // The index for the current epoch (0, 1, ..., m_nepoch-1)

            Callback() :
                m_nbatch(0), m_batch_id(0), m_nepoch(0), m_epoch_id(0)
            {}

            virtual ~Callback() {}

            // Before training a mini-batch
            virtual void pre_training_batch(const Network* net, const Matrix& x,
                                            const Matrix& y) {}
            virtual void pre_training_batch(const Network* net, const Matrix& x,
                                            const IntegerVector& y) {}

            // After a mini-batch is trained
            virtual void post_training_batch(const Network* net, const Matrix& x,
                                            const Matrix& y) {}
            virtual void post_training_batch(const Network* net, const Matrix& x,
                                            const IntegerVector& y) {}
    };
    
    ///
    /// \defgroup Activations Activation Functions
    ///

    ///
    /// \ingroup Activations
    ///
    /// The identity activation function
    ///
    class Identity
    {
        private:
            typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Matrix;

        public:
            // a = activation(z) = z
            // Z = [z1, ..., zn], A = [a1, ..., an], n observations
            static inline void activate(const Matrix& Z, Matrix& A)
            {
                A.noalias() = Z;
            }

            // Apply the Jacobian matrix J to a vector f
            // J = d_a / d_z = I
            // g = J * f = f
            // Z = [z1, ..., zn], G = [g1, ..., gn], F = [f1, ..., fn]
            // Note: When entering this function, Z and G may point to the same matrix
            static inline void apply_jacobian(const Matrix& Z, const Matrix& A,
                                            const Matrix& F, Matrix& G)
            {
                G.noalias() = F;
            }

            static std::string return_type()
            {
                return "Identity";
            }
    };

    ///
    /// \ingroup Activations
    ///
    /// The Mish activation function
    ///
    /// From: https://arxiv.org/abs/1908.08681
    ///
    class Mish
    {
        private:
            typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Matrix;

        public:
            // Mish(x) = x * tanh(softplus(x))
            // softplus(x) = log(1 + exp(x))
            // a = activation(z) = Mish(z)
            // Z = [z1, ..., zn], A = [a1, ..., an], n observations
            static inline void activate(const Matrix& Z, Matrix& A)
            {
                // h(x) = tanh(softplus(x)) = (1 + exp(x))^2 - 1
                //                            ------------------
                //                            (1 + exp(x))^2 + 1
                // Let s = exp(-abs(x)), t = 1 + s
                // If x >= 0, then h(x) = (t^2 - s^2) / (t^2 + s^2)
                // If x <= 0, then h(x) = (t^2 - 1) / (t^2 + 1)
                Matrix S = (-Z.array().abs()).exp();
                A.array() = (S.array() + Scalar(1)).square();  // t^2
                S.noalias() = (Z.array() >= Scalar(0)).select(S.cwiseAbs2(), Scalar(1));  // s^2 or 1
                A.array() = (A.array() - S.array()) / (A.array() + S.array());
                A.array() *= Z.array();
            }

            // Apply the Jacobian matrix J to a vector f
            // J = d_a / d_z = diag(Mish'(z))
            // g = J * f = Mish'(z) .* f
            // Z = [z1, ..., zn], G = [g1, ..., gn], F = [f1, ..., fn]
            // Note: When entering this function, Z and G may point to the same matrix
            static inline void apply_jacobian(const Matrix& Z, const Matrix& A,
                                            const Matrix& F, Matrix& G)
            {
                // Let h(x) = tanh(softplus(x))
                // Mish'(x) = h(x) + x * h'(x)
                // h'(x) = tanh'(softplus(x)) * softplus'(x)
                //       = [1 - h(x)^2] * exp(x) / (1 + exp(x))
                //       = [1 - h(x)^2] / (1 + exp(-x))
                // Mish'(x) = h(x) + [x - Mish(x) * h(x)] / (1 + exp(-x))
                // A = Mish(Z) = Z .* h(Z) => h(Z) = A ./ Z, h(0) = 0.6
                G.noalias() = (Z.array() == Scalar(0)).select(Scalar(0.6), A.cwiseQuotient(Z));
                G.array() += (Z.array() - A.array() * G.array()) / (Scalar(1) + (-Z).array().exp());
                G.array() *= F.array();
            }

            static std::string return_type()
            {
                return "Mish";
            }
        };

        ///
    /// \ingroup Activations
    ///
    /// The ReLU activation function
    ///
    class ReLU
    {
        private:
            typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Matrix;

        public:
            // a = activation(z) = max(z, 0)
            // Z = [z1, ..., zn], A = [a1, ..., an], n observations
            static inline void activate(const Matrix& Z, Matrix& A)
            {
                A.array() = Z.array().cwiseMax(Scalar(0));
            }

            // Apply the Jacobian matrix J to a vector f
            // J = d_a / d_z = diag(sign(a)) = diag(a > 0)
            // g = J * f = (a > 0) .* f
            // Z = [z1, ..., zn], G = [g1, ..., gn], F = [f1, ..., fn]
            // Note: When entering this function, Z and G may point to the same matrix
            static inline void apply_jacobian(const Matrix& Z, const Matrix& A,
                                            const Matrix& F, Matrix& G)
            {
                G.array() = (A.array() > Scalar(0)).select(F, Scalar(0));
            }

            static std::string return_type()
            {
                return "ReLU";
            }
    };

    ///
    /// \ingroup Activations
    ///
    /// The sigmoid activation function
    ///
    class Sigmoid
    {
        private:
            typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Matrix;

        public:
            // a = activation(z) = 1 / (1 + exp(-z))
            // Z = [z1, ..., zn], A = [a1, ..., an], n observations
            static inline void activate(const Matrix& Z, Matrix& A)
            {
                A.array() = Scalar(1) / (Scalar(1) + (-Z.array()).exp());
            }

            // Apply the Jacobian matrix J to a vector f
            // J = d_a / d_z = diag(a .* (1 - a))
            // g = J * f = a .* (1 - a) .* f
            // Z = [z1, ..., zn], G = [g1, ..., gn], F = [f1, ..., fn]
            // Note: When entering this function, Z and G may point to the same matrix
            static inline void apply_jacobian(const Matrix& Z, const Matrix& A,
                                            const Matrix& F, Matrix& G)
            {
                G.array() = A.array() * (Scalar(1) - A.array()) * F.array();
            }

            static std::string return_type()
            {
                return "Sigmoid";
            }
    };

    ///
    /// \ingroup Activations
    ///
    /// The softmax activation function
    ///
    class Softmax
    {
        private:
            typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Matrix;
            typedef Eigen::Array<Scalar, 1, Eigen::Dynamic> RowArray;

        public:
            // a = activation(z) = softmax(z)
            // Z = [z1, ..., zn], A = [a1, ..., an], n observations
            static inline void activate(const Matrix& Z, Matrix& A)
            {
                A.array() = (Z.rowwise() - Z.colwise().maxCoeff()).array().exp();
                RowArray colsums = A.colwise().sum();
                A.array().rowwise() /= colsums;
            }

            // Apply the Jacobian matrix J to a vector f
            // J = d_a / d_z = diag(a) - a * a'
            // g = J * f = a .* f - a * (a' * f) = a .* (f - a'f)
            // Z = [z1, ..., zn], G = [g1, ..., gn], F = [f1, ..., fn]
            // Note: When entering this function, Z and G may point to the same matrix
            static inline void apply_jacobian(const Matrix& Z, const Matrix& A,
                                            const Matrix& F, Matrix& G)
            {
                RowArray a_dot_f = A.cwiseProduct(F).colwise().sum();
                G.array() = A.array() * (F.array().rowwise() - a_dot_f);
            }

            static std::string return_type()
            {
                return "Softmax";
            }
    };


    ///
    /// \ingroup Activations
    ///
    /// The tanh activation function
    ///
    class Tanh
    {
        private:
            typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Matrix;

        public:
            // a = activation(z) = tanh(z)
            // Z = [z1, ..., zn], A = [a1, ..., an], n observations
            static inline void activate(const Matrix& Z, Matrix& A)
            {
                A.array() = Z.array().tanh();
            }

            // Apply the Jacobian matrix J to a vector f
            // tanh'(x) = 1 - tanh(x)^2
            // J = d_a / d_z = diag(1 - a^2)
            // g = J * f = (1 - a^2) .* f
            // Z = [z1, ..., zn], G = [g1, ..., gn], F = [f1, ..., fn]
            // Note: When entering this function, Z and G may point to the same matrix
            static inline void apply_jacobian(const Matrix& Z, const Matrix& A,
                                            const Matrix& F, Matrix& G)
            {
                G.array() = (Scalar(1) - A.array().square()) * F.array();
            }

            static std::string return_type()
            {
                return "Tanh";
            }
    };


    ///
    /// \ingroup Optimizers
    ///
    /// The AdaGrad algorithm
    ///
    class AdaGrad: public Optimizer
    {
        private:
            typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
            typedef Eigen::Array<Scalar, Eigen::Dynamic, 1> Array;
            typedef Vector::ConstAlignedMapType ConstAlignedMapVec;
            typedef Vector::AlignedMapType AlignedMapVec;

            std::map<const Scalar*, Array> m_history;

        public:
            Scalar m_lrate;
            Scalar m_eps;

            AdaGrad(const Scalar& lrate = Scalar(0.001), const Scalar& eps = Scalar(1e-6)) :
                m_lrate(lrate), m_eps(eps)
            {}

            void reset()
            {
                m_history.clear();
            }

            void update(ConstAlignedMapVec& dvec, AlignedMapVec& vec)
            {
                // Get the accumulated squared gradient associated with this gradient
                Array& grad_square = m_history[dvec.data()];

                // If length is zero, initialize it
                if (grad_square.size() == 0)
                {
                    grad_square.resize(dvec.size());
                    grad_square.setZero();
                }

                // Update accumulated squared gradient
                grad_square += dvec.array().square();
                // Update parameters
                vec.array() -= m_lrate * dvec.array() / (grad_square.sqrt() + m_eps);
            }
    };

    ///
    /// \ingroup Optimizers
    ///
    /// The Adam algorithm
    ///
    class Adam: public Optimizer
    {
        private:
            typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
            typedef Eigen::Array<Scalar, Eigen::Dynamic, 1> Array;
            typedef Vector::ConstAlignedMapType ConstAlignedMapVec;
            typedef Vector::AlignedMapType AlignedMapVec;

            std::map<const Scalar*, Array> m_history_m;
            std::map<const Scalar*, Array> m_history_v;
            Scalar m_beta1t;
            Scalar m_beta2t;

        public:
            Scalar m_lrate;
            Scalar m_eps;
            Scalar m_beta1;
            Scalar m_beta2;

            Adam(const Scalar& lrate = Scalar(0.001), const Scalar& eps = Scalar(1e-6),
                const Scalar& beta1 = Scalar(0.9), const Scalar& beta2 = Scalar(0.999)) :
                m_beta1t(beta1), m_beta2t(beta2),
                m_lrate(lrate), m_eps(eps),
                m_beta1(beta1), m_beta2(beta2)
            {}

            void reset()
            {
                m_history_m.clear();
                m_history_v.clear();
                m_beta1t = m_beta1;
                m_beta2t = m_beta2;
            }

            // https://ruder.io/optimizing-gradient-descent/index.html
            void update(ConstAlignedMapVec& dvec, AlignedMapVec& vec)
            {
                using std::sqrt;
                // Get the m and v vectors associated with this gradient
                Array& mvec = m_history_m[dvec.data()];
                Array& vvec = m_history_v[dvec.data()];

                // If length is zero, initialize it
                if (mvec.size() == 0)
                {
                    mvec.resize(dvec.size());
                    mvec.setZero();
                }

                if (vvec.size() == 0)
                {
                    vvec.resize(dvec.size());
                    vvec.setZero();
                }

                // Update m and v vectors
                mvec = m_beta1 * mvec + (Scalar(1) - m_beta1) * dvec.array();
                vvec = m_beta2 * vvec + (Scalar(1) - m_beta2) * dvec.array().square();
                // Correction coefficients
                const Scalar correct1 = Scalar(1) / (Scalar(1) - m_beta1t);
                const Scalar correct2 = Scalar(1) / sqrt(Scalar(1) - m_beta2t);
                // Update parameters
                vec.array() -= (m_lrate * correct1) * mvec / (correct2 * vvec.sqrt() + m_eps);
                m_beta1t *= m_beta1;
                m_beta2t *= m_beta2;
            }
    };

    ///
    /// \ingroup Optimizers
    ///
    /// The RMSProp algorithm
    ///
    class RMSProp: public Optimizer
    {
        private:
            typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
            typedef Eigen::Array<Scalar, Eigen::Dynamic, 1> Array;
            typedef Vector::ConstAlignedMapType ConstAlignedMapVec;
            typedef Vector::AlignedMapType AlignedMapVec;

            std::map<const Scalar*, Array> m_history;

        public:
            Scalar m_lrate;
            Scalar m_eps;
            Scalar m_gamma;

            RMSProp(const Scalar& lrate = Scalar(0.001), const Scalar& eps = Scalar(1e-6),
                    const Scalar& gamma = Scalar(0.9)) :
                m_lrate(lrate), m_eps(eps), m_gamma(gamma)
            {}

            void reset()
            {
                m_history.clear();
            }

            void update(ConstAlignedMapVec& dvec, AlignedMapVec& vec)
            {
                // Get the accumulated squared gradient associated with this gradient
                Array& grad_square = m_history[dvec.data()];

                // If length is zero, initialize it
                if (grad_square.size() == 0)
                {
                    grad_square.resize(dvec.size());
                    grad_square.setZero();
                }

                // Update accumulated squared gradient
                grad_square = m_gamma * grad_square + (Scalar(1) - m_gamma) *
                            dvec.array().square();
                // Update parameters
                vec.array() -= m_lrate * dvec.array() / (grad_square + m_eps).sqrt();
            }
    };

    ///
    /// \ingroup Optimizers
    ///
    /// The Stochastic Gradient Descent (SGD) algorithm
    ///
    class SGD: public Optimizer
    {
        private:
            typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
            typedef Vector::ConstAlignedMapType ConstAlignedMapVec;
            typedef Vector::AlignedMapType AlignedMapVec;

        public:
            Scalar m_lrate;
            Scalar m_decay;

            SGD(const Scalar& lrate = Scalar(0.001), const Scalar& decay = Scalar(0)) :
                m_lrate(lrate), m_decay(decay)
            {}

            void update(ConstAlignedMapVec& dvec, AlignedMapVec& vec)
            {
                vec.noalias() -= m_lrate * (dvec + m_decay * vec);
            }
    };
    
    ///
    /// \ingroup Outputs
    ///
    /// Binary classification output layer using cross-entropy criterion
    ///
    class BinaryClassEntropy: public Output
    {
        private:
            Matrix m_din;  // Derivative of the input of this layer.
            // Note that input of this layer is also the output of previous layer

        public:
            void check_target_data(const Matrix& target)
            {
                // Each element should be either 0 or 1
                const int nelem = target.size();
                const Scalar* target_data = target.data();

                for (int i = 0; i < nelem; i++)
                {
                    if ((target_data[i] != Scalar(0)) && (target_data[i] != Scalar(1)))
                    {
                        throw std::invalid_argument("[class BinaryClassEntropy]: Target data should only contain zero or one");
                    }
                }
            }

            void check_target_data(const IntegerVector& target)
            {
                // Each element should be either 0 or 1
                const int nobs = target.size();

                for (int i = 0; i < nobs; i++)
                {
                    if ((target[i] != 0) && (target[i] != 1))
                    {
                        throw std::invalid_argument("[class BinaryClassEntropy]: Target data should only contain zero or one");
                    }
                }
            }

            void evaluate(const Matrix& prev_layer_data, const Matrix& target)
            {
                // Check dimension
                const int nobs = prev_layer_data.cols();
                const int nvar = prev_layer_data.rows();

                if ((target.cols() != nobs) || (target.rows() != nvar))
                {
                    throw std::invalid_argument("[class BinaryClassEntropy]: Target data have incorrect dimension");
                }

                // Compute the derivative of the input of this layer
                // L = -y * log(phat) - (1 - y) * log(1 - phat)
                // in = phat
                // d（L） / d（in） = -y / phat + (1 - y) / (1 - phat), y is either 0 or 1
                m_din.resize(nvar, nobs);
                m_din.array() = (target.array() < Scalar(0.5)).select((Scalar(
                                    1) - prev_layer_data.array()).cwiseInverse(),
                                -prev_layer_data.cwiseInverse());
            }

            void evaluate(const Matrix& prev_layer_data, const IntegerVector& target)
            {
                // Only when the last hidden layer has only one unit can we use this version
                const int nvar = prev_layer_data.rows();

                if (nvar != 1)
                {
                    throw std::invalid_argument("[class BinaryClassEntropy]: Only one response variable is allowed when class labels are used as target data");
                }

                // Check dimension
                const int nobs = prev_layer_data.cols();

                if (target.size() != nobs)
                {
                    throw std::invalid_argument("[class BinaryClassEntropy]: Target data have incorrect dimension");
                }

                // Same as above
                m_din.resize(1, nobs);
                m_din.array() = (target.array() == 0).select((Scalar(1) -
                                prev_layer_data.array()).cwiseInverse(),
                                -prev_layer_data.cwiseInverse());
            }

            const Matrix& backprop_data() const
            {
                return m_din;
            }

            Scalar loss() const
            {
                // L = -y * log(phat) - (1 - y) * log(1 - phat)
                // y = 0 => L = -log(1 - phat)
                // y = 1 => L = -log(phat)
                // m_din contains 1/(1 - phat) if y = 0, and -1/phat if y = 1, so
                // L = log(abs(m_din)).sum()
                return m_din.array().abs().log().sum() / m_din.cols();
            }

            std::string output_type() const
            {
                return "BinaryClassEntropy";
            }
    };


    ///
    /// \ingroup Outputs
    ///
    /// Multi-class classification output layer using cross-entropy criterion
    ///
    class MultiClassEntropy: public Output
    {
        private:
            Matrix m_din;  // Derivative of the input of this layer.
            // Note that input of this layer is also the output of previous layer

        public:
            void check_target_data(const Matrix& target)
            {
                // Each element should be either 0 or 1
                // Each column has and only has one 1
                const int nobs = target.cols();
                const int nclass = target.rows();

                for (int i = 0; i < nobs; i++)
                {
                    int one = 0;

                    for (int j = 0; j < nclass; j++)
                    {
                        if (target(j, i) == Scalar(1))
                        {
                            one++;
                            continue;
                        }

                        if (target(j, i) != Scalar(0))
                        {
                            throw std::invalid_argument("[class MultiClassEntropy]: Target data should only contain zero or one");
                        }
                    }

                    if (one != 1)
                    {
                        throw std::invalid_argument("[class MultiClassEntropy]: Each column of target data should only contain one \"1\"");
                    }
                }
            }

            void check_target_data(const IntegerVector& target)
            {
                // All elements must be non-negative
                const int nobs = target.size();

                for (int i = 0; i < nobs; i++)
                {
                    if (target[i] < 0)
                    {
                        throw std::invalid_argument("[class MultiClassEntropy]: Target data must be non-negative");
                    }
                }
            }

            // target is a matrix with each column representing an observation
            // Each column is a vector that has a one at some location and has zeros elsewhere
            void evaluate(const Matrix& prev_layer_data, const Matrix& target)
            {
                // Check dimension
                const int nobs = prev_layer_data.cols();
                const int nclass = prev_layer_data.rows();

                if ((target.cols() != nobs) || (target.rows() != nclass))
                {
                    throw std::invalid_argument("[class MultiClassEntropy]: Target data have incorrect dimension");
                }

                // Compute the derivative of the input of this layer
                // L = -sum(log(phat) * y)
                // in = phat
                // d(L) / d(in) = -y / phat
                m_din.resize(nclass, nobs);
                m_din.noalias() = -target.cwiseQuotient(prev_layer_data);
            }

            // target is a vector of class labels that take values from [0, 1, ..., nclass - 1]
            // The i-th element of target is the class label for observation i
            void evaluate(const Matrix& prev_layer_data, const IntegerVector& target)
            {
                // Check dimension
                const int nobs = prev_layer_data.cols();
                const int nclass = prev_layer_data.rows();

                if (target.size() != nobs)
                {
                    throw std::invalid_argument("[class MultiClassEntropy]: Target data have incorrect dimension");
                }

                // Compute the derivative of the input of this layer
                // L = -log(phat[y])
                // in = phat
                // d(L) / d(in) = [0, 0, ..., -1/phat[y], 0, ..., 0]
                m_din.resize(nclass, nobs);
                m_din.setZero();

                for (int i = 0; i < nobs; i++)
                {
                    m_din(target[i], i) = -Scalar(1) / prev_layer_data(target[i], i);
                }
            }

            const Matrix& backprop_data() const
            {
                return m_din;
            }

            Scalar loss() const
            {
                // L = -sum(log(phat) * y)
                // in = phat
                // d(L) / d(in) = -y / phat
                // m_din contains 0 if y = 0, and -1/phat if y = 1
                Scalar res = Scalar(0);
                const int nelem = m_din.size();
                const Scalar* din_data = m_din.data();

                for (int i = 0; i < nelem; i++)
                {
                    if (din_data[i] < Scalar(0))
                    {
                        res += std::log(-din_data[i]);
                    }
                }

                return res / m_din.cols();
            }

            std::string output_type() const
            {
                return "MultiClassEntropy";
            }
    };


    ///
    /// \ingroup Outputs
    ///
    /// Regression output layer using Mean Squared Error (MSE) criterion
    ///
    class RegressionMSE: public Output
    {
        private:
            typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Matrix;
            typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;

            Matrix m_din;  // Derivative of the input of this layer.
            // Note that input of this layer is also the output of previous layer

        public:
            void evaluate(const Matrix& prev_layer_data, const Matrix& target)
            {
                // Check dimension
                const int nobs = prev_layer_data.cols();
                const int nvar = prev_layer_data.rows();

                if ((target.cols() != nobs) || (target.rows() != nvar))
                {
                    throw std::invalid_argument("[class RegressionMSE]: Target data have incorrect dimension");
                }

                // Compute the derivative of the input of this layer
                // L = 0.5 * ||yhat - y||^2
                // in = yhat
                // d(L) / d(in) = yhat - y
                m_din.resize(nvar, nobs);
                m_din.noalias() = prev_layer_data - target;
            }

            const Matrix& backprop_data() const
            {
                return m_din;
            }

            Scalar loss() const
            {
                // L = 0.5 * ||yhat - y||^2
                return m_din.squaredNorm() / m_din.cols() * Scalar(0.5);
            }

            std::string output_type() const
            {
                return "RegressionMSE";
            }
    };
    ///
    /// \ingroup Layers
    ///
    /// Fully connected hidden layer
    ///
    template <typename Activation>
    class FullyConnected: public Layer
    {
        private:
            typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Matrix;
            typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
            typedef Vector::ConstAlignedMapType ConstAlignedMapVec;
            typedef Vector::AlignedMapType AlignedMapVec;
            typedef std::map<std::string, int> MetaInfo;

            Matrix m_weight;  // Weight parameters, W(in_size x out_size)
            Vector m_bias;    // Bias parameters, b(out_size x 1)
            Matrix m_dw;      // Derivative of weights
            Vector m_db;      // Derivative of bias
            Matrix m_z;       // Linear term, z = W' * in + b
            Matrix m_a;       // Output of this layer, a = act(z)
            Matrix m_din;     // Derivative of the input of this layer.
                            // Note that input of this layer is also the output of previous layer

        public:
            ///
            /// Constructor
            ///
            /// \param in_size  Number of input units.
            /// \param out_size Number of output units.
            ///
            FullyConnected(const int in_size, const int out_size) :
                Layer(in_size, out_size)
            {}

            void init(const Scalar& mu, const Scalar& sigma, RNG& rng)
            {
                // Set parameter dimension
                init();
                // Set random coefficients
                internal::set_normal_random(m_weight.data(), m_weight.size(), rng, mu, sigma);
                internal::set_normal_random(m_bias.data(), m_bias.size(), rng, mu, sigma);
            }

            void init()
            {
                // Set parameter dimension
                m_weight.resize(this->m_in_size, this->m_out_size);
                m_bias.resize(this->m_out_size);
                m_dw.resize(this->m_in_size, this->m_out_size);
                m_db.resize(this->m_out_size);
            }

            // prev_layer_data: in_size x nobs
            void forward(const Matrix& prev_layer_data)
            {
                const int nobs = prev_layer_data.cols();
                // Linear term z = W' * in + b
                m_z.resize(this->m_out_size, nobs);
                m_z.noalias() = m_weight.transpose() * prev_layer_data;
                m_z.colwise() += m_bias;
                // Apply activation function
                m_a.resize(this->m_out_size, nobs);
                Activation::activate(m_z, m_a);
            }

            const Matrix& output() const
            {
                return m_a;
            }

            // prev_layer_data: in_size x nobs
            // next_layer_data: out_size x nobs
            void backprop(const Matrix& prev_layer_data, const Matrix& next_layer_data)
            {
                const int nobs = prev_layer_data.cols();
                // After forward stage, m_z contains z = W' * in + b
                // Now we need to calculate d(L) / d(z) = [d(a) / d(z)] * [d(L) / d(a)]
                // d(L) / d(a) is computed in the next layer, contained in next_layer_data
                // The Jacobian matrix J = d(a) / d(z) is determined by the activation function
                Matrix& dLz = m_z;
                Activation::apply_jacobian(m_z, m_a, next_layer_data, dLz);
                // Now dLz contains d(L) / d(z)
                // Derivative for weights, d(L) / d(W) = [d(L) / d(z)] * in'
                m_dw.noalias() = prev_layer_data * dLz.transpose() / nobs;
                // Derivative for bias, d(L) / d(b) = d(L) / d(z)
                m_db.noalias() = dLz.rowwise().mean();
                // Compute d(L) / d_in = W * [d(L) / d(z)]
                m_din.resize(this->m_in_size, nobs);
                m_din.noalias() = m_weight * dLz;
            }

            const Matrix& backprop_data() const
            {
                return m_din;
            }

            void update(Optimizer& opt)
            {
                ConstAlignedMapVec dw(m_dw.data(), m_dw.size());
                ConstAlignedMapVec db(m_db.data(), m_db.size());
                AlignedMapVec      w(m_weight.data(), m_weight.size());
                AlignedMapVec      b(m_bias.data(), m_bias.size());
                opt.update(dw, w);
                opt.update(db, b);
            }

            std::vector<Scalar> get_parameters() const
            {
                std::vector<Scalar> res(m_weight.size() + m_bias.size());
                // Copy the data of weights and bias to a long vector
                std::copy(m_weight.data(), m_weight.data() + m_weight.size(), res.begin());
                std::copy(m_bias.data(), m_bias.data() + m_bias.size(),
                        res.begin() + m_weight.size());
                return res;
            }

            void set_parameters(const std::vector<Scalar>& param)
            {
                if (static_cast<int>(param.size()) != m_weight.size() + m_bias.size())
                {
                    throw std::invalid_argument("[class FullyConnected]: Parameter size does not match");
                }

                std::copy(param.begin(), param.begin() + m_weight.size(), m_weight.data());
                std::copy(param.begin() + m_weight.size(), param.end(), m_bias.data());
            }

            std::vector<Scalar> get_derivatives() const
            {
                std::vector<Scalar> res(m_dw.size() + m_db.size());
                // Copy the data of weights and bias to a long vector
                std::copy(m_dw.data(), m_dw.data() + m_dw.size(), res.begin());
                std::copy(m_db.data(), m_db.data() + m_db.size(), res.begin() + m_dw.size());
                return res;
            }

            std::string layer_type() const
            {
                return "FullyConnected";
            }

            std::string activation_type() const
            {
                return Activation::return_type();
            }

            void fill_meta_info(MetaInfo& map, int index) const
            {
                std::string ind = internal::to_string(index);
                map.insert(std::make_pair("Layer" + ind, internal::layer_id(layer_type())));
                map.insert(std::make_pair("Activation" + ind, internal::activation_id(activation_type())));
                map.insert(std::make_pair("in_size" + ind, in_size()));
                map.insert(std::make_pair("out_size" + ind, out_size()));
            }
    };


    ///
    /// \ingroup Layers
    ///
    /// Convolutional hidden layer
    ///
    /// Currently only supports the "valid" rule of convolution.
    ///
    template <typename Activation>
    class Convolutional: public Layer
    {
        private:
            typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Matrix;
            typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
            typedef Matrix::ConstAlignedMapType ConstAlignedMapMat;
            typedef Vector::ConstAlignedMapType ConstAlignedMapVec;
            typedef Vector::AlignedMapType AlignedMapVec;
            typedef std::map<std::string, int> MetaInfo;

            const internal::ConvDims m_dim; // Various dimensions of convolution


            Vector m_filter_data;  // Filter parameters. Total length is
                                // (in_channels x out_channels x filter_rows x filter_cols)
                                // See Utils/Convolution.h for its layout

            Vector m_df_data;      // Derivative of filters, same dimension as m_filter_data

            Vector m_bias;         // Bias term for the output channels, out_channels x 1. (One bias term per channel)
            Vector m_db;           // Derivative of bias, same dimension as m_bias

            Matrix m_z;            // Linear term, z = conv(in, w) + b. Each column is an observation
            Matrix m_a;            // Output of this layer, a = act(z)
            Matrix m_din;          // Derivative of the input of this layer
                                // Note that input of this layer is also the output of previous layer

        public:
            ///
            /// Constructor
            ///
            /// \param in_width      Width of the input image in each channel.
            /// \param in_height     Height of the input image in each channel.
            /// \param in_channels   Number of input channels.
            /// \param out_channels  Number of output channels.
            /// \param window_width  Width of the filter.
            /// \param window_height Height of the filter.
            ///
            Convolutional(const int in_width, const int in_height,
                        const int in_channels, const int out_channels,
                        const int window_width, const int window_height) :
                Layer(in_width * in_height * in_channels,
                    (in_width - window_width + 1) * (in_height - window_height + 1) * out_channels),
                m_dim(in_channels, out_channels, in_height, in_width, window_height,
                    window_width)
            {}

            void init(const Scalar& mu, const Scalar& sigma, RNG& rng)
            {
                // Set data dimension
                init();
                // Random initialization of filter parameters
                const int filter_data_size = m_dim.in_channels * m_dim.out_channels *
                                            m_dim.filter_rows * m_dim.filter_cols;
                internal::set_normal_random(m_filter_data.data(), filter_data_size, rng, mu,
                                            sigma);
                // Bias term
                internal::set_normal_random(m_bias.data(), m_dim.out_channels, rng, mu, sigma);
            }

            void init()
            {
                // Set parameter dimension
                const int filter_data_size = m_dim.in_channels * m_dim.out_channels *
                                            m_dim.filter_rows * m_dim.filter_cols;
                // Filter parameters
                m_filter_data.resize(filter_data_size);
                m_df_data.resize(filter_data_size);
                // Bias term
                m_bias.resize(m_dim.out_channels);
                m_db.resize(m_dim.out_channels);
            }

            // http://cs231n.github.io/convolutional-networks/
            void forward(const Matrix& prev_layer_data)
            {
                // Each column is an observation
                const int nobs = prev_layer_data.cols();
                // Linear term, z = conv(in, w) + b
                m_z.resize(this->m_out_size, nobs);
                // Convolution
                internal::convolve_valid(m_dim, prev_layer_data.data(), true, nobs,
                                        m_filter_data.data(), m_z.data()
                                        );
                // Add bias terms
                // Each column of m_z contains m_dim.out_channels channels, and each channel has
                // m_dim.conv_rows * m_dim.conv_cols elements
                int channel_start_row = 0;
                const int channel_nelem = m_dim.conv_rows * m_dim.conv_cols;

                for (int i = 0; i < m_dim.out_channels; i++, channel_start_row += channel_nelem)
                {
                    m_z.block(channel_start_row, 0, channel_nelem, nobs).array() += m_bias[i];
                }

                // Apply activation function
                m_a.resize(this->m_out_size, nobs);
                Activation::activate(m_z, m_a);
            }

            const Matrix& output() const
            {
                return m_a;
            }

            // prev_layer_data: in_size x nobs
            // next_layer_data: out_size x nobs
            // https://grzegorzgwardys.wordpress.com/2016/04/22/8/
            void backprop(const Matrix& prev_layer_data, const Matrix& next_layer_data)
            {
                const int nobs = prev_layer_data.cols();
                // After forward stage, m_z contains z = conv(in, w) + b
                // Now we need to calculate d(L) / d(z) = [d(a) / d(z)] * [d(L) / d(a)]
                // d(L) / d(a) is computed in the next layer, contained in next_layer_data
                // The Jacobian matrix J = d(a) / d(z) is determined by the activation function
                Matrix& dLz = m_z;
                Activation::apply_jacobian(m_z, m_a, next_layer_data, dLz);
                // z_j = sum_i(conv(in_i, w_ij)) + b_j
                //
                // d(z_k) / d(w_ij) = 0, if k != j
                // d(L) / d(w_ij) = [d(z_j) / d(w_ij)] * [d(L) / d(z_j)] = sum_i{ [d(z_j) / d(w_ij)] * [d(L) / d(z_j)] }
                // = sum_i(conv(in_i, d(L) / d(z_j)))
                //
                // z_j is an image (matrix), b_j is a scalar
                // d(z_j) / d(b_j) = a matrix of the same size of d(z_j) filled with 1
                // d(L) / d(b_j) = (d(L) / d(z_j)).sum()
                //
                // d(z_j) / d(in_i) = conv_full_op(w_ij_rotate)
                // d(L) / d(in_i) = sum_j((d(z_j) / d(in_i)) * (d(L) / d(z_j))) = sum_j(conv_full(d(L) / d(z_j), w_ij_rotate))
                // Derivative for weights
                internal::ConvDims back_conv_dim(nobs, m_dim.out_channels, m_dim.channel_rows,
                                                m_dim.channel_cols,
                                                m_dim.conv_rows, m_dim.conv_cols);
                internal::convolve_valid(back_conv_dim, prev_layer_data.data(), false,
                                        m_dim.in_channels,
                                        dLz.data(), m_df_data.data()
                                        );
                m_df_data /= nobs;
                // Derivative for bias
                // Aggregate d(L) / d(z) in each output channel
                ConstAlignedMapMat dLz_by_channel(dLz.data(), m_dim.conv_rows * m_dim.conv_cols,
                                                m_dim.out_channels * nobs);
                Vector dLb = dLz_by_channel.colwise().sum();
                // Average over observations
                ConstAlignedMapMat dLb_by_obs(dLb.data(), m_dim.out_channels, nobs);
                m_db.noalias() = dLb_by_obs.rowwise().mean();
                // Compute d(L) / d_in = conv_full(d(L) / d(z), w_rotate)
                m_din.resize(this->m_in_size, nobs);
                internal::ConvDims conv_full_dim(m_dim.out_channels, m_dim.in_channels,
                                                m_dim.conv_rows, m_dim.conv_cols, m_dim.filter_rows, m_dim.filter_cols);
                internal::convolve_full(conv_full_dim, dLz.data(), nobs,
                                        m_filter_data.data(), m_din.data()
                                    );
            }

            const Matrix& backprop_data() const
            {
                return m_din;
            }

            void update(Optimizer& opt)
            {
                ConstAlignedMapVec dw(m_df_data.data(), m_df_data.size());
                ConstAlignedMapVec db(m_db.data(), m_db.size());
                AlignedMapVec      w(m_filter_data.data(), m_filter_data.size());
                AlignedMapVec      b(m_bias.data(), m_bias.size());
                opt.update(dw, w);
                opt.update(db, b);
            }

            std::vector<Scalar> get_parameters() const
            {
                std::vector<Scalar> res(m_filter_data.size() + m_bias.size());
                // Copy the data of filters and bias to a long vector
                std::copy(m_filter_data.data(), m_filter_data.data() + m_filter_data.size(),
                        res.begin());
                std::copy(m_bias.data(), m_bias.data() + m_bias.size(),
                        res.begin() + m_filter_data.size());
                return res;
            }

            void set_parameters(const std::vector<Scalar>& param)
            {
                if (static_cast<int>(param.size()) != m_filter_data.size() + m_bias.size())
                {
                    throw std::invalid_argument("[class Convolutional]: Parameter size does not match");
                }

                std::copy(param.begin(), param.begin() + m_filter_data.size(),
                        m_filter_data.data());
                std::copy(param.begin() + m_filter_data.size(), param.end(), m_bias.data());
            }

            std::vector<Scalar> get_derivatives() const
            {
                std::vector<Scalar> res(m_df_data.size() + m_db.size());
                // Copy the data of filters and bias to a long vector
                std::copy(m_df_data.data(), m_df_data.data() + m_df_data.size(), res.begin());
                std::copy(m_db.data(), m_db.data() + m_db.size(),
                        res.begin() + m_df_data.size());
                return res;
            }

            std::string layer_type() const
            {
                return "Convolutional";
            }

            std::string activation_type() const
            {
                return Activation::return_type();
            }

            void fill_meta_info(MetaInfo& map, int index) const
            {
                std::string ind = internal::to_string(index);
                map.insert(std::make_pair("Layer" + ind, internal::layer_id(layer_type())));
                map.insert(std::make_pair("Activation" + ind, internal::activation_id(activation_type())));
                map.insert(std::make_pair("in_channels" + ind, m_dim.in_channels));
                map.insert(std::make_pair("out_channels" + ind, m_dim.out_channels));
                map.insert(std::make_pair("in_height" + ind, m_dim.channel_rows));
                map.insert(std::make_pair("in_width" + ind, m_dim.channel_cols));
                map.insert(std::make_pair("window_width" + ind, m_dim.filter_cols));
                map.insert(std::make_pair("window_height" + ind, m_dim.filter_rows));
            }
    };
    ///
/// \ingroup Layers
///
/// Max-pooling hidden layer
///
/// Currently only supports the "valid" rule of pooling.
///
template <typename Activation>
class MaxPooling: public Layer
{
    private:
        typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Matrix;
        typedef Eigen::MatrixXi IntMatrix;
        typedef std::map<std::string, int> MetaInfo;

        const int m_channel_rows;
        const int m_channel_cols;
        const int m_in_channels;
        const int m_pool_rows;
        const int m_pool_cols;

        const int m_out_rows;
        const int m_out_cols;

        IntMatrix m_loc;             // Record the locations of maximums
        Matrix m_z;                  // Max pooling results
        Matrix m_a;                  // Output of this layer, a = act(z)
        Matrix m_din;                // Derivative of the input of this layer.
                                     // Note that input of this layer is also the output of previous layer

    public:
        // Currently we only implement the "valid" rule
        // https://stackoverflow.com/q/37674306
        ///
        /// Constructor
        ///
        /// \param in_width       Width of the input image in each channel.
        /// \param in_height      Height of the input image in each channel.
        /// \param in_channels    Number of input channels.
        /// \param pooling_width  Width of the pooling window.
        /// \param pooling_height Height of the pooling window.
        ///
        MaxPooling(const int in_width_, const int in_height_, const int in_channels_,
                   const int pooling_width_, const int pooling_height_) :
            Layer(in_width_ * in_height_ * in_channels_,
                  (in_width_ / pooling_width_) * (in_height_ / pooling_height_) * in_channels_),
            m_channel_rows(in_height_), m_channel_cols(in_width_),
            m_in_channels(in_channels_),
            m_pool_rows(pooling_height_), m_pool_cols(pooling_width_),
            m_out_rows(m_channel_rows / m_pool_rows),
            m_out_cols(m_channel_cols / m_pool_cols)
        {}

        void init(const Scalar& mu, const Scalar& sigma, RNG& rng) {}

        void init() {}

        void forward(const Matrix& prev_layer_data)
        {
            // Each column is an observation
            const int nobs = prev_layer_data.cols();
            m_loc.resize(this->m_out_size, nobs);
            m_z.resize(this->m_out_size, nobs);
            // Use m_loc to store the address of each pooling block relative to the beginning of the data
            int* loc_data = m_loc.data();
            const int channel_end = prev_layer_data.size();
            const int channel_stride = m_channel_rows * m_channel_cols;
            const int col_end_gap = m_channel_rows * m_pool_cols * m_out_cols;
            const int col_stride = m_channel_rows * m_pool_cols;
            const int row_end_gap = m_out_rows * m_pool_rows;

            for (int channel_start = 0; channel_start < channel_end;
                    channel_start += channel_stride)
            {
                const int col_end = channel_start + col_end_gap;

                for (int col_start = channel_start; col_start < col_end;
                        col_start += col_stride)
                {
                    const int row_end = col_start + row_end_gap;

                    for (int row_start = col_start; row_start < row_end;
                            row_start += m_pool_rows, loc_data++)
                    {
                        *loc_data = row_start;
                    }
                }
            }

            // Find the location of the max value in each block
            loc_data = m_loc.data();
            const int* const loc_end = loc_data + m_loc.size();
            Scalar* z_data = m_z.data();
            const Scalar* src = prev_layer_data.data();

            for (; loc_data < loc_end; loc_data++, z_data++)
            {
                const int offset = *loc_data;
                *z_data = internal::find_block_max(src + offset, m_pool_rows, m_pool_cols,
                                                   m_channel_rows, *loc_data);
                *loc_data += offset;
            }

            // Apply activation function
            m_a.resize(this->m_out_size, nobs);
            Activation::activate(m_z, m_a);
        }

        const Matrix& output() const
        {
            return m_a;
        }

        // prev_layer_data: in_size x nobs
        // next_layer_data: out_size x nobs
        void backprop(const Matrix& prev_layer_data, const Matrix& next_layer_data)
        {
            const int nobs = prev_layer_data.cols();
            // After forward stage, m_z contains z = max_pooling(in)
            // Now we need to calculate d(L) / d(z) = [d(a) / d(z)] * [d(L) / d(a)]
            // d(L) / d(z) is computed in the next layer, contained in next_layer_data
            // The Jacobian matrix J = d(a) / d(z) is determined by the activation function
            Matrix& dLz = m_z;
            Activation::apply_jacobian(m_z, m_a, next_layer_data, dLz);
            // d(L) / d(in_i) = sum_j{ [d(z_j) / d(in_i)] * [d(L) / d(z_j)] }
            // d(z_j) / d(in_i) = 1 if in_i is used to compute z_j and is the maximum
            //                  = 0 otherwise
            m_din.resize(this->m_in_size, nobs);
            m_din.setZero();
            const int dLz_size = dLz.size();
            const Scalar* dLz_data = dLz.data();
            const int* loc_data = m_loc.data();
            Scalar* din_data = m_din.data();

            for (int i = 0; i < dLz_size; i++)
            {
                din_data[loc_data[i]] += dLz_data[i];
            }
        }

        const Matrix& backprop_data() const
        {
            return m_din;
        }

        void update(Optimizer& opt) {}

        std::vector<Scalar> get_parameters() const
        {
            return std::vector<Scalar>();
        }

        void set_parameters(const std::vector<Scalar>& param) {}

        std::vector<Scalar> get_derivatives() const
        {
            return std::vector<Scalar>();
        }

        std::string layer_type() const
        {
            return "MaxPooling";
        }

        std::string activation_type() const
        {
            return Activation::return_type();
        }

        void fill_meta_info(MetaInfo& map, int index) const
        {
            std::string ind = internal::to_string(index);
            map.insert(std::make_pair("Layer" + ind, internal::layer_id(layer_type())));
            map.insert(std::make_pair("Activation" + ind, internal::activation_id(activation_type())));
            map.insert(std::make_pair("in_width" + ind, m_channel_cols));
            map.insert(std::make_pair("in_height" + ind, m_channel_rows));
            map.insert(std::make_pair("in_channels" + ind, m_in_channels));
            map.insert(std::make_pair("pooling_width" + ind, m_pool_cols));
            map.insert(std::make_pair("pooling_height" + ind, m_pool_rows));
        }
    };

    namespace internal
    {
            
        // Create a layer from the network meta information and the index of the layer
        inline Layer* create_layer(const std::map<std::string, int>& map, int index)
        {
            std::string ind = internal::to_string(index);
            const int lay_id = map.find("Layer" + ind)->second;
            const int act_id = map.find("Activation" + ind)->second;
            Layer* layer;

            if (lay_id == FULLY_CONNECTED)
            {
                const int in_size = map.find("in_size" + ind)->second;
                const int out_size = map.find("out_size" + ind)->second;

                switch (act_id)
                {
                case IDENTITY:
                    layer = new FullyConnected<Identity>(in_size, out_size);
                    break;
                case RELU:
                    layer = new FullyConnected<ReLU>(in_size, out_size);
                    break;
                case SIGMOID:
                    layer = new FullyConnected<Sigmoid>(in_size, out_size);
                    break;
                case SOFTMAX:
                    layer = new FullyConnected<Softmax>(in_size, out_size);
                    break;
                case TANH:
                    layer = new FullyConnected<Tanh>(in_size, out_size);
                    break;
                case MISH:
                    layer = new FullyConnected<Mish>(in_size, out_size);
                    break;
                default:
                    throw std::invalid_argument("[function create_layer]: Activation is not of a known type");
                }

            } else if (lay_id == CONVOLUTIONAL) {
                const int in_width = map.find("in_width" + ind)->second;
                const int in_height = map.find("in_height" + ind)->second;
                const int in_channels = map.find("in_channels" + ind)->second;
                const int out_channels = map.find("out_channels" + ind)->second;
                const int window_width = map.find("window_width" + ind)->second;
                const int window_height = map.find("window_height" + ind)->second;

                switch(act_id)
                {
                case IDENTITY:
                    layer = new Convolutional<Identity>(in_width, in_height, in_channels,
                                                        out_channels, window_width, window_height);
                    break;
                case RELU:
                    layer = new Convolutional<ReLU>(in_width, in_height, in_channels,
                                                    out_channels, window_width, window_height);
                    break;
                case SIGMOID:
                    layer = new Convolutional<Sigmoid>(in_width, in_height, in_channels,
                                                    out_channels, window_width, window_height);
                    break;
                case SOFTMAX:
                    layer = new Convolutional<Softmax>(in_width, in_height, in_channels,
                                                    out_channels, window_width, window_height);
                    break;
                case TANH:
                    layer = new Convolutional<Tanh>(in_width, in_height, in_channels,
                                                    out_channels, window_width, window_height);
                    break;
                case MISH:
                    layer = new Convolutional<Mish>(in_width, in_height, in_channels,
                                                    out_channels, window_width, window_height);
                    break;
                default:
                    throw std::invalid_argument("[function create_layer]: Activation is not of a known type");
                }

            } else if (lay_id == MAX_POOLING) {
                const int in_width = map.find("in_width" + ind)->second;
                const int in_height = map.find("in_height" + ind)->second;
                const int in_channels = map.find("in_channels" + ind)->second;
                const int pooling_width = map.find("pooling_width" + ind)->second;
                const int pooling_height = map.find("pooling_height" + ind)->second;

                switch (act_id)
                {
                case IDENTITY:
                    layer = new MaxPooling<Identity>(in_width, in_height, in_channels,
                                                    pooling_width, pooling_height);
                    break;
                case RELU:
                    layer = new MaxPooling<ReLU>(in_width, in_height, in_channels,
                                                pooling_width, pooling_height);
                    break;
                case SIGMOID:
                    layer = new MaxPooling<Sigmoid>(in_width, in_height, in_channels,
                                                    pooling_width, pooling_height);
                    break;
                case SOFTMAX:
                    layer = new MaxPooling<Softmax>(in_width, in_height, in_channels,
                                                    pooling_width, pooling_height);
                    break;
                case TANH:
                    layer = new MaxPooling<Tanh>(in_width, in_height, in_channels,
                                                pooling_width, pooling_height);
                    break;
                case MISH:
                    layer = new MaxPooling<Mish>(in_width, in_height, in_channels,
                                                pooling_width, pooling_height);
                    break;
                default:
                    throw std::invalid_argument("[function create_layer]: Activation is not of a known type");
                }

            } else {

                throw std::invalid_argument("[function create_layer]: Layer is not of a known type");
            }

            layer->init();
            return layer;
        }

        // Create an output layer from the network meta information
        inline Output* create_output(const std::map<std::string, int>& map)
        {
            Output* output;
            int out_id = map.find("OutputLayer")->second;

            switch (out_id)
            {
            case REGRESSION_MSE:
                return new RegressionMSE();
            case BINARY_CLASS_ENTROPY:
                return new BinaryClassEntropy();
            case MULTI_CLASS_ENTROPY:
                return new MultiClassEntropy();
            default:
                throw std::invalid_argument("[function create_output]: Output is not of a known type");
            }

            return output;
        }            
    }

    ///
    /// \defgroup Network Neural Network Model
    ///

    ///
    /// \ingroup Network
    ///
    /// This class represents a neural network model that typically consists of a
    /// number of hidden layers and an output layer. It provides functions for
    /// network building, model fitting, and prediction, etc.
    ///
    class Network
    {
        private:
            typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Matrix;
            typedef Eigen::RowVectorXi IntegerVector;
            typedef std::map<std::string, int> MetaInfo;

            RNG                 m_default_rng;      // Built-in RNG
            RNG&                m_rng;              // Reference to the RNG provided by the user,
                                                    // otherwise reference to m_default_rng
            std::vector<Layer*> m_layers;           // Pointers to hidden layers
            Output*             m_output;           // The output layer
            Callback            m_default_callback; // Default callback function
            Callback*           m_callback;         // Points to user-provided callback function,
                                                    // otherwise points to m_default_callback

            // Check dimensions of layers
            void check_unit_sizes() const
            {
                const int nlayer = num_layers();

                if (nlayer <= 1)
                {
                    return;
                }

                for (int i = 1; i < nlayer; i++)
                {
                    if (m_layers[i]->in_size() != m_layers[i - 1]->out_size())
                    {
                        throw std::invalid_argument("[class Network]: Unit sizes do not match");
                    }
                }
            }

            // Let each layer compute its output
            void forward(const Matrix& input)
            {
                const int nlayer = num_layers();

                if (nlayer <= 0)
                {
                    return;
                }

                // First layer
                if (input.rows() != m_layers[0]->in_size())
                {
                    throw std::invalid_argument("[class Network]: Input data have incorrect dimension");
                }

                m_layers[0]->forward(input);

                // The following layers
                for (int i = 1; i < nlayer; i++)
                {
                    m_layers[i]->forward(m_layers[i - 1]->output());
                }
            }

            // Let each layer compute its gradients of the parameters
            // target has two versions: Matrix and RowVectorXi
            // The RowVectorXi version is used in classification problems where each
            // element is a class label
            template <typename TargetType>
            void backprop(const Matrix& input, const TargetType& target)
            {
                const int nlayer = num_layers();

                if (nlayer <= 0)
                {
                    return;
                }

                Layer* first_layer = m_layers[0];
                Layer* last_layer = m_layers[nlayer - 1];
                // Let output layer compute back-propagation data
                m_output->check_target_data(target);
                m_output->evaluate(last_layer->output(), target);

                // If there is only one hidden layer, "prev_layer_data" will be the input data
                if (nlayer == 1)
                {
                    first_layer->backprop(input, m_output->backprop_data());
                    return;
                }

                // Compute gradients for the last hidden layer
                last_layer->backprop(m_layers[nlayer - 2]->output(), m_output->backprop_data());

                // Compute gradients for all the hidden layers except for the first one and the last one
                for (int i = nlayer - 2; i > 0; i--)
                {
                    m_layers[i]->backprop(m_layers[i - 1]->output(),
                                        m_layers[i + 1]->backprop_data());
                }

                // Compute gradients for the first layer
                first_layer->backprop(input, m_layers[1]->backprop_data());
            }

            // Update parameters
            void update(Optimizer& opt)
            {
                const int nlayer = num_layers();

                if (nlayer <= 0)
                {
                    return;
                }

                for (int i = 0; i < nlayer; i++)
                {
                    m_layers[i]->update(opt);
                }
            }

            // Get the meta information of the network, used to export the NN model
            MetaInfo get_meta_info() const
            {
                const int nlayer = num_layers();
                MetaInfo map;
                map.insert(std::make_pair("Nlayers", nlayer));

                for (int i = 0; i < nlayer; i++)
                {
                    m_layers[i]->fill_meta_info(map, i);
                }

                map.insert(std::make_pair("OutputLayer", internal::output_id(m_output->output_type())));
                return map;
            }

        public:
            ///
            /// Default constructor that creates an empty neural network
            ///
            Network() :
                m_default_rng(1),
                m_rng(m_default_rng),
                m_output(NULL),
                m_default_callback(),
                m_callback(&m_default_callback)
            {}

            ///
            /// Constructor with a user-provided random number generator
            ///
            /// \param rng A user-provided random number generator object that inherits
            ///            from the default RNG class.
            ///
            Network(RNG& rng) :
                m_default_rng(1),
                m_rng(rng),
                m_output(NULL),
                m_default_callback(),
                m_callback(&m_default_callback)
            {}

            ///
            /// Destructor that frees the added hidden layers and output layer
            ///
            ~Network()
            {
                const int nlayer = num_layers();

                for (int i = 0; i < nlayer; i++)
                {
                    delete m_layers[i];
                }

                if (m_output)
                {
                    delete m_output;
                }
            }

            ///
            /// Add a hidden layer to the neural network
            ///
            /// \param layer A pointer to a Layer object, typically constructed from
            ///              layer classes such as FullyConnected and Convolutional.
            ///              **NOTE**: the pointer will be handled and freed by the
            ///              network object, so do not delete it manually.
            ///
            void add_layer(Layer* layer)
            {
                m_layers.push_back(layer);
            }

            ///
            /// Set the output layer of the neural network
            ///
            /// \param output A pointer to an Output object, typically constructed from
            ///               output layer classes such as RegressionMSE and MultiClassEntropy.
            ///               **NOTE**: the pointer will be handled and freed by the
            ///               network object, so do not delete it manually.
            ///
            void set_output(Output* output)
            {
                if (m_output)
                {
                    delete m_output;
                }

                m_output = output;
            }

            ///
            /// Number of hidden layers in the network
            ///
            int num_layers() const
            {
                return m_layers.size();
            }

            ///
            /// Get the list of hidden layers of the network
            ///
            std::vector<const Layer*> get_layers() const
            {
                const int nlayer = num_layers();
                std::vector<const Layer*> layers(nlayer);
                std::copy(m_layers.begin(), m_layers.end(), layers.begin());
                return layers;
            }

            ///
            /// Get the output layer
            ///
            const Output* get_output() const
            {
                return m_output;
            }

            ///
            /// Set the callback function that can be called during model fitting
            ///
            /// \param callback A user-provided callback function object that inherits
            ///                 from the default Callback class.
            ///
            void set_callback(Callback& callback)
            {
                m_callback = &callback;
            }
            ///
            /// Set the default silent callback function
            ///
            void set_default_callback()
            {
                m_callback = &m_default_callback;
            }

            ///
            /// Initialize layer parameters in the network using normal distribution
            ///
            /// \param mu    Mean of the normal distribution.
            /// \param sigma Standard deviation of the normal distribution.
            /// \param seed  Set the random seed of the %RNG if `seed > 0`, otherwise
            ///              use the current random state.
            ///
            void init(const Scalar& mu = Scalar(0), const Scalar& sigma = Scalar(0.01),
                    int seed = -1)
            {
                check_unit_sizes();

                if (seed > 0)
                {
                    m_rng.seed(seed);
                }

                const int nlayer = num_layers();

                for (int i = 0; i < nlayer; i++)
                {
                    m_layers[i]->init(mu, sigma, m_rng);
                }
            }

            ///
            /// Get the serialized layer parameters
            ///
            std::vector< std::vector<Scalar> > get_parameters() const
            {
                const int nlayer = num_layers();
                std::vector< std::vector<Scalar> > res;
                res.reserve(nlayer);

                for (int i = 0; i < nlayer; i++)
                {
                    res.push_back(m_layers[i]->get_parameters());
                }

                return res;
            }

            ///
            /// Set the layer parameters
            ///
            /// \param param Serialized layer parameters
            ///
            void set_parameters(const std::vector< std::vector<Scalar> >& param)
            {
                const int nlayer = num_layers();

                if (static_cast<int>(param.size()) != nlayer)
                {
                    throw std::invalid_argument("[class Network]: Parameter size does not match");
                }

                for (int i = 0; i < nlayer; i++)
                {
                    m_layers[i]->set_parameters(param[i]);
                }
            }

            ///
            /// Get the serialized derivatives of layer parameters
            ///
            std::vector< std::vector<Scalar> > get_derivatives() const
            {
                const int nlayer = num_layers();
                std::vector< std::vector<Scalar> > res;
                res.reserve(nlayer);

                for (int i = 0; i < nlayer; i++)
                {
                    res.push_back(m_layers[i]->get_derivatives());
                }

                return res;
            }

            ///
            /// Debugging tool to check parameter gradients
            ///
            template <typename TargetType>
            void check_gradient(const Matrix& input, const TargetType& target, int npoints,
                                int seed = -1)
            {
                if (seed > 0)
                {
                    m_rng.seed(seed);
                }

                this->forward(input);
                this->backprop(input, target);
                std::vector< std::vector<Scalar> > param = this->get_parameters();
                std::vector< std::vector<Scalar> > deriv = this->get_derivatives();
                const Scalar eps = 1e-5;
                const int nlayer = deriv.size();

                for (int i = 0; i < npoints; i++)
                {
                    // Randomly select a layer
                    const int layer_id = int(m_rng.rand() * nlayer);
                    // Randomly pick a parameter, note that some layers may have no parameters
                    const int nparam = deriv[layer_id].size();

                    if (nparam < 1)
                    {
                        continue;
                    }

                    const int param_id = int(m_rng.rand() * nparam);
                    // Turbulate the parameter a little bit
                    const Scalar old = param[layer_id][param_id];
                    param[layer_id][param_id] -= eps;
                    this->set_parameters(param);
                    this->forward(input);
                    this->backprop(input, target);
                    const Scalar loss_pre = m_output->loss();
                    param[layer_id][param_id] += eps * 2;
                    this->set_parameters(param);
                    this->forward(input);
                    this->backprop(input, target);
                    const Scalar loss_post = m_output->loss();
                    const Scalar deriv_est = (loss_post - loss_pre) / eps / 2;
                    std::cout << "[layer " << layer_id << ", param " << param_id <<
                            "] deriv = " << deriv[layer_id][param_id] << ", est = " << deriv_est <<
                            ", diff = " << deriv_est - deriv[layer_id][param_id] << std::endl;
                    param[layer_id][param_id] = old;
                }

                // Restore original parameters
                this->set_parameters(param);
            }

            ///
            /// Fit the model based on the given data
            ///
            /// \param opt        An object that inherits from the Optimizer class, indicating the optimization algorithm to use.
            /// \param x          The predictors. Each column is an observation.
            /// \param y          The response variable. Each column is an observation.
            /// \param batch_size Mini-batch size.
            /// \param epoch      Number of epochs of training.
            /// \param seed       Set the random seed of the %RNG if `seed > 0`, otherwise
            ///                   use the current random state.
            ///
            template <typename DerivedX, typename DerivedY>
            bool fit(Optimizer& opt, const Eigen::MatrixBase<DerivedX>& x,
                    const Eigen::MatrixBase<DerivedY>& y,
                    int batch_size, int epoch, int seed = -1)
            {
                // We do not directly use PlainObjectX since it may be row-majored if x is passed as mat.transpose()
                // We want to force XType and YType to be column-majored
                typedef typename Eigen::MatrixBase<DerivedX>::PlainObject PlainObjectX;
                typedef typename Eigen::MatrixBase<DerivedY>::PlainObject PlainObjectY;
                typedef Eigen::Matrix<typename PlainObjectX::Scalar, PlainObjectX::RowsAtCompileTime, PlainObjectX::ColsAtCompileTime>
                XType;
                typedef Eigen::Matrix<typename PlainObjectY::Scalar, PlainObjectY::RowsAtCompileTime, PlainObjectY::ColsAtCompileTime>
                YType;
                const int nlayer = num_layers();

                if (nlayer <= 0)
                {
                    return false;
                }

                // Reset optimizer
                opt.reset();

                // Create shuffled mini-batches
                if (seed > 0)
                {
                    m_rng.seed(seed);
                }

                std::vector<XType> x_batches;
                std::vector<YType> y_batches;
                const int nbatch = internal::create_shuffled_batches(x, y, batch_size, m_rng,
                                x_batches, y_batches);
                // Set up callback parameters
                m_callback->m_nbatch = nbatch;
                m_callback->m_nepoch = epoch;

                // Iterations on the whole data set
                for (int k = 0; k < epoch; k++)
                {
                    m_callback->m_epoch_id = k;

                    // Train on each mini-batch
                    for (int i = 0; i < nbatch; i++)
                    {
                        m_callback->m_batch_id = i;
                        m_callback->pre_training_batch(this, x_batches[i], y_batches[i]);
                        this->forward(x_batches[i]);
                        this->backprop(x_batches[i], y_batches[i]);
                        this->update(opt);
                        m_callback->post_training_batch(this, x_batches[i], y_batches[i]);
                    }
                }

                return true;
            }

            ///
            /// Use the fitted model to make predictions
            ///
            /// \param x The predictors. Each column is an observation.
            ///
            Matrix predict(const Matrix& x)
            {
                const int nlayer = num_layers();

                if (nlayer <= 0)
                {
                    return Matrix();
                }

                this->forward(x);
                return m_layers[nlayer - 1]->output();
            }

            ///
            /// Export the network to files.
            ///
            /// \param folder   The folder where the network is saved.
            /// \param fileName The filename for the network.
            ///
            void export_net(const std::string& folder, const std::string& filename) const
            {
                bool created = internal::create_directory(folder);
                if (!created)
                    throw std::runtime_error("[class Network]: Folder creation failed");

                MetaInfo map = this->get_meta_info();
                internal::write_map(folder + "/" + filename, map);
                std::vector< std::vector<Scalar> > params = this->get_parameters();
                internal::write_parameters(folder, filename, params);
            }

            ///
            /// Read in a network from files.
            ///
            /// \param folder   The folder where the network is saved.
            /// \param fileName The filename for the network.
            ///
            void read_net(const std::string& folder, const std::string& filename)
            {
                MetaInfo map;
                internal::read_map(folder + "/" + filename, map);
                int nlayer = map.find("Nlayers")->second;
                std::vector< std::vector<Scalar> > params = internal::read_parameters(folder, filename, nlayer);
                m_layers.clear();

                for (int i = 0; i < nlayer; i++)
                {
                    this->add_layer(internal::create_layer(map, i));
                }

                this->set_parameters(params);
                this->set_output(internal::create_output(map));
            }
    };
        ///
    /// \ingroup Callbacks
    ///
    /// Callback function that prints the loss function value in each mini-batch training
    ///
    class VerboseCallback: public Callback
    {
        public:
            void post_training_batch(const Network* net, const Matrix& x, const Matrix& y)
            {
                const Scalar loss = net->get_output()->loss();
                std::cout << "[Epoch " << m_epoch_id << ", batch " << m_batch_id << "] Loss = "
                        << loss << std::endl;
            }

            void post_training_batch(const Network* net, const Matrix& x,
                                    const IntegerVector& y)
            {
                Scalar loss = net->get_output()->loss();
                std::cout << "[Epoch " << m_epoch_id << ", batch " << m_batch_id << "] Loss = "
                        << loss << std::endl;
            }
    };


} // namespace MiniDNN
