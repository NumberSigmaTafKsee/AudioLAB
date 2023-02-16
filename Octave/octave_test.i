%module octave_test 
%{
    #include <vector>
%}

%include "std_math.i" 
%include "std_vector.i"

%template(float_vector) std::vector<float>;
%template(double_vector) std::vector<double>;

%typemap(in) FloatRowVector(float *in, size_t n) %{
    for (size_t i = 0; i < n; i++)
        $1->float_row_vector_value()(i) = in[i];
%}
%typemap(in) RowVector(double *in, size_t n) %{
    for (size_t i = 0; i < n; i++)
        $1->row_vector_value()(i) = in[i];
    %}
%typemap(in) FloatComplexRowVector(std::complex<float> *in, size_t n) %{
    for (size_t i = 0; i < n; i++)
        $1->float_complex_row_vector_value()(i) = in[i];
    %}
%typemap(in) ComplexRowVector(std::complex<double> *in, size_t n) %{
    for (size_t i = 0; i < n; i++)
        $1->complex_row_vector_value()(i) = in[i];
    %}

%typemap(in) FloatColumnVector(float *in, size_t n) %{
    for (size_t i = 0; i < n; i++)
        $1->float_column_vector_value()(i) = in[i];
    %}
%typemap(in) FloatComplexColumnVector(std::complex<float> *in, size_t n) %{
    for (size_t i = 0; i < n; i++)
        $1->float_complex_column_vector_value()(i) = in[i];
    %}
%typemap(in) ComplexColumnVector(std::complex<double> *in, size_t n) %{
    for (size_t i = 0; i < n; i++)
        $1->complex_row_vector_value()(i) = in[i];
    %}
%typemap(in) ColumnVector(double *in, size_t n) %{
    for (size_t i = 0; i < n; i++)
        $1->column_vector_value()(i) = in[i];
    %}

%typemap(in) FloatRowVector(std::vector<float> in) %{
    for (size_t i = 0; i < in.size(); i++)
        $1->float_row_vector_value()(i) = in[i];
    %}
%typemap(in) RowVector(std::vector<double> in) %{
    for (size_t i = 0; i < in.size(); i++)
        $1->row_vector_value()(i) = in[i];
    %}
%typemap(in) FloatComplexRowVector(std::vector<std::complex<float>> in) %{
    for (size_t i = 0; i < in.size(); i++)
        $1->float_complex_row_vector_value()(i) = in[i];
    %}
%typemap(in) ComplexRowVector(std::vector<std::complex<double>> in) %{
    for (size_t i = 0; i < in.size(); i++)
        $1->complex_row_vector_value()(i) = in[i];
    %}

%typemap(in) FloatColumnVector(std::vector<std::complex<float>> in) %{
    for (size_t i = 0; i < in.size(); i++)
        $1->float_complex_column_vector_value()(i) = in[i];
    %}
%typemap(in) ColumnVector(std::vector<std::complex<double>> in) %{
    for (size_t i = 0; i < in.size(); i++)
        $1->column_complex_vector_value()(i) = in[i];
    %}
%typemap(in) FloatComplexColumnVector(std::vector<std::complex<float>> in) %{
    for (size_t i = 0; i < in.size(); i++)
        $1->float_complex_row_vector_value()(i) = in[i];
    %}
%typemap(in) ComplexColumnVector(std::vector<std::complex<double>> in) %{
    for (size_t i = 0; i < in.size(); i++)
        $1->complex_row_vector_value()(i) = in[i];
    %}

%typemap(out) std::vector<float> in %{
    $result = std::vector<float>(in.size());
    for (size_t i = 0; i < in.size(1); i++)
        $result->float_row_vector()(i) = in(i);
    %}
%typemap(out) std::vector<double> in %{
    $result = std::vector<double>(in.size(1));
    for (size_t i = 0; i < in.size(1); i++)
        $result->row_vector()(i) = in(i);
    %}

%typemap(in) FloatMatrix(float *in, size_t m, size_t n) %{
    for (size_t i = 0; i < m; i++)
        for (size_t j = 0; j < n; j++)
            $1->float_matrix_value()(i, j) = in[i * n + j];
    %}
%typemap(in) Matrix(double *in, size_t m, size_t n) %{
    for (size_t i = 0; i < m; i++)
        for (size_t j = 0; j < n; j++)
            $1->matrix_value()(i, j) = in[i * n + j];
    %}
%typemap(in) FloatComplexMatrix(std::complex<float> *in, size_t m, size_t n) %{
    for (size_t i = 0; i < m; i++)
        for (size_t j = 0; j < n; j++)
            $1->float_complex_matrix_value()(i, j) = in[i * n + j];
    %}
%typemap(in) ComplexMatrix(std::complex<double> *in, size_t m, size_t n) %{
    for (size_t i = 0; i < m; i++)
        for (size_t j = 0; j < n; j++)
            $1->complex_matrix_value()(i, j) = in[i * n + j];
    %}

%inline %{

    octave_value vectorize_row(const std::vector<float> &src)
    {
        FloatRowVector dst(src.size());
        for (size_t i = 0; i < src.size(); i++)
            dst(i) = src[i];
        return octave_value(dst);
    }
    octave_value vectorize_row(const float *src, size_t n)
    {
        FloatRowVector dst(src.size());
        for (size_t i = 0; i < src.size(); i++)
            dst(i) = src[i];
        return octave_value(dst);
    }
    std::vector<float> vectorize_row_float_vector(const octave_value &src)
    {
        FloatRowVector tmp = src.float_row_vector_value();
        std::vector<float> dst(tmp.size(1));
        for (size_t i = 0; i < tmp.size(1); i++)
            dst[i] = tmp(i);
        return dst;
    }

    octave_value vectorize_complex_row(const std::vector<std::complex<float>> &src)
    {
        ComplexFloatRowVector dst(src.size());
        for (size_t i = 0; i < src.size(); i++)
            dst(i) = src[i];
        return octave_value(dst);
    }
    octave_value vectorize_complex_row(const std::complex<float> *src, size_t n)
    {
        ComplexFloatRowVector dst(src.size());
        for (size_t i = 0; i < src.size(); i++)
            dst(i) = src[i];
        return octave_value(dst);
    }
    std::vector<std::complex<float>> vectorize_row_complex_float_vector(const octave_value &src)
    {
        ComplexFloatRowVector tmp = src.float_complex_row_vector_value();
        std::vector<float> dst(tmp.size(1));
        for (size_t i = 0; i < tmp.size(1); i++)
            dst[i] = tmp(i);
        return dst;
    }

    octave_value vectorize_row(const std::vector<double> &src)
    {
        RowVector dst(src.size());
        for (size_t i = 0; i < src.size(); i++)
            dst(i) = src[i];
        return octave_value(dst);
    }
    octave_value vectorize_row(const double *src, size_t n)
    {
        RowVector dst(src.size());
        for (size_t i = 0; i < src.size(); i++)
            dst(i) = src[i];
        return octave_value(dst);
    }
    std::vector<double> vectorize_row_double_vector(const octave_value &src)
    {
        RowVector tmp = src.row_vector_value();
        std::vector<double> dst(tmp.size(1));
        for (size_t i = 0; i < tmp.size(1); i++)
            dst[i] = tmp(i);
        return dst;
    }

    octave_value vectorize_complex_row(const std::vector<std::complex<double>> &src)
    {
        ComplexRowVector dst(src.size());
        for (size_t i = 0; i < src.size(); i++)
            dst(i) = src[i];
        return octave_value(dst);
    }
    octave_value vectorize_complex_row(const std::complex<double> *src, size_t n)
    {
        ComplexRowVector dst(src.size());
        for (size_t i = 0; i < src.size(); i++)
            dst(i) = src[i];
        return octave_value(dst);
    }
    std::vector<std::complex<double>> vectorize_complex_row_double_vector(const octave_value &src)
    {
        ComplexRowVector tmp = src.complex_row_vector_value();
        std::vector<double> dst(tmp.size(1));
        for (size_t i = 0; i < tmp.size(1); i++)
            dst[i] = tmp(i);
        return dst;
    }

    octave_value vectorize_col(const std::vector<float> &src)
    {
        FloatColumnVector dst(src.size());
        for (size_t i = 0; i < src.size(); i++)
            dst(i) = src[i];
        return octave_value(dst);
    }
    octave_value vectorize_col(const float *src, size_t n)
    {
        FloatColumnVector dst(src.size());
        for (size_t i = 0; i < src.size(); i++)
            dst(i) = src[i];
        return octave_value(dst);
    }
    std::vector<float> vectorize_col_float_vector(const octave_value &src)
    {
        FloatColumnVector tmp = src.float_column_vector_value();
        std::vector<float> dst(tmp.size(1));
        for (size_t i = 0; i < tmp.size(1); i++)
            dst[i] = tmp(i);
        return dst;
    }

    octave_value vectorize_complex_col(const std::vector<std::complex<float>> &src)
    {
        ComplexFloatColumnVector dst(src.size());
        for (size_t i = 0; i < src.size(); i++)
            dst(i) = src[i];
        return octave_value(dst);
    }
    octave_value vectorize_complex_col(const std::complex<float> *src, size_t n)
    {
        ComplexFloatColumnVector dst(src.size());
        for (size_t i = 0; i < src.size(); i++)
            dst(i) = src[i];
        return octave_value(dst);
    }
    std::vector<std::complex<float>> vectorize_complex_col_float_vector(const octave_value &src)
    {
        ComplexFloatColumnVector tmp = src.float_complex_column_vector_value();
        std::vector<float> dst(tmp.size(1));
        for (size_t i = 0; i < tmp.size(1); i++)
            dst[i] = tmp(i);
        return dst;
    }

    octave_value vectorize_col(const std::vector<double> &src)
    {
        ColumnVector dst(src.size());
        for (size_t i = 0; i < src.size(); i++)
            dst(i) = src[i];
        return octave_value(dst);
    }
    octave_value vectorize_col(const double *src, size_t n)
    {
        ColumnRowVector dst(src.size());
        for (size_t i = 0; i < src.size(); i++)
            dst(i) = src[i];
        return octave_value(dst);
    }
    std::vector<double> vectorize_col_double_vector(const octave_value &src)
    {
        ColumnVector tmp = src.column_vector_value();
        std::vector<double> dst(tmp.size(1));
        for (size_t i = 0; i < tmp.size(1); i++)
            dst[i] = tmp(i);
        return dst;
    }

    octave_value vectorize_complex_col(const std::vector<std::complex<double>> &src)
    {
        ComplexColumnVector dst(src.size());
        for (size_t i = 0; i < src.size(); i++)
            dst(i) = src[i];
        return octave_value(dst);
    }
    octave_value vectorize_complex_col(const std::complex<double> *src, size_t n)
    {
        ColumnRowVector dst(src.size());
        for (size_t i = 0; i < src.size(); i++)
            dst(i) = src[i];
        return octave_value(dst);
    }
    std::vector<std::complex<double>> vectorize_complex_col_double_vector(const octave_value &src)
    {
        ComplexColumnVector tmp = src.complex_column_vector_value();
        std::vector<std::complex<double>> dst(tmp.size(1));
        for (size_t i = 0; i < tmp.size(1); i++)
            dst[i] = tmp(i);
        return dst;
    }

    octave_value matricize(const std::vector<float> &src, size_t m, size_t n)
    {
        FloatMatrix dst(m, n);
        for (size_t i = 0; i < m; i++)
            for (size_t j = 0; j < n; j++)
                dst(i, j) = src[i * n + j];
        return octave_value(dst);
    }
    octave_value matricize(const float *src, size_t m, size_t n)
    {
        FloatMatrix dst(m, n);
        for (size_t i = 0; i < m; i++)
            for (size_t j = 0; j < n; j++)
                dst(i, j) = src[i * n + j];
        return octave_value(dst);
    }
    std::vector<float> matricize_float_vector(const octave_value &src)
    {
        FloatMatrix tmp = src.float_matrix_value();
        std::vector<float> dst(tmp.rows(), tmp.cols());
        for (size_t i = 0; i < m; i++)
            for (size_t j = 0; j < n; j++)
                dst[i * src.cols() + j] = tmp(i, j);
        return dst;
    }

    octave_value matricize(const std::vector<double> &src, size_t m, size_t n)
    {
        Matrix dst(m, n);
        for (size_t i = 0; i < m; i++)
            for (size_t j = 0; j < n; j++)
                dst(i, j) = src[i * n + j];
        return octave_value(dst);
    }
    octave_value matricize(const double *src, size_t m, size_t n)
    {
        Matrix dst(m, n);
        for (size_t i = 0; i < m; i++)
            for (size_t j = 0; j < n; j++)
                dst(i, j) = src[i * n + j];
        return octave_value(dst);
    }
    std::vector<double> matricize_double_vector(const octave_value &src)
    {
        Matrix tmp = src.matrix_value();
        std::vector<double> dst(tmp.rows(), tmp.cols());
        for (size_t i = 0; i < m; i++)
            for (size_t j = 0; j < n; j++)
                dst[i * src.cols() + j] = tmp(i, j);
        return dst;
    }


    octave_value complex_matricize(const std::vector<std::complex<float>> &src, size_t m, size_t n)
    {
        ComplexFloatMatrix dst(m, n);
        for (size_t i = 0; i < m; i++)
            for (size_t j = 0; j < n; j++)
                dst(i, j) = src[i * n + j];
        return octave_value(dst);
    }
    octave_value complex_matricize(const std::complex<float> *src, size_t m, size_t n)
    {
        ComplexFloatMatrix dst(m, n);
        for (size_t i = 0; i < m; i++)
            for (size_t j = 0; j < n; j++)
                dst(i, j) = src[i * n + j];
        return octave_value(dst);
    }
    std::vector<std::complex<floa>t> complex_matricize_float_vector(const octave_value &src)
    {
        ComplexFloatMatrix tmp = src.float_matrix_value();
        std::vector<std::complex<float>> dst(tmp.rows(), tmp.cols());
        for (size_t i = 0; i < m; i++)
            for (size_t j = 0; j < n; j++)
                dst[i * src.cols() + j] = tmp(i, j);
        return dst;
    }

    octave_value complex_matricize(const std::vector < std::complex<double> & src, size_t m, size_t n)
    {
        ComplexMatrix dst(m, n);
        for (size_t i = 0; i < m; i++)
            for (size_t j = 0; j < n; j++)
                dst(i, j) = src[i * n + j];
        return octave_value(dst);
    }
    octave_value complex_matricize(const std::complex<double> *src, size_t m, size_t n)
    {
        ComplexMatrix dst(m, n);
        for (size_t i = 0; i < m; i++)
            for (size_t j = 0; j < n; j++)
                dst(i, j) = src[i * n + j];
        return octave_value(dst);
    }
    std::vector<std::complex<double>> complex_matricize_double_vector(const octave_value &src)
    {
        ComplexMatrix tmp = src.complex_matrix_value();
        std::vector<std::complex<double>> dst(tmp.rows(), tmp.cols());
        for (size_t i = 0; i < m; i++)
            for (size_t j = 0; j < n; j++)
                dst[i * src.cols() + j] = tmp(i, j);
        return dst;
    }
%}