%{
#include <numeric>
%}

namespace std
{
    enum class endian
    {
        little = /*implementation-defined*/,
        big    = /*implementation-defined*/,
        native = /*implementation-defined*/
    };

    template<class T> constexpr T midpoint(T a, T b);
    template<class T> constexpr T midpoint(T *a, T *b);
    constexpr float       lerp( float a, float b, float t ) noexcept;
    constexpr double      lerp( double a, double b, double t ) noexcept;

    template<class T> void iota(std::vector<T> & x, T value) {
        std::iota(x.begin(),x.end(),value);
    }
    template<class T> T accumulate(std::vector<T> & x, T value) {
        return std::accumulate(x.begin(),x.end(),value);
    }
    template<class T> T reduce(std::vector<T> & x, T value) {
        return std::reduce(x.begin(),x.end(),value);
    }
    template<class T> T transform_reduce(std::vector<T> & x, std::vector<T> & y, T value) {
        return std::trandform_reduce(x.begin(),x.end(),y.begin(), value);
    }
    template<class T> void inner_product(std::vector<T> & x,  std::vector<T> & y, T value) {
        std::inner_product(x.begin(),x.end(),y.begin(),value);
    }
    template<class T> void adjacent_difference(std::vector<T> & x, std::vector<T> & y) {
        std::adjacent_difference(x.begin(),x.end(),y.begin());
    }
    template<class T> void partial_sum(std::vector<T> & x, std::vector<T> & y) {
        std::partial_sum(x.begin(),x.end(),y.begin());
    }
    template<class T> void inclusive_scan(std::vector<T> & x, std::vector<T> & y) {
        std::inclusive_scan(x.begin(),x.end(),y.begin());
    }
    template<class T> void exclusive_scan(std::vector<T> & x, std::vector<T> & y, T init) {
        std::exclusive_scan(x.begin(),x.end(),y.begin(), init);
    }

    template< class T > constexpr T byteswap( T n ) noexcept
    template< class T > constexpr bool has_single_bit( T x ) noexcept;
    template< class T > constexpr T bit_ceil( T x );
    template< class T > constexpr T bit_floor( T x ) noexcept;
    template< class T > constexpr T bit_width( T x ) noexcept;    
    template< class T > constexpr T rotl( T x, int s ) noexcept;
    template< class T > constexpr T rotr( T x, int s ) noexcept;
    template< class T > constexpr int countl_zero( T x ) noexcept;
    template< class T > constexpr int countl_one( T x ) noexcept;
    template< class T > constexpr int countr_zero( T x ) noexcept;
    template< class T > constexpr int countr_one( T x ) noexcept;
    template< class T > constexpr int popcount( T x ) noexcept;

}