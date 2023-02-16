#include <vector>
#include <functional>
#include <iostream>

std::vector<float> Sequence(int start, int end, std::function<float (int i)> func)
{
    int cnt = end - start;
    std::vector<float> v(cnt+1);
    for(int i = 0; i <= cnt; i++) v[i] = func(i);
    return v;
}

float foo(int i) {
    return (float)i/256.0;
}
std::ostream& operator << (std::ostream& o, const std::vector<float> & v)
{
    o << "vector[" << v.size() << "]=";
    for(size_t i = 0; i < v.size(); i++) o << v[i] << ",";
    o << std::endl;
    return o;
}
int main() {
    std::vector<float> v = Sequence(1,256,foo);
    std::cout << v << std::endl;
}