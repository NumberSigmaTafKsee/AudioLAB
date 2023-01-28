%module TWaveTable
%{
#include "TWaveTable.hpp"
%}

%include "stdint.i"
%include "std_vector.i"

%template(float_vector) std::vector<float>;
%template(double_vector) std::vector<double>;


%include "TWaveTable.hpp"

%template (FloatWaveTable) SoundWave::TWaveTable<float>;
%template (DoubleWaveTable) SoundWave::TWaveTable<double>;

%template (MakeSinef) SoundWave::MakeSine<float>;
%template (MakeSquaref) SoundWave::MakeSquare<float>;
%template (MakeSawf) SoundWave::MakeSaw<float>;
%template (MakeReverseSawf) SoundWave::MakeReverseSaw<float>;
%template (MakeTrianglef) SoundWave::MakeTriangle<float>;

%template (MakeSined) SoundWave::MakeSine<double>;
%template (MakeSquared) SoundWave::MakeSquare<double>;
%template (MakeSawd) SoundWave::MakeSaw<double>;
%template (MakeReverseSawd) SoundWave::MakeReverseSaw<double>;
%template (MakeTriangled) SoundWave::MakeTriangle<double>;
