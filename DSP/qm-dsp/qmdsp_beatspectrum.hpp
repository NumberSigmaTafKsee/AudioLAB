#pragma once

#include "qmdsp_math.hpp"

namespace qmdsp
{
    /**
    * Given a matrix of "feature values", calculate a self-similarity
    * vector.  The resulting vector will have half as many elements as
    * the number of columns in the matrix.  This is based on the
    * SoundBite rhythmic similarity code.
    */

    class BeatSpectrum
    {
    public:
        BeatSpectrum() { }
        ~BeatSpectrum() { }

        std::vector<double> process(const std::vector<std::vector<double> > &inmatrix);

    };

    std::vector<double> BeatSpectrum::process(const std::vector<std::vector<double> > &m)
    {
        int origin = 0;
        int sz = m.size()/2;

        int i, j, k;

        std::vector<double> v(sz);
        for (i = 0; i < sz; ++i) v[i] = 0.0;

        CosineDistance cd;

        for (i = origin; i < origin + sz; ++i) {

            k = 0;

            for (j = i + 1; j < i + sz + 1; ++j) {

                v[k++] += cd.distance(m[i], m[j]);
            }
        }

        // normalize

        double max = 0.0;

        for (i = 0; i < sz; ++i) {
            if (v[i] > max) max = v[i];
        }

        if (max > 0.0) {
            for (i = 0; i < sz; ++i) {
                v[i] /= max;
            }
        }

        return v;
    }
}