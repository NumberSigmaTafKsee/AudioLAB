#pragma once

#include <complex>
#include <ffts/ffts.h>
#include <stdexcept>
#include <vector>
#include <algorithm>
#include <complex>
#include <numeric>
#include <vector>
#include <cfloat>
#include <cmath>
#include <map>
#include <sstream>
#include <string>
#include <climits>
#include <mlpack/core.hpp>
#include <mlpack/methods/hmm/hmm.hpp>
#include <tuple>


/* ignore me plz */
namespace detail
{
template <typename T>
std::vector<size_t>
bin_pitches(const std::vector<std::pair<T, T>>);

mlpack::hmm::HMM<mlpack::distribution::DiscreteDistribution>
build_hmm();

void
init_pitch_bins();
} // namespace detail

/*
 * The pitch namespace contains the functions:
 *
 * 	pitch::mpm(data, sample_rate)
 * 	pitch::yin(data, sample_rate)
 * 	pitch::pyin(data, sample_rate)
 * 	pitch::pmpm(data, sample_rate)
 *
 * It will auto-allocate any buffers.
 */
namespace pitch
{

template <typename T>
T
yin(const std::vector<T> &, int);

template <typename T>
T
mpm(const std::vector<T> &, int);

template <typename T>
T
swipe(const std::vector<T> &, int);

/*
 * pyin and pmpm emit pairs of pitch/probability
 */
template <typename T>
T
pyin(const std::vector<T> &, int);

template <typename T>
T
pmpm(const std::vector<T> &, int);
} // namespace pitch

/*
 * This namespace is useful for repeated calls to pitch for the same size of
 * buffer.
 *
 * It contains the classes Yin and Mpm which contain the allocated buffers
 * and each implement a `pitch(data, sample_rate)` and
 * `probablistic_pitch(data, sample_rate)` method.
 */
namespace pitch_alloc
{

template <typename T> class BaseAlloc
{
  public:
	long N;
	std::vector<std::complex<float>> out_im;
	std::vector<T> out_real;
	ffts_plan_t *fft_forward;
	ffts_plan_t *fft_backward;
	mlpack::hmm::HMM<mlpack::distribution::DiscreteDistribution> hmm;

	BaseAlloc(long audio_buffer_size)
	    : N(audio_buffer_size), out_im(std::vector<std::complex<float>>(N * 2)),
	      out_real(std::vector<T>(N))
	{
		if (N == 0) {
			throw std::bad_alloc();
		}

		fft_forward = ffts_init_1d(N * 2, FFTS_FORWARD);
		fft_backward = ffts_init_1d(N * 2, FFTS_BACKWARD);
		detail::init_pitch_bins();
		hmm = detail::build_hmm();
	}

	~BaseAlloc()
	{
		ffts_free(fft_forward);
		ffts_free(fft_backward);
	}

  protected:
	void
	clear()
	{
		std::fill(out_im.begin(), out_im.end(), std::complex<T>(0.0, 0.0));
	}
};

/*
 * Allocate the buffers for MPM for re-use.
 * Intended for multiple consistently-sized audio buffers.
 *
 * Usage: pitch_alloc::Mpm ma(1024)
 *
 * It will throw std::bad_alloc for invalid sizes (<1)
 */
template <typename T> class Mpm : public BaseAlloc<T>
{
  public:
	Mpm(long audio_buffer_size) : BaseAlloc<T>(audio_buffer_size){};

	T
	pitch(const std::vector<T> &, int);

	T
	probabilistic_pitch(const std::vector<T> &, int);
};

/*
 * Allocate the buffers for YIN for re-use.
 * Intended for multiple consistently-sized audio buffers.
 *
 * Usage: pitch_alloc::Yin ya(1024)
 *
 * It will throw std::bad_alloc for invalid sizes (<2)
 */
template <typename T> class Yin : public BaseAlloc<T>
{
  public:
	std::vector<T> yin_buffer;

	Yin(long audio_buffer_size)
	    : BaseAlloc<T>(audio_buffer_size),
	      yin_buffer(std::vector<T>(audio_buffer_size / 2))
	{
		if (audio_buffer_size / 2 == 0) {
			throw std::bad_alloc();
		}
	}

	T
	pitch(const std::vector<T> &, int);

	T
	probabilistic_pitch(const std::vector<T> &, int);
};
} // namespace pitch_alloc

namespace util
{
template <typename T>
std::pair<T, T>
parabolic_interpolation(const std::vector<T> &, int);

template <typename T>
void
acorr_r(const std::vector<T> &, pitch_alloc::BaseAlloc<T> *);

template <typename T>
T
pitch_from_hmm(mlpack::hmm::HMM<mlpack::distribution::DiscreteDistribution>,
    const std::vector<std::pair<T, T>>);
} // namespace util


template <typename T>
void
util::acorr_r(const std::vector<T> &audio_buffer, pitch_alloc::BaseAlloc<T> *ba)
{
	if (audio_buffer.size() == 0)
		throw std::invalid_argument("audio_buffer shouldn't be empty");

	std::transform(audio_buffer.begin(), audio_buffer.begin() + ba->N,
	    ba->out_im.begin(), [](T x) -> std::complex<T> {
		    return std::complex<float>(x, static_cast<T>(0.0));
	    });

	ffts_execute(ba->fft_forward, ba->out_im.data(), ba->out_im.data());

	std::complex<float> scale = {
	    1.0f / (float)(ba->N * 2), static_cast<T>(0.0)};
	for (int i = 0; i < ba->N; ++i)
		ba->out_im[i] *= std::conj(ba->out_im[i]) * scale;

	ffts_execute(ba->fft_backward, ba->out_im.data(), ba->out_im.data());

	std::transform(ba->out_im.begin(), ba->out_im.begin() + ba->N,
	    ba->out_real.begin(),
	    [](std::complex<float> cplx) -> T { return std::real(cplx); });
}

template void
util::acorr_r<double>(const std::vector<double> &audio_buffer,
    pitch_alloc::BaseAlloc<double> *ba);

template void
util::acorr_r<float>(
    const std::vector<float> &audio_buffer, pitch_alloc::BaseAlloc<float> *ba);


#define F0 440.0
#define N_BINS 108
#define N_NOTES 12
#define NOTE_OFFSET 57

#define YIN_TRUST 0.5

#define TRANSITION_WIDTH 13
#define SELF_TRANS 0.99

std::vector<double> PITCH_BINS(N_BINS);
std::vector<double> REAL_PITCHES(N_BINS);

const double A = std::pow(2.0, 1.0 / 12.0);

// 108 bins - C0 -> B8
void
detail::init_pitch_bins()
{
	for (int i = 0; i < N_BINS; ++i) {
		auto fi = F0 * std::pow(A, i - NOTE_OFFSET);
		PITCH_BINS[i] = fi;
	}
}

template <typename T>
std::vector<size_t>
detail::bin_pitches(const std::vector<std::pair<T, T>> pitch_candidates)
{
	arma::vec pitch_probs(2 * N_BINS + 1, arma::fill::zeros);
	std::vector<size_t> possible_bins;

	T prob_pitched = 0.0;

	for (auto pitch_candidate : pitch_candidates) {
		// find the most appropriate bin
		T delta = DBL_MAX;
		T prev_delta = DBL_MAX;
		for (int i = 0; i < N_BINS; ++i) {
			delta = std::abs(pitch_candidate.first - PITCH_BINS[i]);
			if (prev_delta < delta) {
				pitch_probs[i - 1] = pitch_candidate.second;
				prob_pitched += pitch_probs[i - 1];
				REAL_PITCHES[i - 1] = pitch_candidate.first;
				break;
			}
			prev_delta = delta;
		}
	}

	T prob_really_pitched = YIN_TRUST * prob_pitched;

	for (int i = 0; i < N_BINS; ++i) {
		if (prob_pitched > 0) {
			pitch_probs[i] *= prob_really_pitched / prob_pitched;
		}
		pitch_probs[i + N_BINS] = (1 - prob_really_pitched) / N_BINS;
	}

	for (size_t i = 0; i < pitch_probs.size(); ++i) {
		auto pitch_probability = pitch_probs[i];
		for (size_t j = 0; j < size_t(100.0 * pitch_probability); ++j)
			possible_bins.push_back(i);
	}

	return possible_bins;
}

mlpack::hmm::HMM<mlpack::distribution::DiscreteDistribution>
detail::build_hmm()
{
	size_t hmm_size = 2 * N_BINS + 1;
	// initial
	arma::vec initial(hmm_size);
	initial.fill(1.0 / double(hmm_size));

	arma::mat transition(hmm_size, hmm_size, arma::fill::zeros);

	// transitions
	for (int i = 0; i < N_BINS; ++i) {
		int half_transition = static_cast<int>(TRANSITION_WIDTH / 2.0);
		int theoretical_min_next_pitch = i - half_transition;
		int min_next_pitch = i > half_transition ? i - half_transition : 0;
		int max_next_pitch =
		    i < N_BINS - half_transition ? i + half_transition : N_BINS - 1;

		double weight_sum = 0.0;
		std::vector<double> weights;

		for (int j = min_next_pitch; j <= max_next_pitch; ++j) {
			if (j <= i) {
				weights.push_back(j - theoretical_min_next_pitch + 1);
			} else {
				weights.push_back(i - theoretical_min_next_pitch + 1 - j + i);
			}
			weight_sum += weights[weights.size() - 1];
		}

		for (int j = min_next_pitch; j <= max_next_pitch; ++j) {
			transition(i, j) =
			    (weights[j - min_next_pitch] / weight_sum * SELF_TRANS);
			transition(i, j + N_BINS) =
			    (weights[j - min_next_pitch] / weight_sum * (1.0 - SELF_TRANS));
			transition(i + N_BINS, j + N_BINS) =
			    (weights[j - min_next_pitch] / weight_sum * SELF_TRANS);
			transition(i + N_BINS, j) =
			    (weights[j - min_next_pitch] / weight_sum * (1.0 - SELF_TRANS));
		}
	}

	// the only valid emissions are exact notes
	// i.e. an identity matrix of emissions
	std::vector<mlpack::distribution::DiscreteDistribution> emissions(hmm_size);

	for (size_t i = 0; i < hmm_size; ++i) {
		emissions[i] = mlpack::distribution::DiscreteDistribution(
		    std::vector<arma::vec>{arma::vec(hmm_size, arma::fill::zeros)});
		emissions[i].Probabilities()[i] = 1.0;
	}

	auto hmm = mlpack::hmm::HMM<mlpack::distribution::DiscreteDistribution>(initial, transition, emissions);
	return hmm;
}

template <typename T>
T
util::pitch_from_hmm(
    mlpack::hmm::HMM<mlpack::distribution::DiscreteDistribution> hmm,
    const std::vector<std::pair<T, T>> pitch_candidates)
{
	if (pitch_candidates.size() == 0) {
		return -1.0;
	}

	std::vector<T> observation_;

	for (auto obs_to_add : detail::bin_pitches(pitch_candidates)) {
		observation_.push_back(obs_to_add);
	}

	if (observation_.size() == 0) {
		return -1.0;
	}

	arma::mat observation(1, observation_.size());
	for (size_t i = 0; i < observation_.size(); ++i) {
		observation(0, i) = observation_[i];
	}

	arma::Row<size_t> state;
	// auto viterbi_out = hmm.Predict(observation, state);
	hmm.Predict(observation, state);

	// count state with most appearances
	std::map<size_t, size_t> counts;
	for (auto state_ : state) {
		counts[state_]++;
	}

	size_t most_frequent;
	size_t max = 0;
	for (auto map_pair : counts) {
		auto state_ = map_pair.first;
		auto count = map_pair.second;
		if (count > max) {
			most_frequent = state_;
			max = count;
		}
	}

	return REAL_PITCHES[most_frequent];
}

template double
util::pitch_from_hmm<double>(
    mlpack::hmm::HMM<mlpack::distribution::DiscreteDistribution> hmm,
    const std::vector<std::pair<double, double>> pitch_candidates);

template float
util::pitch_from_hmm<float>(
    mlpack::hmm::HMM<mlpack::distribution::DiscreteDistribution> hmm,
    const std::vector<std::pair<float, float>> pitch_candidates);


#define MPM_CUTOFF 0.93
#define MPM_SMALL_CUTOFF 0.5
#define MPM_LOWER_PITCH_CUTOFF 80.0
#define PMPM_PA 0.01
#define PMPM_N_CUTOFFS 20
#define PMPM_PROB_DIST 0.05
#define PMPM_CUTOFF_BEGIN 0.8
#define PMPM_CUTOFF_STEP 0.01

template <typename T>
static std::vector<int>
peak_picking(const std::vector<T> &nsdf)
{
	std::vector<int> max_positions{};
	int pos = 0;
	int cur_max_pos = 0;
	ssize_t size = nsdf.size();

	while (pos < (size - 1) / 3 && nsdf[pos] > 0)
		pos++;
	while (pos < size - 1 && nsdf[pos] <= 0.0)
		pos++;

	if (pos == 0)
		pos = 1;

	while (pos < size - 1) {
		if (nsdf[pos] > nsdf[pos - 1] && nsdf[pos] >= nsdf[pos + 1] &&
		    (cur_max_pos == 0 || nsdf[pos] > nsdf[cur_max_pos])) {
			cur_max_pos = pos;
		}
		pos++;
		if (pos < size - 1 && nsdf[pos] <= 0) {
			if (cur_max_pos > 0) {
				max_positions.push_back(cur_max_pos);
				cur_max_pos = 0;
			}
			while (pos < size - 1 && nsdf[pos] <= 0.0) {
				pos++;
			}
		}
	}
	if (cur_max_pos > 0) {
		max_positions.push_back(cur_max_pos);
	}
	return max_positions;
}

template <typename T>
T
pitch_alloc::Mpm<T>::probabilistic_pitch(
    const std::vector<T> &audio_buffer, int sample_rate)
{
	util::acorr_r(audio_buffer, this);

	std::map<T, T> t0_with_probability;
	std::vector<std::pair<T, T>> f0_with_probability;

	T cutoff = PMPM_CUTOFF_BEGIN;

	for (int n = 0; n < PMPM_N_CUTOFFS; ++n) {
		std::vector<int> max_positions = peak_picking(this->out_real);
		std::vector<std::pair<T, T>> estimates;

		T highest_amplitude = -DBL_MAX;

		for (int i : max_positions) {
			highest_amplitude = std::max(highest_amplitude, this->out_real[i]);
			if (this->out_real[i] > MPM_SMALL_CUTOFF) {
				auto x = util::parabolic_interpolation(this->out_real, i);
				estimates.push_back(x);
				highest_amplitude = std::max(highest_amplitude, std::get<1>(x));
			}
		}

		if (estimates.empty())
			continue;

		T actual_cutoff = cutoff * highest_amplitude;
		T period = 0;

		for (auto i : estimates) {
			if (std::get<1>(i) >= actual_cutoff) {
				period = std::get<0>(i);
				break;
			}
		}

		auto a = period != 0 ? 1 : PMPM_PA;

		t0_with_probability[period] += a * PMPM_PROB_DIST;

		cutoff += MPM_CUTOFF;
	}

	for (auto tau_estimate : t0_with_probability) {
		if (tau_estimate.first == 0.0) {
			continue;
		}
		auto f0 = (sample_rate / tau_estimate.first);

		f0 = (f0 > MPM_LOWER_PITCH_CUTOFF) ? f0 : -1;

		if (f0 != -1.0) {
			f0_with_probability.push_back(
			    std::make_pair(f0, tau_estimate.second));
		}
	}
	this->clear();

	return util::pitch_from_hmm(this->hmm, f0_with_probability);
}

template <typename T>
T
pitch_alloc::Mpm<T>::pitch(const std::vector<T> &audio_buffer, int sample_rate)
{
	util::acorr_r(audio_buffer, this);

	std::vector<int> max_positions = peak_picking(this->out_real);
	std::vector<std::pair<T, T>> estimates;

	T highest_amplitude = -DBL_MAX;

	for (int i : max_positions) {
		highest_amplitude = std::max(highest_amplitude, this->out_real[i]);
		if (this->out_real[i] > MPM_SMALL_CUTOFF) {
			auto x = util::parabolic_interpolation(this->out_real, i);
			estimates.push_back(x);
			highest_amplitude = std::max(highest_amplitude, std::get<1>(x));
		}
	}

	if (estimates.empty())
		return -1;

	T actual_cutoff = MPM_CUTOFF * highest_amplitude;
	T period = 0;

	for (auto i : estimates) {
		if (std::get<1>(i) >= actual_cutoff) {
			period = std::get<0>(i);
			break;
		}
	}

	T pitch_estimate = (sample_rate / period);

	this->clear();

	return (pitch_estimate > MPM_LOWER_PITCH_CUTOFF) ? pitch_estimate : -1;
}

template <typename T>
T
pitch::mpm(const std::vector<T> &audio_buffer, int sample_rate)
{
	pitch_alloc::Mpm<T> ma(audio_buffer.size());
	return ma.pitch(audio_buffer, sample_rate);
}

template <typename T>
T
pitch::pmpm(const std::vector<T> &audio_buffer, int sample_rate)
{
	pitch_alloc::Mpm<T> ma(audio_buffer.size());
	return ma.probabilistic_pitch(audio_buffer, sample_rate);
}

template class pitch_alloc::Mpm<double>;
template class pitch_alloc::Mpm<float>;

template double
pitch::mpm<double>(const std::vector<double> &audio_buffer, int sample_rate);

template float
pitch::mpm<float>(const std::vector<float> &audio_buffer, int sample_rate);

template double
pitch::pmpm<double>(const std::vector<double> &audio_buffer, int sample_rate);

template float
pitch::pmpm<float>(const std::vector<float> &audio_buffer, int sample_rate);    

template <typename T>
std::pair<T, T>
util::parabolic_interpolation(const std::vector<T> &array, int x_)
{
	int x_adjusted;
	T x = (T)x_;

	if (x < 1) {
		x_adjusted = (array[x] <= array[x + 1]) ? x : x + 1;
	} else if (x > signed(array.size()) - 1) {
		x_adjusted = (array[x] <= array[x - 1]) ? x : x - 1;
	} else {
		T den = array[x + 1] + array[x - 1] - 2 * array[x];
		T delta = array[x - 1] - array[x + 1];
		return (!den) ? std::make_pair(x, array[x])
		              : std::make_pair(x + delta / (2 * den),
		                    array[x] - delta * delta / (8 * den));
	}
	return std::make_pair(x_adjusted, array[x_adjusted]);
}

template std::pair<double, double>
util::parabolic_interpolation<double>(const std::vector<double> &array, int x);
template std::pair<float, float>
util::parabolic_interpolation<float>(const std::vector<float> &array, int x);


#define SWIPE_DERBS 0.1
#define SWIPE_POLYV 0.0013028
#define SWIPE_DLOG2P 0.0104167
#define SWIPE_ST 0.3
#define SWIPE_MIN 10.0
#define SWIPE_MAX 8000.0
#define SWIPE_YP1 2.0
#define SWIPE_YPN 2.0

extern "C" {
extern int
dgels_(char *trans, int *m, int *n, int *nrhs, double *a, int *lda, double *b,
    int *ldb, double *work, int *lwork, int *info);
}

template <typename T>
static int
bilookv(std::vector<T> &yr_vector, T key, size_t lo)
{
	int md;
	size_t hi = yr_vector.size();
	lo--;
	while (hi - lo > 1) {
		md = (hi + lo) >> 1;
		if (yr_vector[md] > key)
			hi = md;
		else
			lo = md;
	}
	return (hi);
}

template <typename T>
static int
bisectv(std::vector<T> &yr_vector, T key)
{
	return bilookv(yr_vector, key, 2);
}

template <typename T> using matrix = std::vector<std::vector<T>>;

static int
sieve(std::vector<int> &ones)
{
	int k = 0;
	size_t sp = floor(sqrt(ones.size()));
	ones[0] = 0;
	for (size_t i = 1; i < sp; i++) {
		if (ones[i] == 1) {
			for (size_t j = i + i + 1; j < ones.size(); j += i + 1) {
				ones[j] = 0;
			}
			k++;
		}
	}
	for (size_t i = sp; i < ones.size(); ++i) {
		if (ones[i] == 1)
			k++;
	}
	return (k);
}

template <typename T>
static void
spline(std::vector<T> &x, std::vector<T> &y, std::vector<T> &y2)
{
	size_t i, j;
	T p, qn, sig;
	std::vector<T> u((unsigned)(x.size() - 1));
	y2[0] = -.5;
	u[0] = (3. / (x[1] - x[0])) * ((y[1] - y[0]) / (x[1] - x[0]) - SWIPE_YP1);
	for (i = 1; i < x.size() - 1; i++) {
		sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]);
		p = sig * y2[i - 1] + 2.;
		y2[i] = (sig - 1.) / p;
		u[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]) -
		       (y[i] - y[i - 1]) / (x[i] - x[i - 1]);
		u[i] = (6 * u[i] / (x[i + 1] - x[i - 1]) - sig * u[i - 1]) / p;
	}
	qn = .5;
	y2[y2.size() - 1] =
	    ((3. / (x[x.size() - 1] - x[x.size() - 2])) *
	            (SWIPE_YPN - (y[y.size() - 1] - y[y.size() - 2]) /
	                             (x[x.size() - 1] - x[x.size() - 2])) -
	        qn * u[x.size() - 2]) /
	    (qn * y2[y2.size() - 2] + 1.);
	for (j = x.size() - 2; j != (size_t)-1; --j)
		y2[j] = y2[j] * y2[j + 1] + u[j];
	return;
}

template <typename T>
static T
splinv(std::vector<T> &x, std::vector<T> &y, std::vector<T> &y2, T val, int hi)
{
	int lo = hi - 1;
	T h = x[hi] - x[lo];
	T a = (x[hi] - val) / h;
	T b = (val - x[lo]) / h;
	return (
	    a * y[lo] + b * y[hi] +
	    ((a * a * a - a) * y2[lo] * (b * b * b - b) * y2[hi]) * (h * h) / 6.);
}

template <typename T>
static void
polyfit(
    std::vector<T> &A, std::vector<T> &B, std::vector<double> &Bp, int degree)
{
	int info;
	degree++;

	// needs to be a double for LAPACK's DGELS
	std::vector<double> Ap(degree * A.size());

	size_t i, j;
	for (i = 0; i < (size_t)degree; i++)
		for (j = 0; j < A.size(); j++)
			Ap[i * A.size() + j] = pow(A[j], degree - i - 1);
	for (i = 0; i < B.size(); i++)
		Bp[i] = B[i];
	i = 1;
	j = A.size() + degree;
	int a_size = (int)A.size();
	int b_size = (int)B.size();
	int int_i = (int)i;
	int int_j = (int)j;

	// needs to be a double for LAPACK's DGELS
	std::vector<double> work(j);

	dgels_((char *)"N", &a_size, &degree, &int_i, Ap.data(), &b_size, Bp.data(),
	    &degree, work.data(), &int_j, &info);
	if (info < 0) {
		fprintf(stderr, "LAPACK routine dgels() returned error: %d\n", info);
		exit(EXIT_FAILURE);
	}
	return;
}

template <typename T>
static T
polyval(std::vector<double> &coefs, T val)
{
	T sum = 0.;
	for (size_t i = 0; i < coefs.size(); i++)
		sum += coefs[i] * pow(val, coefs.size() - i - 1);
	return (sum);
}

template <typename T>
static T
hz2erb(T hz)
{
	return static_cast<T>(21.4 * log10(1. + hz / 229.));
}

template <typename T>
static T
erb2hz(T erb)
{
	return static_cast<T>((pow(10, erb / 21.4) - 1.) * 229.);
}

template <typename T>
static T
fixnan(T x)
{
	return (std::isnan(x) ? 0. : x);
}

template <typename T>
static void
La(matrix<T> &L, std::vector<T> &f, std::vector<T> &fERBs,
    std::vector<std::complex<float>> &fo, int w2, int hi, int i)
{
	size_t j;
	std::vector<T> a(w2);
	for (j = 0; j < (size_t)w2; j++)
		a[j] = sqrt(std::real(fo[j]) * std::real(fo[j]) +
		            std::imag(fo[j]) * std::imag(fo[j]));
	std::vector<T> a2(f.size());
	spline(f, a, a2);
	L[i][0] = fixnan(sqrt(splinv(f, a, a2, fERBs[0], hi)));
	for (j = 1; j < L[0].size(); j++) {
		hi = bilookv(f, fERBs[j], hi);
		L[i][j] = fixnan(sqrt(splinv(f, a, a2, fERBs[j], hi)));
	}
}

template <typename T>
static matrix<T>
loudness(
    const std::vector<T> &x, std::vector<T> &fERBs, T nyquist, int w, int w2)
{
	size_t i, j;
	int hi;
	int offset = 0;
	T td = nyquist / w2;

	// need to be floats for ffts
	std::vector<std::complex<float>> fi(w);
	std::vector<std::complex<float>> fo(w);
	ffts_plan_t *plan = ffts_init_1d(w, FFTS_FORWARD);
	std::vector<T> hann(w);
	for (i = 0; i < (size_t)w; i++)
		hann[i] = .5 - (.5 * cos(2. * M_PI * ((T)i / w)));
	std::vector<T> f(w2);
	for (i = 0; i < (size_t)w2; i++)
		f[i] = i * td;
	hi = bisectv(f, fERBs[0]);
	matrix<T> L(ceil((T)x.size() / w2) + 1, std::vector<T>(fERBs.size()));
	for (j = 0; j < (size_t)w2; j++)
		fi[j] = {0., 0.};
	for (/* j = w2 */; j < (size_t)w; j++)
		fi[j] = {(float)(x[j - w2] * hann[j]), 0.};
	ffts_execute(plan, fi.data(), fo.data());
	La(L, f, fERBs, fo, w2, hi, 0);
	for (i = 1; i < L.size() - 2; i++) {
		for (j = 0; j < (size_t)w; j++)
			fi[j] = {(float)(x[j + offset] * hann[j]), 0.};
		ffts_execute(plan, fi.data(), fo.data());
		La(L, f, fERBs, fo, w2, hi, i);
		offset += w2;
	}
	for (/* i = L.size() - 2; */; i < L.size(); i++) {
		for (j = 0; j < x.size() - offset; j++)
			fi[j] = {(float)(x[j + offset] * hann[j]), 0.};
		for (/* j = x.size() - offset */; j < (size_t)w; j++)
			fi[j] = {0., 0.};
		ffts_execute(plan, fi.data(), fo.data());
		La(L, f, fERBs, fo, w2, hi, i);
		offset += w2;
	}
	for (i = 0; i < L.size(); i++) {
		td = 0.;
		for (j = 0; j < L[0].size(); j++)
			td += L[i][j] * L[i][j];
		if (td != 0.) {
			td = sqrt(td);
			for (j = 0; j < L[0].size(); j++)
				L[i][j] /= td;
		}
	}
	ffts_free(plan);
	return L;
}

template <typename T>
static void
Sadd(matrix<T> &S, matrix<T> &L, std::vector<T> &fERBs, std::vector<T> &pci,
    std::vector<T> &mu, std::vector<int> &ps, T nyquist2, int lo, int psz,
    int w2)
{
	size_t i, j, k;
	T t = 0.;
	T tp = 0.;
	T td;
	T dtp = w2 / nyquist2;

	matrix<T> Slocal(psz, std::vector<T>(L.size()));
	for (i = 0; i < Slocal.size(); i++) {
		std::vector<T> q(fERBs.size());
		for (j = 0; j < q.size(); j++)
			q[j] = fERBs[j] / pci[i];
		std::vector<T> kernel(fERBs.size());
		for (j = 0; j < ps.size(); j++) {
			if (ps[j] == 1) {
				for (k = 0; k < kernel.size(); k++) {
					td = fabs(q[k] - j - 1.);
					if (td < .25)
						kernel[k] = cos(2. * M_PI * q[k]);
					else if (td < .75)
						kernel[k] += cos(2. * M_PI * q[k]) / 2.;
				}
			}
		}
		td = 0.;
		for (j = 0; j < kernel.size(); j++) {
			kernel[j] *= sqrt(1. / fERBs[j]);
			if (kernel[j] > 0.)
				td += kernel[j] * kernel[j];
		}
		td = sqrt(td);
		for (j = 0; j < kernel.size(); j++)
			kernel[j] /= td;
		for (j = 0; j < L.size(); j++) {
			for (k = 0; k < L[0].size(); k++)
				Slocal[i][j] += kernel[k] * L[j][k];
		}
	}
	k = 0;
	for (j = 0; j < S[0].size(); j++) {
		td = t - tp;
		while (td >= 0.) {
			k++;
			tp += dtp;
			td -= dtp;
		}
		for (int i = 0; i < psz; i++) {
			S[lo + i][j] +=
			    (Slocal[i][k] +
			        (td * (Slocal[i][k] - Slocal[i][k - 1])) / dtp) *
			    mu[i];
		}
	}
}

template <typename T>
static void
Sfirst(matrix<T> &S, const std::vector<T> &x, std::vector<T> &pc,
    std::vector<T> &fERBs, std::vector<T> &d, std::vector<int> &ws,
    std::vector<int> &ps, T nyquist, T nyquist2, int n)
{
	int i;
	int w2 = ws[n] / 2;
	matrix<T> L = loudness(x, fERBs, nyquist, ws[n], w2);
	int lo = 0;
	int hi = bisectv(d, static_cast<T>(2.));
	int psz = hi - lo;
	std::vector<T> mu(psz);
	std::vector<T> pci(psz);
	for (i = 0; i < hi; i++) {
		pci[i] = pc[i];
		mu[i] = 1. - fabs(d[i] - 1.);
	}
	Sadd(S, L, fERBs, pci, mu, ps, nyquist2, lo, psz, w2);
}

template <typename T>
static void
Snth(matrix<T> &S, const std::vector<T> &x, std::vector<T> &pc,
    std::vector<T> &fERBs, std::vector<T> &d, std::vector<int> &ws,
    std::vector<int> &ps, T nyquist, T nyquist2, int n)
{
	int i;
	int w2 = ws[n] / 2;
	matrix<T> L = loudness(x, fERBs, nyquist, ws[n], w2);
	int lo = bisectv(d, static_cast<T>(n));
	int hi = bisectv(d, static_cast<T>(n + 2));
	int psz = hi - lo;
	std::vector<T> mu(psz);
	std::vector<T> pci(psz);
	int ti = 0;
	for (i = lo; i < hi; i++) {
		pci[ti] = pc[i];
		mu[ti] = 1. - fabs(d[i] - (n + 1));
		ti++;
	}
	Sadd(S, L, fERBs, pci, mu, ps, nyquist2, lo, psz, w2);
}

template <typename T>
static void
Slast(matrix<T> &S, const std::vector<T> &x, std::vector<T> &pc,
    std::vector<T> &fERBs, std::vector<T> &d, std::vector<int> &ws,
    std::vector<int> &ps, T nyquist, T nyquist2, int n)
{
	int i;
	int w2 = ws[n] / 2;
	matrix<T> L = loudness(x, fERBs, nyquist, ws[n], w2);
	int lo = bisectv(d, static_cast<T>(n));
	int hi = d.size();
	int psz = hi - lo;
	std::vector<T> mu(psz);
	std::vector<T> pci(psz);
	int ti = 0;
	for (i = lo; i < hi; i++) {
		pci[ti] = pc[i];
		mu[ti] = 1. - fabs(d[i] - (n + 1));
		ti++;
	}
	Sadd(S, L, fERBs, pci, mu, ps, nyquist2, lo, psz, w2);
}

template <typename T>
T
pitch_(matrix<T> &S, std::vector<T> &pc)
{
	size_t i, j;
	size_t maxi = (size_t)-1;
	int search = (int)std::round(
	    (std::log2(pc[2]) - std::log2(pc[0])) / SWIPE_POLYV + 1.);
	T nftc, maxv, log2pc;
	T tc2 = 1. / pc[1];

	std::vector<T> s(3);
	std::vector<T> ntc(3);
	ntc[0] = ((1. / pc[0]) / tc2 - 1.) * 2. * M_PI;
	ntc[1] = (tc2 / tc2 - 1.) * 2. * M_PI;
	ntc[2] = ((1. / pc[2]) / tc2 - 1.) * 2. * M_PI;
	std::vector<T> p;
	for (j = 0; j < S[0].size(); j++) {
		maxv = SHRT_MIN;
		for (i = 0; i < S.size(); i++) {
			if (S[i][j] > maxv) {
				maxv = S[i][j];
				maxi = i;
			}
		}
		if (maxv > SWIPE_ST) {
			if (!(maxi == 0 || maxi == S.size() - 1)) {
				tc2 = 1. / pc[maxi];
				log2pc = std::log2(pc[maxi - 1]);
				s[0] = S[maxi - 1][j];
				s[1] = S[maxi][j];
				s[2] = S[maxi + 1][j];
				// needs to be double for LAPACK's DGELS
				std::vector<double> coefs(2 >= s.size() ? 2 : s.size());
				polyfit(ntc, s, coefs, 2);
				maxv = SHRT_MIN;
				for (i = 0; i < (size_t)search; i++) {
					nftc = polyval(coefs,
					    static_cast<T>(
					        ((1. / pow(2, i * SWIPE_POLYV + log2pc)) / tc2 -
					            1) *
					        (2 * M_PI)));
					if (nftc > maxv) {
						maxv = nftc;
						maxi = i;
					}
				}
				p.push_back(pow(2, log2pc + (maxi * SWIPE_POLYV)));
			}
		}
	}

	return p.size() == 1 ? p[0] : -1.0;
}

template <typename T>
T
pitch::swipe(const std::vector<T> &x, int samplerate)
{
	size_t i;
	T td = 0.;
	T nyquist = samplerate / 2.;
	T nyquist2 = samplerate;
	T nyquist16 = samplerate * 8.;
	std::vector<int> ws(std::round(std::log2((nyquist16) / SWIPE_MIN) -
	                               std::log2((nyquist16) / SWIPE_MAX)) +
	                    1);
	for (i = 0; i < ws.size(); ++i)
		ws[i] =
		    pow(2, std::round(std::log2(nyquist16 / SWIPE_MIN))) / pow(2, i);
	std::vector<T> pc(
	    ceil((std::log2(SWIPE_MAX) - std::log2(SWIPE_MIN)) / SWIPE_DLOG2P));
	std::vector<T> d(pc.size());
	for (i = pc.size() - 1; i != (size_t)-1; i--) {
		td = std::log2(SWIPE_MIN) + (i * SWIPE_DLOG2P);
		pc[i] = pow(2, td);
		d[i] = 1. + td - std::log2(nyquist16 / ws[0]);
	}
	std::vector<T> fERBs(
	    ceil((hz2erb(nyquist) - hz2erb(pow(2, td) / 4)) / SWIPE_DERBS));
	td = hz2erb(SWIPE_MIN / 4.);
	for (i = 0; i < fERBs.size(); i++)
		fERBs[i] = erb2hz(td + (i * SWIPE_DERBS));
	std::vector<int> ps(floor(fERBs[fERBs.size() - 1] / pc[0] - .75), 1);
	sieve(ps);
	ps[0] = 1;
	matrix<T> S(pc.size(), std::vector<T>(1));
	Sfirst(S, x, pc, fERBs, d, ws, ps, nyquist, nyquist2, 0);
	for (i = 1; i < ws.size() - 1; ++i)
		Snth(S, x, pc, fERBs, d, ws, ps, nyquist, nyquist2, i);
	Slast(S, x, pc, fERBs, d, ws, ps, nyquist, nyquist2, i);
	return pitch_(S, pc);
}

template double
pitch::swipe<double>(const std::vector<double> &audio_buffer, int sample_rate);

template float
pitch::swipe<float>(const std::vector<float> &audio_buffer, int sample_rate);

#define YIN_THRESHOLD 0.20
#define PYIN_PA 0.01
#define PYIN_N_THRESHOLDS 100
#define PYIN_MIN_THRESHOLD 0.01

static const float Beta_Distribution[100] = {0.012614, 0.022715, 0.030646,
    0.036712, 0.041184, 0.044301, 0.046277, 0.047298, 0.047528, 0.047110,
    0.046171, 0.044817, 0.043144, 0.041231, 0.039147, 0.036950, 0.034690,
    0.032406, 0.030133, 0.027898, 0.025722, 0.023624, 0.021614, 0.019704,
    0.017900, 0.016205, 0.014621, 0.013148, 0.011785, 0.010530, 0.009377,
    0.008324, 0.007366, 0.006497, 0.005712, 0.005005, 0.004372, 0.003806,
    0.003302, 0.002855, 0.002460, 0.002112, 0.001806, 0.001539, 0.001307,
    0.001105, 0.000931, 0.000781, 0.000652, 0.000542, 0.000449, 0.000370,
    0.000303, 0.000247, 0.000201, 0.000162, 0.000130, 0.000104, 0.000082,
    0.000065, 0.000051, 0.000039, 0.000030, 0.000023, 0.000018, 0.000013,
    0.000010, 0.000007, 0.000005, 0.000004, 0.000003, 0.000002, 0.000001,
    0.000001, 0.000001, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
    0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
    0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
    0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000};

template <typename T>
static int
absolute_threshold(const std::vector<T> &yin_buffer)
{
	ssize_t size = yin_buffer.size();
	int tau;
	for (tau = 2; tau < size; tau++) {
		if (yin_buffer[tau] < YIN_THRESHOLD) {
			while (tau + 1 < size && yin_buffer[tau + 1] < yin_buffer[tau]) {
				tau++;
			}
			break;
		}
	}
	return (tau == size || yin_buffer[tau] >= YIN_THRESHOLD) ? -1 : tau;
}

// pairs of (f0, probability)
template <typename T>
static std::vector<std::pair<T, T>>
probabilistic_threshold(const std::vector<T> &yin_buffer, int sample_rate)
{
	ssize_t size = yin_buffer.size();
	int tau;

	std::map<int, T> t0_with_probability;
	std::vector<std::pair<T, T>> f0_with_probability;

	T threshold = PYIN_MIN_THRESHOLD;

	for (int n = 0; n < PYIN_N_THRESHOLDS; ++n) {
		threshold += n * PYIN_MIN_THRESHOLD;
		for (tau = 2; tau < size; tau++) {
			if (yin_buffer[tau] < threshold) {
				while (
				    tau + 1 < size && yin_buffer[tau + 1] < yin_buffer[tau]) {
					tau++;
				}
				break;
			}
		}
		auto a = yin_buffer[tau] < threshold ? 1 : PYIN_PA;
		t0_with_probability[tau] += a * Beta_Distribution[n];
	}

	for (auto tau_estimate : t0_with_probability) {
		auto f0 = (tau_estimate.first != 0)
		              ? sample_rate / std::get<0>(util::parabolic_interpolation(
		                                  yin_buffer, tau_estimate.first))
		              : -1.0;

		if (f0 != -1.0) {
			f0_with_probability.push_back(
			    std::make_pair(f0, tau_estimate.second));
		}
	}

	return f0_with_probability;
}

template <typename T>
static void
difference(const std::vector<T> &audio_buffer, pitch_alloc::Yin<T> *ya)
{
	util::acorr_r(audio_buffer, ya);

	for (int tau = 0; tau < ya->N / 2; tau++)
		ya->yin_buffer[tau] =
		    ya->out_real[0] + ya->out_real[1] - 2 * ya->out_real[tau];
}

template <typename T>
static void
cumulative_mean_normalized_difference(std::vector<T> &yin_buffer)
{
	double running_sum = 0.0f;

	yin_buffer[0] = 1;

	for (int tau = 1; tau < signed(yin_buffer.size()); tau++) {
		running_sum += yin_buffer[tau];
		yin_buffer[tau] *= tau / running_sum;
	}
}

template <typename T>
T
pitch_alloc::Yin<T>::pitch(const std::vector<T> &audio_buffer, int sample_rate)
{
	int tau_estimate;

	difference(audio_buffer, this);

	cumulative_mean_normalized_difference(this->yin_buffer);
	tau_estimate = absolute_threshold(this->yin_buffer);

	auto ret = (tau_estimate != -1)
	               ? sample_rate / std::get<0>(util::parabolic_interpolation(
	                                   this->yin_buffer, tau_estimate))
	               : -1;

	this->clear();
	return ret;
}

template <typename T>
T
pitch_alloc::Yin<T>::probabilistic_pitch(
    const std::vector<T> &audio_buffer, int sample_rate)
{
	difference(audio_buffer, this);

	cumulative_mean_normalized_difference(this->yin_buffer);

	auto f0_estimates = probabilistic_threshold(this->yin_buffer, sample_rate);

	this->clear();
	return util::pitch_from_hmm(this->hmm, f0_estimates);
}

template <typename T>
T
pitch::yin(const std::vector<T> &audio_buffer, int sample_rate)
{

	pitch_alloc::Yin<T> ya(audio_buffer.size());
	return ya.pitch(audio_buffer, sample_rate);
}

template <typename T>
T
pitch::pyin(const std::vector<T> &audio_buffer, int sample_rate)
{

	pitch_alloc::Yin<T> ya(audio_buffer.size());
	return ya.probabilistic_pitch(audio_buffer, sample_rate);
}

template class pitch_alloc::Yin<double>;
template class pitch_alloc::Yin<float>;

template double
pitch::yin<double>(const std::vector<double> &audio_buffer, int sample_rate);

template float
pitch::yin<float>(const std::vector<float> &audio_buffer, int sample_rate);

template double
pitch::pyin<double>(const std::vector<double> &audio_buffer, int sample_rate);

template float
pitch::pyin<float>(const std::vector<float> &audio_buffer, int sample_rate);