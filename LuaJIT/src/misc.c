#include <math.h>

double const tau = 6.28318530717958647692528676655900576839433879875;
double sinTau(double n) {return sin(tau*n);}
double fractionalPart(double n) {return n - (long)n;}
double lerp(double l, double r, double n) {return l + (r-l)*n;}
double clamp(double l, double h, double n) {
	return n < h ? (n > l ? n : l) : h;
}