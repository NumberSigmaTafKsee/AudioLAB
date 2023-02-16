#include <iostream>
#include <cmath>
#include <iomanip>

/* Chebyshev Filter - Recursion Coefficient Calculator
 *
 * Taken from the book "The Scientist and Engineers Guide to DSP"
 * Chapter 20 - Pages 340 & 341
 * Converted from BASIC to C++
 *
 */

enum class FilterType {LowPass, HighPass};

int main()
{
    //set precision for printout
    std::cout << std::setprecision(15);
    std::cout << std::fixed; //display in fixed decimal, as opposed to exponential notation

    //hold the coefficients we want calculated
    double a[23];
    double b[23];

    //holds internal state used for combining stages
    double ta[23];
    double tb[23];

    for (auto i {0}; i < 23; ++i)
    {
        a[i] = 0;
        b[i] = 0;
    }
    
    a[2] = 1;
    b[2] = 1;

    const double pi {4 * atan(1)};

    double frequencyCutoff;
    FilterType filterType {};
    double percentRipple {};
    int numOfPoles {};
    
    std::cout << "Enter cutoff frequency (0.0 to 0.5): ";
    std::cin >> frequencyCutoff;
    //error checking
    if (frequencyCutoff < 0.0 || frequencyCutoff > 0.5)
    {
        std::cout << "Out of Range" << std::endl;
        return 1; 
    }
    int filterChoice;
    std::cout << "Enter '0' for LP, '1' for HP filter type: ";
    std::cin >> filterChoice;
    if (!(filterChoice == 0 || filterChoice == 1))
    {
        std::cout << "Invalid filter choice" << std::endl;
        return 1;
    }
    filterType = static_cast<FilterType>(filterChoice); 
    std::cout << "Enter percent percentRipple (0.0 to 29.0): ";
    std::cin >> percentRipple;
    if (percentRipple < 0 || percentRipple > 29)
    {
        std::cout << "Out of Range" << std::endl;
        return 1; 

    }
    std::cout << "Enter number of poles (even number 2, 4,...20): ";
    std::cin >> numOfPoles;
    if (numOfPoles < 2 || numOfPoles > 20 || (numOfPoles % 2 != 0))
    {
        std::cout << "Invalid number of poles" << std::endl;
    }


    double a0, a1, a2, b1, b2;

    for (auto p {1}; p <= numOfPoles/2; ++p)
    {
        //std::cout << "p: " << p << std::endl;

        //---call subroutine---//

        //Calculate the pole location on the unit circle
        double rp = -cos(pi/(numOfPoles * 2) + (p - 1) * pi/numOfPoles);
        double ip = -sin(pi/(numOfPoles * 2) + (p - 1) * pi/numOfPoles);

        double es;
        double vx;
        double kx;

        //Warp from a circle to an ellipse
        if (percentRipple > 0)
        {
            es = sqrt(pow((100.0/(100 - percentRipple)),2) - 1);
            vx = (1.0/numOfPoles) * log((1.0/es) + sqrt(1.0/(es * es) + 1));
            kx = (1.0/numOfPoles) * log((1.0/es) + sqrt(1.0/(es * es) - 1));
            kx = (exp(kx) + exp(-kx)) / 2.0;
            rp = rp * ((exp(vx) - exp(-vx)) / 2.0) / kx;
            ip = ip * ((exp(vx) + exp(-vx)) / 2.0) / kx;
        }

        /*
        std::cout << "rp: " << rp << std::endl;
        std::cout << "ip: " << ip << std::endl;
        std::cout << "es: " << es << std::endl;
        std::cout << "vx: " << vx << std::endl;
        std::cout << "kx: " << kx << std::endl;
        */

        //S-domain to Z-domain conversion
        double t = 2.0 * tan(1.0/2.0);
        double w = 2.0 * pi * frequencyCutoff;
        double m = rp * rp + ip * ip;
        double d = 4.0 - 4.0 * rp * t + m * t * t;
        double x0 = (t * t) / d;
        double x1 = (2.0 * t * t) / d;
        double x2 = (t * t) / d;
        double y1 = (8.0 - 2.0 * m * t * t) / d;
        double y2 = (-4.0 - 4.0 * rp * t - m * t * t) / d;  

        /*
        std::cout << "t: " << t << std::endl;
        std::cout << "w: " << w << std::endl;
        std::cout << "m: " << m << std::endl;
        std::cout << "d: " << d << std::endl;
        std::cout << "x0: " << x0 << std::endl;
        std::cout << "x1: " << x1 << std::endl;
        std::cout << "x2: " << x2 << std::endl;
        std::cout << "y1: " << y1 << std::endl;
        std::cout << "y2: " << y2 << std::endl;
        */

        double k;
        //Low pass to Low pass OR Low Pass to High Pass Transform
        if (static_cast<int>(filterType) == 1) // High Pass
        {
            k = -cos(w/2.0 + 1.0/2.0) / cos(w/2.0 - 1.0/2.0);
        }
        else //Low Pass
        {
            k = sin(1.0/2.0 - w/2.0) / sin(1.0/2.0 + w/2.0);
        }

        d = 1 + y1 * k - y2 * k * k;
        a0 = (x0 - x1 * k + x2 * k * k) / d;
        a1 = (-2 * x0 * k + x1 + x1 * k * k - 2 * x2 * k) / d;
        a2 = (x0 * k * k - x1 * k + x2) / d;
        b1 = (2 * k + y1 + y1 * k * k - 2 * y2 * k) / d;
        b2 = (-(k * k) - y1 * k + y2) / d;

        if (static_cast<int>(filterType) == 1) // High Pass
        {
            a1 = -a1;
            b1 = -b1;
        }

        /*
        std::cout << "A0: " << a0 << std::endl;
        std::cout << "A1: " << a1 << std::endl;
        std::cout << "A2: " << a2 << std::endl;
        std::cout << "B1: " << b1 << std::endl;
        std::cout << "B2: " << b2 << std::endl;
        */
 

        //---end subroutine---// 

        //add coefficients to the cascade
        for (auto i {0}; i < 23; ++i)
        {
            ta[i] = a[i];
            tb[i] = b[i];
        }

        for (auto i {2}; i < 23; ++i)
        {
            a[i] = a0 * ta[i] + a1 * ta[i - 1] + a2 * ta[i - 2];
            b[i] = tb[i] - b1 * tb[i - 1] - b2 * tb[i - 2]; 
        }


    } 

    //Finish combining coefficients
    b[2] = 0;
    for (auto i {0}; i < 21; ++i)
    {
        a[i] = a[i + 2];
        b[i] = -b[i + 2];
    }

    //Normalize the gain
    double sa {};
    double sb {};

    for (auto i {0}; i < 21; ++i)
    {
        if (static_cast<int>(filterType) == 1) // High Pass
        {
            sa = sa + a[i] * pow(-1, i);
            sb = sb + b[i] * pow(-1, i);
        }
        else
        {
            sa = sa + a[i];
            sb = sb + b[i];
        }
    }

    double gain = sa / (1 - sb);

    for(auto i {0}; i < 21; ++i)
    {
        a[i] = a[i] / gain;
    }
    
    //print out all the coefficients
    //
    std::cout << "a coefficients: " << std::endl;
    for (auto i {0}; i < 23; ++i)
    {
        std::cout << "a[" << i << "]: " << a[i] << std::endl;
    }

    std::cout << std::endl;

    std::cout << "b coefficients: " << std::endl;
    for (auto i {1}; i < 23; ++i)
    {
        std::cout << "b[" << i << "]: " << b[i] << std::endl;
    }


    return 0;
}