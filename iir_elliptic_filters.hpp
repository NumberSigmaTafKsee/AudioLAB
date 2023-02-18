#pragma once

namespace iir_filters
{
	/////////////////////////////////////////////////////////////////////////////////////////////
	// Elliptic Filter
	/////////////////////////////////////////////////////////////////////////////////////////////    
		BiquadSOS ellipticlp( double Ap,  //maximum pass band loss in dB,
							double As,  //minium stop band loss in dB
							double wp,  //pass band cutoff frequency
							double ws   //stop band cut off frequency
							)
		{
			Filters::CauerFilter filter(Ap,As,wp,ws);
			filter.calculate_coefficients();
			//std::cout << filter.getHo() << "\n";
			//for(size_t i = 0; i < filter.Ai.size(); i++) std::cout << filter.Ai[i] << "\n";
			//for(size_t i = 0; i < filter.Bi.size(); i++) std::cout << filter.Bi[i] << "\n";
			//for(size_t i = 0; i < filter.Ci.size(); i++) std::cout << filter.Ci[i] << "\n";

			double hZero = filter.getHo();
			BiquadSOS sos;
			for(size_t i = 0; i < filter.r; i++)
			{
				BiquadSection c;
				c.z[0] = hZero*filter.Ai[i];
				c.z[1] = 0;
				c.z[2] = hZero;    
				c.p[0] = filter.Ci[i];
				c.p[1] = filter.Bi[i];
				c.p[2] = 1;       
				c.z[0] /= c.p[0];    
				c.z[2] /= c.p[0];    
				c.p[1] /= c.p[0];
				c.p[2] /= c.p[0];
				c.p[0] /= c.p[0];
				sos.push_back(c);
			}
			return sos;
		}       
	}
}
