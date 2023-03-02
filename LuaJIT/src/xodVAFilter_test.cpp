// *===========================================================================* //
//
//  __::((xodVAFilter_test.cpp))::__
//
//  ___::((XODMK Programming Industries))::___
//  ___::((XODMK:CGBW:BarutanBreaks:djoto:2020))::___
//
//
//	Purpose: C++ implementation of Virtual Analog Filters
//			 IIR Toplogy-Preserving Transform (TPT) Filters
//			 "The Art Of VA Filter Design" - Vadim Zavalishin
//
//	Revision History: Feb 08, 2017 - initial
//	Revision History: Mar 10, 2020 - current
//
// *===========================================================================* //
//
//	compiling (GCC): 
//	g++ -Wall -o xodVAFilter xodVAFilter_test.cpp xodVAFilter_base.cpp xodVAFilter.cpp
//
//
//
//
// *===========================================================================* //


#include <iostream>
#include <cmath>
#include <vector>
#include <cstdint>

#include "xodVAFilter_base.h"
#include "xodVAFilter.h"

using namespace std;



// *---------------------------------------------------------------------------* //
// *--- user settings ---* //

// SOURCE_TYPE: 1 = impulse, 2 = step, 3 = rand
//const int SOURCE_TYPE = 3;

//double noise = .001;    // amount of noise [0 <-> 1] - needs more characterization


struct fParamTB_t {
	float K;
	float alpha0;
	float alpha;
	float beta1;
	float beta2;
	float beta3;
	float beta4;
};


struct UserParam {
	string    dataPath;			// path to data directory ; (default ./data)
	string    type;				// filter type: 'LP', 'HP', 'LPHP', 'AP' ;  (default LP)
	uint32_t  numSamples;		// signal test length: n - number of samples of test ; (default 200)
	uint16_t  sampleRate;		// filter sample rate: n ; (default 48000)
	float  cutoff;				// filter cutoff frequency: n ; (default 777)
	float  resonance;		    // filter resonance: n ; (default 1.0)
	uint16_t  srcType;			// SOURCE_TYPE: 1 = impulse, 2 = step, 3 = rand ; (default rand)
};


// *---------------------------------------------------------------------------* //
///// Calculate prewarp & 'big G' (Zavalishin p46) /////////////////////

float onePoleTPT_G(float sampleRate, float cutoff) {

	float G = 0;

	// prewarp the cutoff - bilinear-transform filters
	float wd = 2*pi*cutoff;
	float T = 1/sampleRate;
	float wa = (2/T)*tan(wd*T/2);
	float g = wa*T/2;

	G = g/(1.0 + g);

	return G;
}


// *---------------------------------------------------------------------------* //
///// Display Help & User Parameters /////////////////////


void printParam(UserParam& param) {
    cout << endl<<"Current Parameter Settings:"                       							<< endl
         << "  Source data dir:      " << (param.dataPath.empty() ? "not set" : param.dataPath) << endl
         << "  Filter Type:          " << param.type               								<< endl
         << "  Number of Samples:    " << param.numSamples										<< endl
         << "  Sample Rate:          " << param.sampleRate						                << endl
         << "  Cutoff Freq:          " << param.cutoff   						                << endl
         << "  Resonance:            " << param.resonance   					                << endl
         << "  Test Source Type:     " << param.srcType  						                << endl
         << endl;
}


void help(UserParam& param) {
    cout << "\n__::(( xodVAFilter Base Test ))::__\n"
    	 << "\n  C++ implementation of Virtual Analog Filters\n"
    	 << "  IIR Toplogy-Preserving Transform (TPT) Filters\n"
    	 << "  based on 'The Art Of VA Filter Design' - Vadim Zavalishin\n"
    	 << "\n"
         << "  -h                   Show this help\n"
         << "  -f    <path>         Path to data directory where results are generated\n"
         << "  -t    <string>       Filter Type:\n"
         << "                        - 'LP'   : Lowpass 1-pole TPT Filter\n"
         << "                        - 'HP'   : Highpass Filter\n"
         << "                        - 'LPHP' : Lowpass + Highpass Filter (dual outputs)\n"
         << "                        - 'AP'   : Allpass Filter\n"
         << "                        - 'ML4P' : Moog Ladder 4-Pole Filter\n"
         << "  -n    <uint32_t>     Number of Samples (test length)\n"
         << "  -sr   <uint16_t>     Sample Rate\n"
         << "  -c    <float>        Cutoff Frequency\n"
         << "  -r    <float>        Resonance\n"
         << "  -s    <uint16_t>     Test Source Type\n"
         << endl;
    printParam(param);
    exit(1);
}


// *---------------------------------------------------------------------------* //
/////MAIN /////////////////////

// Main testbench function runs reference filters LP, LPHP, AP
// then runs hardware implementation of selected filter
// The top-level is separate for each filter to isolate each implementation in FPGA

int main(int argc, char *argv[])
{

	// setup default test parameters:
	UserParam param;
    param.dataPath          = "./data/";
    param.type        	 	= "LP";
    param.numSamples		= 1000;
    param.sampleRate		= 48000;
    param.cutoff    		= 777;
    param.resonance    		= 1.0;
    param.srcType 			= 3;

    vector<string> args;
    for ( int i = 1; i < argc; ++i ) {
          args.push_back(argv[i]);
    }

	for ( size_t i = 0; i < args.size(); ++i ) {
		if ( args[i] == "-h" ) {
			help(param);
		}
		if ( args[i] == "-f" && i+1 < args.size() ) {
			param.dataPath = args[++i];
			continue;
        }
        if ( args[i] == "-t" && i+1 < args.size() ) {
            param.type = args[++i];
            continue;
        }
        if ( args[i] == "-n" && i+1 < args.size() ) {
            param.numSamples = atof(args[++i].c_str());
            continue;
        }
        if ( args[i] == "-sr" && i+1 < args.size() ) {
            param.sampleRate = atof(args[++i].c_str());
            continue;
        }
        if ( args[i] == "-c" && i+1 < args.size() ) {
            param.cutoff = atof(args[++i].c_str());
            continue;
        }
        if ( args[i] == "-r" && i+1 < args.size() ) {
            param.resonance = atof(args[++i].c_str());
            continue;
        }
        if ( args[i] == "-s" && i+1 < args.size() ) {
            param.srcType = atof(args[++i].c_str());
            continue;
        }
        cout << endl << "ERROR: Unknown parameter: " << args[i] << endl;
        help(param);
    }

	// arrays for filter data
	float xn[param.numSamples];


	if (param.srcType==1) {
		// *---------------------------------------------------------------------------* //
		///// generate impulse source /////////////////////

		for (uint32_t i = 0; i < param.numSamples; i++) {

			if (i==1) {
				xn[i] = 1;
			} else {
				xn[i] = 0;
			}
		}

	} else if (param.srcType==2) {
		// *---------------------------------------------------------------------------* //
		///// generate step source /////////////////////

		for (uint32_t i = 0; i < param.numSamples; i++) {

			if (i>2) {
				xn[i] = 1;
			} else {
				xn[i] = 0;
			}
		}

	} else if (param.srcType==3) {
		// *---------------------------------------------------------------------------* //
		///// generate random source data /////////////////////

		for (uint32_t i = 0; i < param.numSamples; i++) {

			float r = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);

			xn[i] = 2*r-1;

			//cout<<"xnL["<<i<<"] = "<<xnL[i]<<"\txnR["<<i<<"] = "<<xnR[i]<<endl;

		}
	}

	if(param.type == "LP") {

		// *---------------------------------------------------------------------------* //
		cout << "__(( test Low-Pass filter ))__" << endl;

		printParam(param);

		FILE *f_In, *f_LPOut;

		float ynLP[param.numSamples];

		onePoleTPT_LP vaLPFlt1;

		//float sampleRate = 48000;
		//float cutoff = 5000;

		vaLPFlt1.initialize_LP(param.sampleRate);
		vaLPFlt1.setFc_LP(param.cutoff);

		for (uint32_t i = 0; i < param.numSamples; i++) {
			vaLPFlt1.doFilterStage_LP(xn[i], ynLP[i]);
		}

		string filterIn = "xodVAFilterLP_in.dat";
		string filterOut = "xodVAFilterLP_LPOut.dat";

		string filterIn_fpath = param.dataPath + filterIn;
		string filterOut_fpath = param.dataPath + filterOut;

		// write reference filter results
		f_In=fopen(filterIn_fpath.c_str(),"w");
		f_LPOut=fopen(filterOut_fpath.c_str(),"w");

		for (uint32_t i=0; i<param.numSamples; i++) {

			//cout<<"Ref filter out L: "<<ynL[i]<<"\tRef filter out R: "<<ynR[i]<<endl;

			fprintf(f_In,"%10.7f\n", xn[i]);
			fprintf(f_LPOut,"%10.7f\n", ynLP[i]);
		}
		fclose(f_In);
		fclose(f_LPOut);


		cout<<endl<<"***** Test complete *****"<<endl;
		return 0;
	}


	if(param.type == "HP") {

		// *---------------------------------------------------------------------------* //
		cout << "__(( test High-Pass filter ))__" << endl;

		printParam(param);

		FILE *f_In, *f_HPOut;

		float ynHP[param.numSamples];

		onePoleTPT_HP vaHPFlt1;

		//float sampleRate = 48000;
		//float cutoff = 5000;

		vaHPFlt1.initialize_HP(param.sampleRate);
		vaHPFlt1.setFc_HP(param.cutoff);

		for (uint32_t i = 0; i < param.numSamples; i++) {
			vaHPFlt1.doFilterStage_HP(xn[i], ynHP[i]);
		}

		string filterIn = "xodVAFilterHP_in.dat";
		string filterOut = "xodVAFilterHP_HPOut.dat";

		string filterIn_fpath = param.dataPath + filterIn;
		string filterOut_fpath = param.dataPath + filterOut;

		// write reference filter results
		f_In=fopen(filterIn_fpath.c_str(),"w");
		f_HPOut=fopen(filterOut_fpath.c_str(),"w");

		for (uint32_t i=0; i<param.numSamples; i++) {

			//cout<<"Ref filter out L: "<<ynL[i]<<"\tRef filter out R: "<<ynR[i]<<endl;

			fprintf(f_In,"%10.7f\n", xn[i]);
			fprintf(f_HPOut,"%10.7f\n", ynHP[i]);
		}
		fclose(f_In);
		fclose(f_HPOut);
		

		cout<<endl<<"***** Test complete *****"<<endl;
		return 0;
	}


	if(param.type == "LPHP") {

		// *---------------------------------------------------------------------------* //
		cout << "__(( test Low-Pass + High-Pass filter ))__" << endl;

		printParam(param);

		FILE *f_In, *f_LPOut, *f_HPOut;

		float ynLP[param.numSamples];
		float ynHP[param.numSamples];

		onePoleTPT_LPHP vaLPHPFlt1;

		//float sampleRate = 48000;
		//float cutoff = 3000;

		vaLPHPFlt1.initialize_LPHP(param.sampleRate);
		vaLPHPFlt1.setFc_LPHP(param.cutoff);

		for (uint32_t i = 0; i < param.numSamples; i++) {
			vaLPHPFlt1.doFilterStage_LPHP(xn[i], ynLP[i], ynHP[i]);
		}

		string filterIn = "xodVAFilterLPHP_in.dat";
		string filterLPOut = "xodVAFilterLPHP_LPOut.dat";
		string filterHPOut = "xodVAFilterLPHP_HPOut.dat";

		string filterIn_fpath = param.dataPath + filterLPOut;
		string filterLPOut_fpath = param.dataPath + filterLPOut;
		string filterHPOut_fpath = param.dataPath + filterHPOut;

		// write reference filter results
		f_In=fopen(filterIn_fpath.c_str(),"w");
		f_LPOut=fopen(filterLPOut_fpath.c_str(),"w");
		f_HPOut=fopen(filterHPOut_fpath.c_str(),"w");
		for (uint32_t i=0; i<param.numSamples; i++) {

			//cout<<"filter out L: "<<ynL[i]<<"\tfilter out R: "<<ynR[i]<<endl;

			fprintf(f_In,"%10.7f\n", xn[i]);
			fprintf(f_LPOut,"%10.7f\n", ynLP[i]);
			fprintf(f_HPOut,"%10.7f\n", ynHP[i]);
		}
		fclose(f_In);
		fclose(f_LPOut);
		fclose(f_HPOut);

		cout<<endl<<"***** Test complete *****"<<endl;
		return 0;

	}

	if(param.type == "AP") {

		// *---------------------------------------------------------------------------* //
		cout << "__(( test All-Pass filter ))__" << endl;

		printParam(param);

		FILE *f_In, *f_APOut;

		float ynAP[param.numSamples];

		onePoleTPT_AP vaAPFlt1;

		//float sampleRate = 48000;
		//float cutoff = 3000;

		vaAPFlt1.initialize_AP(param.sampleRate);
		vaAPFlt1.setFc_AP(param.cutoff);

		for (uint32_t i = 0; i < param.numSamples; i++) {
			vaAPFlt1.doFilterStage_AP(xn[i], ynAP[i]);
		}

		string filterIn = "xodVAFilterAP_in.dat";
		string filterOut = "xodVAFilterAP_Out.dat";

		string filterIn_fpath = param.dataPath + filterIn;
		string filterOut_fpath = param.dataPath + filterOut;

		// write reference filter results
		f_In=fopen(filterIn_fpath.c_str(),"w");
		f_APOut=fopen(filterOut_fpath.c_str(),"w");
		for (uint32_t i=0; i<param.numSamples; i++) {

			//cout<<"filter out L: "<<ynL[i]<<"\tfilter out R: "<<ynR[i]<<endl;

			fprintf(f_In,"%10.7f\n", xn[i]);
			fprintf(f_APOut,"%10.7f\n", ynAP[i]);
		}
		fclose(f_In);
		fclose(f_APOut);


		// *---------------------------------------------------------------------------* //
		///// check results /////////////////////

		// cout<<endl<<"***** Simulation done, Check Results *****"<<endl;

		// float error_LP, error_HP, avg_error_LP, avg_error_HP;

		// error_LP = 0;
		// avg_error_LP = 0;
		// error_HP = 0;
		// avg_error_HP = 0;

		// for (int i = 0; i < param.numSamples; i++) {
		// 	error_LP+= abs(ynLP[i].to_float() - ynLP_Ref[i]);
		// 	error_HP+= abs(ynHP[i].to_float() - ynHP_Ref[i]);
		// }


		// cout<<endl<<"total error_LP  = "<<error_LP<<endl;
		// avg_error_LP = error_LP/param.numSamples;
		// cout<<"average error  = "<<avg_error_LP<<endl;

		// cout<<endl<<"total error_HP  = "<<error_HP<<endl;
		// avg_error_HP = error_HP/param.numSamples;
		// cout<<"average error_HP  = "<<avg_error_HP<<endl;


		cout<<endl<<"***** Test complete *****"<<endl;
		return 0;

	}

		if(param.type == "ML4P") {

		// *---------------------------------------------------------------------------* //
		cout << "__(( test Stereo Moog Ladder 4-pole Filter ))__" << endl;

		printParam(param);

		FILE *f_ML4P_In, *f_ML4P_Out;

		float ynML4P[param.numSamples];

		xodMoogLadder4P MoogL4p;

		MoogL4p.initialize(param.sampleRate);
		MoogL4p.setFcAndRes(param.cutoff, param.resonance, param.sampleRate);


		for (uint32_t i = 0; i < param.numSamples; i++) {
			MoogL4p.advance(xn[i], ynML4P[i]);
		}


		// write reference filter results

		string moogL4p_in = "moogL4p_in.dat";
		string moogL4p_inDir = param.dataPath + moogL4p_in;

		string moogL4p_out = "moogL4p_ref_out.dat";
		string moogL4p_outDir = param.dataPath + moogL4p_out;

		f_ML4P_In=fopen(moogL4p_inDir.c_str(),"w");
		f_ML4P_Out=fopen(moogL4p_outDir.c_str(),"w");


		for (uint32_t i=0;i<param.numSamples;i++) {
			fprintf(f_ML4P_In,"%10.7f\n", xn[i]);
			fprintf(f_ML4P_Out,"%10.7f\n", ynML4P[i]);
		}

		fclose(f_ML4P_In);
		fclose(f_ML4P_Out);

		cout<<endl<<"***** Test complete *****"<<endl;
		return 0;

	}

}

