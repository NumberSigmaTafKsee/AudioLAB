// aka Elliptical Filter
#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <fstream>

using namespace std;

struct CauerFilter
{
	DspFloatType Ap; //maximum pass band loss in dB,
	DspFloatType As; //minium stop band loss in dB
	DspFloatType wp; //pass band cutoff frequency
	DspFloatType ws; //stop band cut off frequency


	DspFloatType k, //selectivity factor
	u,
	q, //modular constant
	D, //Discrimination factor
	V,
	po,
	W,
	miu;

	int r;

	DspFloatType n;//order of the normalized equation


	/**Constructor*/
	CauerFilter(DspFloatType Ap, //maximum pass band loss in dB,
					DspFloatType As,//minium stop band loss in dB
					DspFloatType wp,//pass band cutoff frequency
					DspFloatType ws //stop band cut off frequency
					)
	{
		this->Ap = Ap;
		this->As = As;
		this->ws = ws;
		this->wp = wp;

		calculate_constants_of_elliptic_filter();

	}



	CauerFilter(DspFloatType Ap, //maximum pass band loss in dB,
					int n,//order
					DspFloatType wp,//pass band cutoff frequency
					DspFloatType ws, //stop band cut off frequency
					int justThere//just to differentiate the two constructors.
					)

	{
		this->n = n;
		this->Ap =Ap;
		this->ws = ws;
		this->wp = wp;

		k = wp/ws;

		DspFloatType a = (1-(k*k));
		a = pow(a,0.25);

		u = 0.5*(1-a)/(1+a);

		q = u + 2*pow(u,5) + 15*pow(u,9) + 150*pow(u,13);

		int b= (int) n;
		r = (b%2==0?b/2:(b-1)/2);

		//calculate the value of As
		As = pow(10,Ap/10)-1;
		As/=(16*pow(q,n));
		As+=1;
		As = 10*log10(As);
	}


	/**Calculate for normalized transfer functions*/
	void calculate_coefficients();

	DspFloatType getHo();

	void printNumerator(ofstream*);
	void printDenominator(ofstream*);
	
	DspFloatType sum_function_1(int m);
	DspFloatType sum_function_2(int m);
	DspFloatType sum_function_3(int m);
	DspFloatType sum_function_4(int m);
	DspFloatType sum_function(int type,int m);
	/**Calculated the order of the normalized equation*/
	void calculate_constants_of_elliptic_filter();
	DspFloatType summation(int,int);

	vector <DspFloatType> Ai;
	vector <DspFloatType> Bi;
	vector <DspFloatType> Ci;
	
};

void CauerFilter :: calculate_constants_of_elliptic_filter()
{
    k = wp/ws;

    DspFloatType a = (1-(k*k));
    a = pow(a,0.25);

    u = 0.5*(1-a)/(1+a);

    q = u + 2*pow(u,5) + 15*pow(u,9) + 150*pow(u,13);

    D = ((pow(10,As/10))-1)/((pow(10,Ap/10))-1);

    n = ceil(log10(16*D)/log10(1/q));

	int b= (int) n;
	r = (b%2==0?b/2:(b-1)/2);

	As = pow(10,Ap/10)-1;
	As/=(16*pow(q,n));
	As+=1;
	As = 10*log10(As);
}


void CauerFilter::calculate_coefficients()
{
    DspFloatType b = ((pow(10,Ap/20))+1)/((pow(10,Ap/20))-1);

    V = 1/(2*n)*log(b);

	DspFloatType ans = (pow(q,0.25)*summation(1,0))/(0.5+summation(2,1)); 

	po =  abs(ans);


	DspFloatType a = (1 + (po*po/k) ),
		c =(1 + k*po*po);
	W = sqrt(a*c);

	int s= (int) n;

	miu = (s%2==0?0.5:1);
	
	for(DspFloatType i = miu;i<=r;i++)
	{
		miu = i;
		//do everythiung for Xi and all here
		DspFloatType Xi = (2*pow(q,0.25)*summation(3,0))/(1+ 2*summation(4,1));
		DspFloatType Yi = sqrt((1-(Xi*Xi/k))*(1-(k*Xi*Xi)));
		//cout<<i<<"  "<<miu<<"   "<<Xi<<"   "<<Yi<<"    "<<Xi+Yi<<endl;
		
		DspFloatType ai = 1/(Xi*Xi);

		DspFloatType bi = (2*po*Yi)/(1+po*po*Xi*Xi);

		DspFloatType ci = (pow(po*Yi,2) + pow(Xi*W,2))/pow(1+po*po*Xi*Xi,2);

		Ai.push_back(ai);
		Bi.push_back(bi);
		Ci.push_back(ci);

		//cout<<Xi<<"  "<<Yi<<"   ai = "<<ai<<" bi = "<<bi<<" ci = "<<ci<<endl;

	}
	
}

DspFloatType CauerFilter:: getHo()
{
	DspFloatType ans = 1.0;
	for(int i=0;i<r;i++)
	{
		ans*=Ci[i]/Ai[i];
	}
	if((int)n%2==0)
	{
		ans*= pow(10.0,(DspFloatType)-Ap/20);
	}

	else
	{
		ans*=po;
	}
	return ans;
}


void CauerFilter::printNumerator(ofstream* file)
	/**Write numerator to the given file*/
{

	for(int i=0;i<r;i++)
	{
		*file <<"(s^2 + "
			<<Ai[i]
			<<")";
	}
}

void CauerFilter::printDenominator(ofstream* file)
	/**Write denominator to the given file*/
{
	if((int)n%2!=0)
	{
		*file<<"(s + "
			<<po
			<<")";
	}
	for(int i=0;i<r;i++)
	{
		*file<<"(s^2 + "
			<<Bi[i]
			<<"s + "
			<<Ci[i]
			<<")";
	}
}



DspFloatType CauerFilter :: sum_function_1(int m)
	//sum function for the numerator in po
{
	DspFloatType ans = pow(-1.0,m)
		* pow(q,m*(m+1))
		*sinh((1+ 2.0*m)*V);
	return ans;
}

DspFloatType CauerFilter :: sum_function_2(int m)
	//sum function for denominator in po
{
	DspFloatType ans = pow(-1.0,m)
		* pow(q,m*m)
		*cosh(2.0*m*V);
	return ans;
}

DspFloatType CauerFilter :: sum_function_4(int m)
//sum function for denominator in Xi
{
	DspFloatType pi = acos(0.0)*2;
	DspFloatType ans = pow(-1.0,m)
		* pow(q,m*m)
		*cos(2*m*miu*pi/n);
	return ans;
}

DspFloatType CauerFilter :: sum_function_3(int m)
	//sum function for numerator in Xi
{
	DspFloatType pi = acos(0.0)*2;
	DspFloatType ans = pow(-1.0,m)
		* pow(q,m*(m+1))
		*sin((1+ 2*m)*miu*pi/n);
	return ans;
}


DspFloatType CauerFilter :: sum_function(int type,int m)
	//uses swurch statement to select required sum function
{
	DspFloatType answer;
	switch(type)
	{
		case 1:
			answer = sum_function_1(m);
			break;
		
		case 2:
			answer = sum_function_2(m);
			break;
		
		case 3:
			answer = sum_function_3(m);
			break;

		case 4:
			answer = sum_function_4(m);
			break;
	}

	return answer;
}


DspFloatType CauerFilter:: summation(int type,int start)
	//does summation from start to infinity for the chosen sum type
{
	DspFloatType previous = sum_function(type,start);
	
	
	int m=start+1;
	DspFloatType sum = previous;
	while(m<1000)
	{
		DspFloatType current = sum_function(type,m);

		DspFloatType error = abs(current/previous);
		
		if(error<=0.00001)
		{
			break;
		}
		else
		{
			sum+=current;
			previous=current;
			m++;
		}
	}

	return sum;
}