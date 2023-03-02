/* CPSC501: Assignment 4
 * Adam Maidens - 10044293
 *
 * convolve.c
 * takes in a dry recording and an impulse response and produces the convolved signal
 * compile with: g++ concolve.cpp -o convolve
 * run at the cmd line with: ./convolve inputfile IRfile outputfile
 *
 * the time domain version of convolve is commented out in later version of this code
 * a version of four1 given out in class and available on D2L is used instead
 */

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <ctime>

#define SWAP(a,b)  tempr=(a);(a)=(b);(b)=tempr

using namespace std;

/* declare variables to store the parts of a wav file */

// RIFF chunk descriptor
char chunkID[4];
int chunkSize;
char format[2];

// fmt sub-chunk
char subChunk1ID[4];
int subChunk1Size;
int16_t audioFormat;
int16_t numChannels;
int sampleRate;
int byteRate;
int16_t blockAlign;
int16_t bitsPerSample;

// data sub-chunk
char subChunk2ID[4];
int subChunk2Size;
short* fileData;

int size;

float *readWav(char *filename, float *signal);
void writeWav(char *fileName, float *signal, int signalSize);

void four1(float data[], int nn, int isign);
void four1Scaling (float signal[], int N);

//void convolve(float x[], int N, float h[], int M, float y[], int P);
void overlapAdd(float *x, int N, float *h, int M, float *y, int P);

int main(int argc, char* argv[])
{
    clock_t begin = clock();    // used to time program
    
    if (argc != 4) {
        printf("Usage: ./convolve inputFile IRfile outputFile");
        return 0;
    }
    
    char *inputName = argv[1];
    char *impulseName = argv[2];
    char *outputName = argv[3];
    
    int N, M, P;    // input signal size, impulse signal size, output signal size
    
    float *x = NULL, *h = NULL, *y = NULL; // x, h, & y arrays
    
    x = readWav(inputName, x);
    N = size;
    printf("\n");
    
    h = readWav(impulseName, h);
    M = size;
    printf("\n");
    
    P = N + M - 1;
    y = new float[P];
    
    printf("\nConvolving.....\n");
    
    overlapAdd(x, N, h, M, y, P);
    
    printf("Convolution finished\n");
    
    writeWav(outputName, y, P);
    
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    
    printf("Finished with time %f seconds\n", elapsed_secs);
    
    
    return 0;
}

// reads in a wav file into an array
float *readWav (char *filename, float *signal)
{
	
	ifstream inputfile(filename, ios::in | ios::binary);
	
	inputfile.seekg(ios::beg);
	
	inputfile.read(chunkID, 4);
	inputfile.read((char*) &chunkSize, 4);
	inputfile.read(format, 4);
	inputfile.read(subChunk1ID, 4);
	inputfile.read((char*) &subChunk1Size, 4);
	inputfile.read((char*) &audioFormat, 2);
	inputfile.read((char*) &numChannels, 2);
	inputfile.read((char*) &sampleRate, 4);
	inputfile.read((char*) &byteRate, 4);
	inputfile.read((char*) &blockAlign, 2);
	inputfile.read((char*) &bitsPerSample, 2);
	if (subChunk1Size == 18) {
		inputfile.seekg(2, ios::cur);
	}
	inputfile.read(subChunk2ID, 4);
	inputfile.read((char*)&subChunk2Size, 4);
	size = subChunk2Size / 2;
	
	short *data = new short[size];
	for (int i = 0; i < size; i++) {
		inputfile.read((char *) &data[i], 2);
	}
	
	inputfile.close();
	printf("Input file %s \n", filename);
	
	short sample;
	signal = new float[size];
	printf("Size: %d\n", size);
	for (int i = 0; i < size; i++) {
		sample = data[i];
		signal[i] = (sample * 1.0) / (32767);      ///////
		if (signal[i] < -1.0)
			signal[i] = -1.0;
	}
	
	printf("Input file converted to 1 to -1 range\n");
	return signal;
}


// writes a signal to a wav file
void writeWav(char *filename, float *signal, int signalSize)
{
    
	ofstream outputfile(filename, ios::out | ios::binary);
	
	// File corrupted without hardcoded values
	char *ChunkID = "RIFF";
	char *format = "WAVE";
	// PCM = 18 was unnecessary
	subChunk1Size = 16;
	
	subChunk2Size = numChannels * signalSize * (bitsPerSample / 8);
	chunkSize = subChunk2Size + 36;
	
	outputfile.write(chunkID, 4);
	outputfile.write((char*) &chunkSize, 4);
	outputfile.write(format, 4);
	outputfile.write(subChunk1ID, 4);
	outputfile.write((char*) &subChunk1Size, 4);
	outputfile.write((char*) &audioFormat, 2);
	outputfile.write((char*) &numChannels, 2);
	outputfile.write((char*) &sampleRate, 4);
	outputfile.write((char*) &byteRate, 4);
	outputfile.write((char*) &blockAlign, 2);
	outputfile.write((char*) &bitsPerSample, 2);
	outputfile.write(subChunk2ID, 4);
	outputfile.write((char*)&subChunk2Size, 4);
	
	int16_t sample;
	
	// converting float to int between -2^15 to 2^15 - 1
	for(int i = 0; i < signalSize; i++)
	{
		sample = (int16_t)(signal[i] * (32767));
		outputfile.write((char*)&sample, 2);
	}
	outputfile.close();
}

//  The four1 four1 from Numerical Recipes in C,
//  p. 507 - 508.
//  Note:  changed float data types to double.
//  nn must be a power of 2, and use +1 for
//  isign for an four1, and -1 for the Inverse four1.
//  The data is complex, so the array size must be
//  nn*2. This code assumes the array starts
//  at index 1, not 0, so subtract 1 when
//  calling the routine (see main() below).

void four1(float data[], int nn, int isign)
{
    unsigned long n, mmax, m, j, istep, i;
    float wtemp, wr, wpr, wpi, wi, theta;
    float tempr, tempi;
    
    n = nn << 1;
    j = 1;
    
    for (i = 1; i < n; i += 2) {
        if (j > i) {
            SWAP(data[j], data[i]);
            SWAP(data[j+1], data[i+1]);
        }
        m = nn;
        while (m >= 2 && j > m) {
            j -= m;
            m >>= 1;
        }
        j += m;
    }
    
    mmax = 2;
    while (n > mmax) {
        istep = mmax << 1;
        theta = isign * (6.28318530717959 / mmax);
        wtemp = sin(0.5 * theta);
        wpr = -2.0 * wtemp * wtemp;
        wpi = sin(theta);
        wr = 1.0;
        wi = 0.0;
        for (m = 1; m < mmax; m += 2) {
            for (i = m; i <= n; i += istep) {
                j = i + mmax;
                tempr = wr * data[j] - wi * data[j+1];
                tempi = wr * data[j+1] + wi * data[j];
                data[j] = data[i] - tempr;
                data[j+1] = data[i+1] - tempi;
                data[i] += tempr;
                data[i+1] += tempi;
            }
            wr = (wtemp = wr) * wpr - wi * wpi + wr;
            wi = wi * wpr + wtemp * wpi + wi;
        }
        mmax = istep;
    }
}

// scales the numbers in a given array signal[] and stores numbers back in array
void four1Scaling (float signal[], int N)
{
	int k;
	int i;
	for (k = 0, i = 0; k < N; k++, i+=2) {
		signal[i] /= (float)N;
		signal[i+1] /= (float)N;
	}
}

void complexCalculation(float complexInput[],float complexIR[],float complexResult[], int size)
{
	int i = 0;
	int tempI = 0;
	for(i = 0; i < size; i++) {
		tempI = i * 2;
	    complexResult[tempI] = complexInput[tempI] * complexIR[tempI] - 				complexInput[tempI+1] * complexIR[tempI+1];
	    complexResult[tempI+1] = complexInput[tempI+1] * complexIR[tempI] + 				complexInput[tempI] * complexIR[tempI+1];
	}
}

void padZeroes(float toPad[], int size)
{
	memset(toPad, 0, size);
}

void unpadArray(float result[], float complete[], int size)
{
	int i, j;
    
    for(i = 0, j = 0; i < size; i++, j+=2)
    {
	    complete[i] = result[j];
    }
}

void padArray(float output[],float data[], int dataLen, int size)
{
	int i, k;
	for(i = 0, k = 0; i < dataLen; i++, k+=2)
	{
	    output[k] = data[i];
	    output[k + 1] = 0;
	}
	i = k;
    
	memset(output + k, 0, size -1);
}

void scaleSignal(float signal[], int samples)
{
	float min = 0, max = 0;
	int i = 0;
    
	for(i = 0; i < samples; i++)
	{
		if(signal[i] > max)
			max = signal[i];
		if(signal[i] < min)
			min = signal[i];
	}
    
	min = min * -1;
	if(min > max)
		max = min;
    
	for(i = 0; i < samples; i++)
	{
		signal[i] = signal[i] / max;
	}
}

// Uses overlap-add method of four1
void overlapAdd(float *x,int N,float * h,int M, float *y, int P)
{
	int totalSize = 0;
	int paddedTotalSize = 1;
	totalSize = N + M - 1;
    
	int i = 0;
	while (paddedTotalSize < totalSize)
	{
		paddedTotalSize <<= 1;
		i++;
	}
	printf("Padded Size: %i exp: %i\n", paddedTotalSize, i);
	printf("Input size: %i\n",N );
	printf("IR size: %i\n", M);
	printf("Sum IR&Input size: %i\n\n", totalSize);
    
	float *complexResult = new float[2*paddedTotalSize];
	float *input = new float[2*paddedTotalSize];
	float *ir = new float[2*paddedTotalSize];
    
	padArray(input,x, N,2*paddedTotalSize);
	padArray(ir,h, M, 2*paddedTotalSize);
	padZeroes(complexResult, 2*paddedTotalSize);
	four1(input-1, paddedTotalSize, 1);
	four1(ir-1, paddedTotalSize, 1);

	printf("Complex calc\n");
	complexCalculation(input, ir, complexResult, paddedTotalSize);

	printf("Inverse four1\n");
	four1(complexResult-1, paddedTotalSize, -1);
	printf("Scaling\n");
	four1Scaling(complexResult, paddedTotalSize);
	unpadArray(complexResult, y, P);
	scaleSignal(y, P);
}

/* the function convolve is from convolve.c
 * originally given out in class and available on D2l
 */

/*****************************************************************************
 *
 *    Function:     convolve
 *
 *    Description:  Convolves two signals, producing an output signal.
 *                  The convolution is done in the time domain using the
 *                  "Input Side Algorithm" (see Smith, p. 112-115).
 *
 *    Parameters:   x[] is the signal to be convolved
 *                  N is the number of samples in the vector x[]
 *                  h[] is the impulse response, which is convolved with x[]
 *                  M is the number of samples in the vector h[]
 *                  y[] is the output signal, the result of the convolution
 *                  P is the number of samples in the vector y[].  P must
 *                       equal N + M - 1
 *
 *****************************************************************************/

/*
 void convolve(float x[], int N, float h[], int M, float y[], int P)
 {
 int n, m;
 
 //  Make sure the output buffer is the right size: P = N + M - 1
 if (P != (N + M - 1)) {
 printf("Output signal vector is the wrong size\n");
 printf("It is %-d, but should be %-d\n", P, (N + M - 1));
 printf("Aborting convolution\n");
 return;
 }
 
 //  Clear the output buffer y[] to all zero values
 for (n = 0; n < P; n++)
 y[n] = 0.0;
 
 // Do the convolution
 //  Outer loop:  process each input value x[n] in turn
 for (n = 0; n < N; n++) {
 //  Inner loop:  process x[n] with each sample of h[]
 for (m = 0; m < M; m++)
 y[n+m] += x[n] * h[m];
 }
 }
 */