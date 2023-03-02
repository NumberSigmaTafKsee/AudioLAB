#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "halfband.h"

#define BUF_SIZE 2048
#define MAX_FREQS 51

void test_sine(double *buf, int n, double freq)
{   
    int i;
    float ph;
    for(i=0;i<n;i++)    
    {
        ph = i/(double)n;
        buf[i] = sin(ph*M_PI*freq*n);   
    }
}


void test_halfband(void)
{
    int orders[6] = {12,10,8,6,4,2};
    int n_orders = 6;
    int order;
    int steep, i, j, k;
    halfband *half;
    double buf[BUF_SIZE];
        
    
    for(steep=0;steep<2;steep++)
    {
        printf("Steep: %d\n", steep);
        for(i=0;i<n_orders;i++)
        {   
            order = orders[i];
            printf("\tOrder: %d\n", order, steep);
            half = create_halfband(order, steep);
            
            for(k=1;k<MAX_FREQS;k++)
            {
                
                double sum, filt, sine, orig_sum, gain, freq;
                freq = k/((double)MAX_FREQS);
                test_sine(buf, BUF_SIZE, freq);
                
                sum = orig_sum = 0;
                for(j=0;j<BUF_SIZE;j++)                           
                {
                    filt = process_halfband(half, buf[j]);
                    // skip transient
                    if(j>BUF_SIZE/3)
                    {
                        sum += filt * filt;
                        orig_sum += buf[j] * buf[j];
                    }                                    
                }                            
                
                gain = sqrt(sum) / sqrt(orig_sum);                                
                printf("\t\tFreq: %.2f Gain: %6.1fdB\n", freq, 20*log(gain)/log(10));                
            }                        
            destroy_halfband(half);            
        }
        
    }
}

void test_halfband_decimator(void)
{
    half_cascade *cascade;    
    double buf[BUF_SIZE];
    int j,k;
    cascade = create_half_cascade(2, 8, 0);
    printf("\n\nCascade\n");
    for(k=1;k<MAX_FREQS;k++)
        {
            
            double sum, filt, sine, orig_sum, gain, freq;
            freq = k/((double)MAX_FREQS);
            test_sine(buf, BUF_SIZE, freq);            
            sum = orig_sum = 0;
            for(j=0;j<BUF_SIZE;j+=4)                           
            {
                filt = process_half_cascade(cascade, &buf[j]);
                // skip transient
                if(j>BUF_SIZE/3)
                {
                    sum += filt * filt;
                    orig_sum += buf[j] * buf[j];
                }                                    
            }                            
            
            gain = sqrt(sum) / sqrt(orig_sum);                                
            printf("\t\tFreq: %.2f Gain: %6.1fdB\n", freq, 20*log(gain)/log(10));                
        }                        
        destroy_half_cascade(cascade);            
        

}


int main(int argc, char **argv)
{
    test_halfband();
    test_halfband_decimator();
    
    
}