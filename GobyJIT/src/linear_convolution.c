#include<stdio.h>
int main()
{   int x[100],h[100],y[100],i,j,N,M;
    printf("Enter the size of sequence x[n]: ");
    scanf("%d",&N);
    printf("Enter the size of sequence h[n]: ");
    scanf("%d",&M);
    //Input the sequence x[n]
    printf("Enter the sequence x[n]: ");
    for(i=0;i<N;i++)
        scanf("%d",&x[i]);
    for(i=N;i<M+N-1;i++)
        x[i]=0;
    //Input the sequence h[n]
    printf("Enter the sequence h[n]: ");
    for(i=0;i<M;i++)
        scanf("%d",&h[i]);
    for(i=M;i<M+N-1;i++)
        h[i]=0;
    //Compute Convolution
    for(i=0;i<M+N-1;i++)
     {  y[i]=0;
        for(j=0;j<=i;j++)
            y[i]+=x[j]*h[i-j];
     }
    //Display the output sequence
    for(i=0;i<M+N-1;i++)
        printf("%d ",y[i]);

    return 0;
}
