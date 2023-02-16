// program to compute circular convolution of two sequences
#include<stdio.h>
int x[10],h[10],y[10];

main()
{ int i,n,k,l1,l2,j;
  printf("enter the length of the input sequence x(n)");
  scanf("%d",&l1);
  printf("enter the length of the input sequence h(n)");
  scanf("%d",&l2);
  printf("enter the  input sequence x(n)");
  for(i=0;i<l1;i++)
     scanf("%d",&x[i]);
  printf("enter the  input sequence h(n)");
  for(i=0;i<l2;i++)
     scanf("%d",&h[i]);
     if(l1<l2)
     { for(i=l1;i<l2;i++)
          x[i]=0; 
          l1=l2;
     }
  else
     { for(i=l2;i<l1;i++)
       h[i]=0;
     }
  for(n=0;n<l1;n++)
     { y[n]=0;
       for(k=0;k<l1;k++)
          {if(k>n)
           j=n-k+l1;
           else
           j=n-k;
           y[n]=y[n]+x[k]*h[j];
          }
     }
  printf("output sequence y(n)\n");
  for(i=0;i<l1;i++)
     printf("%d\t",y[i]);
}



