#include<stdio.h>
#include<math.h>
#define pi 4*atan(1)
void verify(double,double);
void comquad(double,double,double,double);
void getroots(double,double);// will store the root in Yr , Yi
void cpow(double,double,int);
void comsqrt(double,double);
void comcuberoot(double,double);
void cmul(double,double,double,double);
void cadd(double,double,double,double);
void cdiv(double,double,double,double);//will sore the result in re , im
int flag;
double a1re,a2re,a3re,a4re,a5re,a1im,a2im,a3im,a4im,a5im;
double re1,re2,im1,im2,re3,im3;// to store cuberoots and square roots
double re,im,mod,theta,Yr,Yi,kre,kim;
double a11re,a22re,a33re,yre,a11im,a22im,a33im,yim;
double e1re,e1im,e2re,e2im,e3re,e3im,e4re,e4im,e5re,e5im,e6re,e6im,vre,vim;
void main()
{
 int end = 2;
 double are,bre,cre,dre,aim,bim,cim,dim,pre,qre,rre,pim,qim,rim;
 double x1re,x1im,x2re,x2im,x3re,x3im,x4re,x4im;//roots of actual biquadratic
 double p1re,q1re,p1im,q1im;//coeff of depressed cubic
 double ire1,iim1,jre1,jim1,ire2,iim2,ire3,iim3,jre2,jim2,jre3,jim3; // in solving cubic eqn.
 do{
 flag = 0;
 printf("\nenter the real and imaginary parts of first coefficient\n");
 scanf("%lf %lf",&a1re,&a1im);
 if((a1re!=0)||(a1im!=0))
 {
 printf("\nenter the real and imaginary parts of second coefficient\n");
 scanf("%lf %lf",&a2re,&a2im);
 printf("\nenter the real and imaginary parts of third coefficient\n");
 scanf("%lf %lf",&a3re,&a3im);
 printf("\nenter the real and imaginary parts of fourth coefficient\n");
 scanf("%lf %lf",&a4re,&a4im);
 printf("\nenter the real and imaginary parts of fifth coefficient\n");
 scanf("%lf %lf",&a5re,&a5im);
 cdiv(a2re,a2im,a1re,a1im);
 are = re;
 aim = im;
 cdiv(a3re,a3im,a1re,a1im);
 bre = re;
 bim = im;
 cdiv(a4re,a4im,a1re,a1im);
 cre = re;
 cim = im;
 cdiv(a5re,a5im,a1re,a1im);
 dre = re;
 dim = im;

// p=(8*b-3*a*a)/8;
 e1re = 8*bre;
 e1im = 8*bim;
 cpow(are,aim,2);
 e2re = 3*re;
 e2im = 3*im;
 pre = (e1re - e2re)/8;
 pim = (e1im - e2im)/8;
 //printf("p = %lf +i(%lf)\n",pre,pim);
    // q=(a*a*a-4*a*b+8*c)/8;
 cpow(are,aim,3);
 e1re = re;
 e1im = im;
 cmul(are,aim,bre,bim);
 e2re = 4*re;
 e2im = 4*im;
 e3re = 8*cre;
 e3im = 8*cim;
 qre = (e1re - e2re + e3re)/8;
 qim = (e1im - e2im + e3im)/8;
// printf("q = %lf +i(%lf)\n",qre,qim);

// r=(-3*a*a*a*a+256*d-64*a*c+16*a*a*b)/256;
 cpow(are,aim,4);
 e1re = 3*re;
 e1im = 3*im;
 e2re = 256*dre;
 e2im = 256*dim;
 cmul(are,aim,cre,cim);
 e3re = 64*re;
 e3im = 64*im;
 cpow(are,aim,2);
 cmul(re,im,bre,bim);
 e4re = 16*re;
 e4im = 16*im;
 rre = (-e1re + e2re - e3re + e4re)/256;
 rim = (-e1im + e2im - e3im + e4im)/256;
// printf("r = %lf +i(%lf)\n",rre,rim);

 //a11=2.5*p;
 a11re = 2.5*pre;
 a11im = 2.5*pim;
 //printf("a11 = %lf + i(%lf)\n",a11re,a11im);
// a22=2*p*p-r;
 cpow(pre,pim,2);
 e1re = 2*re;
 e1im = 2*im;
 a22re = e1re - rre;
 a22im = e1im - rim;
// printf("a22 = %lf + i(%lf)\n",a22re,a22im);
// a33=.5*p*p*p-.5*p*r-.125*q*q;
 cpow(pre,pim,3);
 e1re = .5*re;
 e1im = .5*im;
 cmul(pre,pim,rre,rim);
 e2re = .5*re;
 e2im = .5*im;
 cpow(qre,qim,2);
 e3re = .125*re;
 e3im = .125*im;
 a33re = e1re - e2re - e3re;
 a33im = e1im - e2im - e3im;
 //printf("a33 = %lf + i(%lf)\n",a33re,a33im);
//************get the coeff of the cubic to be solved: a11re,a11im ,a22re,a22im, a33re,a33im********************
    // p1=(3*a22-a11*a11)/3;
     cmul(a11re,a11im,a11re,a11im);
     p1re = 3*a22re - re;
     p1im = 3*a22im - im;
     p1re = p1re/3;
     p1im = p1im/3;
 //    printf("p1 = %lf + i(%lf)\n",p1re,p1im);
//     q1=(2*a11*a11*a11-9*a11*a22+27*a33)/27;//t^3+(p1)t+(q1)
     cpow(a11re,a11im,3);
     e1re = 2*re;
     e1im = 2*im;
     cmul(a11re,a11im,a22re,a22im);
     e2re = 9*re;
     e2im = 9*im;
     e3re = 27*a33re;
     e3im = 27*a33im;
     q1re = e1re -e2re +e3re;
     q1im = e1im -e2im +e3im;
     q1re = q1re/27;
     q1im = q1im/27;
    // printf("q1 = %lf + i(%lf)\n",q1re,q1im);

     //not for y but t-a11/3 ****************
   //  e1=.25*q1*q1+(p1*p1*p1/27);
     cpow(p1re,p1im,3);
     e1re = re/27;
     e1im = im/27;
     cpow(q1re,q1im,2);
     e2re = .25*re;
     e2im = .25*im;
     e1re = e1re + e2re;
     e1im = e1im + e2im;
 //    printf("e1 = %lf + i(%lf)\n\n",e1re,e1im);
          comsqrt(e1re,e1im); //re1,im1 has the value
         e1re = re1;
         e1im = im1;
//         i = e1 - .5*q1;
//         j = -e1 - .5*q1;
         ire1 = e1re - .5*q1re;
         iim1 = e1im - .5*q1im;
         jre1 = -e1re - .5*q1re;
         jim1 = -e1im - .5*q1im;
       //  printf("i = %lf + i(%lf)\n",ire1,iim1);
       //  printf("j = %lf + i(%lf)\n",jre1,jim1);
         comcuberoot(ire1,iim1);
         ire1 = re1;
         iim1 = im1;
         ire2 = re2;
         iim2 = im2;
         ire3 = re3;
         iim3 = im3;
         comcuberoot(jre1,jim1);
         jre1 = re1;
         jim1 = im1;
         jre2 = re2;
         jim2 = im2;
         jre3 = re3;
         jim3 = im3;
        yre = ire1 + jre1 - a11re/3;
        yim = iim1 + jim1 - a11im/3;
        getroots(yre,yim);

        yre = ire2 + jre2 - a11re/3;
        yim = iim2 + jim2 - a11im/3;
        getroots(yre,yim);

        yre = ire3 + jre3 - a11re/3;
        yim = iim3 + jim3 - a11im/3;
        getroots(yre,yim);

        yre = ire1 + jre2 - a11re/3;
        yim = iim1 + jim2 - a11im/3;
        getroots(yre,yim);

        yre = ire1 + jre3 - a11re/3;
        yim = iim1 + jim3 - a11im/3;
        getroots(yre,yim);

        yre = ire2 + jre1 - a11re/3;
        yim = iim2 + jim1 - a11im/3;
        getroots(yre,yim);

        yre = ire2 + jre3 - a11re/3;
        yim = iim2 + jim3 - a11im/3;
        getroots(yre,yim);

        yre = ire3 + jre1 - a11re/3;
        yim = iim3 + jim1 - a11im/3;
        getroots(yre,yim);

        yre = ire3 + jre2 - a11re/3;
        yim = iim3 + jim2 - a11im/3;
        getroots(yre,yim);
//*************solved the cubic got the roots as : Yr ,Yi;
     // printf("\nroot: %lf + i%lf\n\n",Yr,Yi);
      kre = pre + 2*Yr;
      kim = pim + 2*Yi;
      comsqrt(kre,kim);
      e1re = re1;
      e1im = im1;
      cdiv(qre,qim,e1re,e1im);
      e2re = pre + Yr - .5*re;
      e2im = pim + Yi - .5*im;
      e3re = pre + Yr + .5*re;
      e3im = pim + Yi + .5*im;
       comquad(e1re,e1im,e2re,e2im);
       x1re = re1 - .25*are;
       x2re = re2 - .25*are;
       x1im = im1 - .25*aim;
       x2im = im2 - .25*aim;
       printf("\n\nROOTS OF THE EQUATION ARE:\n\n");
       if(x1im!=0)
       printf("root1:    %lf + i(%lf)\n\n",x1re,x1im);
       else
       printf("root1:    %lf\n\n",x1re);
       if(x2im!=0)
       printf("root2:    %lf + i(%lf)\n\n",x2re,x2im);
       else
       printf("root2:    %lf\n\n",x2re);
       comquad(-e1re,-e1im,e3re,e3im);
       x3re = re1 - .25*are;
       x4re = re2 - .25*are;
       x3im = im1 - .25*aim;
       x4im = im2 - .25*aim;
       if(x3im!=0)
       printf("root3:    %lf + i(%lf)\n\n",x3re,x3im);
       else
       printf("root3:    %lf\n\n",x3re);
       if(x4im!=0)
       printf("root4:    %lf + i(%lf)\n\n",x4re,x4im);
       else
       printf("root4:    %lf\n\n",x4re);
       verify(x1re,x1im);
       verify(x2re,x2im);
       verify(x3re,x3im);
       verify(x4re,x4im);
 }
     else
    printf("INVALID ENTERY\n\n");
     printf("\nPRESS 1 TO QUIT \n\n");
     scanf("%d",&end);
        }while(end!=1);
     }
     void comsqrt(double x,double y)
{
   // printf("\n\n IN comsqrt() FUNCTION\n");
  //  printf("x + iy = %lf + %lfi\n\n",x,y);
    double mod,re,im;
    mod=x*x+y*y;
    mod=sqrt(mod);
    re1=(mod+x)/2;
    im1=(mod-x)/2;
    re1=sqrt(re1);
    im1=sqrt(im1);
    if(y>0)
    {
        re2=-re1;
        im2=-im1;
    }
    else
    {
        re1=-re1;
        re2=-re1;
        im2=-im1;
    }
       // printf("square roots are:%lf+%lfi ,  %lf+%lfi\n",re1,im1,re2,im2);
       // printf("verify\n");
        re=re1*re1-im1*im1;
        im=2*re1*im1;
      //  printf("%lf+i%lf\n",re,im);
        re=re2*re2-im2*im2;
        im=2*re2*im2;
     //   printf("%lf+i%lf\n",re,im);

}
void cadd(double xre,double xim,double yre,double yim)
{
    re = xre + yre;
    im = xim + yim;
}
void cmul(double xre,double xim,double yre,double yim)
{
    re = xre*yre - xim*yim;
    im = xre*yim + yre*xim;
}
void comcuberoot(double x,double y) //re1,im1, re2,im2, re3,im3
{

mod = x*x + y*y;
mod = pow(mod,1/6.0);
if((x>=0)&&(y>=0))
theta = atan(y/x);
else if((x<=0)&&(y>=0))
{
   x = -x;
   theta = pi - atan(y/x);
}
else if((x<0)&&(y<0))
{
    theta = pi + atan(y/x);
}
else
{
    y = -y;
    theta = 2*pi - atan(y/x);
}
 theta = theta/3;
 re1 = mod*cos(theta);
 im1 = mod*sin(theta);

 theta+=2*pi/3;
 re2 = mod*cos(theta);
 im2 = mod*sin(theta);

theta+=2*pi/3;
re3 = mod*cos(theta);
im3 = mod*sin(theta);

}
void cpow(double x,double y,int i) // changes in re,im;
{
    int k = 1;
    re = 1;
    im = 0;
    for(;k<=i;k++)
     cmul(re,im,x,y);

}

void getroots(double x,double y) //y^3+(a11)y^2+(a22)y+(a33)
{

    float r,i;
//    printf("\nVERIFY FUNCTION\n\n");
     cpow(x,y,3);
        e1re = re;
        e1im = im;
        cpow(x,y,2);
        cmul(a11re,a11im,re,im);
        e2re = re;
        e2im = im;
        cmul(a22re,a22im,x,y);
        e3re = re;
        e3im = im;
        e4re = e1re + e2re + e3re + a33re;
        e4im = e1im + e2im + e3im + a33im;
        r = e4re;
        i = e4im;
        r = r*r;
        i = i*i;
      //  printf("%lf + i(%lf)\n\n",e4re,e4im);
     //   printf("value passed: %lf + i(%lf)\n\n",yre,yim);
        r++;
        i++;

        if((r==1)&&(i==1)&&(flag==0))
        {
            Yr = yre;
            Yi = yim;
            flag = 1;
        }

}
void cdiv(double x1,double y1,double x2,double y2) // will store the result in re,im //divides first by second
{
    double m = x2*x2 + y2*y2;
   // printf("in cdiv function\n");
    cmul(x1,y1,x2,-y2);
    re = re/m;
    im = im/m;

}
void comquad(double xre,double xim,double yre,double yim)
{
    double dre,dim,r1re,r1im,r2re,r2im;
    cpow(xre,xim,2);
    dre = re - 4*yre;
    dim = im - 4*yim;
    comsqrt(dre,dim);
    r1re = -xre + re1;
    r1im = -xim + im1;
    r2re = -xre - re1;
    r2im = -xim - im1;
    re1 = r1re/2;
    im1 = r1im/2;
    re2 = r2re/2;
    im2 = r2im/2;
}
void verify(double x,double y)
{
    cpow(x,y,4);
    cmul(re,im,a1re,a1im);
    e1re = re;
    e1im = im;
    cpow(x,y,3);
    cmul(re,im,a2re,a2im);
    e2re = re;
    e2im = im;
    cpow(x,y,2);
    cmul(re,im,a3re,a3im);
    e3re = re;
    e3im = im;
    cmul(x,y,a4re,a4im);
    e4re = re;
    e4im = im;
    vre = e1re + e2re + e3re + e4re + a5re;
    vim = e1im + e2im + e3im + e4im + a5im;
    printf("verification of %lf + i(%lf):\n",x,y);
    printf("%lf + i(%lf)\n\n\n",vre,vim);
}
