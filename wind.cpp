//
// Created by strikootam on 2-8-2016.
//
#include<math.h>
#include<cmath>
#include "wind.h"
#include<iostream>
using namespace std;

void wind::UeffL(double *y, int gridx, int gridy, double *u, double *du, double **umat, double **dumat) {
    double zo=1e-5, z[gridy] ,k=0.4,ustar=0.3;
    double Bo= 0.015,L, x[gridy],siw[gridy];
    int i,j;

    L=-pow(ustar,3)/(k*Bo);



    for (i=0;i<gridy;++i)
    {
        z[i]= y[i];
        if (z[i]<=0)
        {z[i]=zo;

        }

        x[i]=pow(1 - (16 * z[i] / L),0.25);
/*        siw[i]=2*log((1+x[i])/2)+log((1+pow(x[i],2))/2)-2*atan(x[i])+atan(1)*4/2;*/
        u[i]=abs((ustar/k)*(log(z[i]/zo)));
        du[i]=z[i]*ustar/k;

    }

    for (i=0;i<gridy;++i)
    {for(j=0;j<gridx;j++)
        {
            umat[i][j]=u[i];
            dumat[i][j]=du[i];

        }
    }


    return;
}

