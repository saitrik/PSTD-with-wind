//
// Created by strikootam on 17-8-2016.
//

#include "PML.h"

void PML::PMLlayer(int gridx, int gridy, int PMLcells, double dt, double **PMLpxmat, double **PMLpymat, double **PMLvxmat, double **PMLvymat)
{double amax_PML, sigmaPMLp[PMLcells], sigmaPMLv[PMLcells], PMLp[PMLcells], PMLv[PMLcells];
    double PMLpx[gridy][PMLcells], PMLpy[gridx][PMLcells],PMLvx[gridy][PMLcells], PMLvy[gridx][PMLcells];


    amax_PML=3e4;
    for(int i=0; i<PMLcells;i++)
    {
        sigmaPMLp[i] = amax_PML *pow(double (i+1)/PMLcells,3);
        sigmaPMLv[i] = amax_PML *pow((i+0.5)/PMLcells,3);
        PMLp[i]=exp(-sigmaPMLp[i]*dt);
        PMLv[i]=exp(-sigmaPMLv[i]*dt);
    }

    for(int i=0; i<gridy;i++)
        for (int j=0; j<PMLcells; j++)
        {
            PMLpx[i][j] = PMLp[j];
            PMLvx[i][j]=PMLv[j];

        }

    for(int i=0; i<gridx;i++)
        for (int j=0; j<PMLcells; j++)
        {
            PMLpy[i][j] = PMLp[j];
            PMLvy[i][j]=PMLv[j];
        }

    for(int i=0; i<gridy;i++)
        for (int j=0; j<gridx; j++)
        {   if (j<PMLcells)
                PMLpxmat[i][j]= PMLpx[i][PMLcells-j-1];
            if (j>=PMLcells && j<gridx-PMLcells)
                PMLpxmat[i][j]=1;
            if (j>=gridx-PMLcells)
                PMLpxmat[i][j]=PMLpx[i][-(gridx-j-PMLcells)];

            if (j<PMLcells)
                PMLvxmat[i][j]= PMLvx[i][PMLcells-j-1];
            if (j>=PMLcells && j<gridx-PMLcells)
                PMLvxmat[i][j]=1;
            if (j>=gridx-PMLcells)
                PMLvxmat[i][j]=PMLvx[i][-(gridx-j-PMLcells)];

            if (i<PMLcells)
                PMLpymat[i][j]= PMLpy[j][PMLcells-i-1];
            if (i>=PMLcells && i<gridy-PMLcells)
                PMLpymat[i][j]=1;
            if (i>=gridy-PMLcells)
                PMLpymat[i][j]=PMLpy[j][-(gridy-i-PMLcells)];

            if (i<PMLcells)
                PMLvymat[i][j]= PMLvy[j][PMLcells-i-1];
            if (i>=PMLcells && i<gridy-PMLcells)
                PMLvymat[i][j]=1;
            if (i>=gridy-PMLcells)
                PMLvymat[i][j]=PMLvy[j][-(gridy-i-PMLcells)];
        }




    return;

}
