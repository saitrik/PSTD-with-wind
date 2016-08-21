//
// Created by strikootam on 3-8-2016.
//
#include<iostream>
#include "spatderv.h"
#include<math.h>
#include "fftw3.h"
using namespace std;
int count=1;
void spatderv::derv( double *in_buffer,  double *final_buffer, int gridx, int gridy, double **kmat, double dr,int stag){

    int i,j;
    int rank =1;
    int no[rank]= {gridx};
    fftw_plan plan, plan_inv;
    int howmany = gridy;
    int istride=1;
    int ostride=istride;
    int idist = gridx;
    int odist = (gridx );

    int n2;



    n2 = (gridx ) * gridy;



    fftw_complex *out_buffer;
    fftw_complex *temp;

    out_buffer = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * n2);
    temp = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * n2);

    plan=fftw_plan_many_dft_r2c(rank, no, howmany, in_buffer, NULL, istride, idist,
                                out_buffer, NULL, ostride, odist, FFTW_ESTIMATE);


    fftw_execute(plan);

    for (i=0;i<gridy;i++)
        for(j=0;j<gridx;j++) {
            if (j<(gridx/2))
            {

                out_buffer[j+gridx*i][1]= out_buffer[j + gridx * i][1];
                out_buffer[j+gridx*i][0]= out_buffer[j + gridx * i][0];}
            else{
                out_buffer[j+gridx*i][1]= - out_buffer[(gridx-j) + gridx * i][1];
                out_buffer[j+gridx*i][0]=   out_buffer[(gridx-j) + gridx* i][0];
            }
        }

    if (gridx==6 && stag==1){
        for (i = 0; i < gridx; i++)
            for (j = 0; j < gridy; j++) {
                out_buffer[j+gridy*i][1]=-out_buffer[j+gridy*i][1];
            }}

    for (i=0;i<gridy;i++) {
        for (j = 0; j < gridx; j++) {

            if (stag==0){

                temp[j + gridx * i][0] = -kmat[i][j] * out_buffer[j + gridx * i][1];
                temp[j + gridx * i][1] = kmat[i][j] * out_buffer[j + gridx * i][0];
            }


            else if (stag == 1) {

                temp[j + gridx * i][0] = -kmat[i][j] * ((sin(kmat[i][j] * dr / 2) * out_buffer[j + gridx * i][0]) +
                                                        (cos(kmat[i][j] * dr / 2) * out_buffer[j + gridx * i][1]));
                temp[j + gridx * i][1] = kmat[i][j] * ((cos(kmat[i][j] * dr / 2) * out_buffer[j + gridx * i][0]) -
                                                       (sin(kmat[i][j] * dr / 2) * out_buffer[j + gridx * i][1]));
            }
            else if(stag==2){

                temp[j + gridx * i][0] = ((cos(kmat[i][j] * dr / 2) * out_buffer[j + gridx * i][0]) - (sin(kmat[i][j] * dr / 2) * out_buffer[j + gridx * i][1]));
                temp[j + gridx * i][1] = ((sin(kmat[i][j] * dr / 2) * out_buffer[j + gridx * i][0]) + (cos(kmat[i][j] * dr / 2) * out_buffer[j + gridx * i][1]));


            }

        }

    }


    if (gridx==6 &&stag==1){
        for (i = 0; i < gridx; i++)
            for (j = 0; j < gridy; j++) {
                temp[j+gridy*i][1]=-temp[j+gridy*i][1];
            }}


    plan_inv=fftw_plan_many_dft_c2r(rank, no, howmany, temp, NULL, ostride, odist,
                                    final_buffer, NULL, istride, idist, FFTW_ESTIMATE);
    fftw_execute(plan_inv);
    for(i=0;i<n2;i++){
        final_buffer[i]=final_buffer[i]/gridx;

    }

    fftw_destroy_plan(plan);



    fftw_destroy_plan(plan_inv);

    fftw_free(out_buffer);
    fftw_free(temp);
    count=count+1;
    return ;
}


