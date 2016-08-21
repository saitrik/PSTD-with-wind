
#include<iostream>
#include <fstream>
#include<cmath>
#include "wind.h"
#include "fftw3.h"
#include"spatderv.h"
#include"PML.h"
using namespace std;

int main() {
    spatderv spat;
    wind uwind;
    PML layer;
    int c_air, freq_max, PMLcells, i, j, gridx, gridy, cell_x_r, cell_y_r, cell_x_s, cell_y_s, n, RK, ntotal;
    double rho, fspatial, ppw, dx, dy, dt, tnew, pulse, ttotal, xmin, xmax, ymin, ymax, *x, *y, x_s, y_s, x_r, y_r, **PMLpxmat,**PMLpymat,**PMLvxmat,**PMLvymat;
    double *u, *du, **umat, **dumat, *kx, *ky, **kxmat, **kymat, **p, **px, **py, **v, **vx, **vy, **pxold, **pyold, **vxold, **vyold, Rko6s[6], *P_r;
    double **dx_p_stag, **dy_p_stag, **dx_vx_stag, **dy_vy_stag, **dx_p, **dy_p, **dx_vx, **dy_vy, **dx_vy, **vyold_ystag, **vyold_yxstag, **umat_ystag;

/* Constants and calculation settings */

    c_air = 343;

    rho = 1.21;
    freq_max = 1000;
    ppw = 2;


    fspatial = freq_max * ppw;

    ttotal = 0.1;
    dx = (c_air / fspatial);



    dy = dx;
    dt = 1 / (2 * fspatial);

    ntotal = (int) round(ttotal / dt);
    ntotal = 10;
    PMLcells = 20;

    /*Domain settings */

    xmin = 0;
    xmax = 20;
    ymax = 20;
    ymin = 0;
    x_s = 5;
    y_s = 5;
    x_r = 5;
    y_r = 5;



    xmin = xmin - PMLcells * dx;
    xmax = xmax + PMLcells * dx;
    ymin = ymin - PMLcells * dy;
    ymax = ymax + PMLcells * dy;


    gridx = (int) ((xmax - xmin) / dx) + 1;
    if (gridx % 2 != 0) {
        gridx = gridx - 1;
    }
    gridy = (int) ((ymax - ymin) / dy) + 1;

    if (gridy % 2 != 0) {
        gridy = gridy - 1;
    }

/* Variable initialisation */

    x = new double[gridx];
    y = new double[gridy];
    kx = new double[gridx];
    ky = new double[gridy];
    p = new double *[gridy];
    dx_p_stag = new double *[gridy];
    dy_p_stag = new double *[gridy];
    dx_vx_stag = new double *[gridy];
    dy_vy_stag = new double *[gridy];
    dx_p = new double *[gridy];
    dy_p = new double *[gridy];
    dx_vx = new double *[gridy];
    dy_vy = new double *[gridy];
    dx_vy = new double *[gridy];
    px = new double *[gridy];
    py = new double *[gridy];
    v = new double *[gridy];
    vx = new double *[gridy];
    vy = new double *[gridy];
    u = new double[gridy];
    du = new double[gridy];
    umat = new double *[gridy];
    dumat = new double *[gridy];
    kxmat = new double *[gridy];
    kymat = new double *[gridx];
    vxold = new double *[gridy];
    vyold = new double *[gridy];
    pxold = new double *[gridy];
    pyold = new double *[gridy];
    PMLpxmat = new double *[gridy];
    PMLpymat = new double *[gridy];
    PMLvxmat = new double *[gridy];
    PMLvymat = new double *[gridy];
    vyold_ystag = new double *[gridy];
    vyold_yxstag = new double *[gridy];
    umat_ystag = new double *[gridy];
    P_r = new double[ntotal];

    for (i = 0; i < gridy; i++) {

        p[i] = new double[gridx];
        px[i] = new double[gridx];
        py[i] = new double[gridx];
        v[i] = new double[gridx];
        vx[i] = new double[gridx];
        vy[i] = new double[gridx];
        umat[i] = new double[gridx];
        dumat[i] = new double[gridx];
        kxmat[i] = new double[gridx];
        kymat[i] = new double[gridx];
        dx_p_stag[i] = new double[gridx];
        dy_p_stag[i] = new double[gridx];
        dx_vx_stag[i] = new double[gridx];
        dy_vy_stag[i] = new double[gridx];
        dx_p[i] = new double[gridx];
        dy_p[i] = new double[gridx];
        dx_vx[i] = new double[gridx];
        dy_vy[i] = new double[gridx];
        dx_vy[i] = new double[gridx];
        vxold[i] = new double[gridx];
        vyold[i] = new double[gridx];
        pxold[i] = new double[gridx];
        pyold[i] = new double[gridx];
        vyold_ystag[i] = new double[gridx];
        vyold_yxstag[i] = new double[gridx];
        umat_ystag[i] = new double[gridx];
        PMLpxmat[i] = new double[gridx];
        PMLpymat[i] = new double[gridx];
        PMLvxmat[i] = new double[gridx];
        PMLvymat[i] = new double[gridx];
    }


    for (i = 0; i < gridx; i++) {
        kymat[i] = new double[gridy];
        x[i] = xmin + dx * i;
        if (i <= gridx / 2)
            kx[i] = 2 * atan(1) * 4 * i / (gridx * dx);
        else
            kx[i] = -2 * atan(1) * 4 * (gridx - i) / (gridx * dx);

    }

    for (i = 0; i < gridy; i++) {
        y[i] = ymin + dy * i;
        if (i <= gridy / 2)
            ky[i] = 2 * atan(1) * 4 * i / (gridy * dy);
        else
            ky[i] = -2 * atan(1) * 4 * (gridy - i) / (gridy * dy);

    }


    for (i = 0; i < gridy; i++) {
        for (j = 0; j < gridx; j++) {

            kxmat[i][j] = kx[j];

        }

    }


    for (i = 0; i < gridx; i++) {
        for (j = 0; j < gridy; j++) {

            kymat[i][j] = ky[j];


        }

    }

/* Source and receiver's location */
    cell_x_s = (int) (((gridx) * (x_s - x[0]) / (x[gridx - 1] - x[0])));
    cell_y_s = (int) (((gridy) * (y_s - y[0]) / (y[gridy - 1] - y[0])));


    cell_x_r = (int) (((gridx) * (x_r - x[0]) / (x[gridx - 1] - x[0])));
    cell_y_r = (int) (((gridy) * (y_r - y[0]) / (y[gridy - 1] - y[0])));



    for (i = 0; i < gridy; i++) {
        for (j = 0; j < gridx; j++) {
            p[i][j] = 0;
            px[i][j] = 0;
            py[i][j] = 0;
            v[i][j] = 0;
            vx[i][j] = 0;
            vy[i][j] = 0;
            dx_p_stag[i][j] = 0;
            dy_p_stag[i][j] = 0;
            dx_vx_stag[i][j] = 0;
            dy_vy_stag[i][j] = 0;
            dx_p[i][j] = 0;
            dy_p[i][j] = 0;
            dx_vx[i][j] = 0;
            dy_vy[i][j] = 0;
            vxold[i][j] = 0;
            vyold[i][j] = 0;
            pxold[i][j] = 0;
            pyold[i][j] = 0;
            vyold_ystag[i][j] = 0;
            vyold_yxstag[i][j] = 0;

        }
    }

/* Wind speed calculation */
    uwind.UeffL(y, gridx, gridy, u, du, umat, dumat);


/* Calcuation of PML constants*/

    layer.PMLlayer(gridx, gridy, PMLcells, dt, PMLpxmat,PMLpymat,PMLvxmat,PMLvymat);



    /* Runge kutta coefficients*/

    Rko6s[5] = 1;
    Rko6s[4] = 0.5;
    Rko6s[3] = 0.165919771368 / Rko6s[4];
    Rko6s[2] = 0.040919732041 / (Rko6s[3] * Rko6s[4]);
    Rko6s[1] = 0.007555704391 / (Rko6s[2] * Rko6s[3] * Rko6s[4]);
    Rko6s[0] = 0.000891421261 / (Rko6s[1] * Rko6s[2] * Rko6s[3] * Rko6s[4]);


     double *in_buffer;
     double *final_buffer;
    int n1 = 2 * (gridx / 2) * gridy;
    in_buffer = ( double *) fftw_malloc(sizeof( double) * (n1));
    final_buffer = ( double *) fftw_malloc(sizeof( double) * (n1));


    fstream myfile;
    myfile.open ("C:/Users/strikootam/ClionProjects/untitled/example.txt",ios::out);

    for(n=0;n<ntotal;n++) {
        for (RK = 0; RK < 6; RK++) {

            /* Spatial derivatives */

            for (i = 0; i < gridy; i++)
                for (j = 0; j < gridx; j++) {
                    in_buffer[j + gridx * i] = p[i][j];
                }
        spat.derv(in_buffer, final_buffer, gridx, gridy, kxmat, dx,1);

            for (i = 0; i < gridy; i++) {
                for (j = 0; j < gridx; j++) {
                    dx_p_stag[i][j] = final_buffer[j + gridx * i];

                }
            }

            for (j = 0; j < gridx; j++)
                for (i = 0; i < gridy; i++) {
                    in_buffer[i + gridy * j] = p[i][j];
                }
           spat.derv(in_buffer, final_buffer, gridy, gridx, kymat, dy,1);

            for (j = 0; j < gridx; j++) {
                for (i = 0; i < gridy; i++) {
                    dy_p_stag[i][j] = final_buffer[i + gridy * j];

                }
            }


            for (i = 0; i < gridy; i++)
                for (j = 0; j < gridx; j++) {
                    in_buffer[j + gridx * i] = vx[i][j];
                }

         spat.derv(in_buffer, final_buffer, gridx, gridy, kxmat, -dx,1);


            for (i = 0; i < gridy; i++) {
                for (j = 0; j < gridx; j++) {
                    dx_vx_stag[i][j] = final_buffer[j + gridx * i];
                }
            }




            for (j = 0; j < gridx; j++)
                for (i = 0; i < gridy; i++) {
                    in_buffer[i + gridy * j] = vy[i][j];
                }

         spat.derv(in_buffer, final_buffer, gridy, gridx, kymat, -dy,1);

            for (j = 0; j < gridx; j++) {
                for (i = 0; i < gridy; i++) {
                    dy_vy_stag[i][j] = final_buffer[i + gridy * j];
                }
            }



            for (i = 0; i < gridy; i++)
                for (j = 0; j < gridx; j++) {
                    in_buffer[j + gridx * i] = p[i][j];
                }

         spat.derv(in_buffer, final_buffer, gridx, gridy, kxmat,dx,0);

            for (i = 0; i < gridy; i++) {
                for (j = 0; j < gridx; j++) {
                    dx_p[i][j] = final_buffer[j + gridx * i];

                }
            }


            for (j = 0; j < gridx; j++)
                for (i = 0; i < gridy; i++) {
                    in_buffer[i + gridy * j] = p[i][j];
                }

          spat.derv(in_buffer, final_buffer, gridy, gridx, kymat,dy,0);

            for (j = 0; j < gridx; j++) {
                for (i = 0; i < gridy; i++) {
                    dy_p[i][j] = final_buffer[i + gridy * j];

                }
            }



            for (i = 0; i < gridy; i++)
                for (j = 0; j < gridx; j++) {
                    in_buffer[j + gridx * i] = vx[i][j];
                }

        spat.derv(in_buffer, final_buffer, gridx, gridy, kxmat,dx,0);

            for (i = 0; i < gridy; i++) {
                for (j = 0; j < gridx; j++) {
                    dx_vx[i][j] = final_buffer[j + gridx * i];
                }
            }


            for (j = 0; j < gridx; j++)
                for (i = 0; i < gridy; i++) {
                    in_buffer[i + gridy * j] = vy[i][j];
                }

           spat.derv(in_buffer, final_buffer, gridy, gridx, kymat,dy,0);

            for (j = 0; j < gridx; j++) {
                for (i = 0; i < gridy; i++) {
                    dy_vy[i][j] = final_buffer[i + gridy * j];

                }
            }


            for (i = 0; i < gridy; i++)
                for (j = 0; j < gridx; j++) {
                    in_buffer[j + gridx * i] = vy[i][j];
                }

          spat.derv(in_buffer, final_buffer, gridx, gridy, kxmat,dx,0);

            for (i = 0; i < gridy; i++) {
                for (j = 0; j < gridx; j++) {
                    dx_vy[i][j] = final_buffer[j + gridx * i];
                }
            }


            for (j = 0; j < gridx; j++)
                for (i = 0; i < gridy; i++) {
                    in_buffer[i + gridy * j] = vyold[i][j];
                }

            spat.derv(in_buffer, final_buffer, gridy, gridx, kymat, -dy,2);

            for (j = 0; j < gridx; j++) {
                for (i = 0; i < gridy; i++) {
                    vyold_ystag[i][j] = final_buffer[i + gridy * j];

                }
            }


            for (i = 0; i < gridy; i++)
                for (j = 0; j < gridx; j++) {
                    in_buffer[j + gridx * i] = vyold_ystag[i][j];
                }

          spat.derv(in_buffer, final_buffer, gridx, gridy, kxmat, dx,2);

            for (i = 0; i < gridy; i++) {
                for (j = 0; j < gridx; j++) {

                    vyold_yxstag[i][j] = final_buffer[j + gridx * i];

                }
            }

            for (j = 0; j < gridx; j++) {
                for (i = 0; i < gridy; i++) {

                in_buffer[i + gridy * j] = umat[i][j];

                }
            }

            spat.derv(in_buffer, final_buffer, gridy, gridx, kymat, dy,2);

            for (j = 0; j < gridx; j++) {
                for (i = 0; i < gridy; i++) {
                    umat_ystag[i][j] = final_buffer[i + gridy * j];

                }
            }

                    /* Pressure and velocity calculation */

            for (i = 0; i < gridy; i++){
                for (j = 0; j < gridx; j++) {
                    vx[i][j] = dt * Rko6s[RK] * ((-1 / rho) * dx_p_stag[i][j] - (umat[i][j] * dx_vx[i][j]) - (dumat[i][j] * vyold_yxstag[i][j])) + vxold[i][j];
                    vy[i][j] = dt * Rko6s[RK] * ((-1 / rho) * dy_p_stag[i][j] -(umat_ystag[i][j] * dx_vy[i][j])) + vyold[i][j];
                    px[i][j] = dt * Rko6s[RK] * (-rho * c_air * c_air * dx_vx_stag[i][j] - (umat[i][j] * dx_p[i][j])) + pxold[i][j];
                    py[i][j] = dt * Rko6s[RK] * (-rho * c_air * c_air * dy_vy_stag[i][j] - (umat[i][j] * dy_p[i][j])) + pyold[i][j];
                    p[i][j] = px[i][j] + py[i][j];

                }

        }
               /*Source functions */
            tnew= (n)*dt+Rko6s[RK]*dt;
            double t0=0.045;
            pulse=sin(2*atan(1)*4*fspatial/8.*tnew)*exp(-3e6*pow((fspatial/4000),2)*pow((t0-tnew),2))/2;
              pulse=sin(2*atan(1)*4*tnew*fspatial/8);


                p[cell_y_s][cell_x_s]=p[cell_y_s][cell_x_s]+Rko6s[RK]*dt*pulse;


        }


        for(i=0;i<gridy;i++) {
            for (j = 0; j < gridx; j++) {
                vx[i][j]=vx[i][j]*PMLvxmat[i][j];
                vy[i][j]=vy[i][j]*PMLvymat[i][j];
                px[i][j]=px[i][j]*PMLpxmat[i][j];
                py[i][j]=py[i][j]*PMLpymat[i][j];
                pxold[i][j] = px[i][j];
                pyold[i][j] = py[i][j];
                vxold[i][j] = vx[i][j];
                vyold[i][j] = vy[i][j];

            }
        }

          P_r[n]=p[cell_y_r][cell_x_r];

        cout<<P_r[n]<<endl;

        if(myfile.is_open()) {

            myfile << P_r[n];
            myfile<<"\n";
        }else cout<<"\n Not able to open the file\n";

    }
    myfile.close();





    /*for(i = 0; i <gridy; i++) {
        delete[] p[i];
        delete[] v[i];
        delete[] px[i];
        delete[] py[i];
        delete[] vx[i];
        delete[] vy[i];
        delete[] pxold[i];
        delete[] pyold[i];
        delete[] vxold[i];
        delete[] vyold[i];
        delete[] pxold[i];
        delete[] pyold[i];
        delete[] vxold[i];
        delete[] vyold[i];
        delete[] dx_p_stag[i];
        delete[] dy_p_stag[i];
        delete[] dx_vx_stag[i];
        delete[] dy_vy_stag[i];
        delete[] vyold_ystag[i];
        delete[] vyold_yxstag[i];
        delete[] umat[i];
        delete[] umat_ystag[i];
        delete[] dumat[i];
        delete[] kxmat[i];
        delete[] kymat[i];
    }
    delete [] p;
    delete[] v;
    delete[] px;
    delete[] py;
    delete[] vx;
    delete[] vy;
    delete[] pxold;
    delete[] pyold;
    delete[] vxold;
    delete[] vyold;
    delete[] pxold;
    delete[] pyold;
    delete[] vxold;
    delete[] vyold;
    delete[] dx_p_stag;
    delete[] dy_p_stag;
    delete[] dx_vx_stag;
    delete[] dy_vy_stag;
    delete[] vyold_ystag;
    delete[] vyold_yxstag;
    delete[] x;
    delete[] y;
    delete[] kx;
    delete[] ky;
    delete[] u;
    delete[] du;
    delete[] umat;
    delete[] umat_ystag;
    delete[] dumat;
    delete[] kxmat;
    delete[] kymat;*/


        return 0;



}