//
// Created by strikootam on 17-8-2016.
//

#ifndef UNTITLED_PML_H
#define UNTITLED_PML_H
#include<cmath>

class PML {
public:
   void PMLlayer(int gridx, int gridy, int PMLcells, double dt, double **PMLpxmat, double **PMLpymat, double **PMLvxmat, double **PMLvymat);

};


#endif //UNTITLED_PML_H

