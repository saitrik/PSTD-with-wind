//
// Created by strikootam on 9-8-2016.
//

#ifndef UNTITLED_SPATDERV_H
#define UNTITLED_SPATDERV_H

class spatderv{
public:
    void derv(  double *in_buffer, double *final_buffer, int gridx, int gridy, double **kmat, double dr,int stag);

};
#endif //UNTITLED_SPATDERV_H