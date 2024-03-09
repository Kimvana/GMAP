#include <math.h>calcPot_subbox
#include <stdio.h>
#include <algorithm>

// None of this code is tested, it might not even compile!

extern "C" {
    // Generates a positions list in boxdim space from a normal positions list.
    void tric_positions(float *positions, int Natoms,
                        float inverse_matrix[3][3], float *tric_positions){
        for (int i = 0; i < Natoms; i++){
            for (int row = 0; row < 3; row++){
                for (int col = 0; col < 3; col++){
                    tric_positions[i * 3 + row] += positions[i + col] * inverse_matrix[col][row];
                }
            }
        }
    }


    // Uses vectors in space defined by boxvectors, returns distance 
    // vector in Cartesian coordinates accounting for PBC.
    void tric_PBC_diff(float *tric_vect1, float *tric_vect2,
                       float matrix[3][3], float *vectout){
        float temp_vect[3];
        for (int i = 0; i < 3; i++) {
            vectout[i] = 0;
            temp_vect[i] = tric_vect1[i] - tric_vect2[i];
            if (temp_vect[i] > 0.5){
                temp_vect[i] -= 1;
            }
            if (temp_vect[i] < -0.5){
                temp_vect[i] += 1;
            }
        }
        for (int row = 0; row < 3; row++){
            for (int col = 0; col < 3; col++){
                vectout[row] += temp_vect[col] * matrix[col][row];
                }
            }
    }

    // this function has the same signature as the 'smart' one for aliassing 
    // purposes, but the halfbox and boxdims arguments are not actually used!
    void PBC_diff_dumb(float *vect1, float *vect2,
                       float matrix[3][3], float *vectout) {
        for (int i = 0; i < 3; i++) {
            vectout[i] = vect1[i] - vect2[i];
        }
    }

    float distance2(float *vect){
        return vect[0] * vect[0] + vect[1] * vect[1] + vect[2] * vect[2];
    }

    float distance(float *vect){
        float distance2 = vect[0] * vect[0] + vect[1] * vect[1]
        + vect[2] * vect[2];
        return sqrt(distance2);
    }

    void potential(float *vect, float distance, float charge, float potential){
        potential += (charge / distance);
    }

    void E_vect(float *vect, float distance, float charge, float E_vect[3]){
        float E_vect_base = charge / distance * distance * distance;
        
        E_vect[0] = E_vect_base * vect[0];
        E_vect[1] = E_vect_base * vect[1];
        E_vect[2] = E_vect_base * vect[2];
    }

    void G_vect(){}

}