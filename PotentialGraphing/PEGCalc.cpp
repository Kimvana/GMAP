#include <math.h>calcPot_subbox
#include <stdio.h>
#include <algorithm>

extern "C" {

    float charge_hard(){return 0};
    float charge_linear(){return 0};


    // Generates a positions list in boxdim space from a normal positions list.
    void tric_positions(float *positions, int Natoms, float inverse_matrix[3][3], float *tric_positions){
        for (int i = 0; i < Natoms; i++){
            for (int row = 0; row < 3; row++){
                for (int col = 0; col < 3; col++){
                    tric_positions[i * 3 + row] += positions[i + col] * inverse_matrix[col][row];
                }
            }
        }
    }
    
    // Based on boxdim space input vectors and PBC, returns distance 
    // vector in Cartesian coordinates .
    void tric_PBC_diff(float *tric_vect1, float *tric_vect2, float matrix[3][3], float *vectout){
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

    void calc(int smoothing_type, int distance_type, 
    float smoothing_distance){
        typedef float (*charge_weight)(float, float);
    }

}