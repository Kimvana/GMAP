#include <stdio.h>
#include <string.h>  // memset!


void VM_vect_at_matrix33(float *vect, float *matrix, float *outvect) {
    int a, b;
    memset(outvect, 0, 12);  // 12 bytes need to be set - float is 4 each
    for (a = 0; a < 3; a++) {
        for (b = 0; b < 3; b++) {
            outvect[b] += vect[a] * matrix[3*a + b];
        }
    }
}


float VM_veclen2(float *vec) {
    /*Calculates the square of the length of the input vector.

    Parameters
    ----------
    vec : float[3]
        The input vector.

    Returns
    -------
    len2 : float
        The length of the input vector, squared.
    */

    return (vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
}


void VM_PBC_diff_cubic(
    float *vect1, float *vect2, float *halfbox, float *boxdims,
    float *vectout
) {
    /*Calculates the difference vect1 - vect2, assuming a cubic MD
    system, and assuming vect1 and vect2 are inside the system
    currently.

    Parameters
    ----------
    vect1, vect2 : float[3]
        The positions between which the difference vector should be
        calculated.
    halfbox, boxdims : float[3]
        The size of the PBC (MD system size). Halfbox is assumed to
        equal boxdims/2.
    vectout : float[3]
        The output will be written here. It is the difference vector
        between vect1 and vect2, corrected for the PBC.
    */

    for (int i = 0; i < 3; i++){
        vectout[i] = vect1[i] - vect2[i];
        if (vectout[i] > halfbox[i]) {
            vectout[i] -= boxdims[i];
        }
        else if (vectout[i] < -1 * halfbox[i]) {
            vectout[i] += boxdims[i];
        }
    }
}
