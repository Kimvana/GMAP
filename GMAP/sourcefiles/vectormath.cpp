#include <math.h>
#include <stdio.h>
#include <string.h>  // memset!


void VM_vect_at_scalar(int nentries, float *vect, float scalar) {
    for (int i = 0; i < nentries; i++) {
        vect[i] *= scalar;
    }
}


void VM_vect_at_matrix33(float *vect, float *matrix, float *outvect) {
    int a, b;
    memset(outvect, 0, 12);  // 12 bytes need to be set - float is 4 each
    for (a = 0; a < 3; a++) {
        for (b = 0; b < 3; b++) {
            outvect[b] += vect[a] * matrix[3*a + b];
        }
    }
}


void VM_many_vect_at_matrix33(
    int nvect, float *vects, float *matrix, float *outvects
) {
    for (int i = 0; i < nvect; i++) {
        VM_vect_at_matrix33(&vects[i * 3], matrix, &outvects[i * 3]);
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


float VM_dotprod(float *vec1, float *vec2) {
    /*Calculates the dot product between the two provided vectors.
    Parameters
    ----------
    vec1, vec2 : float[3]
        The vectors to calculate the dot product between.

    Returns
    -------
    dot_product : float
        The dot product of the input vectors
    */

   return (vec1[0] * vec2[0] + vec1[1] * vec2[1] + vec1[2] * vec2[2]);
}


void VM_crossprod(float *vec1, float *vec2, float *vecout) {
    vecout[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1];
    vecout[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2];
    vecout[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0];
}


void VM_project(float *vecbase, float *vecadj) {
    /*Adjusts vecadj in place to be perpendicular to vecbase.*/

    float inprod = VM_dotprod(vecbase, vecadj) / VM_dotprod(vecbase, vecbase);
    for (int i = 0; i < 3; i++) {
        vecadj[i] -= inprod * vecbase[i];
    }
}


void VM_normalize(float *vec) {
    VM_vect_at_scalar(3, vec, 1 / sqrt(VM_veclen2(vec)));
}


// original implementation
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

    // printf("remainder -2, 10 = %f\n", remainderf(-2, 10));
    // printf("remainder -8, 10 = %f\n", remainderf(-8, 10));
    for (int i = 0; i < 3; i++){
        vectout[i] = vect1[i] - vect2[i];
        if (vectout[i] > halfbox[i]) {
            vectout[i] -= boxdims[i];
        }
        else if (vectout[i] < -halfbox[i]) {
            vectout[i] += boxdims[i];
        }
    }
}

inline float posmodf(float divident, float divisor) {
    /*
    Always returns a positive mod. i.e. -2 % 10 = 8.
    */
   return fmodf(fmodf(divident, divisor) + divisor, divisor);
}


// modulo implementation (is noticably slower than orig)
// void VM_PBC_diff_cubic(
//     float *vect1, float *vect2, float *halfbox, float *boxdims, float *vectout
// ) {
//     // printf("hi, im running! -0.2 mod 1 = %f, -0.2 ownmod 1 = %f\n", fmodf(-0.2, 1.0), posmodf(-0.2, 1));
//     for (int i = 0; i < 3; i++) {
//         // printf("((%f - %f + %f) mod %f) - %f = ", vect1[i], vect2[i], halfbox[i], boxdims[i], halfbox[i]);
//         vectout[i] = posmodf(vect1[i] - vect2[i] + halfbox[i], boxdims[i]);
//         // vectout[i] -= fmodf(vectout[i], boxdims[i]);
//         vectout[i] -= halfbox[i];
//         // printf("%f\n", vectout[i]);
//     }
// }


// remainder implementation (is slightly slower than orig)
// void VM_PBC_diff_cubic(
//     float *vect1, float *vect2, float *halfbox, float *boxdims, float *vectout
// ) {
//     // printf("hi, im running! -0.2 mod 1 = %f, -0.2 ownmod 1 = %f\n", fmodf(-0.2, 1.0), posmodf(-0.2, 1));
//     for (int i = 0; i < 3; i++) {
//         // printf("((%f - %f + %f) mod %f) - %f = ", vect1[i], vect2[i], halfbox[i], boxdims[i], halfbox[i]);
//         vectout[i] = remainderf(vect1[i] - vect2[i], boxdims[i]);
//         // vectout[i] -= fmodf(vectout[i], boxdims[i]);
//         // vectout[i] -= halfbox[i];
//         // printf("%f\n", vectout[i]);
//     }
// }

// set divisor == 1 implementation (slightly slower than orig)
// void VM_PBC_diff_mod1(
//     float *vect1, float *vect2, float *halfbox, float *boxdims, float *vectout
// ) {
//     for (int i = 0; i < 3; i++){
//         vectout[i] = remainderf(vect1[i] - vect2[i], 1);
//     }
// }


// set divisor == 1 implementation (slightly slower than orig)
// void VM_PBC_diff_mod1(
//     float *vect1, float *vect2, float *halfbox, float *boxdims, float *vectout
// ) {
//     for (int i = 0; i < 3; i++){
//         vectout[i] = vect1[i] - vect2[i] + 0.5;
//         vectout[i] -= floorf(vectout[i]) + 0.5;
//     }
// }


// a mod 1 version of the original PBCdiff - same speed.
void VM_PBC_diff_mod1(
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

    // printf("remainder -2, 10 = %f\n", remainderf(-2, 10));
    // printf("remainder -8, 10 = %f\n", remainderf(-8, 10));
    for (int i = 0; i < 3; i++){
        vectout[i] = vect1[i] - vect2[i];
        if (vectout[i] > 0.5) {
            vectout[i] -= 1;
        }
        else if (vectout[i] < -0.5) {
            vectout[i] += 1;
        }
    }
}