#include <math.h>
#include <stdio.h>
#include <stdlib.h>  // calloc!

// These files are not from default libraries, or within this map folder.
// Instead, they should be manually included during installation:
// windows (my machine, edit path):
// cl.exe /LD /Fe: TCC_clib_Win64bit /I\github\GEMAIM-dev\GMAP\sourcefiles TCC_clib.cpp
#include "vectormath.cpp"  // in GMAP sourcefiles directory

#ifdef _WIN32
    extern "C" {
        __declspec(dllexport) void prep_coupling(int nosc, int noscats, int *oscixarr, int *osc_used_ats, float *positions_box, float *boxvects, int *dopro, float alpha_gen, float alpha_pro, float *v_gen, float *v_pro, float *tcc_v);
        __declspec(dllexport) void calc_coupling(int npairs, int *allpairs, int noscats, int *oscix_to_ix, int *osc_used_ats, float *positions_box, float *boxvects, int *dopro, float *q_gen, float *q_pro, float *dq_gen, float *dq_pro, float fourPiEps, float *tcc_v, int totosc, float *hamiltonian);
    }
#endif


extern "C" {

    void prep_coupling(
        int nosc, int noscats, int *oscixarr, int *osc_used_ats,
        float *positions_box, float *boxvects, int *dopro, float alpha_gen,
        float alpha_pro, float *v_gen, float *v_pro, float *tcc_v) {
        /*
        nosc - length of oscixlist provided to the map in the
        prep_coupling step)
        noscats - the amount of atoms in a single oscillator
        oscixarr - a np array version of oscixlist.
        osc_used_ats - a np array version of used_atoms of ALL oscillators
        positions_box - the positions of all atoms in the system in box coordinates
        boxvects - the boxvects array
        */

        int ix;
        // as we need all 3 vectors in a single matrix, only make that final
        // matrix. All mathfuncts expect a 3-vector, so by giving the
        // correct position in the matrix, only 3 entries will actually
        // be written/used.
        float tempvec[3], rotmat[9];
        float halfbox[3] = {0.5, 0.5, 0.5};
        float boxdims[3] = {1., 1., 1.};

        // we loop over all oscillators to be coupled by this map.
        for (ix = 0; ix < nosc; ix++) {
            /*----------------------------------------
            Build up rotation matrix of the oscillator
            ----------------------------------------*/
            // CObox = boxpos[osc.used_atoms[1]] - boxpos[osc.used_atoms[0]]
            VM_PBC_diff_cubic(
                &positions_box[osc_used_ats[noscats * ix + 1] * 3],
                &positions_box[osc_used_ats[noscats * ix] * 3],
                halfbox, boxdims, tempvec);
            // COvec = CObox @ boxvects
            VM_vect_at_matrix33(tempvec, boxvects, rotmat);
            VM_normalize(rotmat);  // COvec /= vec3len(COvec)

            // CNbox = boxpos[osc.used_atoms[3]] - boxpos[osc.used_atoms[0]]
            VM_PBC_diff_cubic(
                &positions_box[osc_used_ats[noscats * ix + 3] * 3],
                &positions_box[osc_used_ats[noscats * ix] * 3],
                halfbox, boxdims, tempvec);
            // CNvec = CNbox @ boxvects
            VM_vect_at_matrix33(tempvec, boxvects, &rotmat[3]);
            VM_project(rotmat, &rotmat[3]);  // CNvec = project (COvec, CNvec)
            VM_normalize(&rotmat[3]);  // CNvec /= vec3len(CNvec)

            // zvec = cross(COvec, CNvec)
            VM_crossprod(rotmat, &rotmat[3], &rotmat[6]);
            VM_normalize(&rotmat[6]);  // zvec /= vec3len(zvec)

            /*---------------------------------------
            Calculate the v-term for the TCC coupling
            ---------------------------------------*/
            if (dopro[ix] == 1) {
                VM_many_vect_at_matrix33(6, v_pro, rotmat, &tcc_v[ix * 18]);
                VM_vect_at_scalar(18, &tcc_v[ix * 18], alpha_pro);
            } else {
                VM_many_vect_at_matrix33(6, v_gen, rotmat, &tcc_v[ix * 18]);
                VM_vect_at_scalar(18, &tcc_v[ix * 18], alpha_gen);
            }
        }
    }


    void calc_coupling(int npairs, int *allpairs, int noscats, int *oscix_to_ix, int *osc_used_ats, float *positions_box, float *boxvects, int *dopro, float *q_gen, float *q_pro, float *dq_gen, float *dq_pro, float fourPiEps, float *tcc_v, int totosc, float *hamiltonian) {
        int pairix, oscix1, oscix2, ix1, ix2, a, b;
        float *q1, *q2, *dq1, *dq2, tempvec[3], diff[3], *v1, *v2;
        float r2, ir2, ir, ir3, ir5, J;
        float halfbox[3] = {0.5, 0.5, 0.5};
        float boxdims[3] = {1., 1., 1.};

        for (pairix = 0; pairix < npairs; pairix++) {
            oscix1 = allpairs[pairix * 2];
            ix1 = oscix_to_ix[oscix1];
            oscix2 = allpairs[pairix * 2 + 1];
            ix2 = oscix_to_ix[oscix2];
            if (dopro[ix1] == 0) {
                q1 = q_gen;
                dq1 = dq_gen;
            } else {
                q1 = q_pro;
                dq1 = dq_pro;
            }
            if (dopro[ix2] == 0) {
                q2 = q_gen;
                dq2 = dq_gen;
            } else {
                q2 = q_pro;
                dq2 = dq_pro;
            }

            v1 = &tcc_v[ix1 * 18];
            v2 = &tcc_v[ix2 * 18];

            J = 0;
            
            // loop over all pairs of atoms between these two oscs
            for (a = 0; a < 6; a++) {
                for (b = 0; b < 6; b++) {
                    // diff vector between two atoms
                    VM_PBC_diff_cubic(
                        &positions_box[osc_used_ats[noscats * ix2 + b] * 3],
                        &positions_box[osc_used_ats[noscats * ix1 + a] * 3],
                        halfbox, boxdims, tempvec);
                    VM_vect_at_matrix33(tempvec, boxvects, diff);

                    // get all distance terms
                    r2 = VM_veclen2(diff);
                    if (r2 < 0.01) {r2 = 1.00;}
                    ir2 = 1 / r2;
                    ir = sqrt(ir2);
                    ir3 = ir * ir2;
                    ir5 = ir3 * ir2;

                    J -= 3 * ir5 * q1[a] * q2[b] * VM_dotprod(&v2[b * 3], diff) * VM_dotprod(&v1[a * 3], diff);
                    J -= ir3 * (dq1[a] * q2[b] * VM_dotprod(&v2[b * 3], diff) - q1[a] * dq2[b] * VM_dotprod(&v1[a * 3], diff) - VM_dotprod(&v1[a * 3], &v2[b * 3]) * q1[a] * q2[b]);
                    J += ir * dq1[a] * dq2[b];
                }
            }
            J *= fourPiEps;
            hamiltonian[oscix1 * totosc + oscix2] = J;
            hamiltonian[oscix2 * totosc + oscix1] = J;
        }
    }

}
