#include <math.h>
#include <stdio.h>
#include <stdlib.h>  // calloc!

// These files are not from default libraries, or within this map folder.
// Instead, they should be manually included during installation:

// windows (my machine, edit path):
// cl.exe /LD /Fe: TRESP_clib_Win64bit /I\github\GEMAIM-dev\GMAP\sourcefiles TRESP_clib.cpp

// linux (Kai's cluster, edit path):
// g++ -fPIC -shared -o TRESP_clib_Linux.so -I/scratch/p302934/GMAP_fin/GMAP/GMAP/sourcefiles TRESP_clib.cpp
#include "vectormath.cpp"  // in GMAP sourcefiles directory

#ifdef _WIN32
    extern "C" {
        __declspec(dllexport) void calc_coupling(int npairs, int *allpairs, int *noscats, int *oscstart, float *diff_charges, int *osc_used_ats, float *positions_box, float *boxvects, int totosc, float *hamiltonian);
    }
#endif


extern "C" {
    
    void calc_coupling(
        int npairs, int *allpairs, int *noscats, int *oscstart,
        float *diff_charges, int *osc_used_ats, float *positions_box,
        float *boxvects, int totosc, float *hamiltonian
    ) {
        int pairix, oscix1, oscix2, ix1, ix2, osc1len, osc2len, TRix1, TRix2;
        float tempvec[3], diff[3], r2, ir, J;
        for (pairix = 0; pairix < npairs; pairix++) {
            // oscillator indices
            oscix1 = allpairs[pairix * 2];
            oscix2 = allpairs[pairix * 2 + 1];
            
            J = 0;

            // amount of atoms per oscilator
            osc1len = noscats[oscix1];
            osc2len = noscats[oscix2];

            // index in all-atom TRESP arrays
            TRix1 = oscstart[oscix1];
            for (ix1 = 0; ix1 < osc1len; ix1++) {
                TRix2 = oscstart[oscix2];
                for (ix2 = 0; ix2 < osc2len; ix2++) {
                    // calculate difference vector between the two atoms
                    VM_PBC_diff_mod1(
                        &positions_box[osc_used_ats[TRix1] * 3],
                        &positions_box[osc_used_ats[TRix2] * 3],
                        // 3rd and 4th argument aren't used; last is.
                        diff, diff, tempvec);
                    VM_vect_at_matrix33(tempvec, boxvects, diff);
                    r2 = VM_veclen2(diff);
                    ir = r2 < 0.01 ? 1 : 1/sqrt(r2);

                    J += diff_charges[TRix1] * diff_charges[TRix2] * ir;
                    TRix2++;
                }
                TRix1++;
            }
            J *= 116141.70590152;
            hamiltonian[oscix1 * totosc + oscix2] = J;
            hamiltonian[oscix2 * totosc + oscix1] = J;
        }
    }
}
