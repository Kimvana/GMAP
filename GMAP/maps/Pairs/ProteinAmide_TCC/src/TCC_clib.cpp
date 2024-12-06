#include <math.h>
#include <stdio.h>
#include <stdlib.h>  // calloc!

#ifdef _WIN32
    extern "C" {
        __declspec(dllexport) void testme(
        );
    }
#endif


extern "C" {
    void prep_coupling(int nosc, int *oscixarr) {
        int ix, oscix;
        for (ix = 0; ix < nosc; ix++) {
            oscix = oscixarr[ix];
            
        }
    }
}
