#pragma once
#include <cmath>

//  - taken from hnswlib/.../hnswlib.h
template<typename MTYPE>
using DISTFUNC = MTYPE(*)(const void *, const void *, const void *);

/*
    Borrowing from https://github.com/nmslib/hnswlib
*/
namespace metrics {

    //  - taken from hnswlib/.../space_l2.h
    static float
    EuclideanDistance(const void *pVect1v, const void *pVect2v, const void *qty_ptr) {
        float *pVect1 = (float *) pVect1v;
        float *pVect2 = (float *) pVect2v;
        size_t qty = *((size_t *) qty_ptr);

        float res = 0;
        for (size_t i = 0; i < qty; i++) {
            float t = *pVect1 - *pVect2;
            pVect1++;
            pVect2++;
            res += t * t;
        }
        return std::sqrt(res);
    }
};