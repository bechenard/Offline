#ifndef Mu2eUtilities_polySinCos_hh
#define Mu2eUtilities_polySinCos_hh
//
//
//
// Original author Giani Pezzullo
//

// Modern CLHEP
#include "CLHEP/Vector/ThreeVector.h"


namespace mu2e {

  inline void polySinCos(float a, float& s, float& c)
  {
    //-- Step 1: Range reduction in double to [-pi/4, pi/4]
    constexpr float  TWO_OVER_PI_F = 6.36619772367581382e-1f;
    constexpr double PI_OVER_2_D   = 1.57079632679489661923e+0;

    const float  jf = rintf(a * TWO_OVER_PI_F);
    const int    j  = static_cast<int>(jf);
    const float  ar = static_cast<float>(
                        static_cast<double>(a) - jf * PI_OVER_2_D);

    //-- Step 2: Minimax polynomial evaluation on [-pi/4, pi/4]
    const float a2 = ar * ar;

    // sin(ar) ~ ar * P(ar^2),  polynomial degree 7 in ar,  < 1 ULP
    const float sa = ar * (1.0f
                   + a2 * (-1.6666654611e-1f
                   + a2 * ( 8.3321608736e-3f
                   + a2 * (-1.9515295891e-4f))));

    // cos(ar) ~ Q(ar^2),  polynomial degree 8 in ar,  < 1 ULP
    const float ca = 1.0f
                   + a2 * (-4.9999999630e-1f
                   + a2 * ( 4.1666645680e-2f
                   + a2 * (-1.3887731625e-3f
                   + a2 *   2.4431511800e-5f)));

    //-- Step 3: Branch-free quadrant reconstruction
    const int q = j & 3;

    const bool swap    = (q == 1) || (q == 3);
    const bool sin_neg = (q == 2) || (q == 3);
    const bool cos_neg = (q == 1) || (q == 2);

    const float s_base = swap ? ca : sa;
    const float c_base = swap ? sa : ca;

    s = sin_neg ? -s_base : s_base;
    c = cos_neg ? -c_base : c_base;
  }

}
#endif /* Mu2eUtilities_polySinCos_hh */
