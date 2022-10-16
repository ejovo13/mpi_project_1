#ifndef GEODESY_CONSTANTS_H
#define GEODESY_CONSTANTS_H

/**========================================================================
 * ?                          geodesy.h
 * @brief   : Define special values to be used throughout the 
 *            geodesy project
 * @details : 
 * @author  : Evan Voyles
 * @email   : ejovo13@yahoo.com
 * @date    : 2022-10-15
 *========================================================================**/
// #ifndef PI
// #define PI 3.1415926535897932384626433832795028841971693993751058209749445923078164062
// // const extern int PI;
// #endif



// #ifndef TWO_PI
// #define TWO_PI (2.0 * PI)
// // const extern int two_pi;
// #endif

// #ifndef HALF_PI
// #define HALF_PI (PI / 2.0)
// #endif

extern const double PI;
extern const double TWO_PI;
extern const double HALF_PI;
extern const double EPS;

// #ifndef SPHERICAL_MODEL_CONSTANTS
// #define SPHERICAL_MODEL_CONSTANTS

//     // ranges of th and ph
//     #define TH_0 0
//     #define TH_F PI
//     #define PH_0 0
//     #define PH_F TWO_PI

// #endif

/**========================================================================
 *!                           Looping macros
 *========================================================================**/
#define FOR_I_IN(mat, inner) for (size_t __i = 0, __n = Matrix_size_d(mat), i = 0; __i < __n; __i++, i = mat->data[__i]) { \
    i = mat->data[__i]; \
    inner \
    }
// iterate along the rows
#define FOR_I(mat, inner) for (size_t i = 0, __n = mat->nrows; i < __n; i++) { inner \ }
// #define FOR_N(mat, inner) for (size_t i = 0; i < Matrix_size_d(mat))
#define FOR_IJ(mat, inner) for (size_t i = 0; i < mat->nrows; i++) { \
    for (size_t j = 0; j < mat->ncols; j++) { \
        inner \
    }\
}

// iterate along the cols
#define FOR_J(mat, inner) for (size_t __j = 0, __n = mat->nrows)

#endif // GEODESY_CONSTANTS_H

