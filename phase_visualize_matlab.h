/*
 * phase_visualize_matlab.h
 *
 * Plain-C declarations for use with MATLAB's loadlibrary().
 * MATLAB's loadlibrary cannot handle __declspec(dllimport) or many macros,
 * so this header restates each export with no decorations.
 *
 * Keep these signatures byte-for-byte in sync with phase_visualize.h.
 */
#ifndef PHASE_VISUALIZE_MATLAB_H
#define PHASE_VISUALIZE_MATLAB_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

int phase_unwrap_goldstein(const double *wrapped_phase, double *unwrapped_phase,
                           int32_t rows, int32_t cols);

int phase_unwrap_goldstein_ex(const double *wrapped_phase, double *unwrapped_phase,
                              int32_t rows, int32_t cols, int32_t max_box);

int phase_sincos_decompose(const double *phase, double *cos_out, double *sin_out,
                           int32_t numel);

int phase_gradient(const double *phase,
                   double *grad_row, double *grad_col, double *fringe_rate,
                   int32_t rows, int32_t cols);

int phase_coherence(const double *phase, double *coherence,
                    int32_t rows, int32_t cols, int32_t window);

int phase_to_rgb_hsv(const double *phase, const double *coherence,
                     double *r_out, double *g_out, double *b_out,
                     int32_t rows, int32_t cols);

#ifdef __cplusplus
}
#endif

#endif /* PHASE_VISUALIZE_MATLAB_H */
