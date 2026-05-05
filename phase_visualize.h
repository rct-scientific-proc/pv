#ifndef PHASE_VISUALIZE_H
#define PHASE_VISUALIZE_H

#include <stdint.h>

/* ----- Shared-library export/import macros -----
 *
 * When building the shared library, the build system should define
 *   PHASE_VISUALIZE_BUILD_DLL
 * When consuming the shared library on Windows, no define is required
 * (declarations default to dllimport).  On non-Windows platforms PV_API
 * expands to default visibility on GCC/Clang and nothing otherwise.
 */
#if defined(_WIN32) || defined(__CYGWIN__)
  #ifdef PHASE_VISUALIZE_BUILD_DLL
    #define PV_API __declspec(dllexport)
  #else
    #define PV_API __declspec(dllimport)
  #endif
#else
  #if defined(__GNUC__) || defined(__clang__)
    #define PV_API __attribute__((visibility("default")))
  #else
    #define PV_API
  #endif
#endif

#ifdef __cplusplus
extern "C" {
#endif

// Function declarations for phase visualization
// 1. Phase unwrapping - Goldstein's branch cut method
// 2. Sine/cosine decomposition
// 3. Phase gradient calculation / fringe-rate maps
// 4. Coherence-weighted phase


// 1. Phase unwrapping - Goldstein's branch cut method
// Data should be column-major. Like matlab. Fastest dimension is down the rows.
//
// Inputs:
//   wrapped_phase   : rows*cols array, values in (-pi, pi]
//   rows, cols      : image dimensions
// Output:
//   unwrapped_phase : rows*cols array (caller-allocated)
//
// Returns 0 on success, negative on error:
//   -1 : null pointer
//   -2 : invalid dimensions (rows < 2 or cols < 2)
//   -3 : memory allocation failure
PV_API int phase_unwrap_goldstein(const double *wrapped_phase, double *unwrapped_phase, int32_t rows, int32_t cols);

// Same as phase_unwrap_goldstein, but exposes the maximum branch-cut search
// box size (Goldstein "box growing" cap). This is the dominant performance
// knob: search cost grows with max_box^2 per residue.
//
//   max_box <= 0 : use the built-in default (9, per Ghiglia & Pritt).
//   otherwise    : value is clamped to [3, min(rows,cols)] and forced odd.
//
// Smaller values (e.g. 5-9) are dramatically faster and usually adequate
// for noisy interferograms. Larger values can resolve more residue pairs at
// significant runtime cost.
PV_API int phase_unwrap_goldstein_ex(const double *wrapped_phase, double *unwrapped_phase,
                                     int32_t rows, int32_t cols, int32_t max_box);


// 2. Sine/cosine decomposition.
// Splits a wrapped phase image into its in-phase (cosine) and quadrature
// (sine) components. Each output is continuous in [-1, 1] and free of the
// 2*pi wrap discontinuity, which makes them suitable as two real channels for
// visualization or downstream filtering.
//
// Data is column-major; the fastest dimension is down the rows.
//
// Inputs:
//   phase  : numel-element array of wrapped phase values (any range; typically
//            (-pi, pi]).
//   numel  : total number of elements (e.g. rows*cols).
// Outputs:
//   cos_out : numel-element array, cos(phase[i]).  May be NULL to skip.
//   sin_out : numel-element array, sin(phase[i]).  May be NULL to skip.
//
// At least one of cos_out / sin_out must be non-NULL.
//
// Returns 0 on success, negative on error:
//   -1 : null pointer (phase is NULL, or both outputs are NULL)
PV_API int phase_sincos_decompose(const double *phase, double *cos_out, double *sin_out, int32_t numel);


// 3. Phase gradient / fringe-rate map.
// Computes the spatial gradient of a wrapped phase image with correct 2*pi
// wrap handling (each finite difference is wrapped into (-pi, pi] before
// being scaled), so a wrap discontinuity does not produce a spurious huge
// gradient. Centered differences are used in the interior; one-sided
// (forward/backward) differences are used at the borders.
//
// Data is column-major; the fastest dimension is down the rows.
//
// Conventions:
//   grad_row[i] approximates d(phase) / d(row)  (units: radians per pixel)
//   grad_col[i] approximates d(phase) / d(col)  (units: radians per pixel)
//   fringe_rate[i] = sqrt( grad_row[i]^2 + grad_col[i]^2 )
//
// Any of grad_row, grad_col, fringe_rate may be NULL to skip that output, but
// at least one must be non-NULL.
//
// Inputs:
//   phase       : rows*cols array of wrapped phase values, typically (-pi, pi].
//   rows, cols  : image dimensions (each must be >= 2).
// Outputs:
//   grad_row    : rows*cols array (or NULL).
//   grad_col    : rows*cols array (or NULL).
//   fringe_rate : rows*cols array (or NULL).
//
// Returns 0 on success, negative on error:
//   -1 : null pointer (phase NULL, or all three outputs NULL)
//   -2 : invalid dimensions (rows < 2 or cols < 2)
PV_API int phase_gradient(const double *phase,
                          double *grad_row, double *grad_col, double *fringe_rate,
                          int32_t rows, int32_t cols);


// 4a. Local phase coherence (phase-only estimator).
// For each pixel, computes the magnitude of the mean unit phasor in a
// (2*window+1) x (2*window+1) neighborhood:
//
//      gamma(r,c) = | (1/N) * sum_{i,j in window} exp(j * phase(i,j)) |
//
// gamma is in [0, 1].  gamma == 1 means all phases in the window are equal
// (perfectly coherent); gamma near 0 means random phase (incoherent / noisy).
// This is useful as a quality / weighting map even when no SAR magnitude is
// available.
//
// Data is column-major; the fastest dimension is down the rows. Borders are
// handled by clipping the window to the image (no padding).
//
// Inputs:
//   phase       : rows*cols array of wrapped phase values, typically (-pi, pi].
//   rows, cols  : image dimensions (each must be >= 1).
//   window      : half-width of the window (>= 1). Window size is
//                 (2*window+1)^2 pixels.
// Output:
//   coherence   : rows*cols array (caller-allocated), values in [0, 1].
//
// Returns 0 on success, negative on error:
//   -1 : null pointer
//   -2 : invalid dimensions or window < 1
PV_API int phase_coherence(const double *phase, double *coherence,
                           int32_t rows, int32_t cols, int32_t window);


// 4b. Coherence-weighted phase visualization (HSV -> RGB).
// Builds an RGB image where:
//   - Hue        = phase mapped from (-pi, pi] to [0, 1)  (cyclic)
//   - Saturation = coherence (or constant 1.0 if coherence is NULL)
//   - Value      = coherence (or constant 1.0 if coherence is NULL)
// The result dims low-coherence (noisy) pixels toward black/grey while
// preserving full-color contrast in coherent regions.
//
// Outputs are in [0, 1].  Data is column-major; each of r_out, g_out, b_out
// is a separate rows*cols plane (planar RGB), matching the rest of this API.
//
// Inputs:
//   phase     : rows*cols array of wrapped phase values, typically (-pi, pi].
//   coherence : rows*cols array of values in [0, 1], or NULL for unweighted.
//   rows, cols: image dimensions (each must be >= 1).
// Outputs:
//   r_out, g_out, b_out : rows*cols arrays, values in [0, 1] (caller-allocated).
//
// Returns 0 on success, negative on error:
//   -1 : null pointer (phase or any output channel NULL)
//   -2 : invalid dimensions
PV_API int phase_to_rgb_hsv(const double *phase, const double *coherence,
                            double *r_out, double *g_out, double *b_out,
                            int32_t rows, int32_t cols);


#ifdef __cplusplus
}
#endif

#endif /* PHASE_VISUALIZE_H */
