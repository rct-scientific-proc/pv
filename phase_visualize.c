#include "phase_visualize.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ----- Column-major indexing helper ----- */
#define IDX(r, c, rows) ((size_t)(c) * (size_t)(rows) + (size_t)(r))

/* Wrap a value into (-pi, pi]. */
static inline double wrap_pi(double x)
{
    /* Use the standard "shift, fmod, shift" trick. */
    double y = fmod(x + M_PI, 2.0 * M_PI);
    if (y <= 0.0) {
        y += 2.0 * M_PI;
    }
    return y - M_PI;
}

/* ----- Bresenham line: mark every (r,c) on the line into the cut mask. ----- */
static void draw_cut_line(uint8_t *cut, int32_t rows, int32_t cols,
                          int32_t r0, int32_t c0, int32_t r1, int32_t c1)
{
    int32_t dr = abs(r1 - r0);
    int32_t dc = abs(c1 - c0);
    int32_t sr = (r0 < r1) ? 1 : -1;
    int32_t sc = (c0 < c1) ? 1 : -1;
    int32_t err = (dr > dc ? dr : -dc) / 2;
    int32_t r = r0, c = c0;

    for (;;) {
        if (r >= 0 && r < rows && c >= 0 && c < cols) {
            cut[IDX(r, c, rows)] = 1;
        }
        if (r == r1 && c == c1) {
            break;
        }
        int32_t e2 = err;
        if (e2 > -dc) { err -= dc; r += sr; }
        if (e2 <  dr) { err += dr; c += sc; }
    }
}

/*
 * Goldstein box-growing branch-cut placement.
 *
 * residues: (rows-1) x (cols-1) array, values in {-1, 0, +1}.
 *           Indexed as IDX(rr, cc, rows-1) with rr in [0, rows-2], cc in [0, cols-2].
 *
 * Each unbalanced residue (r,c) on the residue grid corresponds geometrically
 * to the corner shared by the four phase pixels (r,c), (r+1,c), (r,c+1),
 * (r+1,c+1).  For cut-marking we use the phase pixel (r,c) as its location.
 */
static void place_branch_cuts(int8_t *residues, uint8_t *cut,
                              int32_t rows, int32_t cols, int32_t max_box)
{
    const int32_t rrows = rows - 1;
    const int32_t rcols = cols - 1;
    const int32_t total = rrows * rcols;

    /* "active" tracks which residues still need balancing in the current
     * box-growing search; reset between seeds. */
    uint8_t *active = (uint8_t *)calloc((size_t)total, sizeof(uint8_t));
    if (active == NULL) {
        return; /* fall through silently; caller has its own error path */
    }

    for (int32_t rc = 0; rc < rcols; ++rc) {
        for (int32_t rr = 0; rr < rrows; ++rr) {
            const size_t sidx = IDX(rr, rc, rrows);
            if (residues[sidx] == 0 || active[sidx]) {
                continue;
            }

            /* Start a new box-growing search seeded on this residue. */
            int32_t charge = residues[sidx];
            active[sidx] = 1;
            int32_t cr = rr, cc = rc;       /* center of growing box */
            int32_t balanced = 0;
            int32_t connected_to_border = 0;

            for (int32_t bs = 3; bs <= max_box && !balanced; bs += 2) {
                const int32_t half = bs / 2;
                const int32_t r_lo = cr - half;
                const int32_t r_hi = cr + half;
                const int32_t c_lo = cc - half;
                const int32_t c_hi = cc + half;

                /* If the box reaches the residue-grid border (which is one
                 * pixel inside the phase border), connect the seed to the
                 * nearest phase border pixel and stop. */
                if (r_lo <= 0 || r_hi >= rrows - 1 ||
                    c_lo <= 0 || c_hi >= rcols - 1) {
                    int32_t br = (rr < rrows - 1 - rr) ? 0 : rows - 1;
                    int32_t bc = (rc < rcols - 1 - rc) ? 0 : cols - 1;
                    /* Pick whichever border is closer (row vs column). */
                    int32_t drow = (rr < rrows - 1 - rr) ? rr : (rrows - 1 - rr);
                    int32_t dcol = (rc < rcols - 1 - rc) ? rc : (rcols - 1 - rc);
                    if (drow <= dcol) {
                        draw_cut_line(cut, rows, cols, rr, rc, br, rc);
                    } else {
                        draw_cut_line(cut, rows, cols, rr, rc, rr, bc);
                    }
                    connected_to_border = 1;
                    balanced = 1;
                    break;
                }

                /* Scan all residues inside the box. */
                for (int32_t jc = c_lo; jc <= c_hi && !balanced; ++jc) {
                    for (int32_t jr = r_lo; jr <= r_hi && !balanced; ++jr) {
                        if (jr == cr && jc == cc) continue;
                        const size_t jidx = IDX(jr, jc, rrows);
                        const int8_t res = residues[jidx];
                        if (res == 0 || active[jidx]) continue;

                        /* Connect this residue to the seed with a cut. */
                        draw_cut_line(cut, rows, cols, rr, rc, jr, jc);
                        active[jidx] = 1;
                        charge += res;

                        if (charge == 0) {
                            balanced = 1;
                            break;
                        }
                        /* Otherwise grow box centered on the new residue. */
                        cr = jr;
                        cc = jc;
                    }
                }
            }

            /* If we exhausted the maximum box size without balancing,
             * connect the unbalanced charge to the nearest border. */
            if (!balanced && !connected_to_border) {
                int32_t br = (rr < rrows - 1 - rr) ? 0 : rows - 1;
                int32_t bc = (rc < rcols - 1 - rc) ? 0 : cols - 1;
                int32_t drow = (rr < rrows - 1 - rr) ? rr : (rrows - 1 - rr);
                int32_t dcol = (rc < rcols - 1 - rc) ? rc : (rcols - 1 - rc);
                if (drow <= dcol) {
                    draw_cut_line(cut, rows, cols, rr, rc, br, rc);
                } else {
                    draw_cut_line(cut, rows, cols, rr, rc, rr, bc);
                }
            }

            /* Reset active mask for next seed (cuts already drawn). */
            memset(active, 0, (size_t)total * sizeof(uint8_t));
        }
    }

    free(active);
}

int phase_unwrap_goldstein(const double *wrapped_phase, double *unwrapped_phase,
                           int32_t rows, int32_t cols)
{
    /* ---- Validate inputs ---- */
    if (wrapped_phase == NULL || unwrapped_phase == NULL) {
        return -1;
    }
    if (rows < 2 || cols < 2) {
        return -2;
    }

    const size_t N        = (size_t)rows * (size_t)cols;
    const int32_t rrows   = rows - 1;
    const int32_t rcols   = cols - 1;
    const size_t Nres     = (size_t)rrows * (size_t)rcols;

    /* ---- Allocate ALL working buffers up-front ---- */
    int8_t   *residues = (int8_t   *)calloc(Nres, sizeof(int8_t));
    uint8_t  *cut      = (uint8_t  *)calloc(N,    sizeof(uint8_t));
    uint8_t  *visited  = (uint8_t  *)calloc(N,    sizeof(uint8_t));
    int32_t  *queue    = (int32_t  *)calloc(N,    sizeof(int32_t)); /* indices */

    if (residues == NULL || cut == NULL || visited == NULL || queue == NULL) {
        free(residues); free(cut); free(visited); free(queue);
        return -3;
    }

    /* ---- Stage 1: identify residues (parallel over plaquettes) ---- */
    {
        const int32_t total = rrows * rcols;
        int32_t k;
        #pragma omp parallel for
        for (k = 0; k < total; ++k) {
            const int32_t rr = k % rrows;
            const int32_t rc = k / rrows;

            const double p00 = wrapped_phase[IDX(rr,     rc,     rows)];
            const double p10 = wrapped_phase[IDX(rr + 1, rc,     rows)];
            const double p01 = wrapped_phase[IDX(rr,     rc + 1, rows)];
            const double p11 = wrapped_phase[IDX(rr + 1, rc + 1, rows)];

            /* Sum wrapped differences around the 2x2 loop:
             *   (0,0) -> (0,1) -> (1,1) -> (1,0) -> (0,0) */
            const double d1 = wrap_pi(p01 - p00);
            const double d2 = wrap_pi(p11 - p01);
            const double d3 = wrap_pi(p10 - p11);
            const double d4 = wrap_pi(p00 - p10);
            const double s  = d1 + d2 + d3 + d4;

            if      (s >  M_PI) residues[k] = +1;
            else if (s < -M_PI) residues[k] = -1;
            else                residues[k] =  0;
        }
    }

    /* ---- Stage 2: place branch cuts (Goldstein box-growing) ---- */
    {
        int32_t max_box = (rows < cols ? rows : cols) / 2;
        if (max_box < 3)  max_box = 3;
        if ((max_box & 1) == 0) max_box += 1; /* must be odd */
        place_branch_cuts(residues, cut, rows, cols, max_box);
    }

    /* ---- Stage 3: flood-fill integration ---- */
    /* Find first non-cut pixel as seed. */
    int32_t seed = -1;
    for (size_t i = 0; i < N; ++i) {
        if (!cut[i]) { seed = (int32_t)i; break; }
    }

    if (seed < 0) {
        /* Whole image is cut; just copy wrapped phase. */
        memcpy(unwrapped_phase, wrapped_phase, N * sizeof(double));
        free(residues); free(cut); free(visited); free(queue);
        return 0;
    }

    /* Initialize output to NaN so unfilled pixels are obvious during cleanup. */
    {
        const double NaN = 0.0 / 0.0;
        for (size_t i = 0; i < N; ++i) unwrapped_phase[i] = NaN;
    }

    int32_t qhead = 0, qtail = 0;
    unwrapped_phase[seed] = wrapped_phase[seed];
    visited[seed] = 1;
    queue[qtail++] = seed;

    while (qhead < qtail) {
        const int32_t idx = queue[qhead++];
        const int32_t r   = idx % rows;
        const int32_t c   = idx / rows;
        const double  v   = unwrapped_phase[idx];
        const double  pv  = wrapped_phase[idx];

        /* 4-connected neighbors */
        const int32_t nbr_r[4] = { r - 1, r + 1, r,     r     };
        const int32_t nbr_c[4] = { c,     c,     c - 1, c + 1 };

        for (int k = 0; k < 4; ++k) {
            const int32_t nr = nbr_r[k];
            const int32_t nc = nbr_c[k];
            if (nr < 0 || nr >= rows || nc < 0 || nc >= cols) continue;
            const size_t nidx = IDX(nr, nc, rows);
            if (visited[nidx] || cut[nidx]) continue;

            unwrapped_phase[nidx] = v + wrap_pi(wrapped_phase[nidx] - pv);
            visited[nidx] = 1;
            queue[qtail++] = (int32_t)nidx;
        }
    }

    /* ---- Post-process: unwrap cut pixels from any unwrapped neighbor.
     * Iterate until no progress (some cut pixels may be enclosed entirely). */
    int32_t progress = 1;
    while (progress) {
        progress = 0;
        for (int32_t c = 0; c < cols; ++c) {
            for (int32_t r = 0; r < rows; ++r) {
                const size_t idx = IDX(r, c, rows);
                if (visited[idx]) continue;

                const int32_t nbr_r[4] = { r - 1, r + 1, r,     r     };
                const int32_t nbr_c[4] = { c,     c,     c - 1, c + 1 };
                for (int k = 0; k < 4; ++k) {
                    const int32_t nr = nbr_r[k];
                    const int32_t nc = nbr_c[k];
                    if (nr < 0 || nr >= rows || nc < 0 || nc >= cols) continue;
                    const size_t nidx = IDX(nr, nc, rows);
                    if (!visited[nidx]) continue;

                    unwrapped_phase[idx] = unwrapped_phase[nidx]
                        + wrap_pi(wrapped_phase[idx] - wrapped_phase[nidx]);
                    visited[idx] = 1;
                    progress = 1;
                    break;
                }
            }
        }
    }

    /* Any pixel still NaN was completely enclosed by cuts; fall back to the
     * wrapped value so the output contains no NaNs. */
    for (size_t i = 0; i < N; ++i) {
        if (!visited[i]) {
            unwrapped_phase[i] = wrapped_phase[i];
        }
    }

    free(residues);
    free(cut);
    free(visited);
    free(queue);
    return 0;
}


/* ========================================================================
 * 2. Sine/cosine decomposition.
 *
 * Trivially data-parallel: each output element depends only on phase[i].
 * Loop index uses a signed type (int32_t) to satisfy OpenMP 2.0, which
 * requires a signed integer loop counter for `parallel for`.
 * Uses sincos() when available (a single transcendental call instead of two).
 * ======================================================================== */
int phase_sincos_decompose(const double *phase, double *cos_out, double *sin_out, int32_t numel)
{
    if (phase == NULL) {
        return -1;
    }
    if (cos_out == NULL && sin_out == NULL) {
        return -1;
    }
    if (numel <= 0) {
        return 0;
    }

    int32_t i;

    if (cos_out != NULL && sin_out != NULL) {
        /* Both outputs requested: compute sin and cos together. */
        #pragma omp parallel for
        for (i = 0; i < numel; ++i) {
#if defined(__GNUC__) && !defined(__clang__)
            double s, c;
            sincos(phase[i], &s, &c);
            sin_out[i] = s;
            cos_out[i] = c;
#else
            const double p = phase[i];
            cos_out[i] = cos(p);
            sin_out[i] = sin(p);
#endif
        }
    } else if (cos_out != NULL) {
        #pragma omp parallel for
        for (i = 0; i < numel; ++i) {
            cos_out[i] = cos(phase[i]);
        }
    } else {
        #pragma omp parallel for
        for (i = 0; i < numel; ++i) {
            sin_out[i] = sin(phase[i]);
        }
    }

    return 0;
}


/* ========================================================================
 * 3. Phase gradient / fringe-rate map.
 *
 * Centered differences in the interior; one-sided at the borders.  Each
 * finite difference is wrapped into (-pi, pi] BEFORE being scaled by the
 * sample spacing, so a 2*pi jump in the wrapped phase produces a small
 * gradient (the true local rate) instead of a huge spurious one.
 *
 * For column-major data:
 *   - "row gradient" walks down a column: neighbors are (r-1,c) and (r+1,c)
 *     -- contiguous in memory.
 *   - "col gradient" walks across a row:  neighbors are (r,c-1) and (r,c+1)
 *     -- separated by `rows` elements in memory.
 *
 * The outer loop is parallelized over columns (cheap, well-balanced, and
 * gives each thread a contiguous block of memory to write).  Loop counters
 * are signed (int32_t) for OpenMP 2.0 compatibility.
 * ======================================================================== */
int phase_gradient(const double *phase,
                   double *grad_row, double *grad_col, double *fringe_rate,
                   int32_t rows, int32_t cols)
{
    if (phase == NULL) {
        return -1;
    }
    if (grad_row == NULL && grad_col == NULL && fringe_rate == NULL) {
        return -1;
    }
    if (rows < 2 || cols < 2) {
        return -2;
    }

    /* If the caller wants only the magnitude, we still need both partial
     * derivatives.  Use small per-thread scratch when those buffers weren't
     * provided. */
    const int need_gr = (grad_row != NULL) || (fringe_rate != NULL);
    const int need_gc = (grad_col != NULL) || (fringe_rate != NULL);

    int32_t c;

    #pragma omp parallel for
    for (c = 0; c < cols; ++c) {
        int32_t r;
        for (r = 0; r < rows; ++r) {
            const size_t idx = IDX(r, c, rows);
            double gr = 0.0, gc = 0.0;

            if (need_gr) {
                if (rows == 1) {
                    gr = 0.0;
                } else if (r == 0) {
                    /* forward difference */
                    gr = wrap_pi(phase[IDX(1, c, rows)] - phase[IDX(0, c, rows)]);
                } else if (r == rows - 1) {
                    /* backward difference */
                    gr = wrap_pi(phase[IDX(rows - 1, c, rows)] -
                                 phase[IDX(rows - 2, c, rows)]);
                } else {
                    /* centered difference: 0.5 * wrap(phase[r+1] - phase[r-1]) */
                    gr = 0.5 * wrap_pi(phase[IDX(r + 1, c, rows)] -
                                       phase[IDX(r - 1, c, rows)]);
                }
            }

            if (need_gc) {
                if (cols == 1) {
                    gc = 0.0;
                } else if (c == 0) {
                    gc = wrap_pi(phase[IDX(r, 1, rows)] - phase[IDX(r, 0, rows)]);
                } else if (c == cols - 1) {
                    gc = wrap_pi(phase[IDX(r, cols - 1, rows)] -
                                 phase[IDX(r, cols - 2, rows)]);
                } else {
                    gc = 0.5 * wrap_pi(phase[IDX(r, c + 1, rows)] -
                                       phase[IDX(r, c - 1, rows)]);
                }
            }

            if (grad_row    != NULL) grad_row[idx]    = gr;
            if (grad_col    != NULL) grad_col[idx]    = gc;
            if (fringe_rate != NULL) fringe_rate[idx] = sqrt(gr * gr + gc * gc);
        }
    }

    return 0;
}


/* ========================================================================
 * 4a. Local phase coherence.
 *
 *   gamma(r,c) = | (1/N) * sum exp(j * phase(i,j)) |  in a window
 *
 * Computed as sqrt(<cos>^2 + <sin>^2).  Each pixel is independent, so this
 * is a perfect OpenMP target.  Window is clipped at borders (no padding).
 *
 * Loop counters use signed int32_t for OpenMP 2.0 compatibility.
 * ======================================================================== */
int phase_coherence(const double *phase, double *coherence,
                    int32_t rows, int32_t cols, int32_t window)
{
    if (phase == NULL || coherence == NULL) {
        return -1;
    }
    if (rows < 1 || cols < 1 || window < 1) {
        return -2;
    }

    int32_t c;

    #pragma omp parallel for
    for (c = 0; c < cols; ++c) {
        const int32_t c_lo = (c - window <  0   ) ? 0          : c - window;
        const int32_t c_hi = (c + window >= cols) ? cols - 1   : c + window;

        int32_t r;
        for (r = 0; r < rows; ++r) {
            const int32_t r_lo = (r - window <  0   ) ? 0        : r - window;
            const int32_t r_hi = (r + window >= rows) ? rows - 1 : r + window;

            double sum_c = 0.0, sum_s = 0.0;
            int32_t count = 0;
            int32_t jc;
            for (jc = c_lo; jc <= c_hi; ++jc) {
                int32_t jr;
                for (jr = r_lo; jr <= r_hi; ++jr) {
                    const double p = phase[IDX(jr, jc, rows)];
#if defined(__GNUC__) && !defined(__clang__)
                    double s, cc_;
                    sincos(p, &s, &cc_);
                    sum_s += s;
                    sum_c += cc_;
#else
                    sum_c += cos(p);
                    sum_s += sin(p);
#endif
                    ++count;
                }
            }

            const double inv_n = 1.0 / (double)count;
            const double mc = sum_c * inv_n;
            const double ms = sum_s * inv_n;
            coherence[IDX(r, c, rows)] = sqrt(mc * mc + ms * ms);
        }
    }

    return 0;
}


/* ========================================================================
 * 4b. HSV -> RGB composite for coherence-weighted phase visualization.
 *
 * Hue        = (phase + pi) / (2*pi)  in [0, 1)
 * Saturation = coherence (or 1)
 * Value      = coherence (or 1)
 *
 * Standard 6-sector HSV->RGB conversion.  Trivially data-parallel.
 * ======================================================================== */
int phase_to_rgb_hsv(const double *phase, const double *coherence,
                     double *r_out, double *g_out, double *b_out,
                     int32_t rows, int32_t cols)
{
    if (phase == NULL || r_out == NULL || g_out == NULL || b_out == NULL) {
        return -1;
    }
    if (rows < 1 || cols < 1) {
        return -2;
    }

    const int32_t total = rows * cols;
    const double inv_two_pi = 1.0 / (2.0 * M_PI);

    int32_t i;

    #pragma omp parallel for
    for (i = 0; i < total; ++i) {
        /* Hue in [0, 1) from phase */
        double h = (phase[i] + M_PI) * inv_two_pi;
        /* Robustly clamp into [0, 1) for any input range */
        h -= floor(h);

        double s, v;
        if (coherence != NULL) {
            double q = coherence[i];
            if (q < 0.0) q = 0.0;
            else if (q > 1.0) q = 1.0;
            s = q;
            v = q;
        } else {
            s = 1.0;
            v = 1.0;
        }

        /* HSV -> RGB (6-sector form) */
        const double hh = h * 6.0;
        const int    sector = (int)hh;        /* 0..5 */
        const double f  = hh - (double)sector;
        const double p  = v * (1.0 - s);
        const double q  = v * (1.0 - s * f);
        const double t  = v * (1.0 - s * (1.0 - f));

        double rr, gg, bb;
        switch (sector) {
            case 0:  rr = v; gg = t; bb = p; break;
            case 1:  rr = q; gg = v; bb = p; break;
            case 2:  rr = p; gg = v; bb = t; break;
            case 3:  rr = p; gg = q; bb = v; break;
            case 4:  rr = t; gg = p; bb = v; break;
            default: rr = v; gg = p; bb = q; break; /* sector 5 (or 6 from rounding) */
        }

        r_out[i] = rr;
        g_out[i] = gg;
        b_out[i] = bb;
    }

    return 0;
}
