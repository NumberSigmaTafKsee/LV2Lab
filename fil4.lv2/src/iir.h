/* fil4.lv2
 *
 * Copyright (C) 2008, 2015 Robin Gareus <robin@gareus.org>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _FIL4_IIR_H
#define _FIL4_IIR_H

#include <stdint.h>
#include <string.h>
#include <math.h>

typedef struct {
	float a1, a2, b0, b1, b2;
	float y0, y1, y2;

	double rate;
	float gain, freq, q;

	float lpf;
	float f_l, f_u;
} IIRProc;

static void iir_init (IIRProc *f, double rate) {
	memset(f, 0, sizeof(IIRProc));
	f->rate = rate;
	f->gain = 1.0;
	f->freq = 1000;
	f->q    = 0.525;
	f->lpf  = 440.f / rate;
	f->f_l  = 0.0004 * rate;
	f->f_u  = 0.4700 * rate;
}

static int iir_interpolate (IIRProc *f, const float gain, float freq, float q) {
	if (q < .25f) { q = .25f; }
	if (q > 2.0f) { q = 2.0f; }
	if (freq < f->f_l) { freq = f->f_l; }
	if (freq > f->f_u) { freq = f->f_u; }

#ifndef NO_NAN_PROTECTION
	if (isnan(f->y1)) f->y1 = 0;
	if (isnan(f->y2)) f->y2 = 0;
#endif

	if (f->freq == freq && f->gain == gain && f->q == q) {
		return 0;
	}

	f->freq += f->lpf * (freq - f->freq);
	f->gain += f->lpf * (gain - f->gain);
	f->q    += f->lpf * (q    - f->q);

	if ((fabsf(f->gain - gain)) < 1e-4) {
		f->gain = gain;
	}
	if ((fabsf(f->freq - freq)) < 3e-1) {
		f->freq = freq;
	}
	if ((fabsf(f->q - q))       < 1e-3) {
		f->q = q;
	}
	return 1;
}

static void iir_calc_lowshelf (IIRProc *f) {
	const double w0 = 2. * M_PI * (f->freq / f->rate);
	const double _cosW = cos (w0);

	const double A  = sqrt (f->gain);
	const double As = sqrt (A);
	const double a  = sinf (w0) / 2 * (1 / f->q);
	const double b0 =  A *      ((A + 1) - (A - 1) * _cosW + 2 * As * a);
	const double b1 =  2 * A  * ((A - 1) - (A + 1) * _cosW);
	const double b2 =  A *      ((A + 1) - (A - 1) * _cosW - 2 * As * a);
	const double a0 = (A + 1) +  (A - 1) * _cosW + 2 * As * a;
	const double a1 = -2 *      ((A - 1) + (A + 1) * _cosW);
	const double a2 = (A + 1) +  (A - 1) * _cosW - 2 * As * a;

	f->b0 = b0 / a0;
	f->b1 = b1 / a0;
	f->b2 = b2 / a0;
	f->a1 = a1 / a0;
	f->a2 = a2 / a0;
}

static void iir_calc_highshelf (IIRProc *f) {
	const double w0 = 2. * M_PI * (f->freq / f->rate);
	const double _cosW = cos (w0);

	const double A  = sqrt (f->gain);
	const double As = sqrt (A);
	const double a  = sinf (w0) / 2 * (1 / f->q);
	const double b0 =  A *      ((A + 1) + (A - 1) * _cosW + 2 * As * a);
	const double b1 = -2 * A  * ((A - 1) + (A + 1) * _cosW);
	const double b2 =  A *      ((A + 1) + (A - 1) * _cosW - 2 * As * a);
	const double a0 = (A + 1) -  (A - 1) * _cosW + 2 * As * a;
	const double a1 =  2 *      ((A - 1) - (A + 1) * _cosW);
	const double a2 = (A + 1) -  (A - 1) * _cosW - 2 * As * a;

	f->b0 = b0 / a0;
	f->b1 = b1 / a0;
	f->b2 = b2 / a0;
	f->a1 = a1 / a0;
	f->a2 = a2 / a0;
}

static void iir_compute (IIRProc *f, uint32_t n_samples, float *buf) {
	// this depends on prior processors adding denormal protection
	for (uint32_t i = 0; i < n_samples; ++i) {
		const float xn = buf[i];
		const float y = f->b0 * xn + f->y1;
		f->y1         = f->b1 * xn - f->a1 * y + f->y2;
		f->y2         = f->b2 * xn - f->a2 * y;
		buf[i] = y;
	}
}
#endif
