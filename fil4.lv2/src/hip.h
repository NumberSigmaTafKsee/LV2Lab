/* fil4.lv2 - highpass
 *
 * Copyright (C) 2015 Robin Gareus <robin@gareus.org>
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

#ifndef _FIL4_HIP_H
#define _FIL4_HIP_H

#include <math.h>

typedef struct {
	float y2;
	float z1, z2;
	float a, q, g;
	float alpha, omega, q2;
	float freq, qual; // last settings
	float rate;
	bool  en;
} HighPass;

static void hip_setup (HighPass *f, float rate, float freq, float q) {
	memset (f, 0, sizeof(HighPass));
	f->rate = rate;
	f->freq = freq;

	f->qual = q;
	f->q2 = RESHP(q);
	if (f->q2 < 0.f)  f->q2 = 0.f;
	if (f->q2 > 1.6f) f->q2 = 1.6f;

	if (freq > rate / 12.f) freq = rate / 12.f;
	f->alpha = exp (-2.0 * M_PI * freq / rate);

	f->a = 1.0;
	f->q = 0.0; // start bypassed
	f->g = 0.0; // fade in
	f->en = false;
}

static bool hip_interpolate (HighPass *f, bool en, float freq, float q) {
	// called an interval of max 48 samples
	bool changed = f->en != en;
	f->en = en;

	if (freq != f->freq) {
		f->freq = freq;
		if (freq > f->rate / 12.f) {
			freq = f->rate / 12.f;
		}
		if (freq < 5.f) {
			freq = 5.f;
		}
		f->omega = freq / f->rate;
		f->alpha = exp (-2.0 * M_PI * f->omega);
		changed = true;
	}

	if (f->qual != q) {
		f->q2 = RESHP(q);
		if (f->q2 < 0.f)  f->q2 = 0.f;
		if (f->q2 > 1.6f) f->q2 = 1.6f;
		//printf("HI: %f -> %f\n", q, f->q2);
		f->qual = q;
		changed = true;
	}

	const float to = en ? f->alpha : 1.0;
	if (fabsf(to - f->a) < 1e-5) {
		f->a = to;
	} else {
		f->a += .01 * (to - f->a);
		changed = true;
	}

	const float tq = en ? f->q2 : 0;
	if (fabsf(tq - f->q) < 1e-5) {
		f->q = tq;
	} else {
		f->q += .01 * (tq - f->q);
		changed = true;
	}

	//target gain = 1 + (.5 + q) * 2 * w;
	const float tg = en ? (1.f + f->omega + 2.f * f->q * f->omega) : 1.0;
	if (fabsf(tg - f->g) < 1e-5) {
		f->g = tg;
	} else {
		f->g += .01 * (tg - f->g);
		changed = true;
	}

	if (!en) {
		f->z1 += .01 * (f->z2 - f->z1) + 1e-12;
		//f->z2 += .01 * (f->y2 - f->z2) + 1e-12;
	}

#ifndef NO_NAN_PROTECTION
	if (isnan(f->z1)) f->z1 = 0;
	if (isnan(f->z2)) f->z2 = 0;
	if (isnan(f->y2)) f->y2 = 0;
#endif
	return changed;
}

static void hip_compute (HighPass *f, uint32_t n_samples, float *buf) {
	const float a = f->a;
	const float q = f->q;
	const float g = f->g;

	const float m1 = g/a;
	const float m2 = g*q;

	if (a == 1.0 && q == 0.0 && g == 1.0) {
		// might as well save some computing
		// (all values incl state are filtered)
		return;
	}

	float z1 = f->z1;
	float z2 = f->z2;
	float y2 = f->y2;

	for (uint32_t i = 0; i < n_samples; ++i) {
		const float _z1 = z1; // remember previous input
		const float _z2 = z2; // since buf[] is processed in-place

		z1 = m1 * buf[i] - m2 * (y2 - z2); // == g * (buf[i] / a - q * (y2 - z2))
		z2 = a * (z2 + z1 - _z1);
		y2 = a * (y2 + z2 - _z2);
		buf[i] = y2;
	}

	f->y2 = y2;
	f->z1 = z1 + 1e-12;
	f->z2 = z2 + 1e-12;
}
#endif
