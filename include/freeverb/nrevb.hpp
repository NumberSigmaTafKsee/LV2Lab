/**
 *  CCRMA NRev v2 (nrevb)
 *
 *  Copyright (C) 2006-2018 Teru Kamogashira
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#ifndef _FV3_NREVB_HPP
#define _FV3_NREVB_HPP

#include "freeverb/nrev.hpp"

#define FV3_NREVB_NUM_COMB_2    12
#define FV3_NREVB_NUM_ALLPASS_2 3
#define FV3_NREVB_SCALE_WET    0.4

namespace fv3
{

#define _fv3_float_t float
#define _FV3_(name) name ## _f
#include "freeverb/nrevb_t.hpp"
#undef _FV3_
#undef _fv3_float_t

#define _fv3_float_t double
#define _FV3_(name) name ## _
#include "freeverb/nrevb_t.hpp"
#undef _FV3_
#undef _fv3_float_t

#define _fv3_float_t long double
#define _FV3_(name) name ## _l
#include "freeverb/nrevb_t.hpp"
#undef _FV3_
#undef _fv3_float_t

};

#endif