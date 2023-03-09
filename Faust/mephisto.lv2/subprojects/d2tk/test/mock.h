/*
 * SPDX-FileCopyrightText: Hanspeter Portner <dev@open-music-kontrollers.ch>
 * SPDX-License-Identifier: Artistic-2.0
 */

#ifndef _D2TK_MOCK_H
#define _D2TK_MOCK_H

#include <stddef.h>
#include <unistd.h>

#include "src/core_internal.h"

#define DIM_W 640
#define DIM_H 480

#define CLIP_X 10
#define CLIP_Y 10
#define CLIP_W 100
#define CLIP_H 100

typedef void (*d2tk_check_t)(const d2tk_com_t *com, const d2tk_clip_t *clip);
typedef struct _d2tk_mock_ctx_t d2tk_mock_ctx_t;

struct _d2tk_mock_ctx_t {
	d2tk_check_t check;
};

extern const d2tk_core_driver_t d2tk_mock_driver;
extern const d2tk_core_driver_t d2tk_mock_driver_triple;
extern const d2tk_core_driver_t d2tk_mock_driver_lazy;

#endif // _D2TK_MOCK_H
