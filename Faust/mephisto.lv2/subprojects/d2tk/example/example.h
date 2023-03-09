/*
 * SPDX-FileCopyrightText: Hanspeter Portner <dev@open-music-kontrollers.ch>
 * SPDX-License-Identifier: Artistic-2.0
 */

#ifndef _D2TK_EXAMPLE_H
#define _D2TK_EXAMPLE_H

#include <d2tk/base.h>
#include <d2tk/frontend.h>

#ifdef __cplusplus
extern "C" {
#endif

D2TK_API int
d2tk_example_init(void);

D2TK_API void
d2tk_example_deinit();

D2TK_API void
d2tk_example_run(d2tk_frontend_t *frontend, d2tk_base_t *base,
	d2tk_coord_t w, d2tk_coord_t h);

#ifdef __cplusplus
}
#endif

#endif // _D2TK_EXAMPLE_H
