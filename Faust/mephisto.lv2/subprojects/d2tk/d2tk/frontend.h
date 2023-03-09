/*
 * SPDX-FileCopyrightText: Hanspeter Portner <dev@open-music-kontrollers.ch>
 * SPDX-License-Identifier: Artistic-2.0
 */

#ifndef _D2TK_FRONTEND_H
#define _D2TK_FRONTEND_H

#include <signal.h>

#include <d2tk/base.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef int (*d2tk_frontend_expose_t)(void *data, d2tk_coord_t w, d2tk_coord_t h);
typedef struct _d2tk_frontend_t d2tk_frontend_t;

D2TK_API void
d2tk_frontend_free(d2tk_frontend_t *dpugl);

D2TK_API int
d2tk_frontend_step(d2tk_frontend_t *dpugl);

D2TK_API int
d2tk_frontend_poll(d2tk_frontend_t *dpugl, double timeout);

D2TK_API int
d2tk_frontend_get_file_descriptors(d2tk_frontend_t *dpugl, int *fds, int numfds);

D2TK_API void
d2tk_frontend_run(d2tk_frontend_t *dpugl, const sig_atomic_t *done);

D2TK_API void
d2tk_frontend_redisplay(d2tk_frontend_t *dpugl);

D2TK_API int
d2tk_frontend_set_size(d2tk_frontend_t *dpugl, d2tk_coord_t w, d2tk_coord_t h);

D2TK_API int
d2tk_frontend_get_size(d2tk_frontend_t *dpugl, d2tk_coord_t *w, d2tk_coord_t *h);

D2TK_API d2tk_base_t *
d2tk_frontend_get_base(d2tk_frontend_t *dpugl);

D2TK_API int
d2tk_frontend_set_clipboard(d2tk_frontend_t *dpugl, const char *type,
	const void *buf, size_t buf_len);

D2TK_API const void *
d2tk_frontend_get_clipboard(d2tk_frontend_t *dpugl, const char **type,
	size_t *buf_len);

D2TK_API float
d2tk_frontend_get_scale();

#ifdef __cplusplus
}
#endif

#endif // _D2TK_FRONTEND_H
