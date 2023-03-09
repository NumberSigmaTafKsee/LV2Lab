/*
 * SPDX-FileCopyrightText: Hanspeter Portner <dev@open-music-kontrollers.ch>
 * SPDX-License-Identifier: Artistic-2.0
 */

#ifndef _D2TK_FRONTEND_GLFW_H
#define _D2TK_FRONTEND_GLFW_H

#include <signal.h>

#include <d2tk/base.h>
#include <d2tk/frontend.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct _d2tk_glfw_config_t d2tk_glfw_config_t;

struct _d2tk_glfw_config_t {
	d2tk_coord_t w;
	d2tk_coord_t h;
	const char *bundle_path;
	d2tk_frontend_expose_t expose;
	void *data;
};

D2TK_API d2tk_frontend_t *
d2tk_glfw_new(const d2tk_glfw_config_t *config);

#ifdef __cplusplus
}
#endif

#endif // _D2TK_FRONTEND_GLFW_H
