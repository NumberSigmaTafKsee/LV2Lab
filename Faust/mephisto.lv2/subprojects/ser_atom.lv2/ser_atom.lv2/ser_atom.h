/*
 * SPDX-FileCopyrightText: Hanspeter Portner <dev@open-music-kontrollers.ch>
 * SPDX-License-Identifier: Artistic-2.0
 */

#ifndef _SER_ATOM_H
#define _SER_ATOM_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>

#include <lv2/atom/forge.h>

typedef void *(*ser_atom_realloc_t)(void *data, void *buf, size_t size);
typedef void  (*ser_atom_free_t)(void *data, void *buf);

typedef struct _ser_atom_t ser_atom_t;

struct _ser_atom_t {
	ser_atom_realloc_t realloc;
	ser_atom_free_t free;
	void *data;

	size_t size;
	size_t offset;
	union {
		uint8_t *buf;
		LV2_Atom *atom;
	};
};

int
ser_atom_init(ser_atom_t *ser);

int
ser_atom_funcs(ser_atom_t *ser, ser_atom_realloc_t realloc,
	ser_atom_free_t free, void *data);

int
ser_atom_reset(ser_atom_t *ser, LV2_Atom_Forge *forge);

LV2_Atom *
ser_atom_get(ser_atom_t *ser);

int
ser_atom_deinit(ser_atom_t *ser);

#ifdef __cplusplus
}
#endif

#endif //_SER_ATOM_H
