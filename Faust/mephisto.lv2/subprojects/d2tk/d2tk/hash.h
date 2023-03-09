/*
 * SPDX-FileCopyrightText: Hanspeter Portner <dev@open-music-kontrollers.ch>
 * SPDX-License-Identifier: Artistic-2.0
 */

#ifndef _D2TK_HASH_H
#define _D2TK_HASH_H

#include <stdint.h>
#include <unistd.h>

#include <d2tk/d2tk.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct _d2tk_hash_dict_t d2tk_hash_dict_t;

struct _d2tk_hash_dict_t {
	const void *key;
	size_t len;
};

D2TK_API uint64_t
d2tk_hash(const void *data, size_t nbytes);

D2TK_API uint64_t
d2tk_hash_foreach(const void *data, size_t nbytes, ...) __attribute__((sentinel));

D2TK_API uint64_t
d2tk_hash_dict(const d2tk_hash_dict_t *dict);

#ifdef __cplusplus
}
#endif

#endif // _D2TK_HASH_H
