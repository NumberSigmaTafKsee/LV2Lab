/*
 * SPDX-FileCopyrightText: Hanspeter Portner <dev@open-music-kontrollers.ch>
 * SPDX-License-Identifier: Artistic-2.0
 */

#include "base_internal.h"

D2TK_API void
d2tk_base_custom(d2tk_base_t *base, uint64_t dhash, const void *data,
	const d2tk_rect_t *rect, d2tk_core_custom_t custom)
{
	const d2tk_hash_dict_t dict [] = {
		{ rect, sizeof(d2tk_rect_t) } ,
		{ &dhash, sizeof(uint64_t)},
		{ NULL, 0 }
	};
	const uint64_t hash = d2tk_hash_dict(dict);

	d2tk_core_t *core = base->core;;

	D2TK_CORE_WIDGET(core, hash, widget)
	{
		const size_t ref = d2tk_core_bbox_push(core, true, rect);

		d2tk_core_custom(core, rect, dhash, data, custom);

		d2tk_core_bbox_pop(core, ref);
	}
}
