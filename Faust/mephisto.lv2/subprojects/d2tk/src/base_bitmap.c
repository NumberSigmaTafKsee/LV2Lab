/*
 * SPDX-FileCopyrightText: Hanspeter Portner <dev@open-music-kontrollers.ch>
 * SPDX-License-Identifier: Artistic-2.0
 */

#include "base_internal.h"

D2TK_API void
d2tk_base_bitmap(d2tk_base_t *base, uint32_t w, uint32_t h, uint32_t stride,
	const uint32_t *argb, uint64_t rev, const d2tk_rect_t *rect,
	d2tk_align_t align)
{
	const d2tk_hash_dict_t dict [] = {
		{ rect, sizeof(d2tk_rect_t) },
		{ &w, sizeof(uint32_t) },
		{ &h, sizeof(uint32_t) },
		{ &stride, sizeof(uint32_t) },
		{ &rev, sizeof(uint64_t) },
		{ NULL, 0 }
	};
	const uint64_t hash = d2tk_hash_dict(dict);

	d2tk_core_t *core = base->core;;

	D2TK_CORE_WIDGET(core, hash, widget)
	{
		const size_t ref = d2tk_core_bbox_push(core, true, rect);

		d2tk_core_bitmap(core, rect, w, h, stride, argb, rev, align);

		d2tk_core_bbox_pop(core, ref);
	}
}
