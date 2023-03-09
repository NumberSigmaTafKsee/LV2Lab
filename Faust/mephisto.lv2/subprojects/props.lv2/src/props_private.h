/*
 * SPDX-FileCopyrightText: Hanspeter Portner <dev@open-music-kontrollers.ch>
 * SPDX-License-Identifier: Artistic-2.0
 */

#include <props.lv2/props.h>

#ifndef _LV2_PROPS_PRIVATE_H_
#define _LV2_PROPS_PRIVATE_H_

// enumerations
typedef enum _props_state_t {
	PROP_STATE_NONE    = 0,
	PROP_STATE_LOCK    = 1,
	PROP_STATE_RESTORE = 2
} props_state_t;

props_impl_t *
_props_impl_get(props_t *props, LV2_URID property);

#endif // _LV2_PROPS_PRIVATE_H_
