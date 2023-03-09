/*
 * SPDX-FileCopyrightText: Hanspeter Portner <dev@open-music-kontrollers.ch>
 * SPDX-License-Identifier: Artistic-2.0
 */

#ifndef _LV2_PROPS_H_
#define _LV2_PROPS_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>
#include <stdatomic.h>
#include <stdio.h>

#include <lv2/core/lv2.h>
#include <lv2/urid/urid.h>
#include <lv2/atom/atom.h>
#include <lv2/atom/forge.h>
#include <lv2/patch/patch.h>
#include <lv2/state/state.h>

// structures
typedef struct _props_def_t props_def_t;
typedef struct _props_impl_t props_impl_t;
typedef struct _props_dyn_t props_dyn_t;
typedef struct _props_t props_t;

typedef enum _props_dyn_ev_t {
	PROPS_DYN_EV_ADD,
	PROPS_DYN_EV_REM,
	PROPS_DYN_EV_SET
} props_dyn_ev_t;

// function callbacks
typedef void (*props_event_cb_t)(
	void *data,
	int64_t frames,
	props_impl_t *impl);

typedef void (*props_dyn_prop_cb_t)(
	void *data,
	props_dyn_ev_t ev,
	LV2_URID subj,
	LV2_URID prop,
	const LV2_Atom *body);

struct _props_def_t {
	const char *property;
	const char *type;
	const char *access;
	size_t offset;
	bool hidden;

	uint32_t max_size;
	props_event_cb_t event_cb;
};

struct _props_impl_t {
	LV2_URID property;
	LV2_URID type;
	LV2_URID access;

	struct {
		uint32_t size;
		void *body;
	} value;
	struct {
		uint32_t size;
		void *body;
	} stash;

	const props_def_t *def;

	atomic_int state;
	bool stashing;
};

struct _props_dyn_t {
	props_dyn_prop_cb_t prop;
};

struct _props_t {
	struct {
		LV2_URID subject;

		LV2_URID patch_get;
		LV2_URID patch_set;
		LV2_URID patch_put;
		LV2_URID patch_patch;
		LV2_URID patch_wildcard;
		LV2_URID patch_add;
		LV2_URID patch_remove;
		LV2_URID patch_subject;
		LV2_URID patch_body;
		LV2_URID patch_property;
		LV2_URID patch_value;
		LV2_URID patch_writable;
		LV2_URID patch_readable;
		LV2_URID patch_sequence;
		LV2_URID patch_error;
		LV2_URID patch_ack;

		LV2_URID atom_int;
		LV2_URID atom_long;
		LV2_URID atom_float;
		LV2_URID atom_double;
		LV2_URID atom_bool;
		LV2_URID atom_urid;
		LV2_URID atom_path;
		LV2_URID atom_literal;
		LV2_URID atom_vector;
		LV2_URID atom_object;
		LV2_URID atom_sequence;

		LV2_URID state_StateChanged;
	} urid;

	void *data;

	bool stashing;
	atomic_bool restoring;

	uint32_t max_size;

	const props_dyn_t *dyn;

	unsigned nimpls;
	props_impl_t impls [1];
};

#define PROPS_T(PROPS, MAX_NIMPLS) \
	props_t (PROPS); \
	props_impl_t _impls [MAX_NIMPLS]

// rt-safe
int
props_init(props_t *props, const char *subject,
	const props_def_t *defs, int nimpls,
	void *value_base, void *stash_base,
	LV2_URID_Map *map, void *data);

// rt-safe
void
props_dyn(props_t *props, const props_dyn_t *dyn);

// rt-safe
void
props_idle(props_t *props, LV2_Atom_Forge *forge, uint32_t frames,
	LV2_Atom_Forge_Ref *ref);

// rt-safe
int
props_advance(props_t *props, LV2_Atom_Forge *forge, uint32_t frames,
	const LV2_Atom_Object *obj, LV2_Atom_Forge_Ref *ref);

// rt-safe
void
props_set(props_t *props, LV2_Atom_Forge *forge, uint32_t frames,
	LV2_URID property, LV2_Atom_Forge_Ref *ref);

// rt-safe
void
props_get(props_t *props, LV2_Atom_Forge *forge, uint32_t frames,
	LV2_URID property, LV2_Atom_Forge_Ref *ref);

// rt-safe
void
props_stash(props_t *props, LV2_URID property);

// rt-safe
LV2_URID
props_map(props_t *props, const char *property);

// rt-safe
const char *
props_unmap(props_t *props, LV2_URID property);

void
props_update(props_t *props, LV2_URID property, LV2_URID type, uint32_t size,
	const void *body);

// non-rt
LV2_State_Status
props_save(props_t *props, LV2_State_Store_Function store,
	LV2_State_Handle state, uint32_t flags, const LV2_Feature *const *features);

// non-rt
LV2_State_Status
props_restore(props_t *props, LV2_State_Retrieve_Function retrieve,
	LV2_State_Handle state, uint32_t flags, const LV2_Feature *const *features);

#ifdef __cplusplus
}
#endif

#endif // _LV2_PROPS_H_
