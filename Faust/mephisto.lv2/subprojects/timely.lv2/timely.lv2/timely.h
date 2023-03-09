/*
 * SPDX-FileCopyrightText: Hanspeter Portner <dev@open-music-kontrollers.ch>
 * SPDX-License-Identifier: Artistic-2.0
 */

#ifndef _LV2_TIMELY_H_
#define _LV2_TIMELY_H_

#include <lv2/core/lv2.h>
#include <lv2/urid/urid.h>
#include <lv2/atom/atom.h>
#include <lv2/atom/forge.h>
#include <lv2/time/time.h>

typedef struct _timely_t timely_t;
typedef void (*timely_cb_t)(timely_t *timely, int64_t frames, LV2_URID type,
	void *data);

typedef enum _timely_mask_t {
	TIMELY_MASK_BAR_BEAT					= (1 << 0),
	TIMELY_MASK_BAR								= (1 << 1),
	TIMELY_MASK_BEAT_UNIT					= (1 << 2),
	TIMELY_MASK_BEATS_PER_BAR			= (1 << 3),
	TIMELY_MASK_BEATS_PER_MINUTE	= (1 << 4),
	TIMELY_MASK_FRAME							= (1 << 5),
	TIMELY_MASK_FRAMES_PER_SECOND	= (1 << 6),
	TIMELY_MASK_SPEED							= (1 << 7),
	TIMELY_MASK_BAR_BEAT_WHOLE		= (1 << 8),
	TIMELY_MASK_BAR_WHOLE					= (1 << 9)
} timely_mask_t;

struct _timely_t {
	struct {
		LV2_URID atom_object;
		LV2_URID atom_blank;
		LV2_URID atom_resource;

		LV2_URID time_position;
		LV2_URID time_barBeat;
		LV2_URID time_bar;
		LV2_URID time_beatUnit;
		LV2_URID time_beatsPerBar;
		LV2_URID time_beatsPerMinute;
		LV2_URID time_frame;
		LV2_URID time_framesPerSecond;
		LV2_URID time_speed;
	} urid;

	struct {
		float bar_beat;
		int64_t bar;

		int32_t beat_unit;
		float beats_per_bar;
		float beats_per_minute;

		int64_t frame;
		float frames_per_second;

		float speed;
	} pos;

	float multiplier;

	double frames_per_beat;
	double frames_per_bar;

	struct {
		double beat;
		double bar;
	} offset;

	bool first;
	timely_mask_t mask;
	timely_cb_t cb;
	void *data;
};

#define TIMELY_URI_BAR_BEAT(timely)						((timely)->urid.time_barBeat)
#define TIMELY_URI_BAR(timely)								((timely)->urid.time_bar)
#define TIMELY_URI_BEAT_UNIT(timely)					((timely)->urid.time_beatUnit)
#define TIMELY_URI_BEATS_PER_BAR(timely)			((timely)->urid.time_beatsPerBar)
#define TIMELY_URI_BEATS_PER_MINUTE(timely)		((timely)->urid.time_beatsPerMinute)
#define TIMELY_URI_FRAME(timely)							((timely)->urid.time_frame)
#define TIMELY_URI_FRAMES_PER_SECOND(timely)	((timely)->urid.time_framesPerSecond)
#define TIMELY_URI_SPEED(timely)							((timely)->urid.time_speed)

#define TIMELY_BAR_BEAT_RAW(timely)						((timely)->pos.bar_beat)
#define TIMELY_BAR_BEAT(timely)								(floor((timely)->pos.bar_beat) \
	+ (timely)->offset.beat / (timely)->frames_per_beat)
#define TIMELY_BAR(timely)										((timely)->pos.bar)
#define TIMELY_BEAT_UNIT(timely)							((timely)->pos.beat_unit)
#define TIMELY_BEATS_PER_BAR(timely)					((timely)->pos.beats_per_bar)
#define TIMELY_BEATS_PER_MINUTE(timely)				((timely)->pos.beats_per_minute)
#define TIMELY_FRAME(timely)									((timely)->pos.frame)
#define TIMELY_FRAMES_PER_SECOND(timely)			((timely)->pos.frames_per_second)
#define TIMELY_SPEED(timely)									((timely)->pos.speed)

#define TIMELY_FRAMES_PER_BEAT(timely)				((timely)->frames_per_beat)
#define TIMELY_FRAMES_PER_BAR(timely)					((timely)->frames_per_bar)

void
timely_init(timely_t *timely, LV2_URID_Map *map, double rate,
	timely_mask_t mask, timely_cb_t cb, void *data);

void
timely_set_multiplier(timely_t *timely, float multiplier);

int
timely_advance_body(timely_t *timely, uint32_t size, uint32_t type,
	const LV2_Atom_Object_Body *body, uint32_t from, uint32_t to);

int
timely_advance(timely_t *timely, const LV2_Atom_Object *obj,
	uint32_t from, uint32_t to);

#endif // _LV2_TIMELY_H_
