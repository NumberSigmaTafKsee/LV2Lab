/*
 * SPDX-FileCopyrightText: Hanspeter Portner <dev@open-music-kontrollers.ch>
 * SPDX-License-Identifier: Artistic-2.0
 */

#include <ser_atom.lv2/ser_atom.h>

static LV2_Atom_Forge_Ref
_ser_atom_sink(LV2_Atom_Forge_Sink_Handle handle, const void *buf,
	uint32_t size)
{
	ser_atom_t *ser = handle;
	const size_t needed = ser->offset + size;

	while(needed > ser->size)
	{
		const size_t augmented = ser->size
			? ser->size << 1
			: 1024;
		uint8_t *grown = ser->realloc(ser->data, ser->buf, augmented);
		if(!grown) // out-of-memory
		{
			return 0;
		}

		ser->buf = grown;
		ser->size = augmented;
	}

	const LV2_Atom_Forge_Ref ref = ser->offset + 1;
	memcpy(&ser->buf[ser->offset], buf, size);
	ser->offset += size;

	return ref;
}

static LV2_Atom *
_ser_atom_deref(LV2_Atom_Forge_Sink_Handle handle, LV2_Atom_Forge_Ref ref)
{
	ser_atom_t *ser = handle;

	if(!ref) // invalid reference
	{
		return NULL;
	}

	const size_t offset = ref - 1;
	return (LV2_Atom *)&ser->buf[offset];
}

static void *
_ser_atom_realloc(void *data, void *buf, size_t size)
{
	(void)data;

	return realloc(buf, size);
}

static void
_ser_atom_free(void *data, void *buf)
{
	(void)data;

	free(buf);
}

int
ser_atom_funcs(ser_atom_t *ser, ser_atom_realloc_t realloc,
	ser_atom_free_t free, void *data)
{
	if(!ser || !realloc || !free || ser_atom_deinit(ser))
	{
		return -1;
	}

	ser->realloc = realloc;
	ser->free = free;
	ser->data = data;

	return 0;
}

int
ser_atom_init(ser_atom_t *ser)
{
	if(!ser)
	{
		return -1;
	}

	ser->size = 0;
	ser->offset = 0;
	ser->buf = NULL;

	return ser_atom_funcs(ser, _ser_atom_realloc, _ser_atom_free, NULL);
}

int
ser_atom_reset(ser_atom_t *ser, LV2_Atom_Forge *forge)
{
	if(!ser || !forge)
	{
		return -1;
	}

	lv2_atom_forge_set_sink(forge, _ser_atom_sink, _ser_atom_deref, ser);

	ser->offset = 0;

	return 0;
}

LV2_Atom *
ser_atom_get(ser_atom_t *ser)
{
	if(!ser)
	{
		return NULL;
	}

	return ser->atom;
}

int
ser_atom_deinit(ser_atom_t *ser)
{
	if(!ser)
	{
		return -1;
	}

	if(ser->buf)
	{
		ser->free(ser->data, ser->buf);
	}

	ser->size = 0;
	ser->offset = 0;
	ser->buf = NULL;

	return 0;
}
