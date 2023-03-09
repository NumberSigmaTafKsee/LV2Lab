// SPDX-FileCopyrightText: Hanspeter Portner <dev@open-music-kontrollers.ch>
// SPDX-License-Identifier: CC0-1.0

import("stdfaust.lib");

process = _ <: attach(_,abs : ba.linear2db : hbargraph("Level [0]", -60, 0));

// vim: set syntax=faust:
