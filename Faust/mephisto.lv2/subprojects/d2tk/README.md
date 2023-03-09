<!--
  -- (SPDX-FileCopyrightText: Hanspeter Portner <dev@open-music-kontrollers.ch>
  -- (SPDX-License-Identifier: CC0-1.0
  -->
# d2tk

## Data Driven Tool Kit

A performant, dyamic, immediate-mode GUI tool kit in C which partially renders
on-change only by massively hashing-and-cashing of vector drawing instructions
and on-demand rendered sprites.

### Build / test

	git clone https://git.open-music-kontrollers.ch/lad/d2tk
	cd d2tk 
	meson build
	cd build
	ninja -j4
	ninja test

### Dependencies

#### Optional

* [REUSE](https://git.fsfe.org/reuse/tool) (tool for compliance with the REUSE recommendations)

#### Pugl/NanoVG backend

	./d2tk.nanovg

#### Pugl/Cairo backend

	./d2tk.cairo

#### FBdev/Cairo backend

	./d2tk.fbdev

### Screenshots

![Screenshot 1](/screenshots/screenshot_1.png)

![Screenshot 2](/screenshots/screenshot_2.png)

![Screenshot 3](/screenshots/screenshot_3.png)

![Screenshot 4](/screenshots/screenshot_4.png)

![Screenshot 5](/screenshots/screenshot_5.png)

![Screenshot 6](/screenshots/screenshot_6.png)

![Screenshot 7](/screenshots/screenshot_7.png)

![Screenshot 8](/screenshots/screenshot_8.png)

### License

SPDX-FileCopyrightText: Hanspeter Portner <dev@open-music-kontrollers.ch>
SPDX-License-Identifier: Artistic-2.0
