/*
 * Copyright (C) 2020 Hermann Meyer, Andreas Degert
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 * --------------------------------------------------------------------------
 */

/*---------------------------------------------------------------------
-----------------------------------------------------------------------    
                define PortIndex and plugin uri
-----------------------------------------------------------------------
----------------------------------------------------------------------*/

#include "gx_vibe.h"

/*---------------------------------------------------------------------
-----------------------------------------------------------------------    
                define controller numbers
-----------------------------------------------------------------------
----------------------------------------------------------------------*/

#define CONTROLS 8

/*---------------------------------------------------------------------
-----------------------------------------------------------------------    
                include the LV2 interface
-----------------------------------------------------------------------
----------------------------------------------------------------------*/

#include "lv2_plugin.cc"

/*---------------------------------------------------------------------
-----------------------------------------------------------------------    
                define the plugin settings
-----------------------------------------------------------------------
----------------------------------------------------------------------*/

void plugin_value_changed(X11_UI *ui, Widget_t *w, PortIndex index) {
    // do special stuff when needed
}

void plugin_set_window_size(int *w,int *h,const char * plugin_uri) {
    (*w) = 400; //initial width of main window
    (*h) = 340; //initial height of main window
}

const char* plugin_set_name() {
    return "GxVibe Mono"; //plugin name to display on UI
}

void plugin_create_controller_widgets(X11_UI *ui, const char * plugin_uri) {
    ui->widget[0] = add_my_knob(ui->widget[0], WIDTH,"Width", ui,40, 40, 80, 105);
    set_adjustment(ui->widget[0]->adj,0.5, 0.5, 0.0, 1.0, 0.01, CL_CONTINUOS);

    ui->widget[1] = add_my_knob(ui->widget[1], DEPTH,"Depth", ui,160, 40, 80, 105);
    set_adjustment(ui->widget[1]->adj,0.5, 0.5, 0.0, 1.0, 0.01, CL_CONTINUOS);

    ui->widget[2] = add_my_knob(ui->widget[2], WETDRY,"Dry/Wet", ui,280, 40, 80, 105);
    set_adjustment(ui->widget[2]->adj,0.5, 0.5, 0.0, 1.0, 0.01, CL_CONTINUOS);

    if (strcmp(GXPLUGIN_URI "#_vibe_", plugin_uri) == 0) {
        ui->widget[3] = add_my_knob(ui->widget[3], TEMPO,"Tempo", ui,30, 150, 60, 85);
        set_adjustment(ui->widget[3]->adj,4.4, 4.4, 0.1, 10.0, 0.1, CL_CONTINUOS);

        ui->widget[4] = add_my_knob(ui->widget[4], FB,"Feedback", ui,310, 150, 60, 85);
        set_adjustment(ui->widget[4]->adj,-0.6, -0.6, -1.0, 1.0, 0.01, CL_CONTINUOS);

        ui->win->label = "GxVibe"; //plugin name to display on UI
        ui->widget[5] = add_my_knob(ui->widget[5], DF,"St.df", ui,100, 150, 60, 85);
        set_adjustment(ui->widget[5]->adj,0.5, 0.5, 0.0, 1.0, 0.01, CL_CONTINUOS);

        ui->widget[6] = add_my_knob(ui->widget[6], PAN,"Pan", ui,170, 150, 60, 85);
        set_adjustment(ui->widget[6]->adj,0.0, 0.0, -1.0, 1.0, 0.01, CL_CONTINUOS);

        ui->widget[7] = add_my_knob(ui->widget[7], CROSS,"Cross", ui,240, 150, 60, 85);
        set_adjustment(ui->widget[7]->adj,0.0, 0.0, -1.0, 1.0, 0.01, CL_CONTINUOS);
    } else {
        ui->widget[3] = add_my_knob(ui->widget[3], TEMPO,"Tempo", ui,100, 150, 80, 105);
        set_adjustment(ui->widget[3]->adj,4.4, 4.4, 0.1, 10.0, 0.1, CL_CONTINUOS);

        ui->widget[4] = add_my_knob(ui->widget[4], FB,"Feedback", ui,220, 150, 80, 105);
        set_adjustment(ui->widget[4]->adj,-0.6, -0.6, -1.0, 1.0, 0.01, CL_CONTINUOS);
    }

}

void plugin_cleanup(X11_UI *ui) {
}

void plugin_port_event(LV2UI_Handle handle, uint32_t port_index,
                        uint32_t buffer_size, uint32_t format,
                        const void * buffer) {
}
