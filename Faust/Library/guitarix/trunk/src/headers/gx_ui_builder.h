/*
 * Copyright (C) 2009, 2010 Hermann Meyer, James Warden, Andreas Degert
 * Copyright (C) 2011 Pete Shorthose
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

/* ------- This is the GUI namespace ------- */

#pragma once

#ifndef SRC_HEADERS_GX_UI_BUILDER_H_
#define SRC_HEADERS_GX_UI_BUILDER_H_

#include <gxwmm/controlparameter.h>
#include <gxwmm/regler.h>

class PluginUI;
class PluginDict;

namespace gx_gui {

gx_engine::Parameter *check_get_parameter(
    gx_engine::GxMachineBase& machine, const std::string& id, Gtk::Widget *w);

gx_engine::FloatParameter *check_get_float_parameter(
    gx_engine::GxMachineBase& machine, const std::string& id, Gtk::Widget *w);

bool ui_connect_switch(
    gx_engine::GxMachineBase& machine, Gtk::ToggleButton *b,
    const string& id, sigc::signal<void(bool)> *out_ctr, bool disable);

bool ui_connect(gx_engine::GxMachineBase& machine, Gtk::Widget *w, const std::string& id,
                sigc::signal<void(bool)> *out_ctr);

/**************************************************************//**
 * class GxBuilder
 *
 * Attempts to correct the mismatch between GtkBuilder and
 * Gtk::Builder wrt. reference counting.
 *
 * Use get_toplevel or get_toplevel_derived to become the owner of a
 * toplevel widget of the loaded builder widget tree (specifying
 * object id's for the loader means that only part of the defined
 * widget tree is loaded). These pointers must be delete'd to destroy
 * the widget (and its child widgets). If you loaded a window
 * (GtkWindow or derived) you must use one of the get_toplevel
 * functions and delete the instance to avoid a memory leak.
 *
 * Use find_widget for getting a pointer to a widget (but you don't
 * become owner of that widget, don't use delete). If you want to add
 * the loaded widget tree to a container, use this function (the
 * container will ref the widget and manage its lifetime).
 *
 * template function code mostly copied from Gtk::Builder, look
 * there for comments.
 */

class GxBuilder: public Gtk::Builder {
private:
    static bool show_tooltips;
    // only implemented for base class, make inaccessible
    static Glib::RefPtr<GxBuilder> create_from_file(const std::string& filename, const char* object_id);
    static Glib::RefPtr<GxBuilder> create_from_file(const std::string& filename, const Glib::ustring& object_id);
    static Glib::RefPtr<GxBuilder> create_from_file(const std::string& filename, const Glib::StringArrayHandle& object_ids);
    static Glib::RefPtr<GxBuilder> create_from_string(const Glib::ustring& buffer, const char* object_id);
    static Glib::RefPtr<GxBuilder> create_from_string(const Glib::ustring& buffer, const Glib::ustring& object_id);
    static Glib::RefPtr<GxBuilder> create_from_string(const Glib::ustring& buffer, const Glib::StringArrayHandle& object_ids);
    GObject* get_cobject(const Glib::ustring& name);
public:
    static inline Glib::RefPtr<GxBuilder> create() { return Glib::RefPtr<GxBuilder>(new GxBuilder()); }

    static Glib::RefPtr<GxBuilder> create_from_file(
	const std::string& filename, gx_engine::GxMachineBase* pmach = 0, const char* object_id = 0, sigc::signal<void(bool)> *out_ctr = 0);

    static Glib::RefPtr<GxBuilder> create_from_file(
	const std::string& filename, gx_engine::GxMachineBase* pmach, const Glib::StringArrayHandle& object_ids, sigc::signal<void(bool)> *out_ctr = 0);

    static Glib::RefPtr<GxBuilder> create_from_string(
	const Glib::ustring& buffer, gx_engine::GxMachineBase* pmach = 0, const char* object_id = 0, sigc::signal<void(bool)> *out_ctr = 0);

    static Glib::RefPtr<GxBuilder> create_from_string(
	const Glib::ustring& buffer, gx_engine::GxMachineBase* pmach, const Glib::StringArrayHandle& object_ids, sigc::signal<void(bool)> *out_ctr = 0);

    static bool get_show_tooltips() { return show_tooltips; }
    static void set_show_tooltips(bool v) { show_tooltips = v; }
    static void connect_gx_tooltip_handler(GtkWidget *widget);
    static void set_tooltip_text_connect_handler(GtkWidget *widget, const char *text);
    static void set_tooltip_text_connect_handler(Gtk::Widget& widget, const char *text) {
	set_tooltip_text_connect_handler(widget.gobj(), text); }
    static void set_tooltip_text_connect_handler(Gtk::Widget& widget, const Glib::ustring& text) {
	set_tooltip_text_connect_handler(widget.gobj(), text.c_str()); }

    void fixup_controlparameters(gx_engine::GxMachineBase& machine, sigc::signal<void(bool)> *out_ctr);

    template <class T_Widget> inline
    void find_widget(const Glib::ustring& name, T_Widget*& widget) {
	widget = dynamic_cast<T_Widget*>(get_widget_checked(name, T_Widget::get_base_type()));
	assert(widget);
    }
    // pointer set by find_object only have the lifetime of the underlying GxBuilder
    // if you need a longer lifetime, use get_object()!
    template <class T_Object> inline
    void find_object(const Glib::ustring& name, T_Object*& object) {
	object = dynamic_cast<T_Object*>(get_object(name).get());
	assert(object);
    }

    template <class T_Widget, class F> inline
    void find_widget_derived(const Glib::ustring& name, T_Widget*& widget, F f) {
	widget = 0;
	typedef typename T_Widget::BaseObjectType cwidget_type;
	cwidget_type* pCWidget = (cwidget_type*)get_cobject(name);
	if(!pCWidget) {
	    return;
	}
	Glib::ObjectBase* pObjectBase = ObjectBase::_get_current_wrapper((GObject*)pCWidget);
	if (pObjectBase) {
	    widget = dynamic_cast<T_Widget*>( Glib::wrap((GtkWidget*)pCWidget) );
	    if (!widget) {
		g_critical("GxBuilder::get_widget_derived(): dynamic_cast<> failed. An existing C++ instance, of a different type, seems to exist.");      
	    }
	} else {
	    widget = f(pCWidget);
	}
    }

    bool has_object(const Glib::ustring& name) {
	return gtk_builder_get_object(gobj(), name.c_str()) != 0;
    }

    template <class T_Widget> inline
    void get_toplevel(const Glib::ustring& name, T_Widget*& widget) {
	GType type = T_Widget::get_base_type();
	widget = dynamic_cast<T_Widget*>(get_widget_checked(name, type));
	assert(widget);
	assert(!widget->get_parent());
    }

    Gtk::Window *get_first_window();

    template <class T_Widget, class F> inline
    void get_toplevel_derived(const Glib::ustring& name, T_Widget*& widget, F f) {
	widget = 0;
	typedef typename T_Widget::BaseObjectType cwidget_type;
	cwidget_type* pCWidget = (cwidget_type*)get_cobject(name);
	if(!pCWidget) {
	    return;
	}
	if (!g_type_is_a(G_OBJECT_TYPE(pCWidget), GTK_TYPE_WINDOW)) {
	    g_object_ref(pCWidget);
	}
	Glib::ObjectBase* pObjectBase = ObjectBase::_get_current_wrapper((GObject*)pCWidget);
	if (pObjectBase) {
	    widget = dynamic_cast<T_Widget*>( Glib::wrap((GtkWidget*)pCWidget) );
	    if (!widget) {
		g_critical("GxBuilder::get_widget_derived(): dynamic_cast<> failed. An existing C++ instance, of a different type, seems to exist.");      
	    }
	} else {
	    widget = f(pCWidget);
	}
	assert(!widget->get_parent());
    }
};

class StackBoxBuilder;

class UiBuilderImpl: public gx_engine::UiBuilderBase {
protected:
    PluginDict& plugin_dict;
    std::vector<PluginUI*> *pluginlist;
    static StackBoxBuilder *intf;
    static void openTabBox_(const char* label);
    static void openVerticalBox_(const char* label);
    static void openVerticalBox1_(const char* label);
    static void openVerticalBox2_(const char* label);
    static void openHorizontalBox_(const char* label);
    //static void openHorizontalBox_(const char* label, int spacing);
    static void openHorizontalhideBox_(const char* label);
    static void openHorizontalTableBox_(const char* label);
    static void openFrameBox_(const char *label);
    static void openFlipLabelBox_(const char* label);
    static void openpaintampBox_(const char* label);
    static void insertSpacer_();
    static void set_next_flags_(int flags);
    static void create_mid_rackknob_(const char *id, const char *label);
    static void create_small_rackknob_(const char *id, const char *label);
    static void create_small_rackknobr_(const char *id, const char *label);
    static void create_big_rackknob_(const char *id, const char *label);
    static void create_master_slider_(const char *id, const char *label);
    static void create_feedback_slider_(const char *id, const char *label);
    static void create_selector_no_caption_(const char *id);
    static void create_selector_(const char *id, const char *label);
    static void create_simple_meter_(const char *id);
    static void create_simple_c_meter_(const char *id, const char *idl, const char *label);
    static void create_spin_value_(const char *id, const char *label);
    static void create_switch_no_caption_(const char *sw_type,const char * id);
    static void create_feedback_switch_(const char *sw_type,const char * id);
    static void create_fload_switch_(const char *sw_type,const char * id,const char * idf);
    static void create_switch_(const char *sw_type,const char * id, const char *label);
    static void create_wheel_(const char * id, const char *label);
    static void create_port_display_(const char *id, const char *label);
    static void create_p_display_(const char *id, const char *idl, const char *idh);
    static void create_simple_spin_value_(const char *id);
    static void create_eq_rackslider_no_caption_(const char *id);
    static void closeBox_();
    static void load_glade_(const char *data);
    static void load_glade_file_(const char *fname);
    virtual bool load(gx_engine::Plugin *p);
public:
    UiBuilderImpl(PluginDict *i, StackBoxBuilder *b, std::vector<PluginUI*> *pl, PluginUI *pluginui);
    ~UiBuilderImpl();
    bool load_unit(PluginDef *pl);
    friend class gx_engine::GxMachineRemote;
};

GtkWidget *load_toplevel(GtkBuilder *builder, const char* filename, const char* windowname);

} /* end of gx_gui namespace */
#endif  // SRC_HEADERS_GX_UI_BUILDER_H_
