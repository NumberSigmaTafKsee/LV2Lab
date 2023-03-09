// Generated by gmmproc 2.64.2 -- DO NOT MODIFY!
#ifndef _GXWMM_FASTMETER_H
#define _GXWMM_FASTMETER_H


#include <glibmm/ustring.h>
#include <sigc++/sigc++.h>

/*
 * Copyright (C) 2009, 2010 Hermann Meyer, James Warden, Andreas Degert
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
 */

#include <gtkmm/drawingarea.h>
#include <gtkmm/orientable.h>
#include <gxwmm/controlparameter.h>


#ifndef DOXYGEN_SHOULD_SKIP_THIS
using GxFastMeter = struct _GxFastMeter;
using GxFastMeterClass = struct _GxFastMeterClass;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */


#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace Gxw
{ class FastMeter_Class; } // namespace Gxw
#endif //DOXYGEN_SHOULD_SKIP_THIS

namespace Gxw {


class FastMeter: public Gtk::DrawingArea, public Gtk::Orientable, public ControlParameter {
	public:
#ifndef DOXYGEN_SHOULD_SKIP_THIS
  typedef FastMeter CppObjectType;
  typedef FastMeter_Class CppClassType;
  typedef GxFastMeter BaseObjectType;
  typedef GxFastMeterClass BaseClassType;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

  FastMeter(FastMeter&& src) noexcept;
  FastMeter& operator=(FastMeter&& src) noexcept;

  // noncopyable
  FastMeter(const FastMeter&) = delete;
  FastMeter& operator=(const FastMeter&) = delete;

  ~FastMeter() noexcept override;

#ifndef DOXYGEN_SHOULD_SKIP_THIS

private:
  friend class FastMeter_Class;
  static CppClassType fastmeter_class_;

protected:
  explicit FastMeter(const Glib::ConstructParams& construct_params);
  explicit FastMeter(GxFastMeter* castitem);

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

public:

  /** Get the GType for this class, for use with the underlying GObject type system.
   */
  static GType get_type()      G_GNUC_CONST;

#ifndef DOXYGEN_SHOULD_SKIP_THIS


  static GType get_base_type() G_GNUC_CONST;
#endif

  /// Provides access to the underlying C GObject.
  GxFastMeter*       gobj()       { return reinterpret_cast<GxFastMeter*>(gobject_); }

  /// Provides access to the underlying C GObject.
  const GxFastMeter* gobj() const { return reinterpret_cast<GxFastMeter*>(gobject_); }


public:
  //C++ methods used to invoke GTK+ virtual functions:

protected:
  //GTK+ Virtual Functions (override these to change behaviour):

  //Default Signal Handlers::


private:

	public:
	FastMeter();
	
  void set(double lvl);
	
  void set_by_power(double lvl);
	
  void set_c_level(double lvl);
	
  void clear();
	
  void set_hold_count(int val);
	/** Count of cycles for which the peak value is held on display.
   *
   * Default value: 2
   *
   * @return A PropertyProxy that allows you to get or set the value of the property,
   * or receive notification when the value of the property changes.
   */
  Glib::PropertyProxy< int > property_hold() ;

/** Count of cycles for which the peak value is held on display.
   *
   * Default value: 2
   *
   * @return A PropertyProxy_ReadOnly that allows you to get the value of the property,
   * or receive notification when the value of the property changes.
   */
  Glib::PropertyProxy_ReadOnly< int > property_hold() const;

	/** Size of meter.
   *
   * Default value: 2
   *
   * @return A PropertyProxy that allows you to get or set the value of the property,
   * or receive notification when the value of the property changes.
   */
  Glib::PropertyProxy< int > property_dimen() ;

/** Size of meter.
   *
   * Default value: 2
   *
   * @return A PropertyProxy_ReadOnly that allows you to get the value of the property,
   * or receive notification when the value of the property changes.
   */
  Glib::PropertyProxy_ReadOnly< int > property_dimen() const;

	/** Meter peak falloff.
   *
   * Default value: <tt>false</tt>
   *
   * @return A PropertyProxy that allows you to get or set the value of the property,
   * or receive notification when the value of the property changes.
   */
  Glib::PropertyProxy< bool > property_falloff() ;

/** Meter peak falloff.
   *
   * Default value: <tt>false</tt>
   *
   * @return A PropertyProxy_ReadOnly that allows you to get the value of the property,
   * or receive notification when the value of the property changes.
   */
  Glib::PropertyProxy_ReadOnly< bool > property_falloff() const;

	/** Meter is showing signal power (input range: 0 .. 2).
   *
   * Default value: <tt>false</tt>
   *
   * @return A PropertyProxy that allows you to get or set the value of the property,
   * or receive notification when the value of the property changes.
   */
  Glib::PropertyProxy< bool > property_power() ;

/** Meter is showing signal power (input range: 0 .. 2).
   *
   * Default value: <tt>false</tt>
   *
   * @return A PropertyProxy_ReadOnly that allows you to get the value of the property,
   * or receive notification when the value of the property changes.
   */
  Glib::PropertyProxy_ReadOnly< bool > property_power() const;

	/** The id of the linked variable.
   *
   * Default value: ""
   *
   * @return A PropertyProxy that allows you to get or set the value of the property,
   * or receive notification when the value of the property changes.
   */
  Glib::PropertyProxy< Glib::ustring > property_var_id() ;

/** The id of the linked variable.
   *
   * Default value: ""
   *
   * @return A PropertyProxy_ReadOnly that allows you to get the value of the property,
   * or receive notification when the value of the property changes.
   */
  Glib::PropertyProxy_ReadOnly< Glib::ustring > property_var_id() const;


};

} // namespace Gxw


namespace Glib
{
  /** A Glib::wrap() method for this object.
   *
   * @param object The C instance.
   * @param take_copy False if the result should take ownership of the C instance. True if it should take a new copy or ref.
   * @result A C++ instance that wraps this C instance.
   *
   * @relates Gxw::FastMeter
   */
  Gxw::FastMeter* wrap(GxFastMeter* object, bool take_copy = false);
} //namespace Glib


#endif /* _GXWMM_FASTMETER_H */
