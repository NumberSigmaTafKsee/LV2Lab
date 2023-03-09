// Generated by gmmproc 2.64.2 -- DO NOT MODIFY!
#ifndef _GXWMM_TUNER_H
#define _GXWMM_TUNER_H


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


#ifndef DOXYGEN_SHOULD_SKIP_THIS
using GxTuner = struct _GxTuner;
using GxTunerClass = struct _GxTunerClass;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */


#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace Gxw
{ class Tuner_Class; } // namespace Gxw
#endif //DOXYGEN_SHOULD_SKIP_THIS

namespace Gxw {


class Tuner: public Gtk::DrawingArea {
	public:
#ifndef DOXYGEN_SHOULD_SKIP_THIS
  typedef Tuner CppObjectType;
  typedef Tuner_Class CppClassType;
  typedef GxTuner BaseObjectType;
  typedef GxTunerClass BaseClassType;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

  Tuner(Tuner&& src) noexcept;
  Tuner& operator=(Tuner&& src) noexcept;

  // noncopyable
  Tuner(const Tuner&) = delete;
  Tuner& operator=(const Tuner&) = delete;

  ~Tuner() noexcept override;

#ifndef DOXYGEN_SHOULD_SKIP_THIS

private:
  friend class Tuner_Class;
  static CppClassType tuner_class_;

protected:
  explicit Tuner(const Glib::ConstructParams& construct_params);
  explicit Tuner(GxTuner* castitem);

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

public:

  /** Get the GType for this class, for use with the underlying GObject type system.
   */
  static GType get_type()      G_GNUC_CONST;

#ifndef DOXYGEN_SHOULD_SKIP_THIS


  static GType get_base_type() G_GNUC_CONST;
#endif

  /// Provides access to the underlying C GObject.
  GxTuner*       gobj()       { return reinterpret_cast<GxTuner*>(gobject_); }

  /// Provides access to the underlying C GObject.
  const GxTuner* gobj() const { return reinterpret_cast<GxTuner*>(gobject_); }


public:
  //C++ methods used to invoke GTK+ virtual functions:

protected:
  //GTK+ Virtual Functions (override these to change behaviour):

  //Default Signal Handlers::


private:

	public:
	Tuner();
	
  void set_freq(double p1);
	
  double get_freq();
	
  void set_reference_pitch(double p1);
	
  double get_reference_pitch();
	
  void set_scale(double p1);
	
  double get_scale();
        
	/** The frequency for which tuning is displayed.
   *
   * Default value: 0
   *
   * @return A PropertyProxy that allows you to get or set the value of the property,
   * or receive notification when the value of the property changes.
   */
  Glib::PropertyProxy< double > property_freq() ;

/** The frequency for which tuning is displayed.
   *
   * Default value: 0
   *
   * @return A PropertyProxy_ReadOnly that allows you to get the value of the property,
   * or receive notification when the value of the property changes.
   */
  Glib::PropertyProxy_ReadOnly< double > property_freq() const;

	/** The frequency for which tuning is displayed.
   *
   * Default value: 440
   *
   * @return A PropertyProxy that allows you to get or set the value of the property,
   * or receive notification when the value of the property changes.
   */
  Glib::PropertyProxy< double > property_reference_pitch() ;

/** The frequency for which tuning is displayed.
   *
   * Default value: 440
   *
   * @return A PropertyProxy_ReadOnly that allows you to get the value of the property,
   * or receive notification when the value of the property changes.
   */
  Glib::PropertyProxy_ReadOnly< double > property_reference_pitch() const;

	/** scale the tuner area to make it bigger or smaller.
   *
   * Default value: 1
   *
   * @return A PropertyProxy that allows you to get or set the value of the property,
   * or receive notification when the value of the property changes.
   */
  Glib::PropertyProxy< double > property_scale() ;

/** scale the tuner area to make it bigger or smaller.
   *
   * Default value: 1
   *
   * @return A PropertyProxy_ReadOnly that allows you to get the value of the property,
   * or receive notification when the value of the property changes.
   */
  Glib::PropertyProxy_ReadOnly< double > property_scale() const;


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
   * @relates Gxw::Tuner
   */
  Gxw::Tuner* wrap(GxTuner* object, bool take_copy = false);
} //namespace Glib


#endif /* _GXWMM_TUNER_H */
