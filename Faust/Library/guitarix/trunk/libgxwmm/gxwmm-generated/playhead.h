// Generated by gmmproc 2.64.2 -- DO NOT MODIFY!
#ifndef _GXWMM_PLAYHEAD_H
#define _GXWMM_PLAYHEAD_H


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

#include <gxwmm/regler.h>
#include <gtkmm/adjustment.h>


#ifndef DOXYGEN_SHOULD_SKIP_THIS
using GxPlayHead = struct _GxPlayHead;
using GxPlayHeadClass = struct _GxPlayHeadClass;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */


#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace Gxw
{ class PlayHead_Class; } // namespace Gxw
#endif //DOXYGEN_SHOULD_SKIP_THIS

namespace Gxw {


class PlayHead: public Regler {
	public:
#ifndef DOXYGEN_SHOULD_SKIP_THIS
  typedef PlayHead CppObjectType;
  typedef PlayHead_Class CppClassType;
  typedef GxPlayHead BaseObjectType;
  typedef GxPlayHeadClass BaseClassType;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

  PlayHead(PlayHead&& src) noexcept;
  PlayHead& operator=(PlayHead&& src) noexcept;

  // noncopyable
  PlayHead(const PlayHead&) = delete;
  PlayHead& operator=(const PlayHead&) = delete;

  ~PlayHead() noexcept override;

#ifndef DOXYGEN_SHOULD_SKIP_THIS

private:
  friend class PlayHead_Class;
  static CppClassType playhead_class_;

protected:
  explicit PlayHead(const Glib::ConstructParams& construct_params);
  explicit PlayHead(GxPlayHead* castitem);

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

public:

  /** Get the GType for this class, for use with the underlying GObject type system.
   */
  static GType get_type()      G_GNUC_CONST;

#ifndef DOXYGEN_SHOULD_SKIP_THIS


  static GType get_base_type() G_GNUC_CONST;
#endif

  /// Provides access to the underlying C GObject.
  GxPlayHead*       gobj()       { return reinterpret_cast<GxPlayHead*>(gobject_); }

  /// Provides access to the underlying C GObject.
  const GxPlayHead* gobj() const { return reinterpret_cast<GxPlayHead*>(gobject_); }


public:
  //C++ methods used to invoke GTK+ virtual functions:

protected:
  //GTK+ Virtual Functions (override these to change behaviour):

  //Default Signal Handlers::


private:

	public:
	PlayHead();
	explicit PlayHead(const Glib::RefPtr<Gtk::Adjustment>& adjustment);


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
   * @relates Gxw::PlayHead
   */
  Gxw::PlayHead* wrap(GxPlayHead* object, bool take_copy = false);
} //namespace Glib


#endif /* _GXWMM_PLAYHEAD_H */
