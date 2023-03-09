// Generated by gmmproc 2.64.2 -- DO NOT MODIFY!


#include <glibmm.h>

#include <gxwmm/racktuner.h>
#include <gxwmm/private/racktuner_p.h>


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

#include <gxw/GxRackTuner.h>

namespace
{


static const Glib::SignalProxyInfo RackTuner_signal_frequency_poll_info =
{
  "frequency-poll",
  (GCallback) &Glib::SignalProxyNormal::slot0_void_callback,
  (GCallback) &Glib::SignalProxyNormal::slot0_void_callback
};


static void RackTuner_signal_poll_status_changed_callback(GxRackTuner* self, gboolean p0,void* data)
{
  using namespace Gxw;
  using SlotType = sigc::slot< void,bool >;

  auto obj = dynamic_cast<RackTuner*>(Glib::ObjectBase::_get_current_wrapper((GObject*) self));
  // Do not try to call a signal on a disassociated wrapper.
  if(obj)
  {
    try
    {
      if(const auto slot = Glib::SignalProxyNormal::data_to_slot(data))
        (*static_cast<SlotType*>(slot))(p0
);
    }
    catch(...)
    {
       Glib::exception_handlers_invoke();
    }
  }
}

static const Glib::SignalProxyInfo RackTuner_signal_poll_status_changed_info =
{
  "poll-status-changed",
  (GCallback) &RackTuner_signal_poll_status_changed_callback,
  (GCallback) &RackTuner_signal_poll_status_changed_callback
};


} // anonymous namespace


namespace Glib
{

Gxw::RackTuner* wrap(GxRackTuner* object, bool take_copy)
{
  return dynamic_cast<Gxw::RackTuner *> (Glib::wrap_auto ((GObject*)(object), take_copy));
}

} /* namespace Glib */

namespace Gxw
{


/* The *_Class implementation: */

const Glib::Class& RackTuner_Class::init()
{
  if(!gtype_) // create the GType if necessary
  {
    // Glib::Class has to know the class init function to clone custom types.
    class_init_func_ = &RackTuner_Class::class_init_function;

    // This is actually just optimized away, apparently with no harm.
    // Make sure that the parent type has been created.
    //CppClassParent::CppObjectType::get_type();

    // Create the wrapper type, with the same class/instance size as the base type.
    register_derived_type(gx_rack_tuner_get_type());

    // Add derived versions of interfaces, if the C type implements any interfaces:

  }

  return *this;
}


void RackTuner_Class::class_init_function(void* g_class, void* class_data)
{
  const auto klass = static_cast<BaseClassType*>(g_class);
  CppClassParent::class_init_function(klass, class_data);


  klass->frequency_poll = &frequency_poll_callback;
  klass->poll_status_changed = &poll_status_changed_callback;
}


void RackTuner_Class::frequency_poll_callback(GxRackTuner* self)
{
  const auto obj_base = static_cast<Glib::ObjectBase*>(
      Glib::ObjectBase::_get_current_wrapper((GObject*)self));

  // Non-gtkmmproc-generated custom classes implicitly call the default
  // Glib::ObjectBase constructor, which sets is_derived_. But gtkmmproc-
  // generated classes can use this optimisation, which avoids the unnecessary
  // parameter conversions if there is no possibility of the virtual function
  // being overridden:
  if(obj_base && obj_base->is_derived_())
  {
    const auto obj = dynamic_cast<CppObjectType* const>(obj_base);
    if(obj) // This can be NULL during destruction.
    {
      try // Trap C++ exceptions which would normally be lost because this is a C callback.
      {
        // Call the virtual member method, which derived classes might override.
        obj->on_frequency_poll();
        return;
      }
      catch(...)
      {
        Glib::exception_handlers_invoke();
      }
    }
  }

  const auto base = static_cast<BaseClassType*>(
        g_type_class_peek_parent(G_OBJECT_GET_CLASS(self)) // Get the parent class of the object class (The original underlying C class).
    );

  // Call the original underlying C function:
  if(base && base->frequency_poll)
    (*base->frequency_poll)(self);
}
void RackTuner_Class::poll_status_changed_callback(GxRackTuner* self, gboolean p0)
{
  const auto obj_base = static_cast<Glib::ObjectBase*>(
      Glib::ObjectBase::_get_current_wrapper((GObject*)self));

  // Non-gtkmmproc-generated custom classes implicitly call the default
  // Glib::ObjectBase constructor, which sets is_derived_. But gtkmmproc-
  // generated classes can use this optimisation, which avoids the unnecessary
  // parameter conversions if there is no possibility of the virtual function
  // being overridden:
  if(obj_base && obj_base->is_derived_())
  {
    const auto obj = dynamic_cast<CppObjectType* const>(obj_base);
    if(obj) // This can be NULL during destruction.
    {
      try // Trap C++ exceptions which would normally be lost because this is a C callback.
      {
        // Call the virtual member method, which derived classes might override.
        obj->on_poll_status_changed(p0
);
        return;
      }
      catch(...)
      {
        Glib::exception_handlers_invoke();
      }
    }
  }

  const auto base = static_cast<BaseClassType*>(
        g_type_class_peek_parent(G_OBJECT_GET_CLASS(self)) // Get the parent class of the object class (The original underlying C class).
    );

  // Call the original underlying C function:
  if(base && base->poll_status_changed)
    (*base->poll_status_changed)(self, p0);
}


Glib::ObjectBase* RackTuner_Class::wrap_new(GObject* o)
{
  return manage(new RackTuner((GxRackTuner*)(o)));

}


/* The implementation: */

RackTuner::RackTuner(const Glib::ConstructParams& construct_params)
:
  Gxw::Tuner(construct_params)
{
  }

RackTuner::RackTuner(GxRackTuner* castitem)
:
  Gxw::Tuner((GxTuner*)(castitem))
{
  }


RackTuner::RackTuner(RackTuner&& src) noexcept
: Gxw::Tuner(std::move(src))
{}

RackTuner& RackTuner::operator=(RackTuner&& src) noexcept
{
  Gxw::Tuner::operator=(std::move(src));
  return *this;
}

RackTuner::~RackTuner() noexcept
{
  destroy_();
}

RackTuner::CppClassType RackTuner::racktuner_class_; // initialize static member

GType RackTuner::get_type()
{
  return racktuner_class_.init().get_type();
}


GType RackTuner::get_base_type()
{
  return gx_rack_tuner_get_type();
}


RackTuner::RackTuner()
:
  // Mark this class as non-derived to allow C++ vfuncs to be skipped.
  Glib::ObjectBase(nullptr),
  Gxw::Tuner(Glib::ConstructParams(racktuner_class_.init()))
{
  

}

bool RackTuner::get_poll_status()
{
  return gx_rack_tuner_get_poll_status(gobj());
}

void RackTuner::set_freq(double p1)
{
  gx_rack_tuner_set_freq(gobj(), p1);
}

void RackTuner::set_scale_lim(double p1)
{
  gx_rack_tuner_set_scale_lim(gobj(), p1);
}

double RackTuner::get_scale_lim()
{
  return gx_rack_tuner_get_scale_lim(gobj());
}

void RackTuner::set_speed(double p1)
{
  gx_rack_tuner_set_speed(gobj(), p1);
}

double RackTuner::get_speed()
{
  return gx_rack_tuner_get_speed(gobj());
}

void RackTuner::set_streaming(bool p1)
{
  gx_rack_tuner_set_streaming(gobj(), static_cast<int>(p1));
}

bool RackTuner::get_streaming()
{
  return gx_rack_tuner_get_streaming(gobj());
}

void RackTuner::set_display_flat(bool p1)
{
  gx_rack_tuner_set_display_flat(gobj(), static_cast<int>(p1));
}

bool RackTuner::get_display_flat()
{
  return gx_rack_tuner_get_display_flat(gobj());
}

void RackTuner::set_timestep(int p1)
{
  gx_rack_tuner_set_timestep(gobj(), p1);
}

int RackTuner::get_timestep()
{
  return gx_rack_tuner_get_timestep(gobj());
}

void RackTuner::set_limit_timestep(int p1)
{
  gx_rack_tuner_set_limit_timestep(gobj(), p1);
}

int RackTuner::get_limit_timestep()
{
  return gx_rack_tuner_get_limit_timestep(gobj());
}

void RackTuner::set_temperament(int p1)
{
  gx_rack_tuner_set_temperament(gobj(), p1);
}

int RackTuner::get_temperament()
{
  return gx_rack_tuner_get_temperament(gobj());
}

void RackTuner::clear_notes()
{
  gx_rack_tuner_clear_notes(gobj());
}

void RackTuner::push_note(int p1, int p2, int p3)
{
  gx_rack_tuner_push_note(gobj(), p1, p2, p3);
}


Glib::SignalProxy< void > RackTuner::signal_frequency_poll()
{
  return Glib::SignalProxy< void >(this, &RackTuner_signal_frequency_poll_info);
}


Glib::SignalProxy< void,bool > RackTuner::signal_poll_status_changed()
{
  return Glib::SignalProxy< void,bool >(this, &RackTuner_signal_poll_status_changed_info);
}


Glib::PropertyProxy< double > RackTuner::property_freq() 
{
  return Glib::PropertyProxy< double >(this, "freq");
}

Glib::PropertyProxy_ReadOnly< double > RackTuner::property_freq() const
{
  return Glib::PropertyProxy_ReadOnly< double >(this, "freq");
}

Glib::PropertyProxy< double > RackTuner::property_scale_lim() 
{
  return Glib::PropertyProxy< double >(this, "scale-lim");
}

Glib::PropertyProxy_ReadOnly< double > RackTuner::property_scale_lim() const
{
  return Glib::PropertyProxy_ReadOnly< double >(this, "scale-lim");
}

Glib::PropertyProxy< double > RackTuner::property_speed() 
{
  return Glib::PropertyProxy< double >(this, "speed");
}

Glib::PropertyProxy_ReadOnly< double > RackTuner::property_speed() const
{
  return Glib::PropertyProxy_ReadOnly< double >(this, "speed");
}

Glib::PropertyProxy< bool > RackTuner::property_display_flat() 
{
  return Glib::PropertyProxy< bool >(this, "display-flat");
}

Glib::PropertyProxy_ReadOnly< bool > RackTuner::property_display_flat() const
{
  return Glib::PropertyProxy_ReadOnly< bool >(this, "display-flat");
}

Glib::PropertyProxy< bool > RackTuner::property_streaming() 
{
  return Glib::PropertyProxy< bool >(this, "streaming");
}

Glib::PropertyProxy_ReadOnly< bool > RackTuner::property_streaming() const
{
  return Glib::PropertyProxy_ReadOnly< bool >(this, "streaming");
}

Glib::PropertyProxy< int > RackTuner::property_timestep() 
{
  return Glib::PropertyProxy< int >(this, "timestep");
}

Glib::PropertyProxy_ReadOnly< int > RackTuner::property_timestep() const
{
  return Glib::PropertyProxy_ReadOnly< int >(this, "timestep");
}

Glib::PropertyProxy< int > RackTuner::property_limit_timestep() 
{
  return Glib::PropertyProxy< int >(this, "limit-timestep");
}

Glib::PropertyProxy_ReadOnly< int > RackTuner::property_limit_timestep() const
{
  return Glib::PropertyProxy_ReadOnly< int >(this, "limit-timestep");
}

Glib::PropertyProxy< int > RackTuner::property_temperament() 
{
  return Glib::PropertyProxy< int >(this, "temperament");
}

Glib::PropertyProxy_ReadOnly< int > RackTuner::property_temperament() const
{
  return Glib::PropertyProxy_ReadOnly< int >(this, "temperament");
}


void Gxw::RackTuner::on_frequency_poll()
{
  const auto base = static_cast<BaseClassType*>(
      g_type_class_peek_parent(G_OBJECT_GET_CLASS(gobject_)) // Get the parent class of the object class (The original underlying C class).
  );

  if(base && base->frequency_poll)
    (*base->frequency_poll)(gobj());
}
void Gxw::RackTuner::on_poll_status_changed(bool p1)
{
  const auto base = static_cast<BaseClassType*>(
      g_type_class_peek_parent(G_OBJECT_GET_CLASS(gobject_)) // Get the parent class of the object class (The original underlying C class).
  );

  if(base && base->poll_status_changed)
    (*base->poll_status_changed)(gobj(),static_cast<int>(p1));
}


} // namespace Gxw


