/*
  Copyright (C) 2018 - 2020 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.

  ASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/


#ifndef _aspect_time_stepping_composition_time_step_h
#define _aspect_time_stepping_composition_time_step_h

#include <aspect/time_stepping/interface.h>


namespace aspect
{
  namespace TimeStepping
  {
    using namespace dealii;

    /**
     * Compute the maximum time step based on the current solution
     * of the compositional field and
     * return this as the time step.
     *
     * @ingroup TimeStepping
     */
    template <int dim>
    class CompositionTimeStep : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Constructor.
         */
        CompositionTimeStep () = default;


        /**
         * @copydoc aspect::TimeStepping::Interface<dim>::execute()
         */
        double
        execute() override;
      /**
       * Declare the parameters this class takes through input files.
       */
      static void
      declare_parameters(ParameterHandler &prm);

      /**
       * Read the parameters this class declares from the parameter file.
       */
      void
      parse_parameters(ParameterHandler &prm) override;
      /**
       * @}
       */
      
      private:
        double max_change;

    };
  }
}


#endif
