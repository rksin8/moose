//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "SideIntegralVariablePostprocessor.h"
#include "MaterialPropertyInterface.h"
#include "Material.h"

class MassSideFluxIntegral : public SideIntegralVariablePostprocessor
{
  public:
    static InputParameters validParams();

    MassSideFluxIntegral(const InputParameters & parameters);

  protected:

    virtual Real computeQpIntegral() override;

    /// PorousFlowDictator UserObject
  const PorousFlowDictator & _dictator;

  /// Number of phases in the simulation
  const unsigned _num_phases;

  /// Indicates the type of flux that this BC will apply
  const enum class FluxTypeChoiceEnum { FLUID, HEAT } _flux_type;

  /// The fluid component number (only used if _flux_type==FLUID))
  const unsigned int _sp;

  /// Whether to multiply the Darcy flux by the fluid density.  This is automatically set to true if _flux_type==HEAT
  const bool _multiply_by_density;

  /// Whether to multiply the Darcy flux by the relative permeability
  const bool _include_relperm;

  /// Gravitational acceleration
  const RealVectorValue _gravity;

  /// Multiply the flux by this quantity.  This is mainly used for testing purposes, for instance, testing the Jacobian, and should ordinarily set to its default value of 1.0
  // const Real _multiplier;

  /// Gradient of the pore pressure in each phase
  const MaterialProperty<std::vector<RealGradient>> & _grad_p;

  /// Fluid density for each phase (at the qp)
  const MaterialProperty<std::vector<Real>> & _fluid_density_qp;

  /// Permeability of porous material
  const MaterialProperty<RealTensorValue> & _permeability;

  /// Viscosity of each phase
  const MaterialProperty<std::vector<Real>> & _fluid_viscosity;

  /// Whether there is a fluid_phase_density_nodal Material
  const bool _has_density;

  /// Whether there is a mass_frac_nodal Material
  const bool _has_mass_fraction;

  /// Whether there is a relative_permeability_nodal Material
  const bool _has_relperm;

  /// Whether there is an enthalpy Material
  const bool _has_enthalpy;

  /// Whether there is a thermal_conductivity_qp Material
  const bool _has_thermal_conductivity;

  /// Whether there is a grad_temp Material
  const bool _has_t;

  /// Fluid density for each phase (at the node)
  const MaterialProperty<std::vector<Real>> * const _fluid_density_node;

  /// Relative permeability of each phase
  const MaterialProperty<std::vector<Real>> * const _relative_permeability;

  /// Mass fraction of each component in each phase
  const MaterialProperty<std::vector<std::vector<Real>>> * const _mass_fractions;

  /// Enthalpy of each phase
  const MaterialProperty<std::vector<Real>> * const _enthalpy;

  /// Thermal_Conductivity of porous material
  const MaterialProperty<RealTensorValue> * const _thermal_conductivity;

  /// grad(temperature)
  const MaterialProperty<RealGradient> * const _grad_t;

private:

  /// Darcy contribution to the outflowBC
  Real darcy(unsigned ph) const;

  /// Mobility contribution to the outflowBC
  Real mobility(unsigned ph) const;

  /// Either mass_fraction or enthalpy, depending on _flux_type
  Real prefactor(unsigned ph) const;
    
};








