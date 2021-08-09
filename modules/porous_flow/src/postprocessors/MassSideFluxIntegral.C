//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PorousFlowFluidMass.h"

#include "MooseVariable.h"

#include "libmesh/quadrature.h"

#include "MassSideFluxIntegral.h"

registerMooseObject("PorousFlowApp", MassSideFluxIntegral);

InputParameters MassSideFluxIntegral::validParams()
{
  InputParameters params = SideIntegralVariablePostprocessor::validParams();
  params.addRequiredParam<UserObjectName>(
      "PorousFlowDictator", "The UserObject that holds the list of PorousFlow variable names");
  params.addRequiredParam<RealVectorValue>("gravity",
                                           "Gravitational acceleration vector downwards (m/s^2)");
                                           
  MooseEnum flux_type_enum("fluid heat", "fluid");
  params.addParam<MooseEnum>(
      "flux_type",
      flux_type_enum,
      "The type of boundary condition to apply.  'fluid' means this boundary condition will allow "
      "a fluid component to flow freely from the boundary.  'heat' means this boundary condition "
      "will allow heat energy to flow freely from the boundary");
  params.addParam<unsigned int>(
      "mass_fraction_component",
      0,
      "The index corresponding to a fluid "
      "component.  If supplied, the residual contribution will be "
      "multiplied by the mass fraction, corresponding to allowing the "
      "given mass fraction to flow freely from the boundary.  This is ignored if flux_type = heat");
  params.addParam<bool>("multiply_by_density",
                        true,
                        "If true, this BC represents mass flux.  If false, it represents volume "
                        "flux.  User input of this flag is ignored if flux_type = heat as "
                        "multiply_by_density should always be true in that case");
  params.addParam<bool>(
      "include_relperm",
      true,
      "If true, the Darcy flux will include the relative permeability.  If false, the relative "
      "permeability will not be used, which must only be used for fully-saturated situations "
      "where there is no notion of relative permeability");
//  params.addParam<Real>(
//      "multiplier",
//      1.0,
//      "Multiply the flux by this number.  This is mainly used for testing purposes");
//  params.addParamNamesToGroup("multiplier", "Advanced");
  params.addClassDescription(
      "Applies an 'outflow' boundary condition, which allows fluid components or heat energy to "
      "flow freely out of the boundary as if it weren't there.  This is fully upwinded");
  return params;
}


MassSideFluxIntegral::MassSideFluxIntegral(const InputParameters & parameters):
  SideIntegralVariablePostprocessor(parameters),
  _dictator(getUserObject<PorousFlowDictator>("PorousFlowDictator")),
  _num_phases(_dictator.numPhases()),
  _flux_type(getParam<MooseEnum>("flux_type").getEnum<FluxTypeChoiceEnum>()),
  _sp(getParam<unsigned>("mass_fraction_component")),
  _multiply_by_density(
        _flux_type == FluxTypeChoiceEnum::FLUID ? getParam<bool>("multiply_by_density") : true),
  _include_relperm(getParam<bool>("include_relperm")),
  _gravity(getParam<RealVectorValue>("gravity")),
//  _multiplier(getParam<Real>("multiplier")),
  _grad_p(getMaterialProperty<std::vector<RealGradient>>("PorousFlow_grad_porepressure_qp")),
  _fluid_density_qp(getMaterialProperty<std::vector<Real>>("PorousFlow_fluid_phase_density_qp")),
  _permeability(getMaterialProperty<RealTensorValue>("PorousFlow_permeability_qp")),
  _fluid_viscosity(getMaterialProperty<std::vector<Real>>("PorousFlow_viscosity_nodal")),      
  _has_density(hasMaterialProperty<std::vector<Real>>("PorousFlow_fluid_phase_density_qp")),
  _has_mass_fraction(
        hasMaterialProperty<std::vector<std::vector<Real>>>("PorousFlow_mass_frac_nodal")),
  _has_relperm(hasMaterialProperty<std::vector<Real>>("PorousFlow_relative_permeability_qp")),
  _has_enthalpy(hasMaterialProperty<std::vector<Real>>("PorousFlow_fluid_phase_enthalpy_nodal")),
  _has_thermal_conductivity(
        hasMaterialProperty<RealTensorValue>("PorousFlow_thermal_conductivity_qp")),
  _has_t(hasMaterialProperty<RealGradient>("PorousFlow_grad_temperature_qp")),

 // Nodal value why?
  _fluid_density_node(_has_density ? &getMaterialProperty<std::vector<Real>>(
                                           "PorousFlow_fluid_phase_density_qp")  // nodal
                                     : nullptr),
  _relative_permeability(_has_relperm ? &getMaterialProperty<std::vector<Real>>(
                                              "PorousFlow_relative_permeability_qp") // nodal
                                        : nullptr),
  _mass_fractions(_has_mass_fraction ? &getMaterialProperty<std::vector<std::vector<Real>>>(
                                             "PorousFlow_mass_frac_nodal")
                                       : nullptr), // nodal
                                       
  _enthalpy(_has_enthalpy ? &getMaterialPropertyByName<std::vector<Real>>(
                                  "PorousFlow_fluid_phase_enthalpy_nodal")
                            : nullptr), // nodal 
  _thermal_conductivity(_has_thermal_conductivity ? &getMaterialProperty<RealTensorValue>(
                                                          "PorousFlow_thermal_conductivity_qp")
                                                    : nullptr),                                                         
  _grad_t(_has_t ? &getMaterialProperty<RealGradient>("PorousFlow_grad_temperature_qp")
                   : nullptr)
{

  if (_flux_type == FluxTypeChoiceEnum::FLUID && _sp >= _dictator.numComponents())
    paramError("mass_fraction_component",
               "The Dictator declares that the maximum fluid component index is ",
               _dictator.numComponents() - 1,
               ", but you have set mass_fraction_component to ",
               _sp,
               ". Remember that indexing starts at 0. Please be assured that the Dictator has "
               "noted your error.");

  if (_multiply_by_density && !_has_density)
    mooseError("PorousFlowOutflowBC: You requested to multiply_by_density, but you have no nodal "
               "fluid density Material");
  if (_include_relperm && !_has_relperm)
    mooseError("PorousFlowOutflowBC: You requested to include the relative permeability, but you "
               "have no nodal relative permeability Material");
  if (_flux_type == FluxTypeChoiceEnum::FLUID && !_has_mass_fraction)
    mooseError(
        "PorousFlowOutflowBC: For flux_type = fluid, you need a nodal mass-fraction Material");
  if (_flux_type == FluxTypeChoiceEnum::HEAT && !_has_enthalpy)
    mooseError(
        "PorousFlowOutflowBC: For flux_type = heat, you need a nodal fluid enthalpy Material");
  if (_flux_type == FluxTypeChoiceEnum::HEAT && !_has_thermal_conductivity)
    mooseError(
        "PorousFlowOutflowBC: For flux_type = heat, you need a thermal conductivity Material");
  if (_flux_type == FluxTypeChoiceEnum::HEAT && !_has_t)
    mooseError("PorousFlowOutflowBC: For flux_type = heat, you need a temperature Material");
}



Real
MassSideFluxIntegral::computeQpIntegral()
{
  Real flux = 0.0;
  
  for (unsigned ph = 0; ph < _num_phases; ++ph)
       flux-= darcy(ph) * mobility(ph) * prefactor(ph);
       
//  if (_flux_type == FluxTypeChoiceEnum::FLUID)
//    return _test[_i][_qp] * advective_term * _multiplier;

//  const Real conduction_term = -_normals[_qp] * ((*_thermal_conductivity)[_qp] * (*_grad_t)[_qp]);

  return flux;
}

Real
MassSideFluxIntegral::darcy(unsigned ph) const
{
  return _normals[_qp] *
         (_permeability[_qp] * (_grad_p[_qp][ph] - _fluid_density_qp[_qp][ph] * _gravity));
}
Real
MassSideFluxIntegral::mobility(unsigned ph) const
{         
  return (_multiply_by_density ? (_fluid_density_qp)[_qp][ph] : 1.0) *
         (_include_relperm ? (*_relative_permeability)[_qp][ph] : 1.0) / _fluid_viscosity[_qp][ph];
}

Real
MassSideFluxIntegral::prefactor(unsigned ph) const
{
  // FIXME: with qp properies                                            
  return (_flux_type == FluxTypeChoiceEnum::FLUID ? (*_mass_fractions)[_qp][ph][_sp]
                                                  : (*_enthalpy)[_qp][ph]);  
                                                  
 // return (_flux_type == FluxTypeChoiceEnum::FLUID ? (*_mass_fractions)[_i][ph][_sp]
 //                                                 : (*_enthalpy)[_i][ph]);
}


