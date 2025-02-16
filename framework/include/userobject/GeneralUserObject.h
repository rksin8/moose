//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

// MOOSE includes
#include "UserObject.h"
#include "MaterialPropertyInterface.h"
#include "TransientInterface.h"
#include "DependencyResolverInterface.h"
#include "ReporterInterface.h"

// Forward Declarations
class GeneralUserObject;

template <>
InputParameters validParams<GeneralUserObject>();

/* This class is here to combine the Postprocessor interface and the
 * base class Postprocessor object along with adding MooseObject to the inheritance tree*/
class GeneralUserObject : public UserObject,
                          public MaterialPropertyInterface,
                          public TransientInterface,
                          public DependencyResolverInterface
{
public:
  static InputParameters validParams();

  GeneralUserObject(const InputParameters & parameters);

  const std::set<std::string> & getRequestedItems() override;

  const std::set<std::string> & getSuppliedItems() override;

  ///@{
  /**
   * This method is not used and should not be used in a custom GeneralUserObject.
   */
  virtual void threadJoin(const UserObject &) override;
  virtual void subdomainSetup() override;
  ///@}

  ///@{
  /**
   * Store dependency among same object types for proper execution order
   */
  const PostprocessorValue & getPostprocessorValue(const std::string & name,
                                                   unsigned int index = 0) const override final;
  const PostprocessorValue &
  getPostprocessorValueByName(const PostprocessorName & name) const override final;

  const VectorPostprocessorValue &
  getVectorPostprocessorValue(const std::string & name,
                              const std::string & vector_name) const override final;
  const VectorPostprocessorValue &
  getVectorPostprocessorValueByName(const VectorPostprocessorName & name,
                                    const std::string & vector_name) const override final;

  const VectorPostprocessorValue &
  getVectorPostprocessorValue(const std::string & name,
                              const std::string & vector_name,
                              bool use_broadcast) const override final;
  const VectorPostprocessorValue &
  getVectorPostprocessorValueByName(const VectorPostprocessorName & name,
                                    const std::string & vector_name,
                                    bool use_broadcast) const override final;

  ///@}

protected:
  virtual void addReporterDependencyHelper(const ReporterName & state_name) override;

protected:
  mutable std::set<std::string> _depend_vars;
  std::set<std::string> _supplied_vars;
};
