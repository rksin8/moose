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
#include "TransientMultiApp.h"
#include "SamplerInterface.h"

#include "StochasticToolsTypes.h"

class Sampler;
class StochasticToolsTransfer;

class SamplerTransientMultiApp : public TransientMultiApp, public SamplerInterface
{
public:
  static InputParameters validParams();

  SamplerTransientMultiApp(const InputParameters & parameters);

  /**
   * Override solveStep to allow for batch execution.
   */
  virtual bool solveStep(Real dt, Real target_time, bool auto_advance = true) override;

  /**
   * Override to initialize batch backups.
   */
  virtual void initialSetup() override;

protected:
  /// Sampler to utilize for creating MultiApps
  Sampler & _sampler;

  /// The Sup-application solve mode
  const StochasticTools::MultiAppMode _mode;

  /// Override to allow to get correct cli_args
  virtual std::string getCommandLineArgsParamHelper(unsigned int local_app) override;

private:
  /**
   * Helper method for running in mode='batch'
   * */
  bool solveStepBatch(Real dt, Real target_time, bool auto_advance = true);

  /**
   * Helper for getting StochasticToolsTransfer objects.
   *
   * This is a copy from SamplerFullSolveMultiapp, but the alternative is to create an intermediate
   * base. But, given the inheritance difference between these object that gets rather complex, so
   * a few lines of copied code is better for now.
   */
  std::vector<std::shared_ptr<StochasticToolsTransfer>>
  getActiveStochasticToolsTransfers(Transfer::DIRECTION direction);

  /// Storage for batch-restore mode; the outer vector if for the local stochastic data and the
  /// inner vector is for the number of sub-apps. The later is 1 for this object, but it is included
  /// in case that changes in the future or in child classes
  std::vector<std::vector<std::shared_ptr<Backup>>> _batch_backup;

  /// Store the number of rows initialized, if this changes error because it doesn't make sense
  const dof_id_type _number_of_sampler_rows;

  ///@{
  /// PrefGraph timers
  const PerfID _perf_solve_step;
  const PerfID _perf_solve_batch_step;
  const PerfID _perf_initial_setup;
  ///@}
};
