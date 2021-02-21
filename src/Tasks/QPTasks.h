/*
 * Copyright 2012-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#pragma once

// includes
//
#include <array>
#include <set>

// Eigen
#include <Eigen/Core>

// RBDyn, added by Rafa
#include <RBDyn/CoM.h>
#include <RBDyn/Momentum.h>

// Tasks
#include "QPMotionConstr.h"
#include "QPSolver.h"
#include "Tasks.h"

// forward declaration
// RBDyn
namespace rbd
{
class MultiBody;
struct MultiBodyConfig;
} // namespace rbd

namespace tasks
{

namespace qp
{

class TASKS_DLLAPI SetPointTaskCommon : public Task
{
public:
  SetPointTaskCommon(const std::vector<rbd::MultiBody> & mbs, int robotIndex, HighLevelTask * hlTask, double weight);
  SetPointTaskCommon(const std::vector<rbd::MultiBody> & mbs,
                     int robotIndex,
                     HighLevelTask * hlTask,
                     const Eigen::VectorXd & dimWeight,
                     double weight);

  virtual std::string nameTask() const override
  {
    return hlTask_->nameHighLevelTask();
  }
  
  virtual std::pair<int, int> begin() const override
  {
    return std::make_pair(alphaDBegin_, alphaDBegin_);
  }

  void dimWeight(const Eigen::VectorXd & dim);

  const Eigen::VectorXd & dimWeight() const
  {
    return dimWeight_;
  }

  virtual void updateNrVars(const std::vector<rbd::MultiBody> & mbs, const SolverData & data) override;

  virtual const Eigen::MatrixXd & Q() const override;
  virtual const Eigen::VectorXd & C() const override;

protected:
  void computeQC(Eigen::VectorXd & error);

protected:
  HighLevelTask * hlTask_;
  Eigen::VectorXd error_;

private:
  Eigen::VectorXd dimWeight_;
  int robotIndex_, alphaDBegin_;

  Eigen::MatrixXd Q_;
  Eigen::VectorXd C_;
  // cache
  Eigen::MatrixXd preQ_;
  Eigen::VectorXd preC_;
};

class TASKS_DLLAPI SetPointTask : public SetPointTaskCommon
{
public:
  SetPointTask(const std::vector<rbd::MultiBody> & mbs,
               int robotIndex,
               HighLevelTask * hlTask,
               double stiffness,
               double weight);

  SetPointTask(const std::vector<rbd::MultiBody> & mbs,
               int robotIndex,
               HighLevelTask * hlTask,
               double stiffness,
               const Eigen::VectorXd & dimWeight,
               double weight);

  double stiffness() const
  {
    return stiffness_;
  }

  void stiffness(double stiffness);

  virtual void update(const std::vector<rbd::MultiBody> & mbs,
                      const std::vector<rbd::MultiBodyConfig> & mbcs,
                      const SolverData & data) override;

private:
  double stiffness_, stiffnessSqrt_;
};

class TASKS_DLLAPI TrackingTask : public SetPointTaskCommon
{
public:
  TrackingTask(const std::vector<rbd::MultiBody> & mbs,
               int robotIndex,
               HighLevelTask * hlTask,
               double gainPos,
               double gainVel,
               double weight);

  TrackingTask(const std::vector<rbd::MultiBody> & mbs,
               int robotIndex,
               HighLevelTask * hlTask,
               double gainPos,
               double gainVel,
               const Eigen::VectorXd & dimWeight,
               double weight);

  void setGains(double gainPos, double gainVel);

  void errorPos(const Eigen::VectorXd & errorPos);
  void errorVel(const Eigen::VectorXd & errorVel);
  void refAccel(const Eigen::VectorXd & refAccel);

  virtual void update(const std::vector<rbd::MultiBody> & mbs,
                      const std::vector<rbd::MultiBodyConfig> & mbcs,
                      const SolverData & data) override;

private:
  double gainPos_, gainVel_;
  Eigen::VectorXd errorPos_, errorVel_, refAccel_;
};

class TASKS_DLLAPI TrajectoryTask : public SetPointTaskCommon
{
public:
  TrajectoryTask(const std::vector<rbd::MultiBody> & mbs,
                 int robotIndex,
                 HighLevelTask * hlTask,
                 double gainPos,
                 double gainVel,
                 double weight);

  TrajectoryTask(const std::vector<rbd::MultiBody> & mbs,
                 int robotIndex,
                 HighLevelTask * hlTask,
                 double gainPos,
                 double gainVel,
                 const Eigen::VectorXd & dimWeight,
                 double weight);

  void setGains(double gainPos, double gainVel);
  void setGains(const Eigen::VectorXd & stiffness, const Eigen::VectorXd & damping);
  void stiffness(double gainPos);
  void stiffness(const Eigen::VectorXd & stiffness);
  const Eigen::VectorXd & stiffness() const;
  void damping(double gainVel);
  void damping(const Eigen::VectorXd & damping);
  const Eigen::VectorXd & damping() const;

  void refVel(const Eigen::VectorXd & refVel);
  const Eigen::VectorXd & refVel() const;
  void refAccel(const Eigen::VectorXd & refAccel);
  const Eigen::VectorXd & refAccel() const;

  void update(const std::vector<rbd::MultiBody> & mbs,
              const std::vector<rbd::MultiBodyConfig> & mbcs,
              const SolverData & data) override;

private:
  Eigen::VectorXd stiffness_, damping_;
  Eigen::VectorXd refVel_, refAccel_;
};

/// @deprecated Must be replace by TrackingTask
class TASKS_DLLAPI PIDTask : public SetPointTaskCommon
{
public:
  PIDTask(const std::vector<rbd::MultiBody> & mbs,
          int robotIndex,
          HighLevelTask * hlTask,
          double P,
          double I,
          double D,
          double weight);

  PIDTask(const std::vector<rbd::MultiBody> & mbs,
          int robotIndex,
          HighLevelTask * hlTask,
          double P,
          double I,
          double D,
          const Eigen::VectorXd & dimWeight,
          double weight);

  double P() const;
  void P(double p);
  double I() const;
  void I(double i);
  double D() const;
  void D(double d);

  void error(const Eigen::VectorXd & err);
  void errorD(const Eigen::VectorXd & errD);
  void errorI(const Eigen::VectorXd & errI);

  virtual void update(const std::vector<rbd::MultiBody> & mbs,
                      const std::vector<rbd::MultiBodyConfig> & mbcs,
                      const SolverData & data) override;

private:
  double P_, I_, D_;
  Eigen::VectorXd error_, errorD_, errorI_;
};

class TASKS_DLLAPI TargetObjectiveTask : public Task
{
public:
  TargetObjectiveTask(const std::vector<rbd::MultiBody> & mbs,
                      int robotIndex,
                      HighLevelTask * hlTask,
                      double timeStep,
                      double duration,
                      const Eigen::VectorXd & objDot,
                      double weight);

  TargetObjectiveTask(const std::vector<rbd::MultiBody> & mbs,
                      int robotIndex,
                      HighLevelTask * hlTask,
                      double timeStep,
                      double duration,
                      const Eigen::VectorXd & objDot,
                      const Eigen::VectorXd & dimWeight,
                      double weight);

  double duration() const;
  void duration(double d);

  int iter() const
  {
    return iter_;
  }
  void iter(int i)
  {
    iter_ = i;
  }

  int nrIter() const
  {
    return nrIter_;
  }
  void nrIter(int i)
  {
    nrIter_ = i;
  }

  const Eigen::VectorXd & objDot() const
  {
    return objDot_;
  }
  void objDot(const Eigen::VectorXd & o)
  {
    objDot_ = o;
  }

  const Eigen::VectorXd & dimWeight() const
  {
    return dimWeight_;
  }
  void dimWeight(const Eigen::VectorXd & o)
  {
    dimWeight_ = o;
  }

  const Eigen::VectorXd & phi() const
  {
    return phi_;
  }
  const Eigen::VectorXd & psi() const
  {
    return psi_;
  }

  virtual std::pair<int, int> begin() const override
  {
    return std::make_pair(alphaDBegin_, alphaDBegin_);
  }

  virtual void updateNrVars(const std::vector<rbd::MultiBody> & mbs, const SolverData & data) override;
  virtual void update(const std::vector<rbd::MultiBody> & mbs,
                      const std::vector<rbd::MultiBodyConfig> & mbcs,
                      const SolverData & data) override;

  virtual const Eigen::MatrixXd & Q() const override;
  virtual const Eigen::VectorXd & C() const override;

private:
  HighLevelTask * hlTask_;

  int iter_, nrIter_;
  double dt_;
  Eigen::VectorXd objDot_;
  Eigen::VectorXd dimWeight_;
  int robotIndex_, alphaDBegin_;

  Eigen::VectorXd phi_, psi_;

  Eigen::MatrixXd Q_;
  Eigen::VectorXd C_;
  // cache
  Eigen::MatrixXd preQ_;
  Eigen::VectorXd CVecSum_, preC_;
};

class TASKS_DLLAPI JointsSelector : public HighLevelTask
{
public:
  static JointsSelector ActiveJoints(const std::vector<rbd::MultiBody> & mbs,
                                     int robotIndex,
                                     HighLevelTask * hl,
                                     const std::vector<std::string> & activeJointsName,
                                     const std::map<std::string, std::vector<std::array<int, 2>>> & activeDofs = {});
  static JointsSelector UnactiveJoints(const std::vector<rbd::MultiBody> & mbs,
                                       int robotIndex,
                                       HighLevelTask * hl,
                                       const std::vector<std::string> & unactiveJointsName,
                                       const std::map<std::string, std::vector<std::array<int, 2>>> & unactiveDofs = {});

public:
  struct SelectedData
  {
    int posInDof, dof;
  };

public:
  JointsSelector(const std::vector<rbd::MultiBody> & mbs,
                 int robotIndex,
                 HighLevelTask * hl,
                 const std::vector<std::string> & selectedJointsName,
                 const std::map<std::string, std::vector<std::array<int, 2>>> & activeDofs = {});

  const std::vector<SelectedData> selectedJoints() const
  {
    return selectedJoints_;
  }

  virtual int dim() override;
  virtual void update(const std::vector<rbd::MultiBody> & mbs,
                      const std::vector<rbd::MultiBodyConfig> & mbcs,
                      const SolverData & data) override;

  virtual const Eigen::MatrixXd & jac() override;
  virtual const Eigen::MatrixXd & jacDot() override;
  virtual const Eigen::VectorXd & eval() override;
  virtual const Eigen::VectorXd & speed() override;
  virtual const Eigen::VectorXd & normalAcc() override;

private:
  Eigen::MatrixXd jac_;
  Eigen::MatrixXd jacDot_;
  std::vector<SelectedData> selectedJoints_;
  HighLevelTask * hl_;
};

struct JointStiffness
{
  JointStiffness() : jointName(), stiffness() {}
  JointStiffness(const std::string & jName, double stif) : jointName(jName), stiffness(stif) {}

  std::string jointName;
  double stiffness;
};

struct JointGains
{
  JointGains() : jointName(), stiffness(), damping() {}
  JointGains(const std::string & jName, double stif) : jointName(jName), stiffness(stif)
  {
    damping = 2. * std::sqrt(stif);
  }

  JointGains(const std::string & jName, double stif, double damp) : jointName(jName), stiffness(stif), damping(damp) {}

  std::string jointName;
  double stiffness, damping;
};

class TASKS_DLLAPI TorqueTask : public Task
{
public:
  TorqueTask(const std::vector<rbd::MultiBody> & mbs, int robotIndex,
             const std::shared_ptr<rbd::ForwardDynamics> fd,
             const TorqueBound & tb, double weight);

  TorqueTask(const std::vector<rbd::MultiBody> & mbs, int robotIndex,
             const std::shared_ptr<rbd::ForwardDynamics> fd,
             const TorqueBound & tb, const Eigen::VectorXd & jointSelect,
             double weight);

  TorqueTask(const std::vector<rbd::MultiBody> & mbs, int robotIndex,
             const std::shared_ptr<rbd::ForwardDynamics> fd,
             const TorqueBound & tb, const std::string & efName,
             double weight);

  TorqueTask(const std::vector<rbd::MultiBody> & mbs, int robotIndex,
             const std::shared_ptr<rbd::ForwardDynamics> fd,
             const TorqueBound & tb,
             const TorqueDBound & tdb,
             double dt,
             double weight);

  TorqueTask(const std::vector<rbd::MultiBody> & mbs, int robotIndex,
             const std::shared_ptr<rbd::ForwardDynamics> fd,
             const TorqueBound & tb,
             const TorqueDBound & tdb,
             double dt,
             const Eigen::VectorXd & jointSelect,
             double weight);

  TorqueTask(const std::vector<rbd::MultiBody> & mbs, int robotIndex,
             const std::shared_ptr<rbd::ForwardDynamics> fd,
             const TorqueBound & tb,
             const TorqueDBound & tdb,
             double dt,
             const std::string & efName,
             double weight);

  virtual void updateNrVars(const std::vector<rbd::MultiBody> & mbs, const SolverData & data) override;
  virtual void update(const std::vector<rbd::MultiBody> & mbs,
                      const std::vector<rbd::MultiBodyConfig> & mbcs,
                      const SolverData & data) override;

  virtual std::pair<int, int> begin() const override
  {
    return std::make_pair(0, 0);
  }

  virtual const Eigen::MatrixXd & Q() const override
  {
    return Q_;
  }

  virtual const Eigen::VectorXd & C() const override
  {
    return C_;
  }

  virtual const Eigen::VectorXd & jointSelect() const
  {
    return jointSelector_;
  }

private:
  int robotIndex_;
  int alphaDBegin_, lambdaBegin_;
  MotionConstr motionConstr;
  Eigen::VectorXd jointSelector_;
  Eigen::MatrixXd Q_;
  Eigen::VectorXd C_;
};

class TASKS_DLLAPI PostureTask : public Task
{
public:
  PostureTask(const std::vector<rbd::MultiBody> & mbs,
              int robotIndex,
              std::vector<std::vector<double>> q,
              double stiffness,
              double weight);

  virtual std::string nameTask() const override
  {
    return "PostureTask";
  }

  tasks::PostureTask & task()
  {
    return pt_;
  }

  void posture(std::vector<std::vector<double>> q)
  {
    pt_.posture(q);
  }

  const std::vector<std::vector<double>> posture() const
  {
    return pt_.posture();
  }

  double stiffness() const
  {
    return stiffness_;
  }

  double damping() const
  {
    return damping_;
  }

  void stiffness(double stiffness);

  void gains(double stiffness);
  void gains(double stiffness, double damping);

  void jointsStiffness(const std::vector<rbd::MultiBody> & mbs, const std::vector<JointStiffness> & jsv);

  void jointsGains(const std::vector<rbd::MultiBody> & mbs, const std::vector<JointGains> & jgv);

  virtual std::pair<int, int> begin() const override
  {
    return std::make_pair(alphaDBegin_, alphaDBegin_);
  }

  virtual void updateNrVars(const std::vector<rbd::MultiBody> & mbs, const SolverData & data) override;
  virtual void update(const std::vector<rbd::MultiBody> & mbs,
                      const std::vector<rbd::MultiBodyConfig> & mbcs,
                      const SolverData & data) override;

  virtual const Eigen::MatrixXd & Q() const override;
  virtual const Eigen::VectorXd & C() const override;

  const Eigen::VectorXd & eval() const;

private:
  struct JointData
  {
    double stiffness, damping;
    int start, size;
  };

private:
  tasks::PostureTask pt_;

  double stiffness_;
  double damping_;
  int robotIndex_, alphaDBegin_;

  std::vector<JointData> jointDatas_;

  Eigen::MatrixXd Q_;
  Eigen::VectorXd C_;
  Eigen::VectorXd alphaVec_;
};

class TASKS_DLLAPI PositionTask : public HighLevelTask
{
public:
  PositionTask(const std::vector<rbd::MultiBody> & mbs,
               int robotIndex,
               const std::string & bodyName,
               const Eigen::Vector3d & pos,
               const Eigen::Vector3d & bodyPoint = Eigen::Vector3d::Zero());

  virtual std::string nameHighLevelTask() const override
  {
    return "PositionTask";
  }

  tasks::PositionTask & task()
  {
    return pt_;
  }

  void position(const Eigen::Vector3d & pos)
  {
    pt_.position(pos);
  }

  const Eigen::Vector3d & position() const
  {
    return pt_.position();
  }

  void bodyPoint(const Eigen::Vector3d & point)
  {
    pt_.bodyPoint(point);
  }

  const Eigen::Vector3d & bodyPoint() const
  {
    return pt_.bodyPoint();
  }

  virtual int dim() override;
  virtual void update(const std::vector<rbd::MultiBody> & mb,
                      const std::vector<rbd::MultiBodyConfig> & mbc,
                      const SolverData & data) override;

  virtual const Eigen::MatrixXd & jac() override;
  virtual const Eigen::MatrixXd & jacDot() override;
  virtual const Eigen::VectorXd & eval() override;
  virtual const Eigen::VectorXd & speed() override;
  virtual const Eigen::VectorXd & normalAcc() override;

private:
  tasks::PositionTask pt_;
  int robotIndex_;
};

class TASKS_DLLAPI OrientationTask : public HighLevelTask
{
public:
  OrientationTask(const std::vector<rbd::MultiBody> & mbs,
                  int robodIndex,
                  const std::string & bodyName,
                  const Eigen::Quaterniond & ori);
  OrientationTask(const std::vector<rbd::MultiBody> & mbs,
                  int robodIndex,
                  const std::string & bodyName,
                  const Eigen::Matrix3d & ori);

  virtual std::string nameHighLevelTask() const override
  {
    return "OrientationTask";
  }
  
  tasks::OrientationTask & task()
  {
    return ot_;
  }

  void orientation(const Eigen::Quaterniond & ori)
  {
    ot_.orientation(ori);
  }

  void orientation(const Eigen::Matrix3d & ori)
  {
    ot_.orientation(ori);
  }

  const Eigen::Matrix3d & orientation() const
  {
    return ot_.orientation();
  }

  virtual int dim() override;
  virtual void update(const std::vector<rbd::MultiBody> & mbs,
                      const std::vector<rbd::MultiBodyConfig> & mbcs,
                      const SolverData & data) override;

  virtual const Eigen::MatrixXd & jac() override;
  virtual const Eigen::MatrixXd & jacDot() override;
  virtual const Eigen::VectorXd & eval() override;
  virtual const Eigen::VectorXd & speed() override;
  virtual const Eigen::VectorXd & normalAcc() override;

private:
  tasks::OrientationTask ot_;
  int robotIndex_;
};

template<typename transform_task_t>
class TransformTaskCommon : public HighLevelTask
{
public:
  TransformTaskCommon(const std::vector<rbd::MultiBody> & mbs,
                      int robotIndex,
                      const std::string & bodyName,
                      const sva::PTransformd & X_0_t,
                      const sva::PTransformd & X_b_p = sva::PTransformd::Identity())
  : tt_(mbs[robotIndex], bodyName, X_0_t, X_b_p), robotIndex_(robotIndex)
  {
  }

  transform_task_t & task()
  {
    return tt_;
  }

  void target(const sva::PTransformd & X_0_t)
  {
    tt_.target(X_0_t);
  }

  const sva::PTransformd & target() const
  {
    return tt_.target();
  }

  void X_b_p(const sva::PTransformd & X_b_p)
  {
    tt_.X_b_p(X_b_p);
  }

  const sva::PTransformd & X_b_p() const
  {
    return tt_.X_b_p();
  }

  virtual int dim() override
  {
    return 6;
  }

  virtual const Eigen::MatrixXd & jac() override
  {
    return tt_.jac();
  }

  // Rafa's test, not really implemented
  virtual const Eigen::MatrixXd & jacDot() override
  {
    return tt_.jacDot();
  }

  virtual const Eigen::VectorXd & eval() override
  {
    return tt_.eval();
  }

  virtual const Eigen::VectorXd & speed() override
  {
    return tt_.speed();
  }

  virtual const Eigen::VectorXd & normalAcc() override
  {
    return tt_.normalAcc();
  }

protected:
  transform_task_t tt_;
  int robotIndex_;
};

/// TransformTask in surface frame.
class TASKS_DLLAPI SurfaceTransformTask : public TransformTaskCommon<tasks::SurfaceTransformTask>
{
public:
  SurfaceTransformTask(const std::vector<rbd::MultiBody> & mbs,
                       int robotIndex,
                       const std::string & bodyName,
                       const sva::PTransformd & X_0_t,
                       const sva::PTransformd & X_b_p = sva::PTransformd::Identity());

  virtual void update(const std::vector<rbd::MultiBody> & mb,
                      const std::vector<rbd::MultiBodyConfig> & mbc,
                      const SolverData & data) override;
};

/// TransformTask in world or user frame.
class TASKS_DLLAPI TransformTask : public TransformTaskCommon<tasks::TransformTask>
{
public:
  TransformTask(const std::vector<rbd::MultiBody> & mbs,
                int robotIndex,
                const std::string & bodyName,
                const sva::PTransformd & X_0_t,
                const sva::PTransformd & X_b_p = sva::PTransformd::Identity(),
                const Eigen::Matrix3d & E_0_c = Eigen::Matrix3d::Identity());

  void E_0_c(const Eigen::Matrix3d & E_0_c);
  const Eigen::Matrix3d & E_0_c() const;

  virtual void update(const std::vector<rbd::MultiBody> & mb,
                      const std::vector<rbd::MultiBodyConfig> & mbc,
                      const SolverData & data) override;
};

class TASKS_DLLAPI SurfaceOrientationTask : public HighLevelTask
{
public:
  SurfaceOrientationTask(const std::vector<rbd::MultiBody> & mbs,
                         int robodIndex,
                         const std::string & bodyName,
                         const Eigen::Quaterniond & ori,
                         const sva::PTransformd & X_b_s);
  SurfaceOrientationTask(const std::vector<rbd::MultiBody> & mbs,
                         int robodIndex,
                         const std::string & bodyName,
                         const Eigen::Matrix3d & ori,
                         const sva::PTransformd & X_b_s);

  tasks::SurfaceOrientationTask & task()
  {
    return ot_;
  }

  void orientation(const Eigen::Quaterniond & ori)
  {
    ot_.orientation(ori);
  }

  void orientation(const Eigen::Matrix3d & ori)
  {
    ot_.orientation(ori);
  }

  const Eigen::Matrix3d & orientation() const
  {
    return ot_.orientation();
  }

  virtual int dim() override;
  virtual void update(const std::vector<rbd::MultiBody> & mbs,
                      const std::vector<rbd::MultiBodyConfig> & mbcs,
                      const SolverData & data) override;

  virtual const Eigen::MatrixXd & jac() override;
  virtual const Eigen::MatrixXd & jacDot() override;
  virtual const Eigen::VectorXd & eval() override;
  virtual const Eigen::VectorXd & speed() override;
  virtual const Eigen::VectorXd & normalAcc() override;

private:
  tasks::SurfaceOrientationTask ot_;
  int robotIndex_;
};

class TASKS_DLLAPI GazeTask : public HighLevelTask
{
public:
  GazeTask(const std::vector<rbd::MultiBody> & mbs,
           int robotIndex,
           const std::string & bodyName,
           const Eigen::Vector2d & point2d,
           double depthEstimate,
           const sva::PTransformd & X_b_gaze,
           const Eigen::Vector2d & point2d_ref = Eigen::Vector2d::Zero());
  GazeTask(const std::vector<rbd::MultiBody> & mbs,
           int robotIndex,
           const std::string & bodyName,
           const Eigen::Vector3d & point3d,
           const sva::PTransformd & X_b_gaze,
           const Eigen::Vector2d & point2d_ref = Eigen::Vector2d::Zero());

  tasks::GazeTask & task()
  {
    return gazet_;
  }

  void error(const Eigen::Vector2d & point2d, const Eigen::Vector2d & point2d_ref = Eigen::Vector2d::Zero())
  {
    gazet_.error(point2d, point2d_ref);
  }

  void error(const Eigen::Vector3d & point3d, const Eigen::Vector2d & point2d_ref = Eigen::Vector2d::Zero())
  {
    gazet_.error(point3d, point2d_ref);
  }

  virtual int dim() override;
  virtual void update(const std::vector<rbd::MultiBody> & mbs,
                      const std::vector<rbd::MultiBodyConfig> & mbcs,
                      const SolverData & data) override;

  virtual const Eigen::MatrixXd & jac() override;
  virtual const Eigen::MatrixXd & jacDot() override;
  virtual const Eigen::VectorXd & eval() override;
  virtual const Eigen::VectorXd & speed() override;
  virtual const Eigen::VectorXd & normalAcc() override;

private:
  tasks::GazeTask gazet_;
  int robotIndex_;
};

class TASKS_DLLAPI PositionBasedVisServoTask : public HighLevelTask
{
public:
  PositionBasedVisServoTask(const std::vector<rbd::MultiBody> & mbs,
                            int robotIndex,
                            const std::string & bodyName,
                            const sva::PTransformd & X_t_s,
                            const sva::PTransformd & X_b_s = sva::PTransformd::Identity());

  tasks::PositionBasedVisServoTask & task()
  {
    return pbvst_;
  }

  void error(const sva::PTransformd & X_t_s)
  {
    pbvst_.error(X_t_s);
  }

  virtual int dim() override;
  virtual void update(const std::vector<rbd::MultiBody> & mbs,
                      const std::vector<rbd::MultiBodyConfig> & mbcs,
                      const SolverData & data) override;

  virtual const Eigen::MatrixXd & jac() override;
  virtual const Eigen::MatrixXd & jacDot() override;
  virtual const Eigen::VectorXd & eval() override;
  virtual const Eigen::VectorXd & speed() override;
  virtual const Eigen::VectorXd & normalAcc() override;

private:
  tasks::PositionBasedVisServoTask pbvst_;
  int robotIndex_;
};

class TASKS_DLLAPI CoMTask : public HighLevelTask
{
public:
  CoMTask(const std::vector<rbd::MultiBody> & mbs, int robotIndex, const Eigen::Vector3d & com);
  CoMTask(const std::vector<rbd::MultiBody> & mbs,
          int robotIndex,
          const Eigen::Vector3d & com,
          std::vector<double> weight);

  virtual std::string nameHighLevelTask() const override
  {
    return "CoMTask";
  }
  
  tasks::CoMTask & task()
  {
    return ct_;
  }

  void com(const Eigen::Vector3d & com)
  {
    ct_.com(com);
  }

  const Eigen::Vector3d & com() const
  {
    return ct_.com();
  }

  const Eigen::Vector3d & actual() const
  {
    return ct_.actual();
  }

  void updateInertialParameters(const std::vector<rbd::MultiBody> & mbs);

  virtual int dim() override;
  virtual void update(const std::vector<rbd::MultiBody> & mbs,
                      const std::vector<rbd::MultiBodyConfig> & mbcs,
                      const SolverData & data) override;

  virtual const Eigen::MatrixXd & jac() override;
  virtual const Eigen::MatrixXd & jacDot() override;
  virtual const Eigen::VectorXd & eval() override;
  virtual const Eigen::VectorXd & speed() override;
  virtual const Eigen::VectorXd & normalAcc() override;

private:
  tasks::CoMTask ct_;
  int robotIndex_;
};
 
class TASKS_DLLAPI MultiCoMTask : public Task
{
public:
  MultiCoMTask(const std::vector<rbd::MultiBody> & mb,
               std::vector<int> robotIndexes,
               const Eigen::Vector3d & com,
               double stiffness,
               double weight);
  MultiCoMTask(const std::vector<rbd::MultiBody> & mb,
               std::vector<int> robotIndexes,
               const Eigen::Vector3d & com,
               double stiffness,
               const Eigen::Vector3d & dimWeight,
               double weight);

  tasks::MultiCoMTask & task()
  {
    return mct_;
  }

  void com(const Eigen::Vector3d & com)
  {
    mct_.com(com);
  }

  const Eigen::Vector3d com() const
  {
    return mct_.com();
  }

  void updateInertialParameters(const std::vector<rbd::MultiBody> & mbs);

  double stiffness() const
  {
    return stiffness_;
  }

  void stiffness(double stiffness);

  void dimWeight(const Eigen::Vector3d & dim);

  const Eigen::Vector3d & dimWeight() const
  {
    return dimWeight_;
  }

  virtual std::pair<int, int> begin() const override
  {
    return {alphaDBegin_, alphaDBegin_};
  }

  virtual void updateNrVars(const std::vector<rbd::MultiBody> & mbs, const SolverData & data) override;
  virtual void update(const std::vector<rbd::MultiBody> & mbs,
                      const std::vector<rbd::MultiBodyConfig> & mbcs,
                      const SolverData & data) override;

  virtual const Eigen::MatrixXd & Q() const override;
  virtual const Eigen::VectorXd & C() const override;

  const Eigen::VectorXd & eval() const;
  const Eigen::VectorXd & speed() const;

private:
  void init(const std::vector<rbd::MultiBody> & mbs);

private:
  int alphaDBegin_;
  double stiffness_, stiffnessSqrt_;
  Eigen::Vector3d dimWeight_;
  std::vector<int> posInQ_;
  tasks::MultiCoMTask mct_;
  Eigen::MatrixXd Q_;
  Eigen::VectorXd C_;
  Eigen::Vector3d CSum_;
  // cache
  Eigen::MatrixXd preQ_;
};

class TASKS_DLLAPI MultiRobotTransformTask : public Task
{
public:
  MultiRobotTransformTask(const std::vector<rbd::MultiBody> & mbs,
                          int r1Index,
                          int r2Index,
                          const std::string & r1BodyName,
                          const std::string & r2BodyName,
                          const sva::PTransformd & X_r1b_r1s,
                          const sva::PTransformd & X_r2b_r2s,
                          double stiffness,
                          double weight);

  tasks::MultiRobotTransformTask & task()
  {
    return mrtt_;
  }

  void X_r1b_r1s(const sva::PTransformd & X_r1b_r1s);
  const sva::PTransformd & X_r1b_r1s() const;

  void X_r2b_r2s(const sva::PTransformd & X_r2b_r2s);
  const sva::PTransformd & X_r2b_r2s() const;

  double stiffness() const
  {
    return stiffness_;
  }

  void stiffness(double stiffness);

  void dimWeight(const Eigen::Vector6d & dim);

  const Eigen::VectorXd & dimWeight() const
  {
    return dimWeight_;
  }

  virtual std::pair<int, int> begin() const override
  {
    return {alphaDBegin_, alphaDBegin_};
  }

  virtual void updateNrVars(const std::vector<rbd::MultiBody> & mbs, const SolverData & data) override;
  virtual void update(const std::vector<rbd::MultiBody> & mbs,
                      const std::vector<rbd::MultiBodyConfig> & mbcs,
                      const SolverData & data) override;

  virtual const Eigen::MatrixXd & Q() const override;
  virtual const Eigen::VectorXd & C() const override;

  const Eigen::VectorXd & eval() const;
  const Eigen::VectorXd & speed() const;

private:
  int alphaDBegin_;
  double stiffness_, stiffnessSqrt_;
  Eigen::VectorXd dimWeight_;
  std::vector<int> posInQ_, robotIndexes_;
  tasks::MultiRobotTransformTask mrtt_;
  Eigen::MatrixXd Q_;
  Eigen::VectorXd C_;
  Eigen::VectorXd CSum_;
  // cache
  Eigen::MatrixXd preQ_;
};

class TASKS_DLLAPI MomentumTask : public HighLevelTask
{
public:
  MomentumTask(const std::vector<rbd::MultiBody> & mbs, int robotIndex, const sva::ForceVecd & mom);

  tasks::MomentumTask & task()
  {
    return momt_;
  }

  void momentum(const sva::ForceVecd & mom)
  {
    momt_.momentum(mom);
  }

  const sva::ForceVecd momentum()
  {
    return momt_.momentum();
  }

  virtual int dim() override;
  virtual void update(const std::vector<rbd::MultiBody> & mb,
                      const std::vector<rbd::MultiBodyConfig> & mbc,
                      const SolverData & data) override;

  virtual const Eigen::MatrixXd & jac() override;
  virtual const Eigen::MatrixXd & jacDot() override;
  virtual const Eigen::VectorXd & eval() override;
  virtual const Eigen::VectorXd & speed() override;
  virtual const Eigen::VectorXd & normalAcc() override;

private:
  tasks::MomentumTask momt_;
  int robotIndex_;
};

class TASKS_DLLAPI CentroidalAngularMomentumTask : public Task
{
public:
  CentroidalAngularMomentumTask(const std::vector<rbd::MultiBody>& mbs, int robotIndex,
				double gain, const Eigen::Vector3d angMomentum, double weight);

  virtual std::string nameTask() const override
  {
    return "CentroidalAngularMomentumTask";
  }
  
  virtual std::pair<int, int> begin() const
  {
    return std::make_pair(alphaDBegin_, alphaDBegin_);
  }
  
  void angMomentum(const Eigen::Vector3d & angMomentum)
  {
    angMomentum_ = angMomentum;
  }

  const Eigen::Vector3d angMomentum() const
  {
    return angMomentum_;
  }

  void setGain(double gain)
  {
    gain_ = gain;
  }

  void dimWeight(const Eigen::Vector3d & dim)
  {
    dimWeight_ = dim;
  }

  const Eigen::Vector3d & dimWeight() const
  {
    return dimWeight_;
  }

  virtual void updateNrVars(const std::vector<rbd::MultiBody> & mbs,
                            const SolverData & data);
  
  virtual void update(const std::vector<rbd::MultiBody> & mbs,
                      const std::vector<rbd::MultiBodyConfig> & mbcs,
                      const SolverData & data);
  
  virtual const Eigen::MatrixXd & Q() const
  {
    return Q_;
  }
  
  virtual const Eigen::VectorXd & C() const
  {
    return C_;
  }
  
private:
  
  int robotIndex_, alphaDBegin_;
  Eigen::Vector3d dimWeight_;

  Eigen::Vector3d angMomentum_;
  double gain_;

  rbd::CentroidalMomentumMatrix centroidalMomentumMatrix_;

  Eigen::MatrixXd Q_;
  Eigen::VectorXd C_;
  // cache
  Eigen::MatrixXd jacMat_;
  Eigen::MatrixXd preQ_;
  Eigen::Vector3d CSum_;
  Eigen::Vector3d normalAcc_;
};
 
class TASKS_DLLAPI ContactTask : public Task
{
public:
  ContactTask(ContactId contactId, double stiffness, double weight)
  : Task(weight), contactId_(contactId), begin_(0), stiffness_(stiffness), stiffnessSqrt_(2 * std::sqrt(stiffness)),
    conesJac_(), error_(Eigen::Vector3d::Zero()), errorD_(Eigen::Vector3d::Zero()), Q_(), C_()
  {
  }

  virtual std::pair<int, int> begin() const override
  {
    return std::make_pair(begin_, begin_);
  }

  void error(const Eigen::Vector3d & error);
  void errorD(const Eigen::Vector3d & errorD);

  virtual void updateNrVars(const std::vector<rbd::MultiBody> & mbs, const SolverData & data) override;
  virtual void update(const std::vector<rbd::MultiBody> & mbs,
                      const std::vector<rbd::MultiBodyConfig> & mbcs,
                      const SolverData & data) override;

  virtual const Eigen::MatrixXd & Q() const override;
  virtual const Eigen::VectorXd & C() const override;

private:
  ContactId contactId_;
  int begin_;

  double stiffness_, stiffnessSqrt_;
  Eigen::MatrixXd conesJac_;
  Eigen::Vector3d error_, errorD_;

  Eigen::MatrixXd Q_;
  Eigen::VectorXd C_;
};

class TASKS_DLLAPI GripperTorqueTask : public Task
{
public:
  GripperTorqueTask(ContactId contactId, const Eigen::Vector3d & origin, const Eigen::Vector3d & axis, double weight)
  : Task(weight), contactId_(contactId), origin_(origin), axis_(axis), begin_(0), Q_(), C_()
  {
  }

  virtual std::pair<int, int> begin() const override
  {
    return std::make_pair(begin_, begin_);
  }

  virtual void updateNrVars(const std::vector<rbd::MultiBody> & mbs, const SolverData & data) override;
  virtual void update(const std::vector<rbd::MultiBody> & mbs,
                      const std::vector<rbd::MultiBodyConfig> & mbcs,
                      const SolverData & data) override;

  virtual const Eigen::MatrixXd & Q() const override;
  virtual const Eigen::VectorXd & C() const override;

private:
  ContactId contactId_;
  Eigen::Vector3d origin_;
  Eigen::Vector3d axis_;
  int begin_;

  Eigen::MatrixXd Q_;
  Eigen::VectorXd C_;
};

class TASKS_DLLAPI LinVelocityTask : public HighLevelTask
{
public:
  LinVelocityTask(const std::vector<rbd::MultiBody> & mbs,
                  int robotIndex,
                  const std::string & bodyName,
                  const Eigen::Vector3d & vel,
                  const Eigen::Vector3d & bodyPoint = Eigen::Vector3d::Zero());

  tasks::LinVelocityTask & task()
  {
    return pt_;
  }

  void velocity(const Eigen::Vector3d & s)
  {
    pt_.velocity(s);
  }

  const Eigen::Vector3d & velocity() const
  {
    return pt_.velocity();
  }

  void bodyPoint(const Eigen::Vector3d & point)
  {
    pt_.bodyPoint(point);
  }

  const Eigen::Vector3d & bodyPoint() const
  {
    return pt_.bodyPoint();
  }

  virtual int dim() override;
  virtual void update(const std::vector<rbd::MultiBody> & mbs,
                      const std::vector<rbd::MultiBodyConfig> & mbcs,
                      const SolverData & data) override;

  virtual const Eigen::MatrixXd & jac() override;
  virtual const Eigen::MatrixXd & jacDot() override;
  virtual const Eigen::VectorXd & eval() override;
  virtual const Eigen::VectorXd & speed() override;
  virtual const Eigen::VectorXd & normalAcc() override;

private:
  tasks::LinVelocityTask pt_;
  int robotIndex_;
};

class TASKS_DLLAPI OrientationTrackingTask : public HighLevelTask
{
public:
  OrientationTrackingTask(const std::vector<rbd::MultiBody> & mbs,
                          int robotIndex,
                          const std::string & bodyName,
                          const Eigen::Vector3d & bodyPoint,
                          const Eigen::Vector3d & bodyAxis,
                          const std::vector<std::string> & trackingJointsName,
                          const Eigen::Vector3d & trackedPoint);

  tasks::OrientationTrackingTask & task()
  {
    return ott_;
  }

  void trackedPoint(const Eigen::Vector3d & tp)
  {
    ott_.trackedPoint(tp);
  }

  const Eigen::Vector3d & trackedPoint() const
  {
    return ott_.trackedPoint();
  }

  void bodyPoint(const Eigen::Vector3d & bp)
  {
    ott_.bodyPoint(bp);
  }

  const Eigen::Vector3d & bodyPoint() const
  {
    return ott_.bodyPoint();
  }

  void bodyAxis(const Eigen::Vector3d & ba)
  {
    ott_.bodyAxis(ba);
  }

  const Eigen::Vector3d & bodyAxis() const
  {
    return ott_.bodyAxis();
  }

  virtual int dim() override;
  virtual void update(const std::vector<rbd::MultiBody> & mbs,
                      const std::vector<rbd::MultiBodyConfig> & mbcs,
                      const SolverData & data) override;

  virtual const Eigen::MatrixXd & jac() override;
  virtual const Eigen::MatrixXd & jacDot() override;
  virtual const Eigen::VectorXd & eval() override;
  virtual const Eigen::VectorXd & speed() override;
  virtual const Eigen::VectorXd & normalAcc() override;

private:
  int robotIndex_;
  tasks::OrientationTrackingTask ott_;
  Eigen::VectorXd alphaVec_;
  Eigen::VectorXd speed_, normalAcc_;
};

class TASKS_DLLAPI RelativeDistTask : public HighLevelTask
{
public:
  RelativeDistTask(const std::vector<rbd::MultiBody> & mbs,
                   const int rIndex,
                   const double timestep,
                   tasks::RelativeDistTask::rbInfo & rbi1,
                   tasks::RelativeDistTask::rbInfo & rbi2,
                   const Eigen::Vector3d & u1 = Eigen::Vector3d::Zero(),
                   const Eigen::Vector3d & u2 = Eigen::Vector3d::Zero());

  tasks::RelativeDistTask & task()
  {
    return rdt_;
  }

  void robotPoint(const rbd::MultiBody & mb, const std::string & bName, const Eigen::Vector3d & point)
  {
    int bIndex = mb.bodyIndexByName(bName);
    rdt_.robotPoint(bIndex, point);
  }
  void envPoint(const rbd::MultiBody & mb, const std::string & bName, const Eigen::Vector3d & point)
  {
    int bIndex = mb.bodyIndexByName(bName);
    rdt_.envPoint(bIndex, point);
  }
  void vector(const rbd::MultiBody & mb, const std::string & bName, const Eigen::Vector3d & u)
  {
    int bIndex = mb.bodyIndexByName(bName);
    rdt_.vector(bIndex, u);
  }

  virtual int dim() override;
  virtual void update(const std::vector<rbd::MultiBody> & mbs,
                      const std::vector<rbd::MultiBodyConfig> & mbcs,
                      const SolverData & data) override;

  virtual const Eigen::MatrixXd & jac() override;
  virtual const Eigen::VectorXd & eval() override;
  virtual const Eigen::VectorXd & speed() override;
  virtual const Eigen::VectorXd & normalAcc() override;

private:
  int rIndex_;
  tasks::RelativeDistTask rdt_;
};

class TASKS_DLLAPI VectorOrientationTask : public HighLevelTask
{
public:
  VectorOrientationTask(const std::vector<rbd::MultiBody> & mbs,
                        int robotIndex,
                        const std::string & bodyName,
                        const Eigen::Vector3d & bodyVector,
                        const Eigen::Vector3d & targetVector);

  tasks::VectorOrientationTask & task()
  {
    return vot_;
  }
  void bodyVector(const Eigen::Vector3d & vector)
  {
    vot_.bodyVector(vector);
  }
  const Eigen::Vector3d & bodyVector() const
  {
    return vot_.bodyVector();
  }
  void target(const Eigen::Vector3d & vector)
  {
    vot_.target(vector);
  }
  const Eigen::Vector3d & target() const
  {
    return vot_.target();
  }
  const Eigen::Vector3d & actual() const
  {
    return vot_.actual();
  }

  virtual int dim() override;
  virtual void update(const std::vector<rbd::MultiBody> & mbs,
                      const std::vector<rbd::MultiBodyConfig> & mbcs,
                      const SolverData & data) override;

  virtual const Eigen::MatrixXd & jac() override;
  virtual const Eigen::MatrixXd & jacDot() override;  
  virtual const Eigen::VectorXd & eval() override;
  virtual const Eigen::VectorXd & speed() override;
  virtual const Eigen::VectorXd & normalAcc() override;

private:
  tasks::VectorOrientationTask vot_;
  int robotIndex_;
};

class TASKS_DLLAPI WrenchTask : public Task
{
 public:
  WrenchTask(const std::vector<rbd::MultiBody> & mbs, int robotIndex, const std::string & bodyName,
	     const Eigen::Vector3d & bodyPoint, double weight);

  virtual std::string nameTask() const override
  {
    return "WrenchTask";
  }
  
  virtual std::pair<int, int> begin() const
  {
    return std::make_pair(lambdaBegin_, lambdaBegin_);
  }
  
  void bodyPoint(const Eigen::Vector3d & point)
  {
    bodyPoint_ = point;
  }
  
  const Eigen::Vector3d & bodyPoint() const
  {
    return bodyPoint_;
  }
  
  void treatDesWrenchAsLocal(bool local)
  {
    local_ = local;
  }
  
  bool treatSettingsAsLocal()
  {
    return local_;
  }
  
  void wrench(const std::vector<rbd::MultiBodyConfig> & mbcs, const sva::ForceVecd & wrench);
  
  const sva::ForceVecd & wrench() const
  {
    return wrench_;
  }
  
  void dimWeight(const std::vector<rbd::MultiBodyConfig> & mbcs, const Eigen::Vector6d & dim);
  
  const Eigen::Vector6d & dimWeight() const
  {
    return dimWeight_;
  }
  
  virtual void updateNrVars(const std::vector<rbd::MultiBody> & mbs,
			    const SolverData& data);
  
  virtual void update(const std::vector<rbd::MultiBody> & mbs,
		      const std::vector<rbd::MultiBodyConfig> & mbcs,
		      const SolverData & data);
  
  virtual const Eigen::MatrixXd & Q() const
  {
    return Q_;
  }
  
  virtual const Eigen::VectorXd & C() const
  {
    return C_;
  }

 private:
  
  int robotIndex_, bodyIndex_, lambdaBegin_;
  Eigen::Vector6d dimWeight_;
  
  Eigen::Vector3d bodyPoint_;
  sva::ForceVecd wrench_;
  
  Eigen::MatrixXd W_;
  
  // the following flag indicates if the desired values are expressed or not in the local frame
  bool local_;
  
  Eigen::MatrixXd Q_;
  Eigen::VectorXd C_;
  // cache
  Eigen::MatrixXd preQ_;
  Eigen::MatrixXd preC_;
};

class TASKS_DLLAPI LocalCoPTask : public Task
{
 public:
  LocalCoPTask(const std::vector<rbd::MultiBody> & mbs, int robotIndex,
               const std::string & bodyName, const Eigen::Vector3d & localCoP,
               double weight);
  
  virtual std::string nameTask() const override
  {
    return "LocalCoPTask";
  }
  
  virtual std::pair<int, int> begin() const
  {
    return std::make_pair(lambdaBegin_, lambdaBegin_);
  }

  void localCoP(const Eigen::Vector3d & point)
  {
    localCoP_ = point;
  }
  
  const Eigen::Vector3d & localCoP() const
  {
    return localCoP_;
  }

  void dimWeight(const std::vector<rbd::MultiBodyConfig> & mbcs, const Eigen::Vector3d & dim);
  
  const Eigen::Vector3d & dimWeight() const
  {
    return dimWeight_;
  }

  virtual void updateNrVars(const std::vector<rbd::MultiBody> & mbs,
			    const SolverData& data);
  
  virtual void update(const std::vector<rbd::MultiBody> & mbs,
		      const std::vector<rbd::MultiBodyConfig> & mbcs,
		      const SolverData & data);
  
  virtual const Eigen::MatrixXd & Q() const
  {
    return Q_;
  }
  
  virtual const Eigen::VectorXd & C() const
  {
    return C_;
  }

 private:

  int robotIndex_, bodyIndex_, lambdaBegin_;
  Eigen::Vector3d dimWeight_;
  
  Eigen::Vector3d localCoP_;
  
  Eigen::MatrixXd W_;
  
  Eigen::MatrixXd Q_;
  Eigen::VectorXd C_;
  // cache
  Eigen::MatrixXd preQ_;
};

class TASKS_DLLAPI AdmittanceTaskCommon : public Task
{
public:
  AdmittanceTaskCommon(const std::vector<rbd::MultiBody> & mbs, int robotIndex,
		       const std::string & bodyName,
		       const Eigen::Vector3d & bodyPoint,
		       double timeStep, double gainForceP, double gainForceD,
		       double gainCoupleP, double gainCoupleD, double weight);

  virtual std::pair<int, int> begin() const
  {
    return std::make_pair(alphaDBegin_, alphaDBegin_);
  }
  
  void bodyPoint(const Eigen::Vector3d & point)
  {
    jac_.point(point);
  }
  
  const Eigen::Vector3d & bodyPoint() const
  {
    return jac_.point();
  }
  
  void treatSettingsAsLocal(bool local)
  {
    local_ = local;
  }

  void setForceGains(double gainP, double gainD)
  {
    gainForceP_ = gainP;
    gainForceD_ = gainD;
  }
  
  void setCoupleGains(double gainP, double gainD)
  {
    gainCoupleP_ = gainP;
    gainCoupleD_ = gainD;
  }

  void dimWeight(const std::vector<rbd::MultiBodyConfig> & mbcs,
		 const Eigen::Vector6d & dim);
  
  const Eigen::Vector6d & dimWeight() const
  {
    return dimWeight_;
  }
  
  virtual void updateNrVars(const std::vector<rbd::MultiBody> & mbs,
			    const SolverData & data);
  
  virtual void update(const std::vector<rbd::MultiBody> & mbs,
		      const std::vector<rbd::MultiBodyConfig> & mbcs,
		      const SolverData & data) = 0;
  
  virtual const Eigen::MatrixXd & Q() const
  {
    return Q_;
  }
  
  virtual const Eigen::VectorXd & C() const
  {
    return C_;
  }

  struct BodyLambda
  {
    std::string body;
    int begin;
    int size;
  };

protected:

  sva::ForceVecd computeWrench(const rbd::MultiBodyConfig & mbc,
			       const tasks::qp::BilateralContact & contact,
			       Eigen::VectorXd lambdaVec, int pos);

  int robotIndex_, bodyIndex_, alphaDBegin_;
  std::vector<BodyLambda> contactBodies_, contactBodiesPrev_; 
  double dt_;
  double gainForceP_, gainForceD_;
  double gainCoupleP_, gainCoupleD_;

  rbd::Jacobian jac_;
  Eigen::MatrixXd jacMat_;
  Eigen::Vector6d error_;
  Eigen::Vector6d normalAcc_;
  Eigen::Vector6d dimWeight_;
  
  // the following flag indicates if dimWeight values are related or not to the local frame
  bool local_;
  
  Eigen::MatrixXd Q_;
  Eigen::VectorXd C_;
  // cache
  Eigen::MatrixXd preQ_;
  Eigen::Vector6d preC_;
};
 
class TASKS_DLLAPI AdmittanceTask : public AdmittanceTaskCommon
{
public:
  AdmittanceTask(const std::vector<rbd::MultiBody> & mbs, int robotIndex,
		 const std::string & bodyName,
		 const Eigen::Vector3d & bodyPoint,
		 double timeStep, double gainForceP, double gainForceD,
		 double gainCoupleP, double gainCoupleD, double weight);

  virtual std::string nameTask() const override
  {
    return "AdmittanceTask";
  }
  
  void measuredWrench(const sva::ForceVecd wrench)
  {
    measuredWrench_ = wrench;
  }
  
  const sva::ForceVecd & measuredWrench() const
  {
    return measuredWrench_;
  }
  
  const sva::ForceVecd & calculatedWrench() const
  {
    return calculatedWrench_;
  }
  
  virtual void update(const std::vector<rbd::MultiBody> & mbs,
		      const std::vector<rbd::MultiBodyConfig> & mbcs,
		      const SolverData & data);
  
private:
  
  sva::ForceVecd measuredWrench_, measuredWrenchPrev_;
  sva::ForceVecd calculatedWrench_, calculatedWrenchPrev_;
};

class TASKS_DLLAPI NullSpaceAdmittanceTask : public AdmittanceTaskCommon
{
public:
  NullSpaceAdmittanceTask(const std::vector<rbd::MultiBody> & mbs, int robotIndex,
			  const std::string & bodyName,
			  const Eigen::Vector3d & bodyPoint,
			  double timeStep, double gainForceP, double gainForceD,
			  double gainCoupleP, double gainCoupleD, double weight);

  virtual std::string nameTask() const override
  {
    return "NullSpaceAdmittanceTask";
  }
  
  void measuredWrench(const std::string & bodyName, const sva::ForceVecd & wrench);
  
  void measuredWrenches(const std::map<std::string, sva::ForceVecd> & wrenches);

  const sva::ForceVecd & measuredWrench(const std::string & bodyName) const
  {
    return measuredWrenches_.at(bodyName);
  }

  const std::map<std::string, sva::ForceVecd> & measuredWrenches() const
  {
    return measuredWrenches_;
  }

  const Eigen::Vector3d & calculatedForce(const std::string & bodyName) const
  {
    return calculatedForces_.at(bodyName);
  }

  const std::map<std::string, Eigen::Vector3d> & calculatedForces() const
  {
    return calculatedForces_;
  }

  const sva::ForceVecd & calculatedWrench() const
  {
    return calculatedBodyWrench_;
  }
  
  void fdistRatio(const Eigen::Vector3d & ratio)
  {
    fdistRatio_ = ratio;
  }

  const Eigen::Vector3d & fdistRatio() const
  {
    return fdistRatio_;
  }
  
  virtual void updateNrVars(const std::vector<rbd::MultiBody> & mbs,
			    const SolverData & data) override;
  
  virtual void update(const std::vector<rbd::MultiBody> & mbs,
		      const std::vector<rbd::MultiBodyConfig> & mbcs,
		      const SolverData & data);

private:

  int nrBodies_;

  // std::set<std::string> bodies_;
  
  std::map<std::string, sva::ForceVecd> measuredWrenches_;
  std::map<std::string, Eigen::Vector3d> calculatedForces_;

  sva::ForceVecd calculatedBodyWrench_;
  sva::ForceVecd measuredBodyWrenchPrev_;
  sva::ForceVecd calculatedBodyWrenchPrev_;
  
  Eigen::Vector3d fdistRatio_;

  Eigen::Vector3d projForceErr_;
};

class TASKS_DLLAPI ForceDistributionTaskCommon : public Task
{
 public:
  ForceDistributionTaskCommon(const std::vector<rbd::MultiBody> & mbs,
			      int robotIndex, double weight);

  virtual std::pair<int, int> begin() const = 0;

  void fdistRatio(const std::string & bodyName, const Eigen::Vector3d & ratio)
  {
    //fdistRatios_.at(bodyName) = ratio;
    fdistRatios_[bodyName] = ratio;
  }
  
  void fdistRatios(const std::map<std::string, Eigen::Vector3d> & ratios);

  Eigen::Vector3d fdistRatio(const std::string & bodyName) const;
  
  const std::map<std::string, Eigen::Vector3d> & fdistRatios() const
  {
    return fdistRatios_;
  }

  Eigen::Vector3d refForce(const std::string & bodyName) const;
  
  const std::map<std::string, Eigen::Vector3d> & refForces() const
  {
    return refForces_;
  }

  virtual void updateNrVars(const std::vector<rbd::MultiBody> & mbs,
			    const SolverData& data);
  
  virtual void update(const std::vector<rbd::MultiBody> & mbs,
		      const std::vector<rbd::MultiBodyConfig> & mbcs,
		      const SolverData & data);

  virtual const Eigen::MatrixXd & Q() const
  {
    return Q_;
  }
  
  virtual const Eigen::VectorXd & C() const
  {
    return C_;
  }

  struct BodyLambda
  {
    std::string body;
    int begin;
    int size;
  };
  
 protected:

  int robotIndex_, lambdaBegin_, nrBodies_;
  std::vector<BodyLambda> contactBodies_, contactBodiesPrev_; 

  std::map<std::string, Eigen::Vector3d> fdistRatios_;
  std::map<std::string, Eigen::Vector3d> refForces_;

  Eigen::MatrixXd W_;
  Eigen::MatrixXd fdistRatioMat_;
  Eigen::MatrixXd A_;
    
  Eigen::MatrixXd Q_;
  Eigen::VectorXd C_;
};
 
class TASKS_DLLAPI ForceDistributionTaskOriginal : public ForceDistributionTaskCommon
{
public:
  ForceDistributionTaskOriginal(const std::vector<rbd::MultiBody> & mbs,
				int robotIndex, double weight);

  virtual std::string nameTask() const override
  {
    return "ForceDistributionTaskOriginal";
  }
  
  virtual std::pair<int, int> begin() const
  {
    return std::make_pair(0, 0);
  }

  virtual void updateNrVars(const std::vector<rbd::MultiBody> & mbs,
			    const SolverData& data) override;
  
  virtual void update(const std::vector<rbd::MultiBody> & mbs,
		      const std::vector<rbd::MultiBodyConfig> & mbcs,
		      const SolverData & data) override;
  
private:

  int alphaDBegin_;
  
  double gAcc_;
  double totalMass_;
  
  rbd::CoMJacobian comJac_;

  Eigen::MatrixXd comJacMat_;
  Eigen::Vector3d normalAcc_;

  // cache
  Eigen::VectorXd CSum_;
};

class TASKS_DLLAPI ForceDistributionTaskOptimized : public ForceDistributionTaskCommon
{
 public:
  ForceDistributionTaskOptimized(const std::vector<rbd::MultiBody> & mbs,
				 int robotIndex, double weight);
  
  virtual std::string nameTask() const override
  {
    return "ForceDistributionTaskOriginal";
  }

  virtual std::pair<int, int> begin() const
  {
    return std::make_pair(lambdaBegin_, lambdaBegin_);
  }

  virtual void updateNrVars(const std::vector<rbd::MultiBody> & mbs,
			    const SolverData& data) override;
  
  virtual void update(const std::vector<rbd::MultiBody> & mbs,
		      const std::vector<rbd::MultiBodyConfig> & mbcs,
		      const SolverData & data) override;

 protected:

  // cache
  Eigen::MatrixXd SumMat_;
  Eigen::MatrixXd preA_;
};

class TASKS_DLLAPI ZMPBasedCoMTask : public Task
{
 public:
  ZMPBasedCoMTask(const std::vector<rbd::MultiBody> & mbs, int robotIndex,
		  const Eigen::Vector3d & com, const Eigen::Vector3d & zmp,
		  double weight);

  virtual std::string nameTask() const override
  {
    return "ZMPBasedCoMTask";
  }
  
  virtual std::pair<int, int> begin() const
  {
    return std::make_pair(alphaDBegin_, alphaDBegin_);
  }

  void com(const Eigen::Vector3d & com)
  {
    com_ = com;
  }

  const Eigen::Vector3d & com() const
  {
    return com_;
  }
  
  void zmp(const Eigen::Vector3d & zmp)
  {
    zmp_ = zmp;
  }

  const Eigen::Vector3d & zmp() const
  {
    return zmp_;
  }

  void ddcom(const Eigen::Vector3d & ddcom)
  {
    ddcom_ = ddcom;
  }

  const Eigen::Vector3d & ddcom() const
  {
    return ddcom_;
  }

  void dimWeight(const Eigen::Vector3d & dim)
  {
    dimWeight_ = dim;
  }

  const Eigen::Vector3d & dimWeight() const
  {
    return dimWeight_;
  }

  virtual void updateNrVars(const std::vector<rbd::MultiBody> & mbs,
			    const SolverData& data);
  
  virtual void update(const std::vector<rbd::MultiBody> & mbs,
		      const std::vector<rbd::MultiBodyConfig> & mbcs,
		      const SolverData & data);
  
  virtual const Eigen::MatrixXd & Q() const
  {
    return Q_;
  }

  virtual const Eigen::VectorXd & C() const
  {
    return C_;
  }

 private:

  int robotIndex_, alphaDBegin_;

  Eigen::Vector3d com_, ddcom_, zmp_;
  double gAcc_;
  
  Eigen::Vector3d dimWeight_;

  rbd::CoMJacobian jac_;
  
  Eigen::MatrixXd Q_;
  Eigen::VectorXd C_;
  // cache
  Eigen::MatrixXd jacMat_;
  Eigen::MatrixXd preQ_;
  Eigen::Vector3d CSum_;
  Eigen::Vector3d normalAcc_;
};

class TASKS_DLLAPI ZMPTask : public Task
{
 public:
  ZMPTask(const std::vector<rbd::MultiBody> & mbs, int robotIndex,
	  const Eigen::Vector3d & zmp, double weight);

  virtual std::string nameTask() const override
  {
    return "ZMPTask";
  }
  
  virtual std::pair<int, int> begin() const
  {
    return std::make_pair(lambdaBegin_, lambdaBegin_);
  }

  void zmp(const Eigen::Vector3d & zmp)
  {
    zmp_ = zmp;
  }

  const Eigen::Vector3d & zmp() const
  {
    return zmp_;
  }

  const Eigen::Vector3d & totalForce() const
  {
    return totalForce_;
  }

  const Eigen::Vector3d & totalMomentZMP() const
  {
    return totalMomentZMP_;
  }

  void dimWeight(const Eigen::Vector3d & dim)
  {
    dimWeight_ = dim;
  }

  const Eigen::Vector3d & dimWeight() const
  {
    return dimWeight_;
  }

  virtual void updateNrVars(const std::vector<rbd::MultiBody> & mbs,
			    const SolverData& data);
  
  virtual void update(const std::vector<rbd::MultiBody> & mbs,
		      const std::vector<rbd::MultiBodyConfig> & mbcs,
		      const SolverData & data);
  
  virtual const Eigen::MatrixXd & Q() const
  {
    return Q_;
  }

  virtual const Eigen::VectorXd & C() const
  {
    return C_;
  }

  struct BodyLambda
  {
    std::string body;
    int begin;
    int size;
  };

 private:

  int robotIndex_, lambdaBegin_, nrBodies_;
  std::vector<BodyLambda> contactBodies_, contactBodiesPrev_; 
  Eigen::Vector3d dimWeight_;

  Eigen::Vector3d zmp_;
  
  Eigen::Vector3d totalForce_;
  Eigen::Vector3d totalMomentZMP_;

  Eigen::MatrixXd W_;
  Eigen::MatrixXd pW_;

  Eigen::MatrixXd Q_;
  Eigen::VectorXd C_;
  // cache
  Eigen::MatrixXd preQ_;
};

// Pending to finish...
 
class TASKS_DLLAPI YawMomentCompensationTask : public Task
{
 public:
  YawMomentCompensationTask(const std::vector<rbd::MultiBody>& mbs, int robotIndex, double timeStep,
                            double gainKp, double gainKd, double gainKi, double gainKii, double weight);

  virtual std::string nameTask() const override
  {
    return "YawMomentCompensationTask";
  }
  
  virtual std::pair<int, int> begin() const
  {
    return std::make_pair(alphaDBegin_, alphaDBegin_);
  }
  
  void deltaTauP(const double delta_tau_p)
  {
    delta_tau_p_ = delta_tau_p;
  }
  
  double deltaTauP()
  {
    return delta_tau_p_;
  }
  
  void zmp(const Eigen::Vector3d p)
  {
    zmp_ = p;
  }
  
  const Eigen::Vector3d zmp() const
  {
    return zmp_;
  }
  
  void setGains(double gainKp, double gainKd, double gainKi, double gainKii)
  {
    gainKp_  = gainKp;
    gainKd_  = gainKd;
    gainKi_  = gainKi;
    gainKii_ = gainKii;
  }
  
  virtual void updateNrVars(const std::vector<rbd::MultiBody>& mbs,
                            const SolverData& data);
  
  virtual void update(const std::vector<rbd::MultiBody>& mbs,
                      const std::vector<rbd::MultiBodyConfig>& mbcs,
                      const SolverData& data);
  
  virtual const Eigen::MatrixXd& Q() const
  {
    return Q_;
  }
  
  virtual const Eigen::VectorXd& C() const
  {
    return C_;
  }
  
 private:

  double computeTauP();
  sva::ForceVecd computeWrench(const rbd::MultiBodyConfig& mbc,
                               const tasks::qp::BilateralContact& contact, Eigen::VectorXd lambdaVec, int pos);
  Eigen::Matrix3d skewMatrix(const Eigen::Vector3d& v);
        
  int robotIndex_, alphaDBegin_;
  
  double dt_;
  double gainKp_, gainKd_, gainKi_, gainKii_;
        
  Eigen::Vector3d com_prev_;
  Eigen::Vector3d zmp_, zmp_prev_;
  double delta_tau_p_, delta_tau_p_prev_;
  double Lpz_, iLpz_;
  
  rbd::CentroidalMomentumMatrix centroidalMomentumMatrix_;
  Eigen::MatrixXd jacMat_;
  
  Eigen::MatrixXd Q_;
  Eigen::VectorXd C_;
};
 
} // namespace qp

} // namespace tasks
