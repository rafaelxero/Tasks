/*
 * Copyright 2012-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#pragma once

// includes
// std
#include <memory>
#include <vector>

// boost
#include <boost/timer/timer.hpp>

// Eigen
#include <Eigen/Core>

// Tasks
#include <tasks/config.hh>

#include "QPContacts.h"
#include "QPSolverData.h"

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
class Constraint;
class Equality;
class Inequality;
class GenInequality;
class Bound;
class Task;
class GenQPSolver;

// class MotionConstr;
          
class TASKS_DLLAPI QPSolver
{
public:
  QPSolver();
  virtual ~QPSolver();  // declared as virtual for QPSolver to be polymorphic

  /** solve the problem
   *  \param mbs current multibody
   *  \param mbcs current state of the multibody and result of the solved problem
   */
  bool solve(const std::vector<rbd::MultiBody> & mbs, std::vector<rbd::MultiBodyConfig> & mbcs);

  /** solve the problem but don't fill mbc.
   * This is usefull for the python binding.
   *  \param mbs current multibody
   *  \param mbc current state of the multibody
   */
  bool solveNoMbcUpdate(const std::vector<rbd::MultiBody> & mbs, const std::vector<rbd::MultiBodyConfig> & mbcs);

  /** \brief fill mbc with the solution of the problem
   * \param mbc reslut of the solved problem
   * \param robotIndex robot index associated with mbc
   */
  void updateMbc(rbd::MultiBodyConfig & mbc, int robotIndex) const;

  void updateConstrSize();

  void nrVars(const std::vector<rbd::MultiBody> & mbs,
              std::vector<UnilateralContact> uni,
              std::vector<BilateralContact> bi);
  int nrVars() const;

  /// call updateNrVars on all tasks
  void updateTasksNrVars(const std::vector<rbd::MultiBody> & mbs) const;
  /// call updateNrVars on all constraints
  void updateConstrsNrVars(const std::vector<rbd::MultiBody> & mbs) const;
  /// call updateNrVars on all tasks and constraints
  void updateNrVars(const std::vector<rbd::MultiBody> & mbs) const;

  bool hasConstraint(const Constraint* co);
  // std::shared_ptr<tasks::qp::MotionConstr> getMotionConstr();

  void addEqualityConstraint(Equality * co);
  void removeEqualityConstraint(Equality * co);
  int nrEqualityConstraints() const;

  void addInequalityConstraint(Inequality * co);
  void removeInequalityConstraint(Inequality * co);
  int nrInequalityConstraints() const;

  void addGenInequalityConstraint(GenInequality * co);
  void removeGenInequalityConstraint(GenInequality * co);
  int nrGenInequalityConstraints() const;

  void addBoundConstraint(Bound * co);
  void removeBoundConstraint(Bound * co);
  int nrBoundConstraints() const;

  void addConstraint(Constraint * co);
  void addConstraint(const std::vector<rbd::MultiBody> & mbs, Constraint * co);
  void removeConstraint(Constraint * co);
  int nrConstraints() const;

  void addTask(Task * task);
  void addTask(const std::vector<rbd::MultiBody> & mbs, Task * task);
  void removeTask(Task * task);
  void resetTasks();
  int nrTasks() const;

  const std::vector<Task *> & getTasks() const { return tasks_; }
  const std::vector<Equality *> & getEqConstr() const { return eqConstr_; }
  const std::vector<Inequality *> & getInEqConstr() const { return inEqConstr_; }
  const std::vector<GenInequality *> & getGenInEqConstr() const { return genInEqConstr_; }
  const std::vector<Bound *> & getBoundConstr() const { return boundConstr_; }
  
  void solver(const std::string & name);
  std::string solver() const;

  const SolverData & data() const;
  SolverData & data();

  const Eigen::VectorXd & result() const;
  Eigen::VectorXd alphaDVec() const;
  Eigen::VectorXd alphaDVec(int rIndex) const;

  Eigen::VectorXd lambdaVec() const;
  Eigen::VectorXd lambdaVec(int cIndex) const;

  int contactLambdaPosition(const ContactId & cId) const;

  boost::timer::cpu_times solveTime() const;
  boost::timer::cpu_times solveAndBuildTime() const;

protected:
  void preUpdate(const std::vector<rbd::MultiBody> & mbs, const std::vector<rbd::MultiBodyConfig> & mbcs);
  void postUpdate(const std::vector<rbd::MultiBody> & mbs, std::vector<rbd::MultiBodyConfig> & mbcs, bool success);

  std::vector<Constraint *> constr_;
  std::vector<Equality *> eqConstr_;
  std::vector<Inequality *> inEqConstr_;
  std::vector<GenInequality *> genInEqConstr_;
  std::vector<Bound *> boundConstr_;

  std::vector<Task *> tasks_;

  SolverData data_;

  int maxEqLines_, maxInEqLines_, maxGenInEqLines_;

  std::unique_ptr<GenQPSolver> solver_;

  boost::timer::cpu_timer solverTimer_, solverAndBuildTimer_;
};

class TASKS_DLLAPI PassivityPIDTerm_QPSolver : public QPSolver
{
public:
  PassivityPIDTerm_QPSolver();
  
  bool solve(const std::vector<rbd::MultiBody> & mbs,
             std::vector<rbd::MultiBodyConfig> & mbcs_real,
             std::vector<rbd::MultiBodyConfig> & mbcs_calc);
        
  bool solveNoMbcUpdate(const std::vector<rbd::MultiBody> & mbs,
                        const std::vector<rbd::MultiBodyConfig> & mbcs_real,
                        const std::vector<rbd::MultiBodyConfig> & mbcs_calc);

 protected:
  void preUpdate(const std::vector<rbd::MultiBody> & mbs,
                 const std::vector<rbd::MultiBodyConfig> & mbcs_real,
                 const std::vector<rbd::MultiBodyConfig> & mbcs_calc);
};

class TASKS_DLLAPI Constraint
{
public:
  virtual ~Constraint() {}
  virtual void updateNrVars(const std::vector<rbd::MultiBody> & msb, const SolverData & data) = 0;

  virtual void update(const std::vector<rbd::MultiBody> & mbs,
                      const std::vector<rbd::MultiBodyConfig> & mbcs,
                      const SolverData & data) = 0;
};

template<typename... Fun>
class ConstraintFunction : public Constraint, public Fun...
{
public:
  virtual ~ConstraintFunction() override {}

  void addToSolver(QPSolver & sol)
  {
    sol.addConstraint(this);
    addToSolver_p(sol, (&Fun::addToSolver)...);
  }

  void addToSolver(const std::vector<rbd::MultiBody> & mbs, QPSolver & sol)
  {
    sol.addConstraint(mbs, this);
    addToSolver_p(sol, (&Fun::addToSolver)...);
  }

  void removeFromSolver(QPSolver & sol)
  {
    sol.removeConstraint(this);
    removeFromSolver_p(sol, (&Fun::removeFromSolver)...);
  }

private:
  void addToSolver_p(QPSolver & /* sol */) {}

  template<typename F, typename... NextFun>
  void addToSolver_p(QPSolver & sol, F function, NextFun... nFun)
  {
    (this->*function)(sol);
    addToSolver_p(sol, nFun...);
  }

  void removeFromSolver_p(QPSolver & /* sol */) {}

  template<typename F, typename... NextFun>
  void removeFromSolver_p(QPSolver & sol, F function, NextFun... nFun)
  {
    (this->*function)(sol);
    removeFromSolver_p(sol, nFun...);
  }
};

class TASKS_DLLAPI Equality
{
public:
  virtual ~Equality() {}
  virtual int maxEq() const = 0;
  virtual int nrEq() const
  {
    return maxEq();
  }

  virtual const Eigen::MatrixXd & AEq() const = 0;
  virtual const Eigen::VectorXd & bEq() const = 0;

  virtual std::string nameEq() const = 0;
  virtual std::string descEq(const std::vector<rbd::MultiBody> & mbs, int i) = 0;

  void addToSolver(QPSolver & sol)
  {
    sol.addEqualityConstraint(this);
  }

  void removeFromSolver(QPSolver & sol)
  {
    sol.removeEqualityConstraint(this);
  }
};

class TASKS_DLLAPI Inequality
{
public:
  virtual ~Inequality() {}
  virtual int maxInEq() const = 0;
  virtual int nrInEq() const
  {
    return maxInEq();
  }

  virtual const Eigen::MatrixXd & AInEq() const = 0;
  virtual const Eigen::VectorXd & bInEq() const = 0;

  virtual std::string nameInEq() const = 0;
  virtual std::string descInEq(const std::vector<rbd::MultiBody> & mbs, int i) = 0;

  void addToSolver(QPSolver & sol)
  {
    sol.addInequalityConstraint(this);
  }

  void removeFromSolver(QPSolver & sol)
  {
    sol.removeInequalityConstraint(this);
  }
};

class TASKS_DLLAPI GenInequality
{
public:
  virtual ~GenInequality() {}
  virtual int maxGenInEq() const = 0;
  virtual int nrGenInEq() const
  {
    return maxGenInEq();
  }

  virtual const Eigen::MatrixXd & AGenInEq() const = 0;
  virtual const Eigen::VectorXd & LowerGenInEq() const = 0;
  virtual const Eigen::VectorXd & UpperGenInEq() const = 0;

  virtual std::string nameGenInEq() const = 0;
  virtual std::string descGenInEq(const std::vector<rbd::MultiBody> & mbs, int i) = 0;

  void addToSolver(QPSolver & sol)
  {
    sol.addGenInequalityConstraint(this);
  }

  void removeFromSolver(QPSolver & sol)
  {
    sol.removeGenInequalityConstraint(this);
  }
};

class TASKS_DLLAPI Bound
{
public:
  virtual ~Bound() {}
  virtual int beginVar() const = 0;

  virtual const Eigen::VectorXd & Lower() const = 0;
  virtual const Eigen::VectorXd & Upper() const = 0;

  virtual std::string nameBound() const = 0;
  virtual std::string descBound(const std::vector<rbd::MultiBody> & mbs, int i) = 0;

  void addToSolver(QPSolver & sol)
  {
    sol.addBoundConstraint(this);
  }

  void removeFromSolver(QPSolver & sol)
  {
    sol.removeBoundConstraint(this);
  }
};

class TASKS_DLLAPI Task
{
public:
  Task(double weight) : weight_(weight) {}
  virtual ~Task() {}

  virtual std::string nameTask() const
  {
    return "Task";
  }
  
  virtual double weight() const
  {
    return weight_;
  }

  virtual void weight(double w)
  {
    weight_ = w;
  }

  virtual std::pair<int, int> begin() const = 0;

  virtual void updateNrVars(const std::vector<rbd::MultiBody> & mbs, const SolverData & data) = 0;
  virtual void update(const std::vector<rbd::MultiBody> & mbs,
                      const std::vector<rbd::MultiBodyConfig> & mbcs,
                      const SolverData & data) = 0;

  virtual const Eigen::MatrixXd & Q() const = 0;
  virtual const Eigen::VectorXd & C() const = 0;

private:
  double weight_;
};

class TASKS_DLLAPI HighLevelTask
{
public:
  virtual ~HighLevelTask() {}

  virtual std::string nameHighLevelTask() const
  {
    return "HighLevelTask";
  }

  virtual int dim() = 0;

  virtual void update(const std::vector<rbd::MultiBody> & mbs,
                      const std::vector<rbd::MultiBodyConfig> & mbcs,
                      const SolverData & data) = 0;

  virtual const Eigen::MatrixXd & jac() = 0;
  virtual const Eigen::MatrixXd & jacDot() = 0;
  virtual const Eigen::VectorXd & eval() = 0;
  virtual const Eigen::VectorXd & speed() = 0;
  virtual const Eigen::VectorXd & normalAcc() = 0;
};

template<typename T>
struct constr_traits
{
};

template<>
struct constr_traits<Equality>
{
  static int maxLines(const Equality * constr)
  {
    return constr->maxEq();
  }

  static int nrLines(const Equality * constr)
  {
    return constr->nrEq();
  }

  static std::string name(const Equality * constr)
  {
    return constr->nameEq();
  }

  static std::string desc(Equality * constr, const std::vector<rbd::MultiBody> & mbs, int i)
  {
    return constr->descEq(mbs, i);
  }
};

template<>
struct constr_traits<Inequality>
{
  static int maxLines(const Inequality * constr)
  {
    return constr->maxInEq();
  }

  static int nrLines(const Inequality * constr)
  {
    return constr->nrInEq();
  }

  static std::string name(const Inequality * constr)
  {
    return constr->nameInEq();
  }

  static std::string desc(Inequality * constr, const std::vector<rbd::MultiBody> & mbs, int i)
  {
    return constr->descInEq(mbs, i);
  }
};

template<>
struct constr_traits<GenInequality>
{
  static int maxLines(const GenInequality * constr)
  {
    return constr->maxGenInEq();
  }

  static int nrLines(const GenInequality * constr)
  {
    return constr->nrGenInEq();
  }

  static std::string name(const GenInequality * constr)
  {
    return constr->nameGenInEq();
  }

  static std::string desc(GenInequality * constr, const std::vector<rbd::MultiBody> & mbs, int i)
  {
    return constr->descGenInEq(mbs, i);
  }
};

} // namespace qp

} // namespace tasks
