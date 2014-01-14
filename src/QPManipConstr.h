// This file is part of Tasks.
//
// Tasks is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Tasks is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with Tasks.  If not, see <http://www.gnu.org/licenses/>.

#pragma once

// includes
// Eigen
#include <Eigen/Core>

// RBDyn
#include <RBDyn/FD.h>
#include <RBDyn/Jacobian.h>

// SCD
#include <SCD/Matrix/SCD_Types.h>

// Tasks
#include "QPSolver.h"
#include "QPConstr.h"
// forward declaration
// SCD
namespace SCD
{
class S_Object;
class CD_Pair;
}


namespace tasks
{

namespace qp
{

class MotionManipConstr : public ConstraintFunction<Equality, Bound>
{
public:
	MotionManipConstr(const rbd::MultiBody& mb);

	// Constraint
	 void updateNrVars(const rbd::MultiBody& mb,
		const SolverData& data);

	 void update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc);

	// Equality Constraint
	 int maxEq();

	 const Eigen::MatrixXd& AEq() const;
	 const Eigen::VectorXd& BEq() const;

	// Bound Constraint
	 int beginVar();

	 const Eigen::VectorXd& Lower() const;
	 const Eigen::VectorXd& Upper() const;

private:
	struct ContactData
	{
		ContactData() {}
		ContactData(const rbd::MultiBody& mb, int body,
			std::vector<Eigen::Vector3d> points,
			const std::vector<FrictionCone>& cones);


		rbd::Jacobian jac;
		int body;
		std::vector<Eigen::Vector3d> points;
		std::vector<Eigen::Matrix<double, 3, Eigen::Dynamic> > generators;
		// Hold the translated jacobian
		Eigen::MatrixXd jacTrans;
		// Hold the generator in world frame
		std::vector<Eigen::Matrix<double, 3, Eigen::Dynamic> > generatorsComp;
	};

private:
	rbd::ForwardDynamics fd_;
	rbd::ForwardDynamics fdManip_;

	const rbd::MultiBody *mbManip_;
	const rbd::MultiBodyConfig *mbcManip_;
	
	std::vector<ContactData> cont_;
	Eigen::MatrixXd fullJac_;
	Eigen::MatrixXd fullJacRobot_;
	
	std::vector<ContactData> contManip_;
	std::vector<ContactData> contRobot_;

	Eigen::MatrixXd AEq_;
	Eigen::VectorXd BEq_;

	Eigen::VectorXd XL_, XU_;

	int nrDof_, nrFor_, nrTor_;
};

class ContactManipAccConstr : public ConstraintFunction<Equality>,
		public ContactConstrCommon
{
public:
	ContactManipAccConstr(const rbd::MultiBody& mb);

	// Constraint
	 void updateNrVars(const rbd::MultiBody& mb,
		const SolverData& data);

	 void update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc);

	// Equality Constraint
	 int maxEq();

	 const Eigen::MatrixXd& AEq() const;
	 const Eigen::VectorXd& BEq() const;

private:
	struct ContactData
	{
		ContactData(rbd::Jacobian j):
			jac(j)
		{}


		rbd::Jacobian jac;
	};

private:
	std::vector<ContactData> cont_;
	std::vector<ContactData> contManip_;
	std::vector<ContactData> contRobot_;

	const rbd::MultiBody *mbManip_;
	const rbd::MultiBodyConfig *mbcManip_;

	Eigen::MatrixXd fullJacRobot_;
	Eigen::MatrixXd fullJacManip_;
	Eigen::VectorXd alphaVecRobot_;
	Eigen::VectorXd alphaVecManip_;

	Eigen::MatrixXd AEq_;
	Eigen::VectorXd BEq_;

	int nrDof_, nrFor_, nrTor_;
};

} //Namespace QP
} //Namespace Tasks
