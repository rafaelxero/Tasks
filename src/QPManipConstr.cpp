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

// associated header
#include "QPManipConstr.h"

// includes
// std
#include <set>

// RBDyn
#include <RBDyn/MultiBody.h>
#include <RBDyn/MultiBodyConfig.h>

// SCD
#include <SCD/CD/CD_Pair.h>
#include <SCD/S_Object/S_Object.h>

namespace tasks
{

namespace qp
{

MotionManipConstr::ContactData::ContactData(const rbd::MultiBody& mb, int b,
	std::vector<Eigen::Vector3d> p,
	const std::vector<FrictionCone>& cones):
	jac(mb, b),
	body(jac.jointsPath().back()),
	points(std::move(p)),
	generators(cones.size()),
	jacTrans(6, jac.dof()),
	generatorsComp(cones.size())
{
	for(std::size_t i = 0; i < cones.size(); ++i)
	{
		generators[i].resize(3, cones[i].generators.size());
		generatorsComp[i].resize(3, cones[i].generators.size());
		for(std::size_t j = 0; j < cones[i].generators.size(); ++j)
		{
			generators[i].col(j) = cones[i].generators[j];
		}
	}
}


MotionManipConstr::MotionManipConstr(const rbd::MultiBody& mb):
	fd_(mb),
	cont_(),
	fullJac_(6, mb.nrDof()),
	AEq_(),
	BEq_(),
	XL_(),
	XU_(),
	nrDof_(0),
	nrFor_(0),
	nrTor_(0)
{
}


void MotionManipConstr::updateNrVars(const rbd::MultiBody& mb,
	const SolverData& data)
{
	const auto& uniCont = data.unilateralContacts();
	const auto& biCont = data.bilateralContacts();
	cont_.resize(data.nrContacts());

	const auto& manipCont = data.robotToManipBodyContacts();
	const auto& manipRobot = data.ManipBodyToRobotContacts();
	contManip_.resize(data.nrContactsManip());
	contRobot_.resize(data.nrContactsManip());

	nrDof_ = data.alphaD();
	nrFor_ = data.lambda();
	nrTor_ = data.torque();

	std::size_t iCont = 0;
	for(const UnilateralContact& c: uniCont)
	{
		cont_[iCont] = ContactData(mb, c.bodyId, c.points,
			std::vector<FrictionCone>(c.points.size(), c.cone));

		++iCont;
	}

	for(const BilateralContact& c: biCont)
	{
		cont_[iCont] = ContactData(mb, c.bodyId, c.points, c.cones);

		++iCont;
	}

	iCont = 0;
	for(const UnilateralContact& c: manipCont)
	{
		contManip_[iCont] = ContactData(mb, c.bodyId, c.points,
			std::vector<FrictionCone>(c.points.size(), c.cone));
		++iCont;
	}

	iCont = 0;
	for(const UnilateralContact& c: robotCont)
	{
		contRobot_[iCont] = ContactData(mb, c.bodyId, c.points,
			std::vector<FrictionCone>(c.points.size(), c.cone));
		++iCont;
	}

	AEq_.resize(nrDof_, nrDof_ + nrFor_ + nrTor_);
	BEq_.resize(nrDof_);

	AEq_.setZero();
	BEq_.setZero();

	XL_.resize(data.lambda());
	XU_.resize(data.lambda());

	XL_.fill(0.);
	XU_.fill(std::numeric_limits<double>::infinity());
}


void MotionManipConstr::update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc,
				const rbd::MultiBody& mbManip, const rbd::MultiBodyConfig& mbcManip)
{
	using namespace Eigen;

	fd_.computeH(mb, mbc);
	fd_.computeC(mb, mbc);
	
	fdManip_.computeH(mbManip, mbcManip);
	fdManip_.computeC(mbManip, mbcManip);

	// H*alphaD - tau - tau_c = -C

	// AEq
	//         nrDof      nrFor            nrTor
	// nrDof [   H      -Sum J_i^t*ni     [0 ... -1]

	AEq_.block(0, 0, nrDof_, nrDof_) << fd_.H(), MatrixXd::Zero(nrDof_-6,6),
						MatrixXd::Zero(6,nrDof_-6), fdManip_.H();

	int contPos = nrDof_;
	for(std::size_t i = 0; i < cont_.size(); ++i)
	{
		const MatrixXd& jac = cont_[i].jac.jacobian(mb, mbc);

		// for each contact point we compute all the torques
		// due to each generator of the friction cone
		for(std::size_t j = 0; j < cont_[i].points.size(); ++j)
		{
			cont_[i].generatorsComp[j] =
				mbc.bodyPosW[cont_[i].body].rotation().transpose()*cont_[i].generators[j];

			cont_[i].jac.translateJacobian(jac, mbc,
				cont_[i].points[j], cont_[i].jacTrans);
			cont_[i].jac.fullJacobian(mb, cont_[i].jacTrans, fullJac_);

			AEq_.block(0, contPos, nrDof_, cont_[i].generatorsComp[j].cols()) =
				-fullJac_.block(3, 0, 3, fullJac_.cols()).transpose()*
					cont_[i].generatorsComp[j];

			contPos += int(cont_[i].generatorsComp[j].cols());
		}
	}

	for(std::size_t i = 0; i < contManip_.size(); ++i)
	{
		const MatrixXd& jacManip = contManip_[i].jac.jacobian(mbManip, mbcManip);
		const MatrixXd& jacRobot = contRobot_[i].jac.jacobian(mbManip, mbcManip);

		// for each contact point we compute all the torques
		// due to each generator of the friction cone
		for(std::size_t j = 0; j < contManip_[i].points.size(); ++j)
		{
			contManip_[i].generatorsComp[j] =
				mbcManip.bodyPosW[contManip_[i].body].rotation().transpose()*contManip_[i].generators[j];

			contManip_[i].jac.translateJacobian(jacManip, mbcManip,
				contManip_[i].points[j], contManip_[i].jacTrans);
			contManip_[i].jac.fullJacobian(mbManip, contManip_[i].jacTrans, fullJac_);
			
			contRobot_[i].generatorsComp[j] =
				mbcManip.bodyPosW[contRobot_[i].body].rotation().transpose()*contRobot_[i].generators[j];

			contRobot_[i].jac.translateJacobian(jacManip, mbcManip,
				contManip_[i].points[j], contRobot_[i].jacTrans);
			contRobot_[i].jac.fullJacobian(mbManip, contManip_[i].jacTrans, fullJacRobot_);

			AEq_.block(0, contPos, nrDof_-6, contManip_[i].generatorsComp[j].cols()) =
				-fullJac_.block(3, 0, 3, fullJac_.cols()).transpose()*
					contManip_[i].generatorsComp[j];
			AEq_.block(nrDof_-6, contPos, 6, contManip_[i].generatorsComp[j].cols()) =
				-fullJacRobot_.block(3, 0, 3, fullJacRobot_.cols()).transpose()*
					contRobot_[i].generatorsComp[j];

			contPos += int(cont_[i].generatorsComp[j].cols());
		}
	}



	AEq_.block(mb.joint(0).dof(), contPos, nrTor_, nrTor_) =
		-MatrixXd::Identity(nrTor_, nrTor_);

	// BEq = -C
	BEq_ << -fd_.C(),-fdManip_.C();
}


int MotionManipConstr::maxEq()
{
	return nrDof_;
}


const Eigen::MatrixXd& MotionManipConstr::AEq() const
{
	return AEq_;
}


const Eigen::VectorXd& MotionManipConstr::BEq() const
{
	return BEq_;
}


int MotionManipConstr::beginVar()
{
	return nrDof_;
}


const Eigen::VectorXd& MotionManipConstr::Lower() const
{
	return XL_;
}


const Eigen::VectorXd& MotionManipConstr::Upper() const
{
	return XU_;
}

} // Namespace qp

} //Namespace Tasks
