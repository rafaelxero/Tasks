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
	fdManip_(mb),
	cont_(),
	fullJac_(6, mb.nrDof()),
	fullJacRobot_(6, mb.nrDof()),
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
	fdManip_ = rbd::MultiBody(data.manipBody());
	mbManip_ = &data.manipBody();
	mbcManip_ = &data.manipBodyConfig();
	const auto& uniCont = data.unilateralContacts();
	const auto& biCont = data.bilateralContacts();
	cont_.resize(data.nrContacts());

	const auto& manipCont = data.robotToManipBodyContacts();
	const auto& robotCont = data.manipBodyToRobotContacts();
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
		contManip_[iCont] = ContactData(*mbManip_, c.bodyId, c.points,
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


void MotionManipConstr::update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc)
{
	using namespace Eigen;

	fd_.computeH(mb, mbc);
	fd_.computeC(mb, mbc);
	
	fdManip_.computeH(*mbManip_, *mbcManip_);
	fdManip_.computeC(*mbManip_, *mbcManip_);

	// H*alphaD - tau - tau_c = -C

	// AEq
	//         nrDof      nrFor            nrTor
	// nrDof [   H      -Sum J_i^t*ni     [0 ... -1]
	AEq_.block(0, 0, nrDof_-6, nrDof_-6) << fd_.H();
	AEq_.block(nrDof_-6, nrDof_-6, 6, 6) << fdManip_.H();
	
	fullJac_.resize(6, mb.nrDof());
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

			AEq_.block(0, contPos, nrDof_-6, cont_[i].generatorsComp[j].cols()) =
				-fullJac_.block(3, 0, 3, fullJac_.cols()).transpose()*
					cont_[i].generatorsComp[j];

			contPos += int(cont_[i].generatorsComp[j].cols());
		}
	}
	fullJac_.resize(6,6);
	for(std::size_t i = 0; i < contManip_.size(); ++i)
	{
		const MatrixXd& jacRobot = contRobot_[i].jac.jacobian(mb, mbc);
		const MatrixXd& jacManip = contManip_[i].jac.jacobian(*mbManip_, *mbcManip_);
		// for each contact point we compute all the torques
		// due to each generator of the friction cone
		for(std::size_t j = 0; j < contManip_[i].points.size(); ++j)
		{
			contManip_[i].generatorsComp[j] =
				mbcManip_->bodyPosW[contManip_[i].body].rotation().transpose()*contManip_[i].generators[j];

			contManip_[i].jac.translateJacobian(jacManip, *mbcManip_,
				contManip_[i].points[j], contManip_[i].jacTrans);
			contManip_[i].jac.fullJacobian(*mbManip_, contManip_[i].jacTrans, fullJac_);
			
			contRobot_[i].generatorsComp[j] =
				mbc.bodyPosW[contRobot_[i].body].rotation().transpose()*contRobot_[i].generators[j];

			contRobot_[i].jac.translateJacobian(jacRobot, mbc,
				contRobot_[i].points[j], contRobot_[i].jacTrans);
			contRobot_[i].jac.fullJacobian(mb, contRobot_[i].jacTrans, fullJacRobot_);

			AEq_.block(0, contPos, nrDof_-6, contRobot_[i].generatorsComp[j].cols()) =
				-fullJacRobot_.block(3, 0, 3, fullJacRobot_.cols()).transpose()*
					contRobot_[i].generatorsComp[j];
			AEq_.block(nrDof_-6, contPos, 6, contManip_[i].generatorsComp[j].cols()) =
				-fullJac_.block(3, 0, 3, fullJac_.cols()).transpose()*
					-contRobot_[i].generatorsComp[j];
			contPos += int(contManip_[i].generatorsComp[j].cols());
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

ContactManipAccConstr::ContactManipAccConstr(const rbd::MultiBody& mb):
	contManip_(),
	contRobot_(),
	fullJacRobot_(6, mb.nrDof()),
	fullJacManip_(6,6),
	alphaVecRobot_(mb.nrDof()),
	alphaVecManip_(6),
	AEq_(),
	BEq_(),
	nrDof_(0),
	nrFor_(0),
	nrTor_(0)
{}


void ContactManipAccConstr::updateNrVars(const rbd::MultiBody& mb,
	const SolverData& data)
{
	mbManip_ = &data.manipBody();
	mbcManip_ = &data.manipBodyConfig();
	cont_.clear();
	contManip_.clear();
	contRobot_.clear();
	nrDof_ = data.alphaD();
	nrFor_ = data.lambda();
	nrTor_ = data.torque();

	std::set<int> bodyIdSet = bodyIdInContact(mb, data);
	for(int bodyId: bodyIdSet)
	{
		cont_.emplace_back(rbd::Jacobian(mb, bodyId));
	}

	for(const UnilateralContact& c: data.manipBodyToRobotContacts())
	{
		contRobot_.emplace_back(rbd::Jacobian(mb, c.bodyId));
	}

	for(const UnilateralContact& c: data.robotToManipBodyContacts())
	{
		contManip_.emplace_back(rbd::Jacobian(data.manipBody(), c.bodyId));
	}

	auto size = cont_.size() + contRobot_.size();

	AEq_.resize(size*6 , nrDof_ + nrFor_ + nrTor_);
	BEq_.resize(size*6);

	AEq_.setZero();
	BEq_.setZero();
}


void ContactManipAccConstr::update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc)
{
	using namespace Eigen;

	rbd::paramToVector(mbc.alpha, alphaVecRobot_);
	rbd::paramToVector(mbcManip_->alpha, alphaVecManip_);
	int offset = 0;
	// Robot{J_i*alphaD + JD_i*alpha} =  Manip{J_i*alphaD + JD_i*alpha}

	for(std::size_t i = 0; i < cont_.size(); ++i)
	{
		// AEq = [-Robot{J_i}  Manip{J_i} ]
		const MatrixXd& jacRobot = cont_[i].jac.jacobian(mb, mbc);
		cont_[i].jac.fullJacobian(mb, jacRobot, fullJacRobot_);
		AEq_.block(i*6, 0, 6, mb.nrDof()) = -fullJacRobot_;

		// BEq = -JD_i*alpha
		const MatrixXd& jacDotRobot = cont_[i].jac.jacobianDot(mb, mbc);
		cont_[i].jac.fullJacobian(mb, jacDotRobot, fullJacRobot_);
		BEq_.segment(i*6, 6) = fullJacRobot_*alphaVecRobot_;
		++offset;
	}

	for(std::size_t i = 0; i < contRobot_.size(); ++i)
	{
		// AEq = [-Robot{J_i}  Manip{J_i} ]
		const MatrixXd& jacRobot = contRobot_[i].jac.jacobian(mb, mbc);
		const MatrixXd& jacManip = contManip_[i].jac.jacobian(*mbManip_, *mbcManip_);
		contRobot_[i].jac.fullJacobian(mb, jacRobot, fullJacRobot_);
		contManip_[i].jac.fullJacobian(*mbManip_, jacManip, fullJacManip_);
		AEq_.block((offset+i)*6, 0, 6, mb.nrDof()) = -fullJacRobot_;
		AEq_.block((offset+i)*6,mb.nrDof(),6,6) = fullJacManip_;

		// BEq = -JD_i*alpha
		const MatrixXd& jacDotRobot = contRobot_[i].jac.jacobianDot(mb, mbc);
		const MatrixXd& jacDotManip = contManip_[i].jac.jacobianDot(*mbManip_, *mbcManip_);
		contRobot_[i].jac.fullJacobian(mb, jacDotRobot, fullJacRobot_);
		contManip_[i].jac.fullJacobian(mb, jacDotManip, fullJacManip_);
		BEq_.segment((offset+i)*6, 6) = fullJacRobot_*alphaVecRobot_ - fullJacManip_*alphaVecManip_;
	}
}


int ContactManipAccConstr::maxEq()
{
	return int(AEq_.rows());
}


const Eigen::MatrixXd& ContactManipAccConstr::AEq() const
{
	return AEq_;
}


const Eigen::VectorXd& ContactManipAccConstr::BEq() const
{
	return BEq_;
}

ContactManipSpeedConstr::ContactManipSpeedConstr(const rbd::MultiBody& mb, double timeStep):
	contManip_(),
	contRobot_(),
	fullJacRobot_(6, mb.nrDof()),
	fullJacManip_(6,6),
	alphaVecRobot_(mb.nrDof()),
	alphaVecManip_(6),
	AEq_(),
	BEq_(),
	timeStep_(timeStep),
	nrDof_(0),
	nrFor_(0),
	nrTor_(0)
{}


void ContactManipSpeedConstr::updateNrVars(const rbd::MultiBody& mb,
	const SolverData& data)
{
	mbManip_ = &data.manipBody();
	mbcManip_ = &data.manipBodyConfig();
	cont_.clear();
	contManip_.clear();
	contRobot_.clear();
	nrDof_ = data.alphaD();
	nrFor_ = data.lambda();
	nrTor_ = data.torque();

	std::set<int> bodyIdSet = bodyIdInContact(mb, data);
	for(int bodyId: bodyIdSet)
	{
		cont_.emplace_back(rbd::Jacobian(mb, bodyId), sva::PTransformd::Identity());
	}

	for(const ManipContact& c: data.manipBodyToRobotContacts())
	{
		contRobot_.emplace_back(rbd::Jacobian(mb, c.contact.bodyId), c.toSurface);
	}

	for(const ManipContact& c: data.robotToManipBodyContacts())
	{
		contManip_.emplace_back(rbd::Jacobian(data.manipBody(), c.contact.bodyId), c.toSurface);
	}

	auto size = cont_.size() + contRobot_.size();

	AEq_.resize(size*6 , nrDof_ + nrFor_ + nrTor_);
	BEq_.resize(size*6);

	AEq_.setZero();
	BEq_.setZero();
}


void ContactManipSpeedConstr::update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc)
{
	using namespace Eigen;

	rbd::paramToVector(mbc.alpha, alphaVecRobot_);
	rbd::paramToVector(mbcManip_->alpha, alphaVecManip_);

	for(std::size_t i = 0; i < contRobot_.size(); ++i)
	{
		// AEq = [-Robot{J_i}  Manip{J_i} ]
		const MatrixXd& jacRobot = contRobot_[i].jac.jacobian(mb, mbc);
		const MatrixXd& jacManip = contManip_[i].jac.jacobian(*mbManip_, *mbcManip_);
		contRobot_[i].jac.translateJacobian(jacRobot, mbc, contRobot_[i].toSurface.translation(), contRobot_[i].jacTrans);
		contManip_[i].jac.translateJacobian(jacManip, mbc, contManip_[i].toSurface.translation(), contManip_[i].jacTrans);

		contRobot_[i].jac.fullJacobian(mb, contRobot_[i].jacTrans, fullJacRobot_);
		contManip_[i].jac.fullJacobian(mb, contManip_[i].jacTrans, fullJacManip_);
		//contRobot_[i].jac.fullJacobian(mb, jacRobot, fullJacRobot_); //Untranslated jacobians
		//contManip_[i].jac.fullJacobian(*mbManip_, jacManip, fullJacManip_);
		AEq_.block(i*6, 0, 6, mb.nrDof()) = -fullJacRobot_;
		AEq_.block(i*6,mb.nrDof(),6,6) = fullJacManip_;

		// BEq = JD_i*alpha - JD_manip_i * alpha + delta_speed/timeStep
		const MatrixXd& jacDotRobot = contRobot_[i].jac.jacobianDot(mb, mbc);
		const MatrixXd& jacDotManip = contManip_[i].jac.jacobianDot(*mbManip_, *mbcManip_);

		contManip_[i].jac.translateJacobian(jacDotRobot, mbc, contRobot_[i].toSurface.translation(), contRobot_[i].jacTrans);
		contManip_[i].jac.translateJacobian(jacDotManip, mbc, contManip_[i].toSurface.translation(), contManip_[i].jacTrans);

		contRobot_[i].jac.fullJacobian(mb, contRobot_[i].jacTrans, fullJacRobot_);
		contManip_[i].jac.fullJacobian(mb, contManip_[i].jacTrans, fullJacManip_);
		//contRobot_[i].jac.fullJacobian(mb, jacDotRobot, fullJacRobot_); //Untranslated jacobian dot
		//contManip_[i].jac.fullJacobian(mb, jacDotManip, fullJacManip_);

		sva::PTransformd tf1 = sva::PTransformd(mbc.bodyPosW[contRobot_[i].body].rotation());
		sva::PTransformd tf2 = sva::PTransformd(mbcManip_->bodyPosW[contManip_[i].body].rotation());

		auto bodyVel = tf1.invMul(contRobot_[i].toSurface*mbc.bodyVelB[contRobot_[i].body]);
		auto manipVel = tf2.invMul(contManip_[i].toSurface*mbcManip_->bodyVelB[contManip_[i].body]);

		BEq_.segment(+i*6, 6) = fullJacRobot_*alphaVecRobot_ - fullJacManip_*alphaVecManip_
						+ (bodyVel.vector() - manipVel.vector())/timeStep_;
	}
}


int ContactManipSpeedConstr::maxEq()
{
	return int(AEq_.rows());
}


const Eigen::MatrixXd& ContactManipSpeedConstr::AEq() const
{
	return AEq_;
}


const Eigen::VectorXd& ContactManipSpeedConstr::BEq() const
{
	return BEq_;
}

} // Namespace qp
} //Namespace Tasks
