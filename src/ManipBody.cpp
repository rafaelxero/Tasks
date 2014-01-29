#include "ManipBody.h"

rbd::MultiBody tasks::fuseMultiBody(const rbd::MultiBody& mb, const rbd::MultiBody& mbManip,
		int bodyIdContact, const sva::PTransformd& toSurface)
{
	std::vector<rbd::Body> bodies(mb.bodies());
	std::vector<rbd::Joint> joints(mb.joints());
	std::vector<int> pred(mb.predecessors());
	std::vector<int> succ(mb.successors());
	std::vector<int> parent(mb.parents());
	std::vector<sva::PTransformd> Xt(mb.transforms());
	rbd::Body newBody; //Create copy of mbManip with new index
	rbd::Joint newJoint; //Create new fixed joint

	int bodyIndex = bodies.size();

	newBody = rbd::Body(mbManip.body(0).inertia(), 15000, "ManipBody");
	newJoint = rbd::Joint(rbd::Joint::Fixed, true, 42000, "ManipJoint");

	bodies.push_back(newBody);
	joints.push_back(newJoint);

	pred.push_back(mb.bodyIndexById(bodyIdContact));
	succ.push_back(bodyIndex);
	parent.push_back(mb.bodyIndexById(bodyIdContact));

	Xt.push_back(toSurface);

	return rbd::MultiBody(bodies, joints, pred, succ, parent, Xt);
}

void tasks::fillMbcManip(const rbd::MultiBodyConfig& mbc, rbd::MultiBodyConfig& mbcManip)
{
	int i = 0;
	for(auto vec: mbc.q)
	{
		mbcManip.q[i] = vec;
		++i;
	}
}
