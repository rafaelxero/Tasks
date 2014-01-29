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

#include <SpaceVecAlg/SpaceVecAlg>
#include <RBDyn/MultiBody.h>
#include <RBDyn/MultiBodyConfig.h>

namespace tasks
{

rbd::MultiBody fuseMultiBody(const rbd::MultiBody& mb, const rbd::MultiBody& mbManip,
		int bodyIdContact, const sva::PTransformd& toSurface);

void fillMbcManip(const rbd::MultiBodyConfig& mbc, rbd::MultiBodyConfig& mbcManip);
}
