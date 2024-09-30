#pragma once

#include "discamb/MathUtilities/Vector3.h"
#include "discamb/HC_Model/HC_WfnData.h"
#include "discamb/HC_Model/HC_AtomTypeParameters.h"
#include "discamb/HC_Model/HC_ModelParameters.h"
#include "json.hpp"


#include <vector>

namespace discamb{
    /**
* \addtogroup HC_Model
* @{
*/

	namespace hc_model_utilities{
		void rotate_plm(std::vector<std::vector<double> > &plm,
					const std::vector<Vector3d> &oldXYZ, const std::vector<Vector3d> &newXYZ);
        void serialize_hc_model(const HC_ModelParameters &p, nlohmann::json &data);
        void deserialize_hc_model(const HC_ModelParameters &p, nlohmann::json &data);
	}
    /**@}*/
}

