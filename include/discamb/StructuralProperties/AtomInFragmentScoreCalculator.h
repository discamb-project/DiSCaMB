#pragma once
#include "discamb/CrystalStructure/Crystal.h"
#include "json.hpp"

namespace discamb {

    /**
    * \addtogroup StructuralProperties
    * @{
    */


	class AtomInFragmentScoreCalculator
	{
	public:
		virtual ~AtomInFragmentScoreCalculator() = 0;
		virtual void assignScore(
			const Crystal &crystal,
			const std::vector<std::vector<std::pair<std::string, std::string> > > &clusters,
			std::vector < std::vector < double> > &scores) = 0;
		static AtomInFragmentScoreCalculator* create(const std::string& type, const nlohmann::json& data = nlohmann::json());
	};
    /** @}*/
}