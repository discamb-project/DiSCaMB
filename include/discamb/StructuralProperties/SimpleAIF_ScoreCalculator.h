#pragma once

#include "AtomInFragmentScoreCalculator.h"
#include "discamb/CrystalStructure/UnitCellContent.h"

#include "json.hpp"

namespace discamb {

    /**
    * \addtogroup StructuralProperties
    * @{
    */


	class SimpleAIF_ScoreCalculator : public AtomInFragmentScoreCalculator
	{
	public:
		SimpleAIF_ScoreCalculator();
		SimpleAIF_ScoreCalculator(const nlohmann::json &data);
		virtual ~SimpleAIF_ScoreCalculator();
		virtual void assignScore(
			const Crystal& crystal,
			const std::vector<std::vector<std::pair<std::string, std::string> > >& clusters,
			std::vector < std::vector < double> >& scores);
	private:
		void assignScore(
			UnitCellContent &ucContent,
			std::vector<std::vector<UnitCellContent::AtomID> > &ucConnectivity,
			const std::vector<std::pair<std::string, std::string> >& cluster,
			std::vector < double>& scores);

        void assignScore2(
            UnitCellContent& ucContent,
            std::vector<std::vector<UnitCellContent::AtomID> >& ucConnectivity,
            const std::vector<std::pair<std::string, std::string> >& cluster,
            std::vector < double>& scores);

	};
    /** @}*/
}