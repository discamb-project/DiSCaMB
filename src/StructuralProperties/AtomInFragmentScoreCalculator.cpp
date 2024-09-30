#include "discamb/StructuralProperties/AtomInFragmentScoreCalculator.h"
#include "discamb/StructuralProperties/SimpleAIF_ScoreCalculator.h"

using namespace std;

namespace discamb {

	AtomInFragmentScoreCalculator::~AtomInFragmentScoreCalculator() {};

	AtomInFragmentScoreCalculator* AtomInFragmentScoreCalculator::create(
		const std::string& type,
		const nlohmann::json& data)
	{
		if (string("simple") == type)
			return new SimpleAIF_ScoreCalculator(data);
		return NULL;
	}

}