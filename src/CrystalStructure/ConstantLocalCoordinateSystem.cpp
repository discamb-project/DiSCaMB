#include "discamb/CrystalStructure/ConstantLocalCoordinateSystem.h"
#include "discamb/BasicUtilities/string_utilities.h"
#include "discamb/BasicUtilities/on_error.h"

using namespace std;

namespace discamb {

	ConstantLocalCoordinateSystem::ConstantLocalCoordinateSystem()
	{
		set({ 1, 0, 0 },
			{ 0, 1, 0 },
			{ 0, 0, 1 });
	}

	ConstantLocalCoordinateSystem::ConstantLocalCoordinateSystem(
		const Vector3d &x,
		const Vector3d &y,
		const Vector3d &z)
	{
		set(x, y, z);
	}

	ConstantLocalCoordinateSystem::~ConstantLocalCoordinateSystem()
	{

	}

	void ConstantLocalCoordinateSystem::set(
		const std::string &definition,
		const Crystal &c)
	{
		// [1,1,0],[-0.5,0.5,0],[0,0,0.3] - xyz are normalized
		string s = string_utilities::removeChar(definition, '[');
		s = string_utilities::removeChar(s, ']');
		vector<string> words;
		string_utilities::split(s, words, ',');
		if (words.size() != 9)
			on_error::throwException("invalid definition of local coordinate systems vectors", __FILE__, __LINE__);
		Vector3d x, y, z;
		x.set(stod(words[0]), stod(words[1]), stod(words[2]));
		y.set(stod(words[3]), stod(words[4]), stod(words[5]));
		z.set(stod(words[6]), stod(words[7]), stod(words[8]));
		set(x, y, z);
	}

	void ConstantLocalCoordinateSystem::set(
		const Vector3d &x,
		const Vector3d &y,
		const Vector3d &z)
	{
		mX = x / sqrt(x*x);
		mY = y / sqrt(y*y);
		mZ = z / sqrt(z*z);
	}

	void ConstantLocalCoordinateSystem::calculate(
		Vector3d &x,
		Vector3d &y,
		Vector3d &z,
		const Crystal &c)
		const
	{
		x = mX;
		y = mY;
		z = mZ;
	}


}