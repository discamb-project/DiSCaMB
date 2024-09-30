#ifndef _DISCAMB_IO_NATIVEIAM_READER_HPP_
#define _DISCAMB_IO_NATIVEIAM_READER_HPP_


#include "discamb/CrystalStructure/Crystal.h"
#include "discamb/Scattering/SF_CalcDataTypes.h"
#include <vector>
#include <string>

namespace discamb {

    /**
    * \addtogroup IO IO
    * @{
    */


class NativeIAM_Reader
{
public:
    NativeIAM_Reader();
    virtual ~NativeIAM_Reader();

    void read(const std::string &fileName, Crystal &crystal, ScatteringParameters &scatteringParameters) const;
  
};

/** @} */

}

#endif /*_DISCAMB_IO_NATIVEIAM_READER_HPP_*/
