#ifndef _DISCAMBDEV_IO_HKL_IO_H_
#define _DISCAMBDEV_IO_HKL_IO_H_

#include "discamb/MathUtilities/Vector3.h"

#include <string>
#include <vector>
#include <complex>

namespace discamb {
    /**
    * \addtogroup IO IO
    * @{
    */

namespace hkl_io 
{
    
// - discamb

/** Read hkl indices from text file  - for each line in the file it takes 3 subsequent numbers as h, k, l 
starting from column \p hIndexColumn if there is enough columns in given line. */
void readHklIndices(const std::string &fileName, std::vector<Vector3i> &hklIndices, int hIndexColumn = 1);

/** Read hkl indices from text file in DiSCaMB 'compact' format.*/

//void readCompactHklFile(const std::string &fileName, std::vector<Vector3i> &hkl);


/** Write hkl indices to text file in DiSCaMB 'compact' format.*/

//void writeCompactHklFile(const std::string &fileName, const std::vector<Vector3i> &hkl);

/** Read file with hkl and structure fctors of th following format:
- 1-st line header : "number of reflection N" where N is the number of reflections
- subsequent lines "h_index k_index l_index f_r f_i) where f_r and f_i are real
  and inmaginary components of structure factor 
*/
//void readReflectionFile(const std::string &reflection_file,
  //                      std::vector<Vector3i> &hkl,
    //                    std::vector<std::complex<double> > &sf);



// - discamb
    
    void readShelxHkl(
        const std::string &fileName,
        std::vector<Vector3i> &hkl,
        std::vector<double> &intensities,
        std::vector<double> &sigma,
        std::vector<int> &batchNumber,
        bool freeFormat);

    void writeShelxHkl(
        const std::string &fileName,
        const std::vector<Vector3i> &hkl,
        const std::vector<double> &intensities,
        const std::vector<double> &sigma,
        const std::vector<int> &batchNumber,
        bool freeFormat);

//    void readHklSf(
  //      const std::string &fileName,
    //    std::vector<Vector3i> &hkl,
      //  std::vector<std::complex<double> > &sf);

//    void writeHklSf(
  //     const std::string &fileName,
   ///     const std::vector<Vector3i> &hkl,
      //  const std::vector<std::complex<double> > &sf);
}
/**@}*/
}
#endif
