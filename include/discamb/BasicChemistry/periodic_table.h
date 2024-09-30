#ifndef _DISCAMB_BASICCHEMISTRY_PERIODICTABLE_H_
#define _DISCAMB_BASICCHEMISTRY_PERIODICTABLE_H_

#include <string>
#include <vector>
#include <map>

namespace discamb {

    /**
    * \defgroup BasicChemistry BasicChemistry
    \brief Basic chemistry related utilities.
    * @{
    */

    /** \brief Periodic table related information.*/

    namespace periodic_table {

        std::string symbol(int atomicNumber);

        int atomicNumber(const std::string &symbol, bool caseSensitive = true);

        enum class Block {S, P, D, F, NONE};

        /** S, P, D or F (lanthanides and actinides has no block assigned */
        Block block(int atomicNumber);

        char block2str(Block b);

        int period(int atomicNumber);

        int group(int atomicNumber);

    }
    /**@}*/

} // namespace discamb


#endif

