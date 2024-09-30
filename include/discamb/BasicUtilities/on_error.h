#ifndef _DISCAMB_BASICUTILITIES_ONERROR_H_
#define _DISCAMB_BASICUTILITIES_ONERROR_H_

#include <string>
#include <exception>

#include "Exception.h"
//#include "StringUtilities.h"

namespace discamb{
    /**
    * \addtogroup BasicUtilities BasicUtilities
    * @{

    \brief Utilities for error handling
    */

namespace on_error{



	/** \brief throw standard exception with formatted error message*/

	inline void throwException(const std::string &message,const char *fileName,int line)
	{
		throw Exception(message,fileName,line);
	}

    inline void throwException(const std::string& message)
    {
        throw Exception(message);
    }


	/** \brief throws exception noticing about calling not implemented function */
	inline void not_implemented(const char *fileName,int line_number)
	{
		throw Exception("calling not implemented function",fileName,line_number);
	}


} // namespace on_error

  /**@}*/

}// namespace discamb






#endif /*_DISCAMB_BASICUTILITIES_ONERROR_H_*/
