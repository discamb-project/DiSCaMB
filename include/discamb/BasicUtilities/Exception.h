#ifndef _DISCAMB_BASICUTILITIES_EXCEPTION_H_
#define _DISCAMB_BASICUTILITIES_EXCEPTION_H_



#include <exception>
#include <string>

namespace discamb{

/**
* \defgroup BasicUtilities BasicUtilities
\brief Basic utilities for execption handling, operations on strings and performance measurements.
* @{
*/


/** Extends std::exception with custom error message and information on where the error occurred. */

class Exception: public std::exception
{
    public:
    /**
    \param errorMessage information about the error
    \param fileName name of the file where the error occurred
    \param line number of the line where the error occurred
    */
    Exception(const std::string &errorMessage,const std::string &fileName = std::string(), int line=0): mFileName(fileName), mErrorMessage(errorMessage), mLine(line)
    {
        if (mLine == 0)
            mLongMessage = errorMessage;
        else
        {
            mLongMessage = std::string("The following error has occurred:\n") + errorMessage;
            mLongMessage += std::string("\n\nInfo for developer - the error was signalized from:\n  file : ") + std::string(fileName);
            mLongMessage += std::string("\n  line : ") + std::to_string(line) + std::string("\n");
        }
    }
    Exception(const Exception &e): mFileName(e.mFileName), mErrorMessage(e.mErrorMessage), mLine(e.mLine)
    {
    };
    virtual ~Exception() throw(){};
    /** Provides information on where the error occurred / the exception was thrown from. 
    \param fileName name of the file
    \param line number of the line 
    */
    virtual void thrown_from(std::string &fileName,int &line) const {fileName = mFileName, line = mLine;}
    /**
    returns information about the error (errorMessage provided in the constructor)
    */
    virtual std::string errorMessage() const {return mErrorMessage;}
    /**
    returns information about the error (errorMessage provided in the constructor) 
    and where it occurred / where the exception was thrown from
    */

    virtual const char *what() const throw() {return mLongMessage.c_str();}

private:
    Exception();
    Exception& operator=(const Exception&){return *this;};
    std::string mFileName,mErrorMessage, mLongMessage;
    int mLine;
};

/**@}*/ 

}// namespace discamb





#endif /*_DISCAMB_BASICUTILITIES_EXCEPTION_H_*/
