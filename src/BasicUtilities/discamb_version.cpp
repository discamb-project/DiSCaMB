#include "discamb/BasicUtilities/discamb_version.h"

//#include "discamb/revision.h"

using namespace std;

namespace discamb{
	namespace discamb_version{
        std::string version()
        {
            return string("1.002");
        }

        //std::string revision()
        //{
        //    string revisionString;

        //    #ifdef Project_WC_REVISION
        //        revisionString = Project_WC_REVISION;
        //    #else
        //        revisionString = "?";
        //    #endif

        //    #ifdef GIT_BRANCH
        //        revisionString = GIT_BRANCH;
        //    #else
        //        revisionString = "?";
        //    #endif

        //    return revisionString;
        //}
	}
}
