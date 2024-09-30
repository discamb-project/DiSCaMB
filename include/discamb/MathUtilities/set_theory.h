#include <set>
#include <vector>

namespace discamb {

    namespace set_theory {

        template<typename T>
        bool is_subset(const std::vector<T>& set, const std::vector<T>& subset)
        {
            std::multiset<T> s(set.begin(), set.end());

            for (const auto & element : subset)
            {
                auto matchingElement = s.find(element);
                if (matchingElement != s.end())
                    s.erase(matchingElement);
                else
                    return false;
            }
            return true;
        }

    }

}