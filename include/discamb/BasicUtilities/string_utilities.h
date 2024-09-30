#pragma once

#include "discamb/BasicUtilities/on_error.h"

#include <map>
#include <sstream>
#include <iomanip>
#include <string>
#include <algorithm>
#include <functional>
#include <cctype>
#include <locale>
#include <vector>
#include <exception>


namespace  discamb{

    /**
    * \addtogroup BasicUtilities BasicUtilities
    * @{
    */


namespace string_utilities {

enum class CharacterType {WHITE_SPACE, DIGIT};

/**
\brief Converts the argument to std::string
*/
template<typename T>
std::string convertToString(T const &);

/**
\brief Converts real number to std::string

\param real real number to be converted

\param precision precision of the number as a string

\param setToFixedOrScientific specifies if the last argument (fixed) will be used to set scientific or fixed format
                                (if false default formatting will be used)
\param fixed sets formating to fixed if true and to scientific if false, the argument is used only if 
               setToFixedOrScientific is set to true
*/

template<typename T>
std::string realToString(T const &real,int precision=6,bool setToFixedOrScientific=false,bool fixed=true);
    
/** 
\brief Converts from std::string to type T
*/
template<typename T>
T convertFromString(std::string const &arg);

/**
\brief Converts from std::string to type T
*/
template<typename T>
void convertFromString(std::string const &,T &);


/**
\brief Trims left side of the string
*/

// $todo should it not be const std::strign &s

static inline std::string &trim_left(std::string &s)
{
	std::string chars = "\t\n\v\f\r ";
	s.erase(0, s.find_first_not_of(chars));
    //s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
    return s;
}

/**
\brief Trims right side of the string
*/

// $todo const std::strign &s?

static inline std::string &trim_right(std::string &s)
{
	std::string chars = "\t\n\v\f\r ";
	s.erase(s.find_last_not_of(chars) + 1);
    //s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
    return s;
}

/**
\brief Trims the string s
*/

// $todo const std::strign &s ?

static inline std::string &trim(std::string &s)
{
    return trim_left(trim_right(s));
}

static inline int maxWidth(const std::vector<std::string> & words)
{
    int result = 0;
    for (auto const& s : words)
        result = std::max<int>(result, s.size());
    return result;
}


static inline void findAllOccurences(
    const std::string& s,
    char c,
    std::vector<int>& pos)
{
    pos.clear();
    int idx = 0;
    while ((idx = s.find(c, idx)) != std::string::npos) {
        pos.push_back(idx);
        idx++;
    }
}

/**
\brief Splits string into words
\param str string to be splitted
\param vecStr result of the split
\param delimiter char being a delimiter or a function accepting char as an argument an returning value convertible to bool 
                  (e.g. isdigit or isspace)

*/

template < typename T >
void split(const std::string &str,std::vector<std::string> &vecStr,T delimiter);

/**
\brief Splits string into words - specialization for delimiter of type char
\param str string to be splitted
\param vecStr result of the split
\param c char being a delimiter
(e.g. isdigit or isspace)

*/


template <>
inline
void split<char>(const std::string &str,std::vector<std::string> &vecStr,char delimiter)
{
    //split(str,vecStr,std::bind1st(std::equal_to<char>(),c));
	vecStr.clear();
	bool inWord=false;
	for (char c : str)
		if (c == delimiter)
			inWord = false;
		else
		{
			if (!inWord)
			{
				vecStr.resize(vecStr.size() + 1);
				inWord = true;
			}
			vecStr.back() += c;
		}

}

/**
\brief Splits string into words - specialization for delimiter of type enum CharacterType
\param str string to be splitted
\param vecStr result of the split
\param cType character type being a delimiter
(e.g. isdigit or isspace)

*/

template <>
inline
void split<CharacterType>(const std::string &str,std::vector<std::string> &vecStr,CharacterType cType)
{
    if(cType == string_utilities::CharacterType::WHITE_SPACE)
        split<int (int)>(str,vecStr,std::isspace);
}

/**
\brief Splits string into words using space as delimiter
\param str string to be splitted
\param vecStr result of the split

*/

inline
void split(const std::string &str,std::vector<std::string> &vecStr)
{
    split<char>(str,vecStr,' ');
}

/**
\brief Splits string into words
\param str string to be splitted
\param vecStr result of the split
\param delimiters vector of characters treated as delimiters
*/


inline
void split(const std::string &str,std::vector<std::string> &vecStr,const std::vector<char> &delimiters)
{
    std::vector<std::string> words,toSplit;
    int n,delimiterIndex,j,nDelimiters=delimiters.size();

    vecStr.assign(1,str);
    
     for( delimiterIndex = 0 ; delimiterIndex < nDelimiters ; delimiterIndex++ )
    {
        vecStr.swap(toSplit);
        n = toSplit.size();
        for(j=0;j<n;j++)
        {
            split<char>(toSplit[j], words, delimiters[delimiterIndex]);
            vecStr.insert(vecStr.end(),words.begin(),words.end());
        }
        toSplit.clear();
    }
}


/**
\brief portable version of getline - recognizes LF , CR+LF , LF+CR and CR as end of line
*/

inline
std::istream &portable_getline(std::istream &is,std::string &line)
{
    line.clear();
    char c;
    while(is.get(c))
    {
        if(c==10)
        {
            is.get(c);
            if(c==13)
                return is;
            is.unget();
            return is;
        }
        if(c==13)
        {
            is.get(c);
            if(c==10)
                return is;
            is.unget();
            return is;
        }
        line += c;
    }

    return is;
}

/** \brief Convers all characters in string to lower case. */

inline
void toLower(const std::string &inputString,std::string &transformedString)
{
    transformedString = inputString;
    std::transform(inputString.begin(),inputString.end(),transformedString.begin(),::tolower);
}


inline
std::string toLower(const std::string& inputString)
{
    std::string transformedString;
    transformedString = inputString;
    std::transform(inputString.begin(), inputString.end(), transformedString.begin(), ::tolower);
    return transformedString;
}


/** \brief Convers all characters in string to upper case. */

inline
void toUpper(const std::string &inputString,std::string &transformedString)
{
    transformedString = inputString;
    std::transform(inputString.begin(),inputString.end(),transformedString.begin(),::toupper);
}

inline
std::string toUpper(const std::string& inputString)
{
    std::string transformedString;
    transformedString = inputString;
    std::transform(inputString.begin(), inputString.end(), transformedString.begin(), ::toupper);
    return transformedString;
}

inline 
std::string replace(const std::string& s, char replaceThis, char withThis)
{
    std::string result;
    for (auto c : s)
        if (c == replaceThis)
            result += withThis;
        else
            result += c;
    return result;
}

/** \brief Removes all occurences of character from string
\param s string from which the character will be removed
\param c the character to remove
*/

inline std::string removeChar(const std::string &s,char c)
{
    std::string result;
    int i,n=s.size();
    for(i=0;i<n;i++)
        if(s[i]!=c)
            result+=s[i];
    return result;
}

/**  \brief Merges strings into one string separating them with delimiter
\param words strings to be merged
\param separator string separator

*/

inline std::string merge(const std::vector<std::string> &words,char separator)
{
    std::string result;
    int i,n=words.size();
    for(i=0;i<n;i++)
    {
        if(i!=0)
            result+=separator;
        result+=words[i];
    }
    return result;
}

//############# TEMPLATES IMPLEMENTATION ########################


template<typename T>
std::string convertToString(
T const &anything)
{
    std::stringstream ss;
    std::string result;	
    ss<< anything;
    ss>> result;
    return result;
}

    

template<typename T>
T convertFromString(
std::string const &s)
{
    std::stringstream ss;
    T result;
    ss<< s;
    ss>> result;
    
    return result;   
}



template<typename T>
void convertFromString(
std::string const &s,
T &result)
{
    result=convertFromString<T>(s);
}


template<typename T>
std::string realToString(
T const &number,
int precision,
bool setToFixedOrScientific,
bool fixed)
{
    std::stringstream ss;
    std::string result;

    if(setToFixedOrScientific)
        fixed ? ss.flags(std::ios_base::fixed) : ss.flags(std::ios_base::scientific);

    ss.precision(precision);
    
    ss<< number;
    ss>> result;
    
    return result;   
}


template < typename T >
void split(const std::string &str,std::vector<std::string> &vecStr,T t)
{	
    int i,nCharacters=str.size();
    bool previousWasDelimiter = true;
    bool currentIsDelimiter;
        
    vecStr.clear();
        
    for(i=0;i<nCharacters;i++)
    {
        currentIsDelimiter = ( t(int(str[i])) != 0 );
        if( !currentIsDelimiter )
        {
            if(previousWasDelimiter)
                vecStr.resize(vecStr.size()+1);		
            vecStr.back() += str[i]; 
        }
        previousWasDelimiter = 	currentIsDelimiter;
    }
}

    
inline void fill_template(
    const std::map<std::string, std::string> &dictionary,
    const std::string &input,
    std::string &output,
    char c)
{
    int start = input.find(c);
    int end = 0;
    output.clear();
    std::string key;
    while (start != std::string::npos)
    {
        output += input.substr(end, start - end);
        end = input.find(c, start+1);
        
        if (end == std::string::npos)
            on_error::throwException("missing matching character when processing template", __FILE__, __LINE__);
        key = input.substr(start + 1, end - start - 1);
        if (dictionary.find(key) == dictionary.end())
            on_error::throwException(std::string("encoutered unknown key '") + key + std::string("' when processing template"), __FILE__, __LINE__);
        else
            output += dictionary.find(key)->second;

        end += 1;
        start = input.find(c, end);
    }
    
    output += input.substr(end);
}

/**
converts from 1.234(56) type notation
*/

inline void string_to_number_and_uncertainty(
	const std::string &s,
	double &number,
	double &uncertainty)
{
	// is there any uncertainty reproted - looks for bracket
	std::vector<std::string> words;
	split(s, words, '(');
	if (words.size() == 0)
		on_error::throwException(std::string("an attempt to convert string not representing number: '") + s + std::string("' to number"), __FILE__, __LINE__);
	if (words.size() == 1)
	{
		number = std::stod(s);
		uncertainty = 0.0;
		return;
	}
	if (words.size() == 2)
	{
		if(words[1].back() != ')')
			on_error::throwException(std::string("cannot convert string: '") + s + std::string("' to number and uncertaintity"), __FILE__, __LINE__);
		words[1].pop_back();
		std::string uStringReversed, uString;
		auto uIt = words[1].rbegin();
		for (auto it = words[0].rbegin(); it != words[0].rend(); ++it)
		{
			if (isdigit(*it))
			{
				if (uIt != words[1].rend())
					uStringReversed += *(uIt++);
				else
					uStringReversed += '0';
			}
			else
				uStringReversed += *it;
		}
		while(uIt != words[1].rend())
			uStringReversed += *(uIt++);
		for (uIt = uStringReversed.rbegin(); uIt != uStringReversed.rend(); ++uIt)
			uString += *uIt;
		number = stod(words[0]);
		uncertainty = stod(uString);
		return;
	}

	on_error::throwException(std::string("an attempt to convert string not representing number: '") + s + std::string("' to number"), __FILE__, __LINE__);

}

inline std::string  fixed_precision_number_with_uncertainty(
    double number,
    double uncertainty,
    int max_precision = 15)
{
    std::stringstream ss, ss2;
    std::string s, uncertaintyString;

    if (uncertainty <= 0)
        ss << std::setprecision(max_precision) << std::fixed << number;

    if (uncertainty >= 10)
        ss << std::setprecision(0) << std::fixed << number << "(" << std::setprecision(0) << std::fixed << uncertainty << ")";

    if (uncertainty < 10 && uncertainty >= 1)
    {
        if (max_precision == 0)
            ss << std::setprecision(0) << std::fixed << number << "(" << std::setprecision(0) << std::fixed << uncertainty << ")";
        else
            ss << std::setprecision(1) << std::fixed << number << "(" << std::setprecision(0) << std::fixed << 10 * uncertainty << ")";

    }
    


    if (uncertainty < 1 && uncertainty > 0 && max_precision>0)
    {
        ss << std::setprecision(max_precision) << std::fixed << uncertainty;
        ss >> uncertaintyString;

        // u = 0.0034 n= 0.123456, n is cat at second digit, which is not 0 in u
        // precision of n
        int uncertainty_n_zeros_after_dot = 0;
        if (uncertaintyString.find('.') != std::string::npos)
            for (int i = 2; i < uncertaintyString.size(); i++)
                if (uncertaintyString[i] != '0')
                    break;
                else
                    uncertainty_n_zeros_after_dot++;

        ss.clear();

        if (uncertainty_n_zeros_after_dot == max_precision)
            ss << std::setprecision(max_precision) << std::fixed << number;

        if (uncertainty_n_zeros_after_dot == max_precision - 1)
            ss << std::setprecision(max_precision) << std::fixed << number << "(" << uncertaintyString.back() << ")";

        if (uncertainty_n_zeros_after_dot < max_precision - 1)
            ss << std::setprecision(uncertainty_n_zeros_after_dot + 2) << std::fixed << number << "(" << uncertaintyString[2 + uncertainty_n_zeros_after_dot]
            << uncertaintyString[3 + uncertainty_n_zeros_after_dot] << ")";

    }
    ss >> s;
    return s;
}


} //namespace string_utilities
 /**@}*/
} //namespace  discamb


