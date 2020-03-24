/**
 * \file
 * \author Thomas Fischer
 * \date   2010-06-16
 * \brief  Implementation of string helper functions.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "StringTools.h"

#include <algorithm>
#include <cctype>
#include <cstdarg>
#include <cstdio>
#include <iomanip>

#include <boost/algorithm/string/replace.hpp>
#include "Logging.h"

namespace BaseLib
{
std::vector<std::string> splitString(std::string const& str)
{
    std::istringstream str_stream(str);
    std::vector<std::string> items;
    std::copy(std::istream_iterator<std::string>(str_stream),
        std::istream_iterator<std::string>(),
        std::back_inserter(items));
    return items;
}

std::list<std::string> splitString(const std::string &str, char delim)
{
    std::list<std::string> strList;
    std::stringstream ss(str);
    std::string item;
    while (getline(ss, item, delim))
    {
        strList.push_back(item);
    }
    return strList;
}

std::string replaceString(const std::string &searchString,
                          const std::string &replaceString,
                          std::string stringToReplace)
{
    boost::replace_all(stringToReplace, searchString, replaceString);
    return stringToReplace;
}

void trim(std::string &str, char ch)
{
    std::string::size_type pos = str.find_last_not_of(ch);
    if(pos != std::string::npos)
    {
        str.erase(pos + 1);
        pos = str.find_first_not_of(ch);
        if (pos != std::string::npos)
        {
            str.erase(0, pos);
        }
    }
    else
    {
        str.erase(str.begin(), str.end());
    }
}

void simplify(std::string &str)
{
    trim (str);
    str.erase(
        std::unique(str.begin(), str.end(), [](char a, char b) { return a == ' ' && b == ' '; }),
        str.end()
    );
}

std::string const& tostring(std::string const& value)
{
    return value;
}

}  // end namespace BaseLib
