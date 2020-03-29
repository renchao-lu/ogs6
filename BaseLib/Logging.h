/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <spdlog/spdlog.h>

namespace BaseLib
{
extern std::shared_ptr<spdlog::logger> console;
void setConsoleLogLevel(std::string const& level_string);

}  // namespace BaseLib

template <typename... Args>
void DBUG(char const* fmt, Args const&... args)
{
    BaseLib::console->debug(fmt, args...);
}
template <typename... Args>
void INFO(char const* fmt, Args const&... args)
{
    BaseLib::console->info(fmt, args...);
}
template <typename... Args>
void WARN(char const* fmt, Args const&... args)
{
    BaseLib::console->warn(fmt, args...);
}
template <typename... Args>
void ERR(char const* fmt, Args const&... args)
{
    BaseLib::console->error(fmt, args...);
}
template <typename... Args>
void CRITICAL(char const* fmt, Args const&... args)
{
    BaseLib::console->critical(fmt, args...);
}
