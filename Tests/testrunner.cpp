/**
 * \file
 * \author Lars Bilke
 * \date   2012-04-29
 * \brief  GTest test executables main function.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// ** INCLUDES **
#include <clocale>

#include "gtest/gtest.h"

#include "Applications/ApplicationsLib/LinearSolverLibrarySetup.h"
#include "BaseLib/Logging.h"
#include "NumLib/DOF/GlobalMatrixProviders.h"

#ifdef OGS_BUILD_GUI
#include <QCoreApplication>
#endif

/// Implementation of the googletest testrunner
int main(int argc, char* argv[])
{
#ifdef NDEBUG
    std::string logLevel("info");
#else
    std::string logLevel("all");
#endif
    for (int i = 1; i < argc; i++)
    {
        if (i + 1 == argc)
        {
            break;
        }
        if (std::strcmp(argv[i], "-l") == 0)
        {
            logLevel = argv[i + 1];
        }
    }

    setlocale(LC_ALL, "C");
#ifdef OGS_BUILD_GUI
    QCoreApplication app(argc, argv, false);
#endif

    ApplicationsLib::LinearSolverLibrarySetup linear_solver_library_setup(
                argc, argv);

    /* TODO (naumov) BEFORE MERGING SPDLOG
    logog_setup.setFormatter(
        std::make_unique<BaseLib::TemplateLogogFormatterSuppressedGCC<
            TOPIC_LEVEL_FLAG | TOPIC_FILE_NAME_FLAG |
            TOPIC_LINE_NUMBER_FLAG>>());
            */
    BaseLib::setConsoleLogLevel(logLevel);

    try
    {
        // start google test
        testing::InitGoogleTest ( &argc, argv );
        return RUN_ALL_TESTS();
    }
    catch (char* e)
    {
        ERR(e);
        return 1;
    }
    catch (std::exception& e)
    {
        ERR(e.what());
        return 1;
    }
    catch (...)
    {
        ERR("Unknown exception occurred!");
        return 1;
    }

    return 0;
}
