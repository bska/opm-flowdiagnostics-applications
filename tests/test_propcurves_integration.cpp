/*
  Copyright 2025 Equinor ASA.

  This file is part of the Open Porous Media Project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#define NVERBOSE

#define BOOST_TEST_MODULE Interation_Test_Property_Curves

#include <boost/test/unit_test.hpp>

#include <opm/utility/ECLGraph.hpp>
#include <opm/utility/ECLPhaseIndex.hpp>
#include <opm/utility/ECLResultData.hpp>
#include <opm/utility/ECLSaturationFunc.hpp>

#include <boost/filesystem.hpp>

#include <array>
#include <cstddef>
#include <string>
#include <vector>

namespace {

    boost::filesystem::path norneDataDir()
    {
        return boost::filesystem::path{ "TestInput" }
            / "Integration" / "Models" / "Norne";
    }

    struct NorneFixture
    {
        NorneFixture()
            : init { norneDataDir() / "NORNE_ATW2013.INIT" }
            , G    { Opm::ECLGraph::load(norneDataDir() / "NORNE_ATW2013.EGRID", init) }
        {}

        Opm::ECLInitFileData init;
        Opm::ECLGraph G;
    };

    struct NorneSatfuncFixture : public NorneFixture
    {
        NorneSatfuncFixture()
            : NorneFixture()
            , satFunc { G, init }
        {}

        Opm::ECLSaturationFunc satFunc;
    };

} // Anonymous namespace

BOOST_FIXTURE_TEST_SUITE(Saturation_Functions, NorneSatfuncFixture)

namespace {
    std::array<int, 3> cell_16_31_10()
    {
        return { 16 - 1, 31 - 1, 10 - 1 };
    }
} // Anonymous namespace

BOOST_AUTO_TEST_SUITE(Unscaled)

BOOST_AUTO_TEST_CASE(Relperm_Gas)
{
    auto curves = std::vector<Opm::ECLSaturationFunc::RawCurve>{};

    auto& curve = curves.emplace_back();
    curve.curve = Opm::ECLSaturationFunc::RawCurve::Function::RelPerm;
    curve.subsys = Opm::ECLSaturationFunc::RawCurve::SubSystem::OilGas;
    curve.thisPh = Opm::ECLPhaseIndex::Vapour;

    auto scaling = Opm::ECLSaturationFunc::SatFuncScaling{};
    scaling.enable = static_cast<unsigned char>(0);

    const auto graphs = this->satFunc
        .getSatFuncCurve(curves, this->G.activeCell(cell_16_31_10()), scaling);

    BOOST_REQUIRE_EQUAL(graphs.size(), std::size_t{ 1 });

    // Gas saturations
    {
        const auto expect_sg = std::vector<double> {
            0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45,
            0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95,
            0.9999, 1.0,
        };

        const auto& sg = graphs.front().first;
        BOOST_REQUIRE_EQUAL(sg.size(), expect_sg.size());

        for (auto i = 0 * sg.size(); i < sg.size(); ++i) {
            BOOST_TEST_MESSAGE("Compare Sg[" << i << "] with Expected Value");
            BOOST_CHECK_CLOSE(sg[i], expect_sg[i], 1.0e-6);
        }
    }

    // Relative permeability values
    {
        const auto expect_kr = std::vector<double> {
            0.0, 0.001655, 0.006913, 0.016213, 0.02999, 0.048655, 0.072573,
            0.102046, 0.137287, 0.178402, 0.225368, 0.27803, 0.336093,
            0.399135, 0.466631, 0.538, 0.612665, 0.690169, 0.770395,
            0.854218, 0.9499, 0.95,
        };

        const auto& kr = graphs.front().second;
        BOOST_REQUIRE_EQUAL(kr.size(), expect_kr.size());

        for (auto i = 0 * kr.size(); i < kr.size(); ++i) {
            BOOST_TEST_MESSAGE("Compare Krg[" << i << "] with Expected Value");
            BOOST_CHECK_CLOSE(kr[i], expect_kr[i], 1.0e-6);
        }
    }
}

BOOST_AUTO_TEST_CASE(Relperm_Oil_in_Oil_Gas)
{
    auto curves = std::vector<Opm::ECLSaturationFunc::RawCurve>{};

    auto& curve = curves.emplace_back();
    curve.curve = Opm::ECLSaturationFunc::RawCurve::Function::RelPerm;
    curve.subsys = Opm::ECLSaturationFunc::RawCurve::SubSystem::OilGas;
    curve.thisPh = Opm::ECLPhaseIndex::Liquid;

    auto scaling = Opm::ECLSaturationFunc::SatFuncScaling{};
    scaling.enable = static_cast<unsigned char>(0);

    const auto graphs = this->satFunc
        .getSatFuncCurve(curves, this->G.activeCell(cell_16_31_10()), scaling);

    BOOST_REQUIRE_EQUAL(graphs.size(), std::size_t{ 1 });

    // Gas saturations
    {
        const auto expect_sg = std::vector<double> {
            0.0, 1.0E-4, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4,
            0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9,
            0.95, 0.9999, 1.0,
        };

        const auto& sg = graphs.front().first;
        BOOST_REQUIRE_EQUAL(sg.size(), expect_sg.size());

        for (auto i = 0 * sg.size(); i < sg.size(); ++i) {
            BOOST_TEST_MESSAGE("Compare Sg[" << i << "] with Expected Value");
            BOOST_CHECK_CLOSE(sg[i], expect_sg[i], 1.0e-6);
        }
    }

    // Relative permeability values
    {
        const auto expect_kr = std::vector<double> {
            1.0, 0.999613776, 0.806888, 0.633562, 0.485506, 0.364043,
            0.267589, 0.192992, 0.136554, 0.094671, 0.064151,
            0.042324, 0.027035, 0.016586, 0.009662, 0.005254,
            0.002597, 0.001117, 0.000384, 8.8E-05, 7.0E-06, 0.0, 0.0,
        };

        const auto& kr = graphs.front().second;
        BOOST_REQUIRE_EQUAL(kr.size(), expect_kr.size());

        for (auto i = 0 * kr.size(); i < kr.size(); ++i) {
            BOOST_TEST_MESSAGE("Compare Krog[" << i << "] with Expected Value");
            BOOST_CHECK_CLOSE(kr[i], expect_kr[i], 1.0e-6);
        }
    }
}

BOOST_AUTO_TEST_CASE(Relperm_Oil_in_Oil_Water)
{
    auto curves = std::vector<Opm::ECLSaturationFunc::RawCurve>{};

    auto& curve = curves.emplace_back();
    curve.curve = Opm::ECLSaturationFunc::RawCurve::Function::RelPerm;
    curve.subsys = Opm::ECLSaturationFunc::RawCurve::SubSystem::OilWater;
    curve.thisPh = Opm::ECLPhaseIndex::Liquid;

    auto scaling = Opm::ECLSaturationFunc::SatFuncScaling{};
    scaling.enable = static_cast<unsigned char>(0);

    const auto graphs = this->satFunc
        .getSatFuncCurve(curves, this->G.activeCell(cell_16_31_10()), scaling);

    BOOST_REQUIRE_EQUAL(graphs.size(), std::size_t{ 1 });

    // Water saturations
    {
        const auto expect_sw = std::vector<double>{
            0.0, 1.0E-4, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4,
            0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9,
            0.95, 0.9999, 1.0,
        };

        const auto& sw = graphs.front().first;
        BOOST_REQUIRE_EQUAL(sw.size(), expect_sw.size());

        for (auto i = 0 * sw.size(); i < sw.size(); ++i) {
            BOOST_TEST_MESSAGE("Compare Sw[" << i << "] with Expected Value");
            BOOST_CHECK_CLOSE(sw[i], expect_sw[i], 1.0e-6);
        }
    }

    // Relative permeability values
    {
        const auto expect_kr = std::vector<double> {
            1.0, 0.999, 0.84782, 0.69746, 0.55717,
            0.43286, 0.32757, 0.24177, 0.17415, 0.12237,
            0.08374, 0.05565, 0.03572, 0.02199, 0.01284,
            0.00699, 0.00346, 0.00149, 0.00051, 0.00012,
            1.E-05, 0.0, 0.0,
        };

        const auto& kr = graphs.front().second;
        BOOST_REQUIRE_EQUAL(kr.size(), expect_kr.size());

        for (auto i = 0 * kr.size(); i < kr.size(); ++i) {
            BOOST_TEST_MESSAGE("Compare Krow[" << i << "] with Expected Value");
            BOOST_CHECK_CLOSE(kr[i], expect_kr[i], 1.0e-6);
        }
    }
}

BOOST_AUTO_TEST_CASE(Relperm_Water)
{
    auto curves = std::vector<Opm::ECLSaturationFunc::RawCurve>{};

    auto& curve = curves.emplace_back();
    curve.curve = Opm::ECLSaturationFunc::RawCurve::Function::RelPerm;
    curve.subsys = Opm::ECLSaturationFunc::RawCurve::SubSystem::OilWater;
    curve.thisPh = Opm::ECLPhaseIndex::Aqua;

    auto scaling = Opm::ECLSaturationFunc::SatFuncScaling{};
    scaling.enable = static_cast<unsigned char>(0);

    const auto graphs = this->satFunc
        .getSatFuncCurve(curves, this->G.activeCell(cell_16_31_10()), scaling);

    BOOST_REQUIRE_EQUAL(graphs.size(), std::size_t{ 1 });

    // Water saturations
    {
        const auto expect_sw = std::vector<double> {
            0.0, 0.0001, 0.05, 0.1, 0.15, 0.2, 0.25,
            0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65,
            0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0,
        };

        const auto& sw = graphs.front().first;
        BOOST_REQUIRE_EQUAL(sw.size(), expect_sw.size());

        for (auto i = 0 * sw.size(); i < sw.size(); ++i) {
            BOOST_TEST_MESSAGE("Compare Sw[" << i << "] with Expected Value");
            BOOST_CHECK_CLOSE(sw[i], expect_sw[i], 1.0e-6);
        }
    }

    // Relative permeability values
    {
        const auto expect_kr = std::vector<double> {
            0.0, 0.0, 0.00086, 0.00263, 0.00524, 0.00877, 0.01338,
            0.01927, 0.02672, 0.03608, 0.04781, 0.0625, 0.0809, 0.10394,
            0.13277, 0.16869, 0.21302, 0.26667, 0.32918, 0.39706,
            0.46103, 0.5,
        };

        const auto& kr = graphs.front().second;
        BOOST_REQUIRE_EQUAL(kr.size(), expect_kr.size());

        for (auto i = 0 * kr.size(); i < kr.size(); ++i) {
            BOOST_TEST_MESSAGE("Compare Krw[" << i << "] with Expected Value");
            BOOST_CHECK_CLOSE(kr[i], expect_kr[i], 1.0e-6);
        }
    }
}

BOOST_AUTO_TEST_CASE(CapPress_Gas_Oil)
{
    auto curves = std::vector<Opm::ECLSaturationFunc::RawCurve>{};

    auto& curve = curves.emplace_back();
    curve.curve = Opm::ECLSaturationFunc::RawCurve::Function::CapPress;
    curve.subsys = Opm::ECLSaturationFunc::RawCurve::SubSystem::OilGas;
    curve.thisPh = Opm::ECLPhaseIndex::Vapour;

    auto scaling = Opm::ECLSaturationFunc::SatFuncScaling{};
    scaling.enable = static_cast<unsigned char>(0);

    const auto graphs = this->satFunc
        .getSatFuncCurve(curves, this->G.activeCell(cell_16_31_10()), scaling);

    BOOST_REQUIRE_EQUAL(graphs.size(), std::size_t{ 1 });

    // Gas saturations
    {
        const auto expect_sg = std::vector<double> {
            0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45,
            0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95,
            0.9999, 1.0,
        };

        const auto& sg = graphs.front().first;
        BOOST_REQUIRE_EQUAL(sg.size(), expect_sg.size());

        for (auto i = 0 * sg.size(); i < sg.size(); ++i) {
            BOOST_TEST_MESSAGE("Compare Sg[" << i << "] with Expected Value");
            BOOST_CHECK_CLOSE(sg[i], expect_sg[i], 1.0e-6);
        }
    }

    // Gas/oil capillary pressure (Pg - Po) values
    {
        const auto expect_pc = std::vector<double> {
            0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        };

        const auto& pc = graphs.front().second;
        BOOST_REQUIRE_EQUAL(pc.size(), expect_pc.size());

        for (auto i = 0 * pc.size(); i < pc.size(); ++i) {
            BOOST_TEST_MESSAGE("Compare Pcgo[" << i << "] with Expected Value");
            BOOST_CHECK_CLOSE(pc[i], expect_pc[i], 1.0e-6);
        }
    }
}

BOOST_AUTO_TEST_CASE(CapPress_Oil_Water)
{
    auto curves = std::vector<Opm::ECLSaturationFunc::RawCurve>{};

    auto& curve = curves.emplace_back();
    curve.curve = Opm::ECLSaturationFunc::RawCurve::Function::CapPress;
    curve.subsys = Opm::ECLSaturationFunc::RawCurve::SubSystem::OilWater;
    curve.thisPh = Opm::ECLPhaseIndex::Aqua;

    auto scaling = Opm::ECLSaturationFunc::SatFuncScaling{};
    scaling.enable = static_cast<unsigned char>(0);

    const auto graphs = this->satFunc
        .getSatFuncCurve(curves, this->G.activeCell(cell_16_31_10()), scaling);

    BOOST_REQUIRE_EQUAL(graphs.size(), std::size_t{ 1 });

    // Water saturations
    {
        const auto expect_sw = std::vector<double> {
            0.0, 0.0001, 0.05, 0.1, 0.15, 0.2, 0.25,
            0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65,
            0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0,
        };

        const auto& sw = graphs.front().first;
        BOOST_REQUIRE_EQUAL(sw.size(), expect_sw.size());

        for (auto i = 0 * sw.size(); i < sw.size(); ++i) {
            BOOST_TEST_MESSAGE("Compare Sw[" << i << "] with Expected Value");
            BOOST_CHECK_CLOSE(sw[i], expect_sw[i], 1.0e-6);
        }
    }

    // Oil/water capillary pressure (Po - Pw) values
    {
        const auto expect_pc = std::vector<double> {
            375633.0, 375632.0, 186981.0, 123731.0, 91821.0,
            72451.0, 59341.0, 49811.0, 42511.0, 36691.0,
            31911.0, 27881.0, 24401.0, 21351.0, 18631.0,
            16161.0, 13901.0, 11801.0, 9831.0, 7961.0,
            6161.0, 4408.0,
        };

        const auto& pc = graphs.front().second;
        BOOST_REQUIRE_EQUAL(pc.size(), expect_pc.size());

        for (auto i = 0 * pc.size(); i < pc.size(); ++i) {
            BOOST_TEST_MESSAGE("Compare Pcow[" << i << "] with Expected Value");
            BOOST_CHECK_CLOSE(pc[i], expect_pc[i], 1.0e-6);
        }
    }
}

BOOST_AUTO_TEST_SUITE_END() // Unscaled

// --------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(Scaled)

BOOST_AUTO_TEST_CASE(Relperm_Gas)
{
    auto curves = std::vector<Opm::ECLSaturationFunc::RawCurve>{};

    auto& curve = curves.emplace_back();
    curve.curve = Opm::ECLSaturationFunc::RawCurve::Function::RelPerm;
    curve.subsys = Opm::ECLSaturationFunc::RawCurve::SubSystem::OilGas;
    curve.thisPh = Opm::ECLPhaseIndex::Vapour;

    // Default scaling mode is horizontal + vertical.
    const auto scaling = Opm::ECLSaturationFunc::SatFuncScaling{};

    const auto graphs = this->satFunc
        .getSatFuncCurve(curves, this->G.activeCell(cell_16_31_10()), scaling);

    BOOST_REQUIRE_EQUAL(graphs.size(), std::size_t{ 1 });

    // Gas saturations
    {
        const auto expect_sg = std::vector<double> {
            0.029999999329447746,
            0.072504249735861948,
            0.11500850014227613,
            0.15751275054869032,
            0.20001700095510452,
            0.24252125136151867,
            0.28502550176793290,
            0.32752975217434704,
            0.37003400258076130,
            0.41253825298717545,
            0.45504250339358959,
            0.49754675380000385,
            0.54005100420641805,
            0.58255525461283231,
            0.62505950501924634,
            0.66756375542566060,
            0.71006800583207486,
            0.75257225623848889,
            0.79507650664490315,
            0.83758075705131729,
            0.87999999895691872,
            0.94999998807907104,
        };

        const auto& sg = graphs.front().first;
        BOOST_REQUIRE_EQUAL(sg.size(), expect_sg.size());

        for (auto i = 0 * sg.size(); i < sg.size(); ++i) {
            BOOST_TEST_MESSAGE("Compare Sg[" << i << "] with Expected Value");
            BOOST_CHECK_CLOSE(sg[i], expect_sg[i], 1.0e-6);
        }
    }

    // Relative permeability values
    {
        const auto expect_kr = std::vector<double> {
            0.0, 0.001655, 0.006913, 0.016213, 0.02999, 0.048655, 0.072573,
            0.102046, 0.137287, 0.178402, 0.225368, 0.27803, 0.336093,
            0.399135, 0.466631, 0.538, 0.612665, 0.690169, 0.770395,
            0.854218, 0.9499, 0.95,
        };

        const auto& kr = graphs.front().second;
        BOOST_REQUIRE_EQUAL(kr.size(), expect_kr.size());

        for (auto i = 0 * kr.size(); i < kr.size(); ++i) {
            BOOST_TEST_MESSAGE("Compare Krg[" << i << "] with Expected Value");
            BOOST_CHECK_CLOSE(kr[i], expect_kr[i], 1.0e-6);
        }
    }
}

BOOST_AUTO_TEST_CASE(Relperm_Oil_in_Oil_Gas)
{
    auto curves = std::vector<Opm::ECLSaturationFunc::RawCurve>{};

    auto& curve = curves.emplace_back();
    curve.curve = Opm::ECLSaturationFunc::RawCurve::Function::RelPerm;
    curve.subsys = Opm::ECLSaturationFunc::RawCurve::SubSystem::OilGas;
    curve.thisPh = Opm::ECLPhaseIndex::Liquid;

    // Default scaling mode is horizontal + vertical.
    const auto scaling = Opm::ECLSaturationFunc::SatFuncScaling{};

    const auto graphs = this->satFunc
        .getSatFuncCurve(curves, this->G.activeCell(cell_16_31_10()), scaling);

    BOOST_REQUIRE_EQUAL(graphs.size(), std::size_t{ 1 });

    // Gas saturations
    {
        const auto expect_sg = std::vector<double> {
            0.029999999329447746,
            0.030085007830260579,
            0.072504249735862003,
            0.11500850014227604,
            0.15751275054869029,
            0.20001700095510444,
            0.24252125136151870,
            0.28502550176793284,
            0.32752975217434699,
            0.37003400258076125,
            0.41253825298717550,
            0.45504250339358965,
            0.49754675380000385,
            0.54005100420641794,
            0.58255525461283231,
            0.62505950501924645,
            0.66756375542566060,
            0.71006800583207486,
            0.75257225623848889,
            0.79507650664490315,
            0.83758075705131729,
            0.87999999895691872,
        };

        const auto& sg = graphs.front().first;
        BOOST_REQUIRE_EQUAL(sg.size(), expect_sg.size());

        for (auto i = 0 * sg.size(); i < sg.size(); ++i) {
            BOOST_TEST_MESSAGE("Compare Sg[" << i << "] with Expected Value");
            BOOST_CHECK_CLOSE(sg[i], expect_sg[i], 1.0e-6);
        }
    }

    // Relative permeability values
    {
        const auto expect_kr = std::vector<double>{
            1.0, 0.999613776, 0.806888, 0.633562, 0.485506, 0.364043,
            0.267589, 0.192992, 0.136554, 0.094671, 0.064151,
            0.042324, 0.027035, 0.016586, 0.009662, 0.005254,
            0.002597, 0.001117, 0.000384, 8.8E-05, 7.0E-06, 0.0,
        };

        const auto& kr = graphs.front().second;
        BOOST_REQUIRE_EQUAL(kr.size(), expect_kr.size());

        for (auto i = 0 * kr.size(); i < kr.size(); ++i) {
            BOOST_TEST_MESSAGE("Compare Krog[" << i << "] with Expected Value");
            BOOST_CHECK_CLOSE(kr[i], expect_kr[i], 1.0e-6);
        }
    }
}

BOOST_AUTO_TEST_CASE(Relperm_Oil_in_Oil_Water)
{
    auto curves = std::vector<Opm::ECLSaturationFunc::RawCurve>{};

    auto& curve = curves.emplace_back();
    curve.curve = Opm::ECLSaturationFunc::RawCurve::Function::RelPerm;
    curve.subsys = Opm::ECLSaturationFunc::RawCurve::SubSystem::OilWater;
    curve.thisPh = Opm::ECLPhaseIndex::Liquid;

    // Default scaling mode is horizontal + vertical
    const auto scaling = Opm::ECLSaturationFunc::SatFuncScaling{};

    const auto graphs = this->satFunc
        .getSatFuncCurve(curves, this->G.activeCell(cell_16_31_10()), scaling);

    BOOST_REQUIRE_EQUAL(graphs.size(), std::size_t{ 1 });

    // Water saturations
    {
        const auto expect_sw = std::vector<double> {
            0.050000000745058060,
            0.12999999523162842,
            0.16693338238494260,
            0.20394078434217111,
            0.24094818629939974,
            0.27795558825662836,
            0.31496299021385699,
            0.35197039217108561,
            0.38897779412831412,
            0.42598519608554275,
            0.46299259804277137,
            0.50000000000000000,
            0.53700740195722863,
            0.57401480391445725,
            0.61102220587168588,
            0.64802960782891439,
            0.68503700978614301,
            0.72204441174337164,
            0.75905181370060015,
            0.79605921565782878,
            0.83306661761505740,
            0.87000000476837158,
        };

        const auto& sw = graphs.front().first;
        BOOST_REQUIRE_EQUAL(sw.size(), expect_sw.size());

        for (auto i = 0 * sw.size(); i < sw.size(); ++i) {
            BOOST_TEST_MESSAGE("Compare Sw[" << i << "] with Expected Value");
            BOOST_CHECK_CLOSE(sw[i], expect_sw[i], 1.0e-6);
        }
    }

    // Relative permeability values
    {
        const auto expect_kr = std::vector<double> {
            1.0, 0.999, 0.84782, 0.69746, 0.55717,
            0.43286, 0.32757, 0.24177, 0.17415, 0.12237,
            0.08374, 0.05565, 0.03572, 0.02199, 0.01284,
            0.00699, 0.00346, 0.00149, 0.00051, 0.00012,
            1.E-05, 0.0,
        };

        const auto& kr = graphs.front().second;
        BOOST_REQUIRE_EQUAL(kr.size(), expect_kr.size());

        for (auto i = 0 * kr.size(); i < kr.size(); ++i) {
            BOOST_TEST_MESSAGE("Compare Krow[" << i << "] with Expected Value");
            BOOST_CHECK_CLOSE(kr[i], expect_kr[i], 1.0e-6);
        }
    }
}

BOOST_AUTO_TEST_CASE(Relperm_Water)
{
    auto curves = std::vector<Opm::ECLSaturationFunc::RawCurve>{};

    auto& curve = curves.emplace_back();
    curve.curve = Opm::ECLSaturationFunc::RawCurve::Function::RelPerm;
    curve.subsys = Opm::ECLSaturationFunc::RawCurve::SubSystem::OilWater;
    curve.thisPh = Opm::ECLPhaseIndex::Aqua;

    // Default scaling mode is horizontal + vertical.
    const auto scaling = Opm::ECLSaturationFunc::SatFuncScaling{};

    const auto graphs = this->satFunc
        .getSatFuncCurve(curves, this->G.activeCell(cell_16_31_10()), scaling);

    BOOST_REQUIRE_EQUAL(graphs.size(), std::size_t{ 1 });

    // Water saturations
    {
        const auto expect_sw = std::vector<double> {
            0.12999999523162842,
            0.16693338238494257,
            0.20394078434217117,
            0.24094818629939979,
            0.27795558825662836,
            0.31496299021385699,
            0.35197039217108561,
            0.38897779412831418,
            0.42598519608554281,
            0.46299259804277143,
            0.50000000000000000,
            0.53700740195722863,
            0.57401480391445725,
            0.61102220587168588,
            0.64802960782891439,
            0.68503700978614301,
            0.72204441174337164,
            0.75905181370060026,
            0.79605921565782889,
            0.83306661761505740,
            1.0000000000000000,
        };

        const auto& sw = graphs.front().first;
        BOOST_REQUIRE_EQUAL(sw.size(), expect_sw.size());

        for (auto i = 0 * sw.size(); i < sw.size(); ++i) {
            BOOST_TEST_MESSAGE("Compare Sw[" << i << "] with Expected Value");
            BOOST_CHECK_CLOSE(sw[i], expect_sw[i], 1.0e-6);
        }
    }

    // Relative permeability values
    {
        const auto expect_kr = std::vector<double> {
            0.0, 0.00086, 0.00263, 0.00524, 0.00877, 0.01338,
            0.01927, 0.02672, 0.03608, 0.04781, 0.0625, 0.0809, 0.10394,
            0.13277, 0.16869, 0.21302, 0.26667, 0.32918, 0.39706,
            0.46103, 0.5,
        };

        const auto& kr = graphs.front().second;
        BOOST_REQUIRE_EQUAL(kr.size(), expect_kr.size());

        for (auto i = 0 * kr.size(); i < kr.size(); ++i) {
            BOOST_TEST_MESSAGE("Compare Krw[" << i << "] with Expected Value");
            BOOST_CHECK_CLOSE(kr[i], expect_kr[i], 1.0e-6);
        }
    }
}

BOOST_AUTO_TEST_CASE(CapPress_Gas_Oil)
{
    auto curves = std::vector<Opm::ECLSaturationFunc::RawCurve>{};

    auto& curve = curves.emplace_back();
    curve.curve = Opm::ECLSaturationFunc::RawCurve::Function::CapPress;
    curve.subsys = Opm::ECLSaturationFunc::RawCurve::SubSystem::OilGas;
    curve.thisPh = Opm::ECLPhaseIndex::Vapour;

    // Default scaling mode is horizontal + vertical.
    const auto scaling = Opm::ECLSaturationFunc::SatFuncScaling{};

    const auto graphs = this->satFunc
        .getSatFuncCurve(curves, this->G.activeCell(cell_16_31_10()), scaling);

    BOOST_REQUIRE_EQUAL(graphs.size(), std::size_t{ 1 });

    // Gas saturations
    {
        const auto expect_sg = std::vector<double> {
            0.0000000000000000,
            0.047499999403953552,
            0.094999998807907104,
            0.14249999821186066,
            0.18999999761581421,
            0.23749999701976776,
            0.28499999642372131,
            0.33249999582767487,
            0.37999999523162842,
            0.42749999463558197,
            0.47499999403953552,
            0.52249999344348907,
            0.56999999284744263,
            0.61749999225139618,
            0.66499999165534973,
            0.71249999105930328,
            0.75999999046325684,
            0.80749998986721039,
            0.85499998927116394,
            0.90249998867511749,
            0.94990498808026314,
            0.94999998807907104,
        };

        const auto& sg = graphs.front().first;
        BOOST_REQUIRE_EQUAL(sg.size(), expect_sg.size());

        for (auto i = 0 * sg.size(); i < sg.size(); ++i) {
            BOOST_TEST_MESSAGE("Compare Sg[" << i << "] with Expected Value");
            BOOST_CHECK_CLOSE(sg[i], expect_sg[i], 1.0e-6);
        }
    }

    // Gas/oil capillary pressure (Pg - Po) values
    {
        const auto expect_pc = std::vector<double> {
            0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        };

        const auto& pc = graphs.front().second;
        BOOST_REQUIRE_EQUAL(pc.size(), expect_pc.size());

        for (auto i = 0 * pc.size(); i < pc.size(); ++i) {
            BOOST_TEST_MESSAGE("Compare Pcgo[" << i << "] with Expected Value");
            BOOST_CHECK_CLOSE(pc[i], expect_pc[i], 1.0e-6);
        }
    }
}

BOOST_AUTO_TEST_CASE(CapPress_Oil_Water)
{
    auto curves = std::vector<Opm::ECLSaturationFunc::RawCurve>{};

    auto& curve = curves.emplace_back();
    curve.curve = Opm::ECLSaturationFunc::RawCurve::Function::CapPress;
    curve.subsys = Opm::ECLSaturationFunc::RawCurve::SubSystem::OilWater;
    curve.thisPh = Opm::ECLPhaseIndex::Aqua;

    // Default scaling mode is horizontal + vertical.
    const auto scaling = Opm::ECLSaturationFunc::SatFuncScaling{};

    const auto graphs = this->satFunc
        .getSatFuncCurve(curves, this->G.activeCell(cell_16_31_10()), scaling);

    BOOST_REQUIRE_EQUAL(graphs.size(), std::size_t{ 1 });

    // Water saturations
    {
        const auto expect_sw = std::vector<double> {
            0.050000000745058060,
            0.050095000744983555,
            0.097500000707805151,
            0.14500000067055224,
            0.19250000063329933,
            0.24000000059604645,
            0.28750000055879354,
            0.33500000052154061,
            0.38250000048428773,
            0.43000000044703485,
            0.47750000040978197,
            0.52500000037252903,
            0.57250000033527615,
            0.62000000029802316,
            0.66750000026077039,
            0.71500000022351740,
            0.76250000018626451,
            0.81000000014901163,
            0.85750000011175864,
            0.90500000007450587,
            0.95250000003725288,
            1.0000000000000000,
        };

        const auto& sw = graphs.front().first;
        BOOST_REQUIRE_EQUAL(sw.size(), expect_sw.size());

        for (auto i = 0 * sw.size(); i < sw.size(); ++i) {
            BOOST_TEST_MESSAGE("Compare Sw[" << i << "] with Expected Value");
            BOOST_CHECK_CLOSE(sw[i], expect_sw[i], 1.0e-6);
        }
    }

    // Oil/water capillary pressure (Po - Pw) values
    {
        const auto expect_pc = std::vector<double> {
            375633.0, 375632.0, 186981.0, 123731.0, 91821.0,
            72451.0, 59341.0, 49811.0, 42511.0, 36691.0,
            31911.0, 27881.0, 24401.0, 21351.0, 18631.0,
            16161.0, 13901.0, 11801.0, 9831.0, 7961.0,
            6161.0, 4408.0,
        };

        const auto& pc = graphs.front().second;
        BOOST_REQUIRE_EQUAL(pc.size(), expect_pc.size());

        for (auto i = 0 * pc.size(); i < pc.size(); ++i) {
            BOOST_TEST_MESSAGE("Compare Pcow[" << i << "] with Expected Value");
            BOOST_CHECK_CLOSE(pc[i], expect_pc[i], 1.0e-6);
        }
    }
}

BOOST_AUTO_TEST_SUITE_END() // Scaled

BOOST_AUTO_TEST_SUITE_END() // Saturation_Functions
