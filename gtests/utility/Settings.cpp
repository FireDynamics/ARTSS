/// \file       Settings.cpp
/// \brief      
/// \date       Dec 10, 2021
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2021> Forschungszentrum Juelich All rights reserved.
//

#include <gtest/gtest.h>

#include "src/utility/settings/Settings.h"

TEST(SettingsTest, goodCase) {
    std::string xml = R"(
<ARTSS>
  <physical_parameters>
    <t_end> 1.0 </t_end>
    <dt> 0.1 </dt>
  </physical_parameters>
</ARTSS>)";
    tinyxml2::XMLDocument doc;
    doc.Parse(xml.c_str());
    Settings::physical_parameters physical_parameters = Settings::parse_physical_parameters(doc);
    EXPECT_EQ(physical_parameters.t_end, 1.0);
    EXPECT_EQ(physical_parameters.dt, 0.1);
}

TEST(SettingsTest, goodCaseOptional) {
    std::string xml = R"(
<ARTSS>
  <physical_parameters>
    <t_end> 1.0 </t_end>
    <dt> 0.1 </dt>
    <beta> 0.3 </beta>
  </physical_parameters>
</ARTSS>)";
    tinyxml2::XMLDocument doc;
    doc.Parse(xml.c_str());
    Settings::physical_parameters physical_parameters = Settings::parse_physical_parameters(doc);
    EXPECT_EQ(physical_parameters.t_end, 1.0);
    EXPECT_EQ(physical_parameters.dt, 0.1);
    EXPECT_EQ(physical_parameters.beta, 0.3);
}

TEST(SettingsTest, unknownAttribute) {
    std::string xml = R"(
<ARTSS>
  <physical_parameters>
    <t_end> 1.0 </t_end>
    <dt> 0.1 </dt>
    <bs> 1 </bs>
  </physical_parameters>
</ARTSS>)";
    tinyxml2::XMLDocument doc;
    doc.Parse(xml.c_str());
    Settings::physical_parameters physical_parameters = Settings::parse_physical_parameters(doc);
    EXPECT_EQ(physical_parameters.t_end, 1.0);
    EXPECT_EQ(physical_parameters.dt, 0.1);
//TODO expect what?
}

TEST(SettingsTest, requiredAttributeMissing) {
    std::string xml = R"(
<ARTSS>
  <physical_parameters>
    <t_end> 1.0 </t_end>
  </physical_parameters>
</ARTSS>)";
    tinyxml2::XMLDocument doc;
    doc.Parse(xml.c_str());
    EXPECT_THROW(Settings::parse_physical_parameters(doc), Settings::config_error);
}

TEST(SettingsTest, requiredDomainParameters) {
    std::string xml = R"(
<ARTSS>
  <domain_parameters enable_computational_domain="Yes">
    <X1> 1.0 </X1>
    <Y1> -12.0 </Y1>
    <Z1> -1.0 </Z1>
    <X2> 3.0 </X2>
    <Y2> -4.0 </Y2>
    <Z2> 0 </Z2>
    <nx> 20 </nx>
    <ny> 21 </ny>
    <nz> 1 </nz>
  </domain_parameters>
</ARTSS>)";
    tinyxml2::XMLDocument doc;
    doc.Parse(xml.c_str());
    Settings::domain_parameters domain_parameters = Settings::parse_domain_parameters(doc);
    EXPECT_EQ(domain_parameters.enable_computational_domain, 1.0);
    EXPECT_EQ(domain_parameters.X1, 1.0);
    EXPECT_EQ(domain_parameters.Y1, -12.0);
    EXPECT_EQ(domain_parameters.Z1, -1.0);
    EXPECT_EQ(domain_parameters.X2, 3.0);
    EXPECT_EQ(domain_parameters.Y2, -4.0);
    EXPECT_EQ(domain_parameters.Z2, 0);
    EXPECT_EQ(domain_parameters.nx, 20);
    EXPECT_EQ(domain_parameters.ny, 21);
    EXPECT_EQ(domain_parameters.nz, 1);
}

TEST(SettingsTest, requiredLoggingParameters) {
    std::string xml = R"(
<ARTSS>
    <logging file="tmp.log" level="debug">
    </logging>
</ARTSS>)";
    tinyxml2::XMLDocument doc;
    doc.Parse(xml.c_str());
    Settings::logging_parameters logging_parameters = Settings::parse_logging_parameters(doc);
    EXPECT_EQ(logging_parameters.file, "tmp.log");
    EXPECT_EQ(logging_parameters.level, "debug");
}
TEST(SettingsTest, requiredVisualisationParameters) {
    std::string xml = R"(
<ARTSS>
    <visualisation save_vtk="Yes" save_csv="Yes">
        <vtk_nth_plot> 10 </vtk_nth_plot>
        <csv_nth_plot> 21 </csv_nth_plot>
    </visualisation>
</ARTSS>)";
    tinyxml2::XMLDocument doc;
    doc.Parse(xml.c_str());
    Settings::visualisation_parameters visualisation_parameters = Settings::parse_visualisation_parameters(doc);
    EXPECT_TRUE(visualisation_parameters.save_csv);
    EXPECT_TRUE(visualisation_parameters.save_vtk);
    EXPECT_EQ(visualisation_parameters.vtk_nth_plot, 10);
    EXPECT_EQ(visualisation_parameters.csv_nth_plot, 21);
}
TEST(SettingsTest, optionalVisualisationParameters) {
    std::string xml = R"(
<ARTSS>
    <visualisation save_vtk="No" save_csv="No">
    </visualisation>
</ARTSS>)";
    tinyxml2::XMLDocument doc;
    doc.Parse(xml.c_str());
    Settings::visualisation_parameters visualisation_parameters = Settings::parse_visualisation_parameters(doc);
    EXPECT_FALSE(visualisation_parameters.save_csv);
    EXPECT_FALSE(visualisation_parameters.save_vtk);
    EXPECT_EQ(visualisation_parameters.vtk_nth_plot, 1);
    EXPECT_EQ(visualisation_parameters.csv_nth_plot, 1);
}
TEST(SettingsTest, requiredInitialConditionsParameters) {
    std::string xml = R"(
<ARTSS>
    <initial_conditions usr_fct="Uniform" random="No">
    </initial_conditions>
</ARTSS>)";
    tinyxml2::XMLDocument doc;
    doc.Parse(xml.c_str());
    Settings::initial_conditions_parameters initial_conditions_parameters = Settings::parse_initial_conditions_parameters(doc);
    EXPECT_FALSE(initial_conditions_parameters.random);
    EXPECT_EQ(initial_conditions_parameters.usr_fct, "Uniform");
}
TEST(SettingsTest, requiredObstaclesParameters) {
    std::string xml = R"(
<ARTSS>
    <obstacles enabled="No" />
</ARTSS>)";
    tinyxml2::XMLDocument doc;
    doc.Parse(xml.c_str());
    Settings::obstacles_parameters obstacles_parameters = Settings::parse_obstacles_parameters(doc);
    EXPECT_FALSE(obstacles_parameters.enabled);
}
TEST(SettingsTest, requiredSurfacesParameters) {
    std::string xml = R"(
<ARTSS>
    <surfaces enabled="No" />
</ARTSS>)";
    tinyxml2::XMLDocument doc;
    doc.Parse(xml.c_str());
    Settings::surfaces_parameters surfaces_parameters = Settings::parse_surfaces_parameters(doc);
    EXPECT_FALSE(surfaces_parameters.enabled);
}
