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
    Settings::physical_parameters physical_parameters = Settings::parse_physical_parameters(doc.RootElement());
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
    Settings::physical_parameters physical_parameters = Settings::parse_physical_parameters(doc.RootElement());
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
    Settings::physical_parameters physical_parameters = Settings::parse_physical_parameters(doc.RootElement());
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
    EXPECT_THROW(Settings::parse_physical_parameters(doc.RootElement()), Settings::config_error);
}

TEST(SettingsTest, requiredDomainParameters) {
    std::string xml = R"(
<ARTSS>
  <domain_parameters enable_computational_domain="No">
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
    Settings::domain_parameters domain_parameters = Settings::parse_domain_parameters(doc.RootElement());
    EXPECT_FALSE(domain_parameters.enable_computational_domain);
    EXPECT_EQ(domain_parameters.start_coords_PD[CoordinateAxis::X], 1.0);
    EXPECT_EQ(domain_parameters.start_coords_PD[CoordinateAxis::Y], -12.0);
    EXPECT_EQ(domain_parameters.start_coords_PD[CoordinateAxis::Z], -1.0);
    EXPECT_EQ(domain_parameters.end_coords_PD[CoordinateAxis::X], 3.0);
    EXPECT_EQ(domain_parameters.end_coords_PD[CoordinateAxis::Y], -4.0);
    EXPECT_EQ(domain_parameters.end_coords_PD[CoordinateAxis::Z], 0);
    EXPECT_EQ(domain_parameters.number_of_inner_cells[CoordinateAxis::X], 20);
    EXPECT_EQ(domain_parameters.number_of_inner_cells[CoordinateAxis::Y], 21);
    EXPECT_EQ(domain_parameters.number_of_inner_cells[CoordinateAxis::Z], 1);
}

TEST(SettingsTest, requiredLoggingParameters) {
    std::string xml = R"(
<ARTSS>
    <logging file="tmp.log" level="debug">
    </logging>
</ARTSS>)";
    tinyxml2::XMLDocument doc;
    doc.Parse(xml.c_str());
    Settings::logging_parameters logging_parameters = Settings::parse_logging_parameters(doc.RootElement());
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
    Settings::visualisation_parameters visualisation_parameters = Settings::parse_visualisation_parameters(
            doc.RootElement());
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
    Settings::visualisation_parameters visualisation_parameters = Settings::parse_visualisation_parameters(
            doc.RootElement());
    EXPECT_FALSE(visualisation_parameters.save_csv);
    EXPECT_FALSE(visualisation_parameters.save_vtk);
    EXPECT_EQ(visualisation_parameters.vtk_nth_plot, 0);
    EXPECT_EQ(visualisation_parameters.csv_nth_plot, 0);
}

TEST(SettingsTest, requiredInitialConditionsParameters) {
    std::string xml = R"(
<ARTSS>
    <initial_conditions usr_fct="Uniform" random="No">
        <val> 10 </val>
    </initial_conditions>
</ARTSS>)";
    tinyxml2::XMLDocument doc;
    doc.Parse(xml.c_str());
    Settings::initial_conditions_parameters initial_conditions_parameters = Settings::parse_initial_conditions_parameters(
            doc.RootElement());
    EXPECT_FALSE(initial_conditions_parameters.random);
    EXPECT_EQ(initial_conditions_parameters.usr_fct, "Uniform");
    EXPECT_EQ(std::get<Settings::uniform>(initial_conditions_parameters.ic).value, 10);
}

TEST(SettingsTest, requiredInitialConditionsParameters2) {
    std::string xml = R"(
<ARTSS>
    <initial_conditions usr_fct="Uniform" random="Yes">
        <val> 1 </val>
        <test> 1.1 </test>
        <random absolute="Yes" custom_seed="Yes" custom_steps="Yes">
            <seed> 10 </seed>
            <step_size> 0.1 </step_size>
            <range> 1 </range>
        </random>
    </initial_conditions>
</ARTSS>)";
    tinyxml2::XMLDocument doc;
    doc.Parse(xml.c_str());
    Settings::initial_conditions_parameters initial_conditions_parameters = Settings::parse_initial_conditions_parameters(
            doc.RootElement());
    EXPECT_TRUE(initial_conditions_parameters.random);
    EXPECT_EQ(initial_conditions_parameters.usr_fct, "Uniform");
    EXPECT_EQ(std::get<Settings::uniform>(initial_conditions_parameters.ic).value, 1);
    EXPECT_TRUE(initial_conditions_parameters.random_parameters.absolute);
    EXPECT_TRUE(initial_conditions_parameters.random_parameters.custom_seed);
    EXPECT_TRUE(initial_conditions_parameters.random_parameters.custom_steps);
    EXPECT_EQ(initial_conditions_parameters.random_parameters.range, 1);
    EXPECT_EQ(initial_conditions_parameters.random_parameters.seed, 10);
    EXPECT_EQ(initial_conditions_parameters.random_parameters.step_size, 0.1);
}

TEST(SettingsTest, requiredObstaclesParameters) {
    std::string xml = R"(
<ARTSS>
    <obstacles enabled="No" />
</ARTSS>)";
    tinyxml2::XMLDocument doc;
    doc.Parse(xml.c_str());
    Settings::obstacles_parameters obstacles_parameters = Settings::parse_obstacles_parameters(doc.RootElement());
    EXPECT_FALSE(obstacles_parameters.enabled);
}

TEST(SettingsTest, requiredSurfacesParameters) {
    std::string xml = R"(
<ARTSS>
    <surfaces enabled="No" />
</ARTSS>)";
    tinyxml2::XMLDocument doc;
    doc.Parse(xml.c_str());
    Settings::surfaces_parameters surfaces_parameters = Settings::parse_surfaces_parameters(doc.RootElement());
    EXPECT_FALSE(surfaces_parameters.enabled);
}

TEST(SettingsTest, boundaries) {
    std::string xml = R"(
<ARTSS>
    <boundaries>
        <boundary field="u,v,w" patch="front,back,bottom,top,left,right" type="dirichlet" value="0.0" />
        <boundary field="p" patch="front,back,bottom,top,left,right" type="neumann" value="1.0" />
        <boundary field="T" patch="front,back,top,left,right" type="dirichlet" value="299.14" />
        <boundary field="T" patch="bottom" type="periodic" value="-10.0" />
    </boundaries>
</ARTSS>)";
    tinyxml2::XMLDocument doc;
    doc.Parse(xml.c_str());
    Settings::boundaries_parameters boundaries_parameters = Settings::parse_boundaries_parameters(doc.RootElement());

    // field types
    EXPECT_NE(std::find(boundaries_parameters.boundaries[0].field_type.begin(),
                        boundaries_parameters.boundaries[0].field_type.end(),
                        FieldType::U),
              boundaries_parameters.boundaries[0].field_type.end());
    EXPECT_NE(std::find(boundaries_parameters.boundaries[0].field_type.begin(),
                        boundaries_parameters.boundaries[0].field_type.end(),
                        FieldType::V),
              boundaries_parameters.boundaries[0].field_type.end());
    EXPECT_NE(std::find(boundaries_parameters.boundaries[0].field_type.begin(),
                        boundaries_parameters.boundaries[0].field_type.end(),
                        FieldType::W),
              boundaries_parameters.boundaries[0].field_type.end());
    EXPECT_NE(std::find(boundaries_parameters.boundaries[1].field_type.begin(),
                        boundaries_parameters.boundaries[1].field_type.end(),
                        FieldType::P),
              boundaries_parameters.boundaries[1].field_type.end());
    EXPECT_NE(std::find(boundaries_parameters.boundaries[2].field_type.begin(),
                        boundaries_parameters.boundaries[2].field_type.end(),
                        FieldType::T),
              boundaries_parameters.boundaries[2].field_type.end());
    EXPECT_NE(std::find(boundaries_parameters.boundaries[3].field_type.begin(),
                        boundaries_parameters.boundaries[3].field_type.end(),
                        FieldType::T),
              boundaries_parameters.boundaries[3].field_type.end());

    // patches
    for (size_t p = 0; p < number_of_patches; p++) {
        auto patch = Patch(p);
        EXPECT_NE(std::find(boundaries_parameters.boundaries[0].patch.begin(),
                            boundaries_parameters.boundaries[0].patch.end(),
                            patch),
                  boundaries_parameters.boundaries[0].patch.end());
        EXPECT_NE(std::find(boundaries_parameters.boundaries[1].patch.begin(),
                            boundaries_parameters.boundaries[1].patch.end(),
                            patch),
                  boundaries_parameters.boundaries[1].patch.end());
    }
    for (Patch patch: {Patch::LEFT, Patch::RIGHT, Patch::TOP, Patch::FRONT, Patch::BACK}) {
        EXPECT_NE(std::find(boundaries_parameters.boundaries[2].patch.begin(),
                            boundaries_parameters.boundaries[2].patch.end(),
                            patch),
                  boundaries_parameters.boundaries[2].patch.end());
    }
    EXPECT_NE(std::find(boundaries_parameters.boundaries[3].patch.begin(),
                        boundaries_parameters.boundaries[3].patch.end(),
                        Patch::BOTTOM),
              boundaries_parameters.boundaries[3].patch.end());

    // boundary condition
    EXPECT_EQ(boundaries_parameters.boundaries[0].boundary_condition, BoundaryCondition::DIRICHLET);
    EXPECT_EQ(boundaries_parameters.boundaries[1].boundary_condition, BoundaryCondition::NEUMANN);
    EXPECT_EQ(boundaries_parameters.boundaries[2].boundary_condition, BoundaryCondition::DIRICHLET);
    EXPECT_EQ(boundaries_parameters.boundaries[3].boundary_condition, BoundaryCondition::PERIODIC);

    // value
    EXPECT_DOUBLE_EQ(boundaries_parameters.boundaries[0].value, 0);
    EXPECT_DOUBLE_EQ(boundaries_parameters.boundaries[1].value, 1);
    EXPECT_DOUBLE_EQ(boundaries_parameters.boundaries[2].value, 299.14);
    EXPECT_DOUBLE_EQ(boundaries_parameters.boundaries[3].value, -10);
}
