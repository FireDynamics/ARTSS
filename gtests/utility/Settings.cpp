/// \file       Settings.cpp
/// \brief      
/// \date       Dec 10, 2021
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2021> Forschungszentrum Juelich All rights reserved.
//

#include <gtest/gtest.h>

#include "src/utility/settings/Settings.h"
#include "src/solver/SolverSelection.h"

using namespace Settings::initial_conditions;

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
    EXPECT_DOUBLE_EQ(physical_parameters.t_end, 1.0);
    EXPECT_DOUBLE_EQ(physical_parameters.dt, 0.1);
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
    EXPECT_DOUBLE_EQ(physical_parameters.t_end, 1.0);
    EXPECT_DOUBLE_EQ(physical_parameters.dt, 0.1);
    EXPECT_DOUBLE_EQ(physical_parameters.beta, 0.3);
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
    EXPECT_DOUBLE_EQ(physical_parameters.t_end, 1.0);
    EXPECT_DOUBLE_EQ(physical_parameters.dt, 0.1);
//TODO maybe warning?
}
// TODO wrong casts -> string instead of numeral for example
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
}

TEST(SettingsTest, requiredInitialConditionsParameters) {
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
    EXPECT_EQ(std::get<uniform>(initial_conditions_parameters.ic.value()).value, 1);
    EXPECT_TRUE(initial_conditions_parameters.random_parameters.value().absolute);
    EXPECT_TRUE(initial_conditions_parameters.random_parameters.value().custom_seed);
    EXPECT_TRUE(initial_conditions_parameters.random_parameters.value().custom_steps);
    EXPECT_EQ(initial_conditions_parameters.random_parameters.value().range, 1);
    EXPECT_EQ(initial_conditions_parameters.random_parameters.value().seed, 10);
    EXPECT_EQ(initial_conditions_parameters.random_parameters.value().step_size, 0.1);
}
namespace initial_conditions {
    TEST(SettingsTest, uniform) {
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
        EXPECT_EQ(std::get<uniform>(initial_conditions_parameters.ic.value()).value, 10);
    }

    TEST(SettingsTest, expSinusProd) {
        std::string xml = R"(
<ARTSS>
    <initial_conditions usr_fct="ExpSinusProd"  random="No">
        <l> 111.0 </l>
    </initial_conditions>
</ARTSS>)";
        tinyxml2::XMLDocument doc;
        doc.Parse(xml.c_str());
        Settings::initial_conditions_parameters initial_conditions_parameters = Settings::parse_initial_conditions_parameters(
                doc.RootElement());
        auto esp = std::get<exp_sinus_prod>(initial_conditions_parameters.ic.value());
        EXPECT_DOUBLE_EQ(esp.l, 111);
    }
    TEST(SettingsTest, hat) {
        std::string xml = R"(
<ARTSS>
    <initial_conditions usr_fct="Hat"  random="No">
        <x1> 0.5 </x1>
        <x2> 1.0 </x2>
        <y1> -0.5 </y1>
        <y2> 1.1 </y2>
        <z1> 0.05 </z1>
        <z2> 10 </z2>
        <val_in> 2.0 </val_in>
        <val_out> 1.0 </val_out>
    </initial_conditions>
</ARTSS>)";
        tinyxml2::XMLDocument doc;
        doc.Parse(xml.c_str());
        Settings::initial_conditions_parameters initial_conditions_parameters = Settings::parse_initial_conditions_parameters(
                doc.RootElement());
        auto h = std::get<hat>(initial_conditions_parameters.ic.value());
        EXPECT_DOUBLE_EQ(h.start_coords[CoordinateAxis::X], 0.5);
        EXPECT_DOUBLE_EQ(h.start_coords[CoordinateAxis::Y], -0.5);
        EXPECT_DOUBLE_EQ(h.start_coords[CoordinateAxis::Z], 0.05);
        EXPECT_DOUBLE_EQ(h.end_coords[CoordinateAxis::X], 1);
        EXPECT_DOUBLE_EQ(h.end_coords[CoordinateAxis::Y], 1.1);
        EXPECT_DOUBLE_EQ(h.end_coords[CoordinateAxis::Z], 10);
        EXPECT_EQ(h.val_in, 2.);
        EXPECT_EQ(h.val_out, 1);
    }
    TEST(SettingsTest, gaussBubble) {
        std::string xml = R"(
<ARTSS>
    <initial_conditions usr_fct="GaussBubble"  random="No">  <!-- Gaussian function  -->
        <u_lin> 0.05 </u_lin>        <!-- x-velocity in linear case  -->
        <v_lin> 0.5 </v_lin>        <!-- y-velocity in linear case  -->
        <w_lin> 0.25 </w_lin>       <!-- z-velocity in linear case  -->
        <x_shift> 1.125 </x_shift>  <!-- x_shift of Gauss Bubble in domain  -->
        <y_shift> 1.025 </y_shift>  <!-- y_shift of Gauss Bubble in domain  -->
        <z_shift> 0.5 </z_shift>    <!-- z_shift of Gauss Bubble in domain  -->
        <l> 0.03125 </l>            <!-- sigma in Gaussian -->
    </initial_conditions>
</ARTSS>)";
        tinyxml2::XMLDocument doc;
        doc.Parse(xml.c_str());
        Settings::initial_conditions_parameters initial_conditions_parameters = Settings::parse_initial_conditions_parameters(
                doc.RootElement());
        auto gb = std::get<gauss_bubble>(initial_conditions_parameters.ic.value());
        EXPECT_DOUBLE_EQ(gb.l, 0.03125);
        EXPECT_DOUBLE_EQ(gb.velocity_lin[CoordinateAxis::X], 0.05);
        EXPECT_DOUBLE_EQ(gb.velocity_lin[CoordinateAxis::Y], 0.5);
        EXPECT_DOUBLE_EQ(gb.velocity_lin[CoordinateAxis::Z], 0.25);
        EXPECT_DOUBLE_EQ(gb.shift[CoordinateAxis::X], 1.125);
        EXPECT_DOUBLE_EQ(gb.shift[CoordinateAxis::Y], 1.025);
        EXPECT_DOUBLE_EQ(gb.shift[CoordinateAxis::Z], 0.5);
    }
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

TEST(SettingsTest, advectionSolverSemiLagrangian) {
    std::string xml = R"(
<ARTSS>
<solver description="AdvectionSolver" >
    <advection type="SemiLagrangian" field="u,v,w" />
    <solution available="No"/>
</solver>
</ARTSS>)";
    tinyxml2::XMLDocument doc;
    doc.Parse(xml.c_str());
    Settings::solver_parameters solver_parameters = Settings::parse_solver_parameters(doc.RootElement());
    EXPECT_EQ(solver_parameters.advection.type, AdvectionMethods::SemiLagrangian);

    EXPECT_NE(std::find(solver_parameters.advection.fields.begin(),
                        solver_parameters.advection.fields.end(),
                        FieldType::U),
              solver_parameters.advection.fields.end());
    EXPECT_NE(std::find(solver_parameters.advection.fields.begin(),
                        solver_parameters.advection.fields.end(),
                        FieldType::V),
              solver_parameters.advection.fields.end());
    EXPECT_NE(std::find(solver_parameters.advection.fields.begin(),
                        solver_parameters.advection.fields.end(),
                        FieldType::W),
              solver_parameters.advection.fields.end());

    EXPECT_FALSE(solver_parameters.solution.analytical_solution);
}

namespace diffusion_solver {
    TEST(SettingsTest, diffusionSolverJacobi) {
        std::string xml = R"(
<ARTSS>
<solver description="DiffusionSolver" >
    <diffusion type="Jacobi" field="u,v,w" >
        <max_iter> 77 </max_iter>
        <tol_res> 1e-07 </tol_res>
        <w> 1 </w>
    </diffusion>
    <solution available="Yes">
        <tol> 1e-4 </tol>
    </solution>
</solver>
</ARTSS>)";
        tinyxml2::XMLDocument doc;
        doc.Parse(xml.c_str());
        Settings::solver_parameters solver_parameters = Settings::parse_solver_parameters(doc.RootElement());
        EXPECT_EQ(solver_parameters.diffusion.type, DiffusionMethods::Jacobi);

        EXPECT_NE(std::find(solver_parameters.diffusion.fields.begin(),
                            solver_parameters.diffusion.fields.end(),
                            FieldType::U),
                  solver_parameters.diffusion.fields.end());
        EXPECT_NE(std::find(solver_parameters.diffusion.fields.begin(),
                            solver_parameters.diffusion.fields.end(),
                            FieldType::V),
                  solver_parameters.diffusion.fields.end());
        EXPECT_NE(std::find(solver_parameters.diffusion.fields.begin(),
                            solver_parameters.diffusion.fields.end(),
                            FieldType::W),
                  solver_parameters.diffusion.fields.end());
        auto jacobi = std::get<Settings::solver::diffusion_solvers::jacobi>(solver_parameters.diffusion.solver.value());
        EXPECT_EQ(jacobi.max_iter, 77);
        EXPECT_DOUBLE_EQ(jacobi.tol_res, 1e-7);
        EXPECT_DOUBLE_EQ(jacobi.w, 1);

        EXPECT_TRUE(solver_parameters.solution.analytical_solution);
        EXPECT_DOUBLE_EQ(solver_parameters.solution.solution_tolerance.value(), 1e-4);
    }
    TEST(SettingsTest, diffusionSolverColoredGaussSeidel) {
        std::string xml = R"(
<ARTSS>
<solver description="DiffusionSolver" >
    <diffusion type="ColoredGaussSeidel" field="u,v,w" >
        <max_iter> 767 </max_iter>
        <tol_res> 1e-01 </tol_res>
        <w> 0.666 </w>
    </diffusion>
    <solution available="No"/>
</solver>
</ARTSS>)";
        tinyxml2::XMLDocument doc;
        doc.Parse(xml.c_str());
        Settings::solver_parameters solver_parameters = Settings::parse_solver_parameters(doc.RootElement());
        EXPECT_EQ(solver_parameters.diffusion.type, DiffusionMethods::ColoredGaussSeidel);

        EXPECT_NE(std::find(solver_parameters.diffusion.fields.begin(),
                            solver_parameters.diffusion.fields.end(),
                            FieldType::U),
                  solver_parameters.diffusion.fields.end());
        EXPECT_NE(std::find(solver_parameters.diffusion.fields.begin(),
                            solver_parameters.diffusion.fields.end(),
                            FieldType::V),
                  solver_parameters.diffusion.fields.end());
        EXPECT_NE(std::find(solver_parameters.diffusion.fields.begin(),
                            solver_parameters.diffusion.fields.end(),
                            FieldType::W),
                  solver_parameters.diffusion.fields.end());
        auto colored_gauss_seidel = std::get<Settings::solver::diffusion_solvers::colored_gauss_seidel>(solver_parameters.diffusion.solver.value());
        EXPECT_EQ(colored_gauss_seidel.max_iter, 767);
        EXPECT_DOUBLE_EQ(colored_gauss_seidel.tol_res, 1e-1);
        EXPECT_DOUBLE_EQ(colored_gauss_seidel.w, 0.666);
    }

    TEST(SettingsTest, diffusionSolverExplicit) {
        std::string xml = R"(
<ARTSS>
<solver description="DiffusionSolver" >
    <diffusion type="Explicit" field="u,v,w" >
    </diffusion>
    <solution available="No"/>
</solver>
</ARTSS>)";
        tinyxml2::XMLDocument doc;
        doc.Parse(xml.c_str());
        Settings::solver_parameters solver_parameters = Settings::parse_solver_parameters(doc.RootElement());
        EXPECT_EQ(solver_parameters.diffusion.type, DiffusionMethods::Explicit);

        EXPECT_NE(std::find(solver_parameters.diffusion.fields.begin(),
                            solver_parameters.diffusion.fields.end(),
                            FieldType::U),
                  solver_parameters.diffusion.fields.end());
        EXPECT_NE(std::find(solver_parameters.diffusion.fields.begin(),
                            solver_parameters.diffusion.fields.end(),
                            FieldType::V),
                  solver_parameters.diffusion.fields.end());
        EXPECT_NE(std::find(solver_parameters.diffusion.fields.begin(),
                            solver_parameters.diffusion.fields.end(),
                            FieldType::W),
                  solver_parameters.diffusion.fields.end());
    }
}

namespace turbulence_solvers {
    TEST(SettingsTest, turbulenceSolverConstSmagorinsky) {
        std::string xml = R"(
<ARTSS>
<solver description="DiffusionTurbSolver" >
    <diffusion type="Explicit" field="u,v,w" >
    </diffusion>
    <turbulence type="ConstSmagorinsky">
        <Cs> 0.2 </Cs>
    </turbulence>
    <solution available="No"/>
</solver>
</ARTSS>)";
        tinyxml2::XMLDocument doc;
        doc.Parse(xml.c_str());
        Settings::solver_parameters solver_parameters = Settings::parse_solver_parameters(doc.RootElement());

        EXPECT_EQ(solver_parameters.diffusion.type, DiffusionMethods::Explicit);

        EXPECT_NE(std::find(solver_parameters.diffusion.fields.begin(),
                            solver_parameters.diffusion.fields.end(),
                            FieldType::U),
                  solver_parameters.diffusion.fields.end());
        EXPECT_NE(std::find(solver_parameters.diffusion.fields.begin(),
                            solver_parameters.diffusion.fields.end(),
                            FieldType::V),
                  solver_parameters.diffusion.fields.end());
        EXPECT_NE(std::find(solver_parameters.diffusion.fields.begin(),
                            solver_parameters.diffusion.fields.end(),
                            FieldType::W),
                  solver_parameters.diffusion.fields.end());


        EXPECT_EQ(solver_parameters.turbulence.type, TurbulenceMethods::ConstSmagorinsky);
        EXPECT_DOUBLE_EQ(solver_parameters.turbulence.solver.value().cs, 0.2);
    }
    TEST(SettingsTest, turbulenceSolverDynamicSmagorinsky) {
        std::string xml = R"(
<ARTSS>
<solver description="DiffusionTurbSolver" >
    <diffusion type="Explicit" field="u,v,w" >
    </diffusion>
    <turbulence type="DynamicSmagorinsky">
    </turbulence>
    <solution available="No"/>
</solver>
</ARTSS>)";
        tinyxml2::XMLDocument doc;
        doc.Parse(xml.c_str());
        Settings::solver_parameters solver_parameters = Settings::parse_solver_parameters(doc.RootElement());

        EXPECT_EQ(solver_parameters.diffusion.type, DiffusionMethods::Explicit);

        EXPECT_NE(std::find(solver_parameters.diffusion.fields.begin(),
                            solver_parameters.diffusion.fields.end(),
                            FieldType::U),
                  solver_parameters.diffusion.fields.end());
        EXPECT_NE(std::find(solver_parameters.diffusion.fields.begin(),
                            solver_parameters.diffusion.fields.end(),
                            FieldType::V),
                  solver_parameters.diffusion.fields.end());
        EXPECT_NE(std::find(solver_parameters.diffusion.fields.begin(),
                            solver_parameters.diffusion.fields.end(),
                            FieldType::W),
                  solver_parameters.diffusion.fields.end());


        EXPECT_EQ(solver_parameters.turbulence.type, TurbulenceMethods::DynamicSmagorinsky);
    }
}

namespace pressure_solvers {
    TEST(SettingsTest, pressureSolverVCycleMG) {
        std::string xml = R"(
<ARTSS>
    <solver description="PressureSolver" >
        <pressure type="VCycleMG" field="p">
            <n_level> 5 </n_level>  <!-- number of restriction levels -->
            <n_cycle> 2 </n_cycle> <!-- number of cycles -->
            <max_cycle> 100 </max_cycle>  <!-- maximal number of cycles in first time step -->
            <tol_res> 1e-10 </tol_res>  <!-- tolerance for residuum/ convergence -->
            <n_relax> 4 </n_relax>  <!-- number of iterations -->
            <diffusion type="Jacobi" field="p">
                <max_iter> 100 </max_iter>  <!-- maximal number of iterations in solving at lowest level -->
                <tol_res> 1e-07 </tol_res>  <!-- tolerance for residuum/ convergence -->
                <w> 0.6666666667 </w>  <!-- relaxation parameter  -->
            </diffusion>
        </pressure>
        <solution available="Yes">
            <tol> 1e-03 </tol>  <!-- tolerance for further tests -->
        </solution>
    </solver>
</ARTSS>)";
        tinyxml2::XMLDocument doc;
        doc.Parse(xml.c_str());
        Settings::solver_parameters solver_parameters = Settings::parse_solver_parameters(doc.RootElement());

        EXPECT_EQ(solver_parameters.pressure.type, PressureMethods::VCycleMG);
        EXPECT_EQ(solver_parameters.pressure.field, FieldType::P);
        auto vcycle = solver_parameters.pressure.solver;
        EXPECT_EQ(vcycle.n_level, 5);
        EXPECT_EQ(vcycle.n_cycle, 2);
        EXPECT_EQ(vcycle.max_cycle, 100);
        EXPECT_EQ(vcycle.n_relax, 4);
        EXPECT_DOUBLE_EQ(vcycle.tol_res, 1e-10);

        EXPECT_EQ(vcycle.smoother.type, DiffusionMethods::Jacobi);

        EXPECT_NE(std::find(vcycle.smoother.fields.begin(),
                            vcycle.smoother.fields.end(),
                            FieldType::P),
                  vcycle.smoother.fields.end());
        auto jacobi = std::get<Settings::solver::diffusion_solvers::jacobi>(vcycle.smoother.solver.value());
        EXPECT_DOUBLE_EQ(jacobi.tol_res, 1e-7);
        EXPECT_EQ(jacobi.max_iter, 100);
        EXPECT_DOUBLE_EQ(jacobi.w, 0.6666666667);

        EXPECT_TRUE(solver_parameters.solution.analytical_solution);
        EXPECT_DOUBLE_EQ(solver_parameters.solution.solution_tolerance.value(), 1e-3);
    }
}