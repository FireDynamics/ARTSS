
#include <iostream>
#include <filesystem>

#include <gtest/gtest.h>

#include "src/Domain.h"
#include "src/field/Field.h"
#include "src/utility/Parameters.h"
#include "src/utility/GlobalMacrosTypes.h"


class FieldTest : public testing::Test {
 protected:
    void SetUp() {
        std::FILE* paramTmp = std::tmpfile();
        std::fputs("\
<?xml version=\"1.0\" encoding=\"UTF-8\" ?>\
<ARTSS>\
  <physical_parameters>\
    <t_end> 2.0 </t_end>  <!-- simulation end time -->\
    <dt> 0.001 </dt>  <!-- time stepping, caution: CFL-condition dt < 0.5*dx^2/nu -->\
  </physical_parameters>\
\
  <solver description=\"AdvectionSolver\" >\
    <advection type=\"SemiLagrangian\" field=\"u,v,w\">\
    </advection>\
    <solution available=\"Yes\">\
      <tol> 1e-03 </tol>  <!-- tolerance for further tests -->\
    </solution>\
  </solver>\
\
  <domain_parameters>\
    <X1> 0. </X1>  <!-- physical domain -->\
    <X2> 2.0 </X2>\
    <Y1> 0. </Y1>\
    <Y2> 2.0 </Y2>\
    <Z1> 0. </Z1>\
    <Z2> 1.0 </Z2>\
    <x1> 0. </x1>  <!-- computational domain -->\
    <x2> 2.0 </x2>\
    <y1> 0. </y1>\
    <y2> 2.0 </y2>\
    <z1> 0. </z1>\
    <z2> 1.0 </z2>\
    <nx> 40 </nx>  <!-- grid resolution (number of cells excl. ghost cells) -->\
    <ny> 40 </ny>\
    <nz> 1 </nz>\
  </domain_parameters>\
\
  <adaption dynamic=\"No\" data_extraction=\"No\"> </adaption>\
\
  <boundaries>\
    <boundary field=\"u,v,w\" patch=\"front,back,left,right,bottom,top\" type=\"dirichlet\" value=\"0.0\" />\
  </boundaries>\
\
  <obstacles enabled=\"No\"/>\
\
  <surfaces enabled=\"No\"/>\
\
  <initial_conditions usr_fct=\"GaussBubble\"  random=\"No\">  <!-- Gaussian function  -->\
    <u_lin> 0.5 </u_lin>      <!-- x-velocity in linear case  -->\
    <v_lin> 0.5 </v_lin>      <!-- y-velocity in linear case  -->\
    <w_lin> 0.25 </w_lin>     <!-- z-velocity in linear case  -->\
    <xshift> 1.025 </xshift>  <!-- xshift of Gauss Bubble in domain  -->\
    <yshift> 1.025 </yshift>  <!-- yshift of Gauss Bubble in domain  -->\
    <zshift> 0.5 </zshift>    <!-- zshift of Gauss Bubble in domain  -->\
    <l> 0.03125 </l>          <!-- sigma in Gaussian -->\
  </initial_conditions>\
\
  <visualisation save_vtk=\"Yes\" save_csv=\"No\">\
    <vtk_nth_plot> 200 </vtk_nth_plot>\
  </visualisation>\
\
  <logging file=\"output_test_advection.log\" level=\"info\">\
  </logging>\
</ARTSS>\
                """, paramTmp);
        std::rewind(paramTmp);

        auto params = Parameters::getInstance();
        params->parse(paramTmp, std::to_string(fileno(paramTmp)));
    }
};

TEST_F(FieldTest, swap_field) {
    Field a(UNKNOWN_FIELD, 0.0);
    Field b(UNKNOWN_FIELD, 0.0);

    real x = 0.0;
    for (auto i=0; i < Domain::getInstance()->get_size(); ++i) {
        a.data[i] = x + 0.0;
        b.data[i] = x + 0.5;
        x += 1.0;
    }

    Field::swap(&a, &b);

    x = 0.0;
    for (auto i=0; i < Domain::getInstance()->get_size(); ++i) {
        ASSERT_EQ(a.data[i], x + 0.5);
        ASSERT_EQ(b.data[i], x + 0.0);
        x += 1.0;
    }
}
