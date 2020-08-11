#!/usr/bin/env bash

#----------------------------Default values--------------------------------------
XML=Test_\%d.xml

DOMAINSURFACES=""
# leave empty if not necessary - do not forget to escape "
DOMAINOBSTACLES=""
# leave empty if not necessary - do not forget to escape "
BOUNDARYCONDITIONS="
    <boundaries>
       <boundary field=\"u,v,w\" patch=\"front,back,top,bottom\" type=\"dirichlet\" value=\"0.0\" />
       <boundary field=\"u,v,w\" patch=\"left,right\" type=\"neumann\" value=\"0.0\" />
       <boundary field=\"p\" patch=\"front,back,top,bottom\" type=\"neumann\" value=\"0.0\" />
       <boundary field=\"p\" patch=\"left,right\" type=\"dirichlet\" value=\"0.0\" />
       <boundary field=\"T\" patch=\"front,back,left,right,top\" type=\"dirichlet\" value=\"304.64\" />
       <boundary field=\"T\" patch=\"bottom\" type=\"neumann\" value=\"0.0\" />
    </boundaries>"
INITIALCONDITION="
    <initial_conditions usr_fct = \"LayersT\" dir=\"y\">  <!-- Layers  -->
	<n_layers> 5 </n_layers>
        <border_1> 1. </border_1>  <!-- at cell face -->
        <border_2> 2. </border_2>  <!-- at cell face -->
        <border_3> 3. </border_3>  <!-- at cell face -->
        <border_4> 4. </border_4>  <!-- at cell face -->
        <value_1> 303.64 </value_1>
        <value_2> 304.04 </value_2>
        <value_3> 305.24 </value_3>
        <value_4> 308.84 </value_4>
        <value_5> 310.54 </value_5>
    </initial_conditions>"

SOURCE="
"

TEND=0.1 #simulation end time
NU=0. #physical diffusion
BETA=3.4e-3 #buoancy force
G=-9.81 #gravitational constant
KAPPA=4.2e-5 #thermal diffusion

SOLVER=NSTempTurbSolver

ADVECTIONTYPE=SemiLagrangian
DIFFUSIONTYPE=Jacobi
MAXITER=100
TURBULENCETYPE=ConstSmagorinsky
CS=0.2
SOURCETYPE=ExplicitEuler
FORCEFCT=Buoyancy
FORCEDIR=y
PRESSURETYPE=VCycleMG
PRESSUREDIFFUSIONTYPE=Jacobi
TEMPADVECTIONTYPE=SemiLagrangian
TEMPDIFFUSIONTYPE=Jacobi
PR=0.2
PRT=0.9
TEMPSOURCETYPE=ExplicitEuler
TEMPSOURCEFCT=Zero
HRR=100.
CP=1.
GAUSSX0=0.
GAUSSY0=0.
GAUSSZ0=0.
SIGMAX=0.25
SIGMAY=0.6
SIGMAZ=0.25
TAU=5.
CONADVECTIONTYPE=SemiLagrangian
CONDIFFUSIONTYPE=Jacobi
SCT=0.5
CONSOURCETYPE=ExplicitEuler
CONFORCEFCT=Zero
HC=13100.
YS=0.01
NLEVEL=4
NCYCLE=2
MAXCYCLE=100
MAXSOLVE=100
CALCNLEVEL=1 #use function calc_nlevel

#solution and visualisation
SOLAVAIL=No # solution available
CSV=No
CSVPLOTS=10 #number of csv plots
VTK=Yes
VTKPLOTS=10 #number of vtk plots

#domain parameters
XSTART=-0.1556 #physical domain parameter X1
XEND=0.1556 #physical domain parameter X2
YSTART=-0.1556 #physical domain parameter Y1
YEND=0.1556 #physical domain parameter Y2
ZSTART=-0.1556 #physical domain parameter Z1
ZEND=0.1556 #physical domain parameter Z2

xSTART=-0.1556 #computational domain parameter x1
xEND=0.1556 #computational domain parameter x2
ySTART=-0.1556 #computational domain parameter y1
yEND=0.1556 #computational domain parameter y2
zSTART=-0.1556 #computational domain parameter z1
zEND=0.1556 #computational domain parameter z2

NX=6,10,18,34,66 #Nx
NY=6,10,18,34,66 #Ny
NZ=3,3,3,3,3 #Nz
DT=0.1

#dynamic adaption
ADAPTION=No
ADAPT_PAR=1
ADAPT_CLASS=Layers
ADAPT_DATA=No
ADAPT_BEFORE=No
ADAPT_AFTER=No
ADAPT_TIME=No
ADAPT_END=No

OBSTACLE=1
SURFACE=1

LEADINGNO=2

LOGGINGFILE=tmp
LOGGINGLEVEL=INFO
#------------------------------------------------------------------------------
YELLOW='\033[1;33m'
NC='\033[0;m'

INIC=1
BOUC=1
SOUC=1
TEMPSOUC=1
DOBST=1
DSURF=1

POSITIONAL=()
#----Help text----
DESCRIPTION="Description:
Script to build a XML file for the different cases of ARTSS. Built in: Advection, Burgers (Advection Diffusion), Diffusion, Diffusion Turbulence, Navier Stokes, Navier Stokes Temperature, Navier Stokes Temperatur Turbulence, Navier Stokes Turbulence, Navier Stokes Temperature Turbulence Concentration, Pressure.

There is no validation of the given parameter.\n"
OPTIONSTEXT="Available options:

Basic options:"
OPTIONSA="
${YELLOW}-o${NC}\tname of xml file
${YELLOW}-p${NC}\tnumber of vtk plots
${YELLOW}-t${NC}\tdt parameter, Multiple parameter possible. Delimeter: ','
${YELLOW}-v${NC}\tverbose - print output to screen
${YELLOW}-x${NC}\tdomain parameter Nx values. delimeter: ','
${YELLOW}-y${NC}\tdomain parameter Ny values. delimeter: ','
${YELLOW}-z${NC}\tdomain parameter Nz values. delimeter: ','

Advanced Options:"

OPTIONSB="
${YELLOW}--adaption${NC}\tenable dynamic domain adaption

${YELLOW}--adaptionafter${NC}\theight for after option of domain adaption (default: $ADAPT_AFTER)

${YELLOW}--adaptionbefore${NC}\theight for before option of domain adaption (default: $ADAPT_BEFORE)

${YELLOW}--adaptionendresult${NC}\tenable endresult of domain adaption

${YELLOW}--adaptionfield${NC}\tenable field option of domain adaption

${YELLOW}--adaptionparameter${NC}\tget adaption parameters from file

${YELLOW}--adaptiontime${NC}\tenable time measuring for domain adaption

${YELLOW}--adv${NC}
${YELLOW}--advection${NC}              \tcreate xml file with parameter for advection only

${YELLOW}--advectiontype${NC} \tset advection discretization method (default: $ADVECTIONTYPE)

${YELLOW}--beta${NC}          \tset physical parameter thernal expansion coefficient for buoancy force (default: $BETA)

${YELLOW}--boundaryconditions${NC} \ttake boundary conditions from file

${YELLOW}--bur${NC}
${YELLOW}--burgers${NC}                \tcreate xml file with parameter for Burgers (Advection Diffusion)

${YELLOW}--cn${NC}
${YELLOW}--calcnlevel${NC}\tcalculate nlevel

${YELLOW}--conadvtype${NC}
${YELLOW}--conadvectiontype${NC}       \tset concentration advection discretization method (default: $CONADVECTIONTYPE)

${YELLOW}--condifftype${NC}
${YELLOW}--condiffusiontype${NC}       \tset concentration diffusion discretization method (default: $CONDIFFUSIONTYPE)

${YELLOW}--consourcetype${NC} \tset concentration source discretization method (default: $CONSOURCETYPE)

${YELLOW}--consourcefct${NC} \tset concentration source function (default: $CONFORCEFCT)

${YELLOW}--cp${NC}            \tset specific heat capacity (in kJ/kgK) in case of volumetric heat source (default: $CP)

${YELLOW}--cs${NC}            \tset parameter for ConstSmagorinsky (default: $CS)

${YELLOW}--csv${NC}           \tenable write out of csv files (default: $CSV)

${YELLOW}--csvplots${NC}      \twrite out a csv file every nth timestep (default: $CSVPLOTS)

${YELLOW}--dataextraction${NC}\t enable data extraction

${YELLOW}--diff${NC}
${YELLOW}--diffusion${NC}              \tcreate xml file with parameter for diffusion only

${YELLOW}--difft${NC}
${YELLOW}--diffturb${NC}               \tcreate xml file with parameter for difussion turbluence

${YELLOW}--diffusiontype${NC} \tset diffusion discretization method (default: $DIFFUSIONTYPE)

${YELLOW}--domainobstacles${NC}\ttake domain obstacles from file

${YELLOW}--domainsurfaces${NC}\ttake domain surfaces from file

${YELLOW}--dt${NC}            \tset dt, multiple parameter possible. Delimeter: ',' (default: $DT)

${YELLOW}--forcefct${NC}      \tset force function of source (default: $FORCEFCT)

${YELLOW}--forcedir${NC}      \tset direction of force (default: $FORCEDIR)

${YELLOW}--g${NC}             \tset physical parameter gravitational constant (default: $G)

${YELLOW}--gaussx0${NC}             \tset x-center of Gaussian in case of volumetric heat source (default: $GAUSSX0)
${YELLOW}--gaussy0${NC}             \tset y-center of Gaussian in case of volumetric heat source (default: $GAUSSY0)
${YELLOW}--gaussz0${NC}             \tset z-center of Gaussian in case of volumetric heat source (default: $GAUSSZ0)

${YELLOW}--help${NC}          \tdisplay help

${YELLOW}--hrr${NC}          \tset total heat release rate (less radiative fraction) in kW in case of volumetric heat source (default: $HRR)

${YELLOW}--hc${NC}          \tset energy release per unit mass (in kJ/kg) (default: $HC)

${YELLOW}--initialconditions${NC} \tgive initial conditions in file

${YELLOW}--ys${NC}          \tset soot yield (in %/100)

${YELLOW}--kappa${NC}         \tset physical parameter thermal diffusion (default: $KAPPA)

${YELLOW}--logginglevel${NC}        \tset the logging level \{TRACE, DEBUG, INFO, WARN, ERROR\}
${YELLOW}--loggingfile${NC}         \tset the output file

${YELLOW}--maxcycle${NC}\tset number of maximal cycles for V-Cycle multigrid pre-conditioning (default: $MAXCYCLE)

${YELLOW}--maxiter${NC}\tset number of iteration for diffusion, multiple parameter possible. Delimiter: ',' (default: $MAXITER)

${YELLOW}--maxsolve${NC}\tset number of maximal iterations in solver in multigrid (default: $MAXSOLVE)

${YELLOW}--ncycle${NC}        \tset ncycle for multigrid (default: $NCYCLE)

${YELLOW}--nlevel${NC}        \tset nlevel for multigrid, set dependent on Nx,Ny,Nz (default: $NLEVEL)

${YELLOW}--ns${NC}
${YELLOW}--navierstokes${NC}           \tcreate xml file with parameter for navier stokes

${YELLOW}--nste${NC}
${YELLOW}--nstemp${NC}                 \tcreate xml file with parameter for navier stokes temperature

${YELLOW}--nstt${NC}
${YELLOW}--nstempturb${NC}          \tcreate xml file with parameter for navier stokes temperature with turbulence model

${YELLOW}--nstu${NC}
${YELLOW}--nsturb${NC}                 \tcreate xml file with parameter for navier stokes with turbulence model

${YELLOW}--nstc${NC}
${YELLOW}--nstempcon${NC}              \tcreate xml file with parameter for navier stokes temperature concentration

${YELLOW}--nsttc${NC}
${YELLOW}--nstempturbcon${NC}          \tcreate xml file with parameter for navier stokes temperature concentration with turbulence model

${YELLOW}--nu${NC}            \tset physical parameter nu (default: $NU)

${YELLOW}--nx${NC}            \tset domain parameter nx, multiple parameter possible. Delimiter: ',' (default: $NX)

${YELLOW}--ny${NC}            \tset domain parameter ny, multiple parameter possible. Delimiter: ',' (default: $NY)

${YELLOW}--nz${NC}            \tset domain parameter nz, multiple parameter possible. Delimiter: ',' (default: $NZ)

${YELLOW}--out
${YELLOW}--xml${NC}                    \tset name of xml file. Use %d for multiple files (default: $XML)

${YELLOW}--pre${NC}
${YELLOW}--pressure${NC}               \tcreate XML file for pressure only

${YELLOW}--pressuretype${NC}  \tset pressure discretization method (default: $PRESSURETYPE)

${YELLOW}--predifftype
${YELLOW}--pressurediffusiontype${NC}  \tset pressure diffusion discretization method (default: $PRESSUREDIFFUSIONTYPE)

${YELLOW}--prt${NC}           \tset turbulent Prandtl number Pr_T (default: $PRT)

${YELLOW}--sigmax${NC}             \tset x-spread (FWHM) of Gaussian in case of volumetric heat source (default: $SIGMAX)
${YELLOW}--sigmay${NC}             \tset y-spread (FWHM) of Gaussian in case of volumetric heat source (default: $SIGMAY)
${YELLOW}--sigmaz${NC}             \tset z-spread (FWHM) of Gaussian in case of volumetric heat source (default: $SIGMAZ)

${YELLOW}--solavail${NC}        \tset availability of an analytical solution (in Yes/No) (default: $SOLAVAIL)

${YELLOW}--solver${NC}        \tset solver description (default: $SOLVER)

${YELLOW}--sct${NC}           \tset turbulent Schmidt number Sc_T (default: $SCT)

${YELLOW}--sourcecondition${NC}\ttake source condition from file

${YELLOW}--sourcetype${NC}    \tset source discretization method (default: $SOURCETYPE)

${YELLOW}--tau${NC}             \tset ramp-up time in case of volumetric heat source (default: $TAU)

${YELLOW}--tempadvtype${NC}   \tset temperature advection discretization method (default: $TEMPADVECTIONTYPE)

${YELLOW}--tempdifftype${NC}  \tset temperature diffusion discretization method (default: $TEMPDIFFUSIONTYPE)

${YELLOW}--tempsourcefunction${NC}\ttake temperature source function from file

${YELLOW}--tempsourcetype${NC}\tset temperature source discretization method (default: $TEMPSOURCETYPE)

${YELLOW}--tempsourcefct${NC}\tset temperature source function (default: $TEMPFORCEFCT)

${YELLOW}--tend${NC}          \tset t_end, multiple parameter possible. Delimiter: ',' (default: ${TEND})

${YELLOW}--turbulencetype${NC}\tset turbulence model (default: $TURBULENCETYPE)

${YELLOW}--verbose${NC}       \tprint output to screen (std out)

${YELLOW}--vtk${NC}           \tenable write out of vtk files (default: $VTK)

${YELLOW}--vtkplots${NC}      \twrite out a vtk file every nth timestep (default: $VTKPLOTS)

${YELLOW}--xstart${NC}        \tset both domain parameter (computational [x1] and physical domain [X1])

${YELLOW}--xstartp${NC}       \tset domain parameter of physical domain X1 (default: $XSTART)

${YELLOW}--xstartc${NC}       \tset domain parameter of computational domain x1 (default: $xSTART)

${YELLOW}--xend${NC}        \tset both domain parameter (computational [x2] and physical domain [X2])

${YELLOW}--xendp${NC}         \tset domain parameter of physical domain X2 (default: $XEND)

${YELLOW}--xendc${NC}         \tset domain parameter of computational domain x2 (default: $xEND)

${YELLOW}--ystart${NC}        \tset both domain parameter (computational [y1] and physical domain [Y1])

${YELLOW}--ystartp${NC}       \tset domain parameter of physical domain Y1 (default: $YSTART)

${YELLOW}--ystartc${NC}       \tset domain parameter of computational domain y1 (default: $ySTART)

${YELLOW}--yend${NC}          \tset both domain parameter (computational [y2] and physical domain [Y2])

${YELLOW}--yendp${NC}         \tset domain parameter of physical domain Y2 (default: $YEND)

${YELLOW}--yendc${NC}         \tset domain parameter of computational domain y2 (default: $yEND)

${YELLOW}--zstart${NC}        \tset both domain parameter (computational [z1] and physical domain [Z1])

${YELLOW}--zstartp${NC}       \tset domain parameter of physical domain Z1 (default: $ZSTART)

${YELLOW}--zstartc${NC}       \tset domain parameter of computational domain z1 (default: $zSTART)

${YELLOW}--zend${NC}        \tset both domain parameter (computational [z2] and physical domain [Z2])

${YELLOW}--zendp${NC}         \tset domain parameter of physical domain Z2 (default: $ZEND)

${YELLOW}--zendc${NC}         \tset domain parameter of computational domain z2 (default: $zEND)"
#----help text end----

#parse help text
COUNTER=0
TMP=0
while IFS= read -r LINE
do
  if [ ! -z "$LINE" ]
  then
    tmpa=$(echo -e "$LINE" | cut -f1 | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//' )
    options[$COUNTER]="$tmpa"
    tmpb="$(echo -e "$LINE" | cut -f2 )"
    if [[ "$tmpa" != "$tmpb" ]]
    then
      text[$COUNTER]="$tmpb"
    fi
    ((COUNTER++))
  fi
done < <(printf '%s\n' "$OPTIONSB")

for ((i=0;i<$COUNTER;i++))
do
  OPTIONS=$OPTIONS$'\n'$(printf "%-37s%s\n" ${options[$i]} "${text[$i]}")
done

HELP="${DESCRIPTION}${OPTIONSTEXT}${OPTIONSA}${OPTIONS}"

#different test cases
ADV=1 #advection
BUR=1 #burgers
DIFF=1 #diffusion
DIFFT=1 #diffusion turb
NS=1 #navier stokes
NSTu=1 #navier stokes turb
NSTe=1 #navier stokes temp
NSTT=1 #navier stokes temp turb
NSTC=1 #navier stokes temp con
NSTTC=1 #navier stokes temp turb con
PRE=1 #pressure

#--------------------------------functions---------------------------------------
#parse default values
function parse_default {
  DTVALUES=( $(echo "$DT" | tr "," "\n") )
  TVALUES=( $(echo "$TEND" | tr "," "\n") )
  NXVALUES=( $(echo "$NX" | tr "," "\n") )
  NYVALUES=( $(echo "$NY" | tr "," "\n") )
  NZVALUES=( $(echo "$NZ" | tr "," "\n") )
  MAXITERVALUES=$MAXITER
}
#calculate nlevel
function calc_nlevel {
  NXm2=$(bc <<< "scale=0;$NX-2")
  NYm2=$(bc <<< "scale=0;$NY-2")
  NZm2=$(bc <<< "scale=0;$NZ-2")
  LOG2NXF=$(echo "l($NXm2)/l(2)" | bc -l)
  LOG2NX=$(echo "$LOG2NXF/1" | bc)
  LOG2NYF=$(echo "l($NYm2)/l(2)" | bc -l)
  LOG2NY=$(echo "$LOG2NYF/1" | bc)
  LOG2NZF=$(echo "l($NZm2)/l(2)" | bc -l)
  LOG2NZ=$(echo "$LOG2NZF/1" | bc)
  LOG2NX1=$(bc <<< "scale=0;$LOG2NX-1")
  LOG2NY1=$(bc <<< "scale=0;$LOG2NY-1")
  LOG2NZ1=$(bc <<< "scale=0;$LOG2NZ-1")
  MAXARR=( $LOG2NX1 $LOG2NY1 $LOG2NZ1 )
  max=${MAXARR[0]}
  for n in "${MAXARR[@]}" ; do
    ((n > max)) && max=$n
  done
  NLEVEL=$max
}

#write xml file
function writeXML {
if [ -e $NAME ]
then
  echo "$NAME already exist. Skipping file"
else
  WRITETO="<?xml version=\"1.0\" encoding=\"UTF-8\" ?>
<ARTSS>
  <physical_parameters>
    <t_end> $TEND </t_end>  <!-- simulation end time -->
    <dt> $DT </dt>  <!-- time stepping, caution: CFL-condition dt < 0.5*dx^2/nu -->"
  if [ $BUR -eq 0 -o $DIFF -eq 0 -o $DIFFT -eq 0 -o $NS -eq 0 -o $NSTe -eq 0 -o $NSTC -eq 0 -o $NSTT -eq 0 -o $NSTTC -eq 0 -o $NSTu -eq 0 ]
  then
    WRITETO="$WRITETO
    <nu> $NU </nu>  <!-- kinematic viscosity -->"
  fi
  if [ $NSTe -eq 0 -o $NSTC -eq 0 -o $NSTT -eq 0 -o $NSTTC -eq 0 ]
  then
    if [ \"$FORCEFCT\" == \"Buoyancy\" ]
    then
    WRITETO="$WRITETO
    <beta> $BETA </beta>  <!-- thermal expansion coefficient -->
    <g> $G </g>  <!-- gravitational constant -->"
    fi
  fi
  if [ $NSTe -eq 0 -o $NSTC -eq 0 -o $NSTT -eq 0 -o $NSTTC -eq 0 ]
  then
    WRITETO="$WRITETO
    <kappa> $KAPPA </kappa>  <!-- thermal diffusion -->"
  fi
  WRITETO="$WRITETO
  </physical_parameters>
"
  if [ $ADV -eq 0 ]
  then
    WRITETO="$WRITETO
  <solver description = \"AdvectionSolver\" >"
  elif [ $BUR -eq 0 ]
  then
    WRITETO="$WRITETO
  <solver description = \"AdvectionDiffusionSolver\" >"
  elif [ $DIFF -eq 0 ]
  then
    WRITETO="$WRITETO
  <solver description = \"DiffusionSolver\" >"
  elif [ $DIFFT -eq 0 ]
  then
    WRITETO="$WRITETO
  <solver description = \"DiffusionTurbSolver\" >"
  elif [ $NS -eq 0 ]
  then
    WRITETO="$WRITETO
  <solver description = \"NSSolver\" >"
  elif [ $NSTu -eq 0 ]
  then
    WRITETO="$WRITETO
  <solver description = \"NSTurbSolver\" >"
  elif [ $NSTe -eq 0 ]
  then
    WRITETO="$WRITETO
  <solver description = \"NSTempSolver\" >"
  elif [ $NSTT -eq 0 ]
  then
    WRITETO="$WRITETO
  <solver description = \"NSTempTurbSolver\" >"
  elif [ $NSTC -eq 0 ]
  then
    WRITETO="$WRITETO
  <solver description = \"NSTempConSolver\" >"
  elif [ $NSTTC -eq 0 ]
  then
    WRITETO="$WRITETO
  <solver description = \"NSTempTurbConSolver\" >"
  elif [ $PRE -eq 0 ]
  then
    WRITETO="$WRITETO
  <solver description = \"PressureSolver\" >"
  else
    WRITETO="$WRITETO
  <solver description = \"$SOLVER\" >"
  fi
#-------advection------
  if [ $ADV -eq 0 -o $BUR -eq 0 -o $NS -eq 0 -o $NSTu -eq 0 -o $NSTe -eq 0 -o $NSTT -eq 0 -o $NSTC -eq 0 -o $NSTTC -eq 0 ]
  then
    WRITETO="$WRITETO
    <advection type = \"$ADVECTIONTYPE\" field = \"u,v,w\">
    </advection>"
  fi
#-------diffusion-------
  if [ $DIFF -eq 0 -o $DIFFT -eq 0 -o $BUR -eq 0 -o $NS -eq 0 -o $NSTu -eq 0 -o $NSTe -eq 0 -o $NSTT -eq 0 -o $NSTC -eq 0 -o $NSTTC -eq 0 ]
  then
    WRITETO="$WRITETO
    <diffusion type = \"$DIFFUSIONTYPE\" field = \"u,v,w\">
      <max_iter> $MAXITER </max_iter>  <!-- max number of iterations -->
      <tol_res> 1e-07 </tol_res>  <!-- tolerance for residuum/ convergence -->
      <w> 1 </w>  <!-- relaxation parameter -->
    </diffusion>"
  fi
#-------turbulence-------
  if [ $DIFFT -eq 0 -o $NSTu -eq 0 -o $NSTT -eq 0 -o $NSTTC -eq 0 ]
  then
    WRITETO="$WRITETO
    <turbulence type = \"$TURBULENCETYPE\">
      <Cs> $CS </Cs>
    </turbulence>"
  fi
#--------source-------
  if [ $NS -eq 0 -o $NSTu -eq 0 -o $NSTe -eq 0 -o $NSTT -eq 0 -o $NSTC -eq 0 -o $NSTTC -eq 0 ]
  then
    if [ $SOUC -eq 0 ]
    then
      WRITETO=${WRITETO}"\n"$(cat $SFILE)
    else
      WRITETO="$WRITETO
    <source type = \"$SOURCETYPE\" force_fct = \"$FORCEFCT\" dir = \"$FORCEDIR\">  <!-- Direction of force (x,y,z or combinations xy,xz,yz,xyz) -->
    </source>"
    fi
  fi
#------pressure-------
  if [ $PRE -eq 0 -o $NS -eq 0 -o $NSTu -eq 0 -o $NSTe -eq 0 -o $NSTT -eq 0 -o $NSTC -eq 0 -o $NSTTC -eq 0 ]
  then
    WRITETO="$WRITETO
    <pressure type = \"$PRESSURETYPE\" field = \"p\">
      <n_level> $NLEVEL </n_level>  <!-- number of restriction levels -->
      <n_cycle> $NCYCLE </n_cycle> <!-- number of cycles -->
      <max_cycle> $MAXCYCLE </max_cycle>  <!-- maximal number of cycles in first time step -->
      <tol_res> 1e-07 </tol_res>  <!-- tolerance for residuum/ convergence -->
      <diffusion type = \"$PRESSUREDIFFUSIONTYPE\" field = \"p\">
        <n_relax> 4 </n_relax>  <!-- number of iterations -->
        <max_solve> $MAXSOLVE </max_solve>  <!-- maximal number of iterations in solving at lowest level -->
        <tol_res> 1e-07 </tol_res>  <!-- tolerance for residuum/ convergence -->
        <w> 0.6666666667 </w>  <!-- relaxation parameter  -->
      </diffusion>
    </pressure>"
  fi
#------temperature---------
  if [ $NSTe -eq 0 -o $NSTT -eq 0 -o $NSTC -eq 0 -o $NSTTC -eq 0 ]
  then
    WRITETO="$WRITETO
    <temperature>
      <advection type = \"$TEMPADVECTIONTYPE\" field = \"T\">
      </advection>
      <diffusion type = \"$TEMPDIFFUSIONTYPE\" field = \"T\">
        <max_iter> $MAXITER </max_iter>
        <tol_res> 1e-07 </tol_res>
        <w> 1 </w>
      </diffusion>"
  fi
  if [ $NSTT -eq 0 -o $NSTTC -eq 0 ]
  then
    WRITETO="$WRITETO
      <turbulence include = \"Yes\">
        <Pr_T> $PRT </Pr_T>
      </turbulence>"
  fi
  if [  $NSTe -eq 0 -o $NSTT -eq 0 -o $NSTC -eq 0 -o $NSTTC -eq 0 ]
  then
    if [ $TEMPSOUC -eq 0 ]
    then
      WRITETO="$WRITETO\n$(cat $TSFILE)"
    else
      if [ "$TEMPSOURCEFCT" == "Zero" ]
      then
        WRITETO="$WRITETO
      <source type = \"$TEMPSOURCETYPE\" temp_fct = \"$TEMPSOURCEFCT\" dissipation = \"No\">
      </source>"
      fi
      if [ \"$TEMPSOURCEFCT\" == \"GaussST\" ]
      then
        WRITETO="$WRITETO
      <source type = \"$TEMPSOURCETYPE\" temp_fct = \"GaussST\" ramp_fct = \"RampTanh\" dissipation = \"No\">
        <HRR> $HRR </HRR>  <!-- total heat release rate (in kW) -->
        <cp> $CP </cp>  <!-- specific heat capacity (in kJ/kgK)-->
        <x0> $GAUSSX0  </x0>
        <y0> $GAUSSY0 </y0>
        <z0> $GAUSSZ0 </z0>
        <sigmax> $SIGMAX </sigmax>
        <sigmay> $SIGMAY </sigmay>
        <sigmaz> $SIGMAZ </sigmaz>
        <tau> $TAU </tau>
      </source>"
      fi
    fi
    WRITETO="$WRITETO
    </temperature>"
  fi
#-------concentration---------
  if [ $NSTTC -eq 0 -o $NSTC -eq 0 ]
  then
    WRITETO="$WRITETO
    <concentration>
      <advection type = \"$CONADVECTIONTYPE\" field = \"rho\">
      </advection>
      <diffusion type = \"$CONDIFFUSIONTYPE\" field = \"rho\">
        <gamma> $GAMMA </gamma>  <!-- concentration diffussion -->
        <max_iter> $MAXITER </max_iter>  <!-- max number of iterations -->
        <tol_res> 1e-07 </tol_res>  <!-- tolerance for residuum/ convergence -->
        <w> 1 </w>  <!-- relaxation parameter -->
      </diffusion>"
    if [ $NSTTC -eq 0 ]
    then
      WRITETO="$WRITETO
        <turbulence include = \"Yes\">
          <Sc_T> $SCT </Sc_T>
        </turbulence>"
    fi
    if [ $CONFORCEFCT == "Zero" ]
    then
      WRITETO="$WRITETO
        <source type = \"$CONSOURCETYPE\" con_fct = \"Zero\">
        </source>"
    elif [ $CONFORCEFCT == "GaussSC" ]
    then
      WRITETO="$WRITETO
        <source type = \"$CONSOURCETYPE\" con_fct = \"GaussSC\" ramp_fct = \"RampTanh\">
          <HRR> $HRR </HRR>  <!-- total heat release rate (in kW) -->
          <Hc> $HC </Hc>  <!-- Heating value (in kJ/kg) -->
          <Ys> $YS </Ys>  <!-- Soot yield -->
          <x0> $GAUSSX0 </x0>
          <y0> $GAUSSY0 </y0>
          <z0> $GAUSSZ0 </z0>
          <sigmax> $SIGMAX </sigmax>
          <sigmay> $SIGMAY </sigmay>
          <sigmaz> $SIGMAZ </sigmaz>
          <tau> $TAU </tau>
        </source>"
    fi
  WRITETO=${WRITETO}"
    </concentration>"

  fi
  WRITETO="$WRITETO
    <solution available = \"$SOLAVAIL\">"
  if [ "$SOLAVAIL" == "Yes" ]
  then
    WRITETO="$WRITETO
      <tol> 1e-03 </tol>  <!-- tolerance for further tests -->"
  fi
  WRITETO="$WRITETO
    </solution>
  </solver>"
  WRITETO="$WRITETO

  <domain_parameters>
    <X1> $XSTART </X1>  <!-- physical domain -->
    <X2> $XEND </X2>
    <Y1> $YSTART </Y1>
    <Y2> $YEND </Y2>
    <Z1> $ZSTART </Z1>
    <Z2> $ZEND </Z2>
    <x1> $xSTART </x1>  <!-- computational domain -->
    <x2> $xEND </x2>
    <y1> $ySTART </y1>
    <y2> $yEND </y2>
    <z1> $zSTART </z1>
    <z2> $zEND </z2>
    <nx> $NX </nx>  <!-- grid resolution (number of cells excl. ghost cells) -->
    <ny> $NY </ny>
    <nz> $NZ </nz>
  </domain_parameters>
"

    if [ $ADAPT_PAR -eq 0 ]
    then
      ATEMP="\n$(cat $AFILE)\n"
    else
      WRITETO="${WRITETO}
  <adaption dynamic="
      if [[ $ADAPTION == "No" && $ADAPT_DATA == "No" ]]
      then
        WRITETO="$WRITETO\"No\" data_extraction=\"No\"> </adaption>"
      else
        unset ATEMP
        if [[ $ADAPTION == "Yes" ]]
        then
          WRITETO="$WRITETO\"Yes\""
        else
          WRITETO="$WRITETO\"No\""
        fi
        if [[ $ADAPT_DATA == "Yes" ]]
        then
          WRITETO="$WRITETO data_extraction=\"Yes\">
        $ATEMP
    <data_extraction>"
          if [[ $ADAPT_BEFORE == "No" ]]
          then
            WRITETO="$WRITETO
      <before enabled=\"No\"> </before>"
          else
            WRITETO="$WRITETO
      <before enabled=\"Yes\">$ADAPT_BEFORE</before>"
          fi
          if [[ $ADAPT_AFTER == "No" ]]
          then
            WRITETO="$WRITETO
      <after enabled=\"No\"> </after>"
          else
            WRITETO="$WRITETO
      <after enabled=\"Yes\">$ADAPT_AFTER</after>"
          fi
          if [[ $ADAPT_END == "No" ]]
          then
            WRITETO="$WRITETO
      <endresult enabled=\"No\"> </endresult>"
          else
            WRITETO="$WRITETO
      <endresult enabled=\"Yes\"> </endresult>"
          fi
          if [[ $ADAPT_TIME == "No" ]]
          then
            WRITETO="$WRITETO
      <time_measuring enabled=\"No\"> </time_measuring>"
          else
            WRITETO="$WRITETO
      <time_measuring enabled=\"Yes\"> </time_measuring>"
          fi
          WRITETO="$WRITETO
    </data_extraction>
  </adaption>"
        else
          WRITETO="$WRITETO data_extraction=\"No\">$ATEMP </adaption>"
        fi
      fi
    fi
    if [ $BOUC -eq 0 ]
    then
      WRITETO=${WRITETO}"\n\n"$(cat $BFILE)
    else
      WRITETO=${WRITETO}${BOUNDARYCONDITIONS}
    fi
    if [ $DOBST -eq 0 ]
    then
      WRITETO="${WRITETO}\n$(cat $DOBSTFILE)"
    else 
      if [ $OBSTACLE -eq 1 ]
      then
        WRITETO="${WRITETO}\n
  <obstacles enabled=\"No\"/>"
      else
        WRITETO="${WRITETO}\n
  <obstacles enabled=\"Yes\">"
        WRITETO="${WRITETO}${DOMAINOBSTACLES}
  </obstacles>"
      fi
    fi
    if [ $SURFACE -eq 1 ]
    then
      WRITETO="$WRITETO\n
  <surfaces enabled=\"No\"/>"
    else
      WRITETO="$WRITETO\n
  <surfaces enabled=\"Yes\">"
      if [ $DSURF -eq 0 ]
      then
        WRITETO="${WRITETO}\n$(cat $DSURFFILE)"
      else
        WRITETO="${WRITETO}${DOMAINSURFACES}"
      fi
      WRITETO="${WRITETO}
  </surfaces>"
    fi
    if [ $INIC -eq 0 ]
    then
      WRITETO=${WRITETO}"\n\n"$(cat $IFILE)
    else
      WRITETO=${WRITETO}${INITIALCONDITION}
    fi
    
    WRITETO="${WRITETO}\n
  <visualisation save_vtk=\"$VTK\" save_csv=\"$CSV\">"
    if [ $CSV == "Yes" ]
    then
      WRITETO=${WRITETO}"
    <csv_nth_plot> $CSVPLOTS </csv_nth_plot>"
    fi
    if [ $VTK == "Yes" ]
    then
      WRITETO=${WRITETO}"
    <vtk_nth_plot> $VTKPLOTS </vtk_nth_plot>"
    fi
    WRITETO=${WRITETO}"
  </visualisation>"
    WRITETO=${WRITETO}"
  <logging>
    <level>${LOGGINGLEVEL}</level>
    <file>${LOGGINGFILE}</file>
  </logging>
</ARTSS>"
  echo -e "$WRITETO" >> $NAME
  fi
}

#---------------------------------PROGRAMM--------------------------
parse_default
#parse command line arguments
while [[ $# -gt 0 ]]
do
  key="$1"
  case $key in
    --adaption)
      ADAPTION="Yes"
      shift
      ;;
    --adaptionafter)
      ADAPT_AFTER=$2
      shift
      shift
      ;;
    --adaptionbefore)
      ADAPT_BEFORE=$2
      shift
      shift
      ;;
    --adaptionendresult)
      ADAPT_END=0
      shift
      ;;
    --adaptionparameter)
      ADAPT_PAR=0
      AFILE=$2
      shift
      shift
      ;;
    --adaptiontime)
      ADAPT_TIME=0
      shift
      ;;
    --adv|--advection)
      ADV=0
      SOLVER=AdvectionSolver
      shift
      ;;
    --advectiontype)
      ADVECTIONTYPE=$2
      shift
      shift
      ;;
    --beta)
      BETA="$2"
      shift
      shift
      ;;
    --boundaryconditions)
      BOUC=0
      BFILE=$2
      shift
      shift
      ;;
    --bur|--burgers)
      BUR=0
      SOLVER=AdvectionDiffusionSolver
      shift
      ;;
    --cn|--calcnlevel)
      CALCNLEVEL=0
      shift
      ;;
    --conadvtype|--conadvectiontype)
      CONADVECTIONTYPE=$2
      shift
      shift
      ;;
    --condifftype|--condiffusiontype)
      CONDIFFUSIONTYPE=$2
      shift
      shift
      ;;
    --consourcetype)
      CONSOURCETYPE=$2
      shift
      shift
      ;;
    --consourcefct)
      CONSOURCEFCT=$2
      shift
      shift
      ;;
    --cp)
      CP=$2
      shift
      shift
      ;;
    --cs)
      CS=$2
      shift
      shift
      ;;
    --csv)
      CSV=Yes
      shift
      ;;
    --csvplots)
      CSVPLOTS=$2
      shift
      shift
      ;;
    --dataextraction)
      ADAPT_DATA="Yes"
      shift
      ;;
    --diff|--diffusion)
      DIFF=0
      SOLVER=DiffusionSolver
      shift
      ;;
    --difft|--diffturb|--diffusionturbulence)
      DIFFT=0
      SOLVER=DiffusionTurbSolver
      shift
      ;;
    --diffusiontype)
      DIFFUSIONTYPE=$2
      shift
      shift
      ;;
    --domainobstacles)
      OBSTACLE=0
      DOBST=0
      DOBSTFILE=$2
      shift
      shift
      ;;
    --domainsurfaces)
      SURFACE=0
      DSURF=0
      DSURFFILE=$2
      shift
      shift
      ;;
    --forcefct)
      FORCEFCT=$2
      shift
      shift
      ;;
    --forcedir)
      FORCEDIR=$2
      shift
      shift
      ;;
    --g)
      G=$2
      shift
      shift
      ;;
    --gamma)
      GAMMA=$2
      shift
      shift
      ;;
    --gaussx0)
      GAUSSX0=$2
      shift
      shift
      ;;
    --gaussy0)
      GAUSSY0=$2
      shift
      shift
      ;;
    --gaussz0)
      GAUSSZ0=$2
      shift
      shift
      ;;
    -h|--help)
      echo -e "$HELP"
      exit
      ;;
    --hrr)
      HRR=$2
      shift
      shift
      ;;
    --hc)
      HC=$2
      shift
      shift
      ;;
    --initialconditions)
      INIC=0
      IFILE=$2
      shift
      shift
      ;;
    --kappa)
      KAPPA=$2
      shift
      shift
      ;;
    --logginglevel)
      LOGGINGLEVEL=$2
      shift
      shift
      ;;
    --loggingfile)
      LOGGINGFILE=$2
      shift
      shift
      ;;
    --maxcycle)
      MAXCYCLE=$2
      shift
      shift
      ;;
    --maxiter)
      MIVALUESTMP=$(echo $2 | tr "," "\n")
      MAXITERVALUES=( $MIVALUESTMP )
      shift
      shift
      ;;
    --maxsolve)
      MAXSOLVE=$2
      shift
      shift
      ;;
    --ncycle)
      NCYCLE=$2
      shift
      shift
      ;;
    --nlevel)
      NLEVEL=$2
      shift
      shift
      ;;
    --ns|--navierstokes)
      NS=0
      SOLVER=NSSolver
      shift
      ;;
    --nste|--nstemp|--navierstokestemperature)
      NSTe=0
      SOLVER=NSTempSolver
      shift
      ;;
    --nstu|--nsturb|--navierstokesturbulence)
      NSTu=0
      SOLVER=NSTurbSolver
      shift
      ;;
    --nstc|--nstempcon|--navierstokestemperatureconcentration)
      NSTC=0
      SOLVER=NSTempConSolver
      shift
      ;;
    --nstt|--nstempturb|--navierstokestemperatureturbulence)
      NSTT=0
      SOLVER=NSTempTurbSolver
      shift
      ;;
    --nsttc|--nstempturbcon|--navierstokestemperatureturbulenceconcentration)
      NSTTC=0
      SOLVER=NSTempTurbConSolver
      shift
      ;;
    --nu)
      NU="$2"
      shift
      shift
      ;;
    -o|--out|--xml)
      OUT=$2
      XML=${OUT}
      shift
      shift
      ;;
    --pre|--pressure)
      PRE=0
      SOLVER=PressureSolver
      shift
      ;;
    --pressuretype)
      PRESSURETYPE=$2
      shift
      shift
      ;;
    --prediftype|--pressurediffusiontype)
      PRESSUREDIFFUSIONTYPE=$2
      shift
      shift
      ;;
    --pr)
      PR=$2
      shift
      shift
      ;;
    --prt)
      PRT=$2
      shift
      shift
      ;;
    --solavail)
      SOLAVAIL=$2
      shift
      shift
      ;;
    -s|--solver)
      SOLVER=$2
      shift
      shift
      ;;
    --sct)
      SCT=$2
      shift
      shift
      ;;
    --sigmax)
      SIGMAX=$2
      shift
      shift
      ;;
    --sigmay)
      SIGMAY=$2
      shift
      shift
      ;;
    --sigmaz)
      SIGMAZ=$2
      shift
      shift
      ;;
    --sourceconditions)
      SOUC=0
      SFILE=$2
      shift
      shift
      ;;
    --sourcetype)
      SOURCETYPE=$2
      shift
      shift
      ;;
    -t|--dt)
      DTVALUES=$(echo "$2" | tr "," "\n")
      shift
      shift
      ;;
    --tau)
      TAU=$2
      shift
      shift
      ;;
    --tempadvtype|--tempadvectiontype)
      TEMPADVECTIONTYPE=$2
      shift
      shift
      ;;
    --tempdifftype|--tempdiffusiontype)
      TEMPDIFFUSIONTYPE=$2
      shift
      shift
      ;;
    --tempsourceconditions)
      TEMPSOUC=0
      TSFILE=$2
      shift
      shift
      ;;
    --tempsourcetype)
      TEMPSOURCETYPE=$2
      shift
      shift
      ;;
    --tempsourcefct)
      TEMPSOURCEFCT=$2
      shift
      shift
      ;;
    --tend)
      TVALUES=$(echo "$2" | tr "," "\n")
      shift
      shift
      ;;
    --turbulencetype)
      TURBULENCETYPE=$2
      shift
      shift
      ;;
    -v|--verbose)
      VERBOSE=0
      shift
      ;;
    --vtk)
      VTK=Yes
      shift
      ;;
    --vtkplots)
      VTK=Yes
      VTKPLOTS=$2
      shift
      shift
      ;;
    -x|--nx)
      NXVALUESTMP=$(echo $2 | tr "," "\n")
      NXVALUES=( $NXVALUESTMP )
      shift
      shift
      ;;
    --xstart)
      XSTART=$2
      xSTART=$2
      shift
      shift
      ;;
    --xstartp|--xstartphysical)
      XSTART=$2
      shift
      shift
      ;;
    --xstartc|--xstartcomputational)
      xSTART=$2
      shift
      shift
      ;;
    --xend)
      XEND=$2
      xEND=$2
      shift
      shift
      ;;
    --xendp|--xendphysical)
      XEND=$2
      shift
      shift
      ;;
    --xendc|--xendcomputational)
      xEND=$2
      shift
      shift
      ;;
    -y|--ny)
      NYVALUESTMP=$(echo $2 | tr "," "\n")
      NYVALUES=( $NYVALUESTMP )
      shift
      shift
      ;;
    --ys)
      YS=$2
      shift
      shift
      ;;
    --ystart)
      YSTART=$2
      ySTART=$2
      shift
      shift
      ;;
    --ystartp|--ystartphysical)
      YSTART=$2
      shift
      shift
      ;;
    --ystartc|--ystartcomputational)
      ySTART=$2
      shift
      shift
      ;;
    --yend)
      YEND=$2
      yEND=$2
      shift
      shift
      ;;
    --yendp|--yendphysical)
      YEND=$2
      shift
      shift
      ;;
    --yendc|--yendcomputational)
      yEND=$2
      shift
      shift
      ;;
    -z|--nz)
      NZVALUESTMP=$(echo $2 | tr "," "\n")
      NZVALUES=( $NZVALUESTMP )
      shift
      shift
      ;;
    --zstart)
      ZSTART=$2
      zSTART=$2
      shift
      shift
      ;;
    --zstartp|--zstartphysical)
      ZSTART=$2
      shift
      shift
      ;;
      --zend)
      ZEND=$2
      zEND=$2
      shift
      shift
      ;;
    --zstartc|--zstartcomputational)
      zSTART=$2
      shift
      shift
      ;;
    --zendp|--zendphysical)
      ZEND=$2
      shift
      shift
      ;;
    --zendc|--zendcomputational)
      zEND=$2
      shift
      shift
      ;;
    *)
      POSITIONAL+=("$1")
      echo "unknown option: $1"
      shift
      ;;
  esac
done

COUNTER=1
COUNTERTMP="$(printf "%0${LEADINGNO}d" $COUNTER)"

for MAXITER in ${MAXITERVALUES[@]}
do
  for TEND in ${TVALUES[@]}
  do
    for DT in ${DTVALUES[@]}
    do
      for ((i=0; i<${#NXVALUES[@]}; i++))
        #for NX in ${NXVALUES[@]}
      do
        NX=${NXVALUES[$i]}
        NY=${NYVALUES[$i]}
        NZ=${NZVALUES[$i]}
        if [ $CALCNLEVEL -eq 0 ]; then calc_nlevel ; fi
        NAME=${XML/\%d/$COUNTERTMP}
        writeXML #write xml file
        ((COUNTER++))
        COUNTERTMP="$(printf "%0${LEADINGNO}d" $COUNTER)"
      done
    done
  done
done
