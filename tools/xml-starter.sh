#!/usr/bin/env bash
COMPILE=artss_serial
OUT=convergence_data.csv
FILES=-1

APPEND=1
OVER=1
VERBOSE=1

DELIM="\t"

YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0;m'

OPTIONS="Available options:\n
Executables:
  Production (with data output, visualization and analysis):
   ${YELLOW}-s${NC}
  ${YELLOW}--serial${NC}
  ${YELLOW}--artss_serial${NC}        \t Executable: artss_serial

  ${YELLOW} -m${NC}
  ${YELLOW}--multicore${NC}
  ${YELLOW}--artss_multicore_cpu${NC} \t Executable: artss_multicore_cpu

   ${YELLOW}-g${NC}
  ${YELLOW}--gpu${NC}
  ${YELLOW}--artss_gpu${NC}           \t Executable: artss_gpu
----
  Benchmarking (without output, visualization, analysis but with tracing for profiling):
  ${YELLOW}--sb${NC}
  ${YELLOW}--sbenchmark${NC}
  ${YELLOW}--artss_benchmark${NC}           \t Executable artss_serial_benchmark

  ${YELLOW}--mb${NC}
  ${YELLOW}--multicore_benchmark${NC}
  ${YELLOW}--artss_multicore_cpu_benchmark${NC} \t Executable artss_multicore_cpu_benchmark

  ${YELLOW}--gb${NC}
  ${YELLOW}--gpu_benchmark${NC}
  ${YELLOW}--artss_gpu_benchmark${NC}        \t Executable artss_gpu_benchmark

Other:
   ${YELLOW}-c${NC}
  ${YELLOW}--comment${NC}   \tadd comment to result file (use quotation marks)

  ${YELLOW}--cf${NC}
  ${YELLOW}--commentfile${NC}\tadd comment from specified file to result file

   ${YELLOW}-d${NC}
  ${YELLOW}--dir${NC}       \tspecify a parent directory of JuROr or the absolute path to it (default: '\$HOME')

  ${YELLOW} -e${NC}          \tif output file already exist, then:
           \t  [a]ppend
           \t  [c]ancel
           \t  [o]verwrite
           \t  [q]uestion me (default)

   ${YELLOW}-f${NC}
  ${YELLOW}--files${NC}    \tspecify files to execute (default: all xml in directory -  *.xml)

  ${YELLOW}--hoster${NC}
  ${YELLOW}--hostname${NC}  \tset host name for output file

   ${YELLOW}-o${NC}
  ${YELLOW}--output${NC}    \tspecify output file

   ${YELLOW}-v${NC}
  ${YELLOW}--verbose${NC}   \tverbose - print output to screen
\n"

INIT="Initialised values: output file=$OUT, verbose=false, compile version=$COMPILE, hoster name=$HOSTNAME\n"
HELP="$OPTIONS$INIT$EXAMPLE"

HEADER="Date${DELIM}Server${DELIM}Version${DELIM}Hardware${DELIM}t_end${DELIM}Nt${DELIM}dt${DELIM}nx${DELIM}ny${DELIM}nz${DELIM}dx${DELIM}dy${DELIM}dz${DELIM}RMS Error u${DELIM}RMS Error p${DELIM}RMS Error T${DELIM}Absolute Error u${DELIM}Absolute Error v${DELIM}Absolute Error w${DELIM}Absolute Error p${DELIM}Absolute Error T${DELIM}calctime[s]${DELIM}runtime[s]${DELIM}CUPS${DELIM}XML${DELIM}Output${DELIM}Obstacles${DELIM}Surfaces${DELIM}Dynamic${DELIM}Comment"

HOMEPATH=$HOME
HOSTER=$(hostname)

JUHYDRA=1
JURECA=1

function xml_parser {
   SOLUTION=$(xmllint -xpath "string(//ARTSS/solver/solution/@available)" $1)
   DYNAMIC=$(xmllint -xpath "string(//ARTSS/adaption/@dynamic)" $1)
   OBSTACLES=$(xmllint -xpath "string(//ARTSS/obstacles/@enabled)" $1)
   SURFACES=$(xmllint -xpath "string(//ARTSS/surfaces/@enabled)" $1)
   SOLVER=$(xmllint -xpath "string(//ARTSS/solver/@description)" $1)
   #echo "solver=$SOLVER"
   TEND=$(xmllint -xpath "string(//t_end)" $1 | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')
   DT=$(xmllint -xpath "string(//dt)" $1 | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')
   NT=$(bc <<< "scale=0; $TEND / $DT")
   # echo "t_end=$TEND, DT=$DT, NT=$NT"
   XSTART=$(xmllint -xpath "string(//x1)" $1 | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')
   XEND=$(xmllint -xpath "string(//x2)" $1 | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')
   NX=$(xmllint -xpath "string(//nx)" $1 | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')
   DX=$(bc <<< "scale=5; ($XEND - $XSTART) / ($NX )")
   # echo "x1=$XSTART, x2=$XEND, NX=$NX, DX=$DX"
   YSTART=$(xmllint -xpath "string(//y1)" $1 | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')
   YEND=$(xmllint -xpath "string(//y2)" $1 | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')
   NY=$(xmllint -xpath "string(//ny)" $1 | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')
   DY=$( echo "scale=5; ($YEND - $YSTART) / ($NY)" | bc)
   # echo "y1=$YSTART, y2=$YEND, NY=$NY, DY=$DY"
   ZSTART=$(xmllint -xpath "string(//z1)" $1 | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')
   ZEND=$(xmllint -xpath "string(//z2)" $1 | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')
   NZ=$(xmllint -xpath "string(//nz)" $1 | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')
   DZ=$(bc <<< "scale=5; ( $ZEND - $ZSTART ) / ( $NZ )")
   # echo "z1=$ZSTART, z2=$ZEND, NZ=$NZ, DZ=$DZ"

   MAXITER=$(xmllint -xpath "string(//max_iter)" $1 | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')
}

POSITIONAL=()
#parse command line arguments
while [[ $# -gt 0 ]]
do
  key="$1"
  case $key in
      -g|--gpu|--artss_gpu)
      COMPILE="artss_gpu"
      GPU=0
      shift
      ;;
    --gb|--gpu_benchmark---artss_gpu_benchmark)
      COMPILE="artss_gpu_benchmark "
      GPU=0
      shift
      ;;
    -c|--comment)
      COMMENT=$2
      shift
      shift
      ;;
    --cf|--commentfile)
      COMMENTFILE=$2
      shift
      shift
      ;;
    -h|--help)
      echo -e "$HELP"
      exit
      ;;
    -s|--serial|--artss_serial)
      COMPILE="artss_serial "
      shift
      ;;
    --sb|--serial_benchmark|--artss_serial_benchmark)
      COMPILE="artss_serial_benchmark "
      shift
      ;;
    --manual)
      COMPILE=$2
      shift
      shift
      ;;
    -m|--multicore|--artss_multicore_cpu)
      COMPILE="artss_multicore_cpu "
      GPU=0
      shift
      ;;
    --mb|--multicore_benchmark|--artss_multicore_cpu_benchmark)
      COMPILE="artss_multicore_cpu_benchmark "
      GPU=0
      shift
      ;;
    -d|--dir)
      HOMEPATH=$2
      shift
      shift
      ;;
    -e)
      case $2 in
        [aA] | [Aa][pP][pP][eE][nN][dD] )
          APPEND=0
          ;;
        [cC] | [cC][aA][nN][cC][eE][lL] )
          exit 0;
          ;;
        [oO] | [oO][vV][eE][rR][wW][rR][iI][tT][eE] )
          OVER=0
          ;;
        [qQ] | [qQ][uU][eE][sS][tT][iI][oO][nN] )
          APPEND=1
          ;;
      esac
      shift
      shift
      ;;
    -f|--files)
      FILES=$2
      shift
      shift
      ;;
    --hoster|--hostname)
      HOSTERNAME=$2
      shift
      shift
      ;;
    -o|--output)
      OLD=$OUT
      OUT=$2
      while [ -d $OUT ]
      do
        unset OUT
        read -t 60 -p "Given file is a directory. Please choose another name: " $OUT
        if [ -z ${OUT+x} ]
        then
          OUT=$OLD
          echo "There was no answer within 60 seconds. The script will continue with default name: $OUT"
        fi
      done
      shift
      shift
      ;;
    -v|--verbose)
      VERBOSE=0
      shift
      ;;
    *)
      POSITIONAL+=("$1")
      echo "unknown option: $1"
      shift
      ;;
  esac
done

if [ ${#COMMENTFILE} -gt 0 ]
then
  if [ -f $COMMENTFILE ]
  then
    if [ ${#COMMENT} -gt 0 ]
    then
      COMMENT=$COMMENT" "$(cat $COMMENTFILE)
    else
      COMMENT=$(cat $COMMENTFILE)
    fi
  else
    echo "$COMMENTFILE does not exist"
    unset $COMMENTFILE
  fi
fi

if [ "$FILES" = "-1" ]
then
  if [ $(ls *.xml |wc -l ) -eq 0 ]
  then
    (>&2 echo -e "${RED}XML files missing${NC}")
    exit 1;
  else
    FILES=$(ls *.xml)
  fi
else
  echo "given files: ${FILES}"
fi

HOST=$(hostname)
HARDWARE=$(cat /proc/cpuinfo | grep 'model name' | head -1 | awk '{split($0, a, ":"); print (a[2])}' | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')
SOCKET=$(lscpu | grep 'Socket(s)' | head -1 | awk '{split($0, a, ":"); print (a[2])}' | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')
CORESPERSOCKET=$(lscpu | grep 'Core(s) per socket'| head -1 | awk '{split($0, a, ":"); print (a[2])}' | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')
THREADS=$(lscpu | grep 'Thread(s) per core'| head -1 | awk '{split($0, a, ":"); print (a[2])}'| sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')
THREADSTIMESCORESPERSOCKET=$(bc <<< "scale=0; ${THREADS}*${CORESPERSOCKET}")
CORES=${SOCKET}x${THREADSTIMESCORESPERSOCKET}
HARDWARE="${HARDWARE}, ${CORES} cores"
SERVER=$(hostname)

#configuration for JURECA
if [[ $HOST = jrc* ]]
then
    SERVER="JURECA"
    JURECA=0
    module load CMake
    module load PGI
    module load CUDA
    export CUDA_LIB=${CUDA_ROOT}/lib64/
    export CUDA_INC=${CUDA_ROOT}/include/

    if [[ $COMPILE = *_multicore* ]]
    then export OMP_NUM_CORES=12 export ACC_BIND=yes export MP_BIND=yes export MP_BLIST=0,1,2,3,4,5,6,7,8,9,10,11; fi
    if [[ $COMPILE = *_gpu* ]]
    then HARDWARE="NVIDIA Kepler K80 with 2 GPUs, 562 MHz, 2x240 GB/s, 2x12 GB, 2x2496 cores"; fi
fi
#configuration for ias7139
if [[ "$HOST" = "ias7139" ]]
then
    SERVER="ias7139"
    IAS7139=0
    export CUDA_LIB=$CUDA_ROOT/lib64
    export CUDA_INC=$CUDA_ROOT/include
    if [[ $COMPILE = *_multicore* ]];
    then export OMP_NUM_CORES=8 export ACC_BIND=yes export OMP_PROC_BIND=yes export OMP_PLACES=0,1,2,3,4,5,6,7; fi
    if [[ $COMPILE = *_gpu* ]]
    then HARDWARE="NVIDIA Pascal P100, PCIe, 1382 MHz, 720 GB/s, 16 GB, 3584 cores"; fi
fi
if [[ $COMPILE = *_gpu* && $JURECA != 0 && $IAS7139 != 0 ]]; then HARDWARE=-; fi

FILEPATH=$(find $HOMEPATH -xdev -name $COMPILE)
LINES=$(find $HOMEPATH -xdev -name $COMPILE | wc -l)

if [ $LINES -eq 0 ]
then
   if [ -e $HOMEPATH$COMPILE ]
   then
       FILEPATH=$HOMEPATH$COMPILE
       LINES=$FILEPATH
   elif [ -e $HOMEPATH/$COMPILE ]
   then
       FILEPATH=$HOMEPATH/$COMPILE
       LINES=$FILEPATH
   fi
fi
if [ $LINES -ne 1 ]
then
   echo "There were zero or multiple version of $COMPILE in $HOMEPATH:"
   echo $FILEPATH
   echo "Please specify another location by using -d and check the given command line arguments. Choose option -h for more information."
   exit 1;
fi

if [[ -f $OUT && $OVER -eq 0 ]]; then rm $OUT; fi
if [[ -e $OUT && $APPEND -eq 1 ]]
then
    read -t 60 -p "$OUT already exist. Do you want to append the new output? (yes/overwrite/no) " APPEND
    APPEND=${APPEND}1
    case $APPEND in
       yes1)
           APPEND=0
           ;;
       overwrite1)
           APPEND=1
           rm $OUT
           ;;
       1)
           IFS='.' read -r -a TEMPOUT <<< "$OUT"
           OUT=${TEMPOUT[0]}"_"$TIME"."${TEMPOUT[-1]}
           echo "There was no answer within 60 seconds. The script will continue. Results will be written in: $OUT"
           ;;
       *)
           exit 1
           ;;
    esac
fi
if [ ! -e $OUT ]; then echo -e "$HEADER" >> $OUT; fi

for FILE in $FILES
do
    TIMEANDDATE=$(date)
    TIME=$(date +%s)
    echo -e "\n$TIMEANDDATE starting file $FILE"
    TMP=out_${TIME}.tmp
    xml_parser $FILE
    TMPNAME=$(echo $FILE | tr -dc '[:alnum:]' | tr '[:upper:]' '[:lower:]')
    TMPNAME=${TMPNAME}_${TIME}_$(hostname)

    START=$(date +%s)

    if [ $VERBOSE -eq 0 ]
    then
      $FILEPATH $FILE | tee $TMPNAME
    else
      $FILEPATH $FILE 2>> ${TMP}.error >  $TMPNAME
    fi

    if [ -f ${TMP}.error ]
    then
        ERRCOUNT=$( (wc -l) < ${TMP}.error)
        if [ $ERRCOUNT -eq 0 ]
        then
            rm ${TMP}.error
        else
            ERRORMES="An error occured, please check ${TMP}.error for more information"
            echo -e "${RED}An error occured, please check ${TMP}.error for more information${NC}" 
        fi
    fi

    END=$(date +%s)

    ERRU="-"
    ERRV="-"
    ERRW="-"
    ERRP="-"
    ERRT="-"
    RMSU="-"
    RMSP="-"
    RMST="-"
    if [ "$SOLUTION" = "Yes" ]
    then
      if [[ "$SOLVER" == *Advection* || "$SOLVER" == *Diffusion* ]] 
      then
        TIME=`cat $TMPNAME | tail -n 12 | head -1 | awk '{split($0, a, " "); print(a[3]/1000);}'`
        ERRU=`cat $TMPNAME | tail -n 8  | head -1 | awk '{split($0, a, "="); print(a[2]);}'`
        ERRV=`cat $TMPNAME | tail -n 5  | head -1 | awk '{split($0, a, "="); print(a[2]);}'`
        ERRW=`cat $TMPNAME | tail -n 2  | head -1 | awk '{split($0, a, "="); print(a[2]);}'`
        RMSU=`cat $TMPNAME | tail -n 18 | head -1 | awk '{split($0, a, "="); print(a[2]);}'`
        RMSP=`cat $TMPNAME | tail -n 17 | head -1 | awk '{split($0, a, "="); print(a[2]);}'`
        RMST=`cat $TMPNAME | tail -n 16 | head -1 | awk '{split($0, a, "="); print(a[2]);}'`
      elif [[ "$SOLVER" == NS* ]]
      then
        if [[ "$SOLVER" == *Temp* ]]
        then
          TIME=`cat $TMPNAME | tail -n 18 | head -1 | awk '{split($0, a, " "); print(a[3]/1000);}'`
          ERRU=`cat $TMPNAME | tail -n 14 | head -1 | awk '{split($0, a, "="); print(a[2]);}'`
          ERRV=`cat $TMPNAME | tail -n 11 | head -1 | awk '{split($0, a, "="); print(a[2]);}'`
          ERRW=`cat $TMPNAME | tail -n 8  | head -1 | awk '{split($0, a, "="); print(a[2]);}'`
          ERRP=`cat $TMPNAME | tail -n 5  | head -1 | awk '{split($0, a, "="); print(a[2]);}'`
          ERRT=`cat $TMPNAME | tail -n 2  | head -1 | awk '{split($0, a, "="); print(a[2]);}'`
          RMSU=`cat $TMPNAME | tail -n 24 | head -1 | awk '{split($0, a, "="); print(a[2]);}'`
          RMSP=`cat $TMPNAME | tail -n 23 | head -1 | awk '{split($0, a, "="); print(a[2]);}'`
          RMST=`cat $TMPNAME | tail -n 22 | head -1 | awk '{split($0, a, "="); print(a[2]);}'`
        else
          TIME=`cat $TMPNAME | tail -n 15 | head -1 | awk '{split($0, a, " "); print(a[3]/1000);}'`
          ERRU=`cat $TMPNAME | tail -n 11 | head -1 | awk '{split($0, a, "="); print(a[2]);}'`
          ERRV=`cat $TMPNAME | tail -n 8  | head -1 | awk '{split($0, a, "="); print(a[2]);}'`
          ERRW=`cat $TMPNAME | tail -n 5  | head -1 | awk '{split($0, a, "="); print(a[2]);}'`
          ERRP=`cat $TMPNAME | tail -n 2  | head -1 | awk '{split($0, a, "="); print(a[2]);}'`
          RMSU=`cat $TMPNAME | tail -n 21 | head -1 | awk '{split($0, a, "="); print(a[2]);}'`
          RMSP=`cat $TMPNAME | tail -n 20 | head -1 | awk '{split($0, a, "="); print(a[2]);}'`
          RMST=`cat $TMPNAME | tail -n 19 | head -1 | awk '{split($0, a, "="); print(a[2]);}'`
        fi
      elif [[ "$SOLVER" == "PressureSolver" ]]
      then
        TIME=`cat $TMPNAME | tail -n 5 | head -1 | awk '{split($0, a, " "); print(a[3]/1000);}'`
        ERRP=`cat $TMPNAME | tail -n 1  | head -1 | awk '{split($0, a, "="); print(a[2]);}'`
        RMSU=`cat $TMPNAME | tail -n 11 | head -1 | awk '{split($0, a, "="); print(a[2]);}'`
        RMSP=`cat $TMPNAME | tail -n 10 | head -1 | awk '{split($0, a, "="); print(a[2]);}'`
        RMST=`cat $TMPNAME | tail -n 9  | head -1 | awk '{split($0, a, "="); print(a[2]);}'`
      fi
    else
      TIME=`cat $TMPNAME | tail -n 1  | head -1 | awk '{split($0, a, " "); print(a[3]/1000);}'`
    fi

    NUMBEROFLINES=$(cat $TMPNAME | wc -l)
    if [ "$TIME" = "0" ]
    then
      CUPS="-"
    else
      CUPS=$(bc <<< "scale=2; ($NX+2) * ($NY+2) * ($NZ+2) * $NT / $TIME")
    fi
    TTIME=$(bc <<< "scale=2; $END - $START")
    NEWNAME=$TMP
    if [ -e $TMP ]
    then
      COUNTER=1
      while [ -e ${TMP}_${COUNTER} ]
      do
        ((COUNTER++))
      done
      NEWNAME=${TMP}_${COUNTER}
    fi
    mv $TMPNAME $NEWNAME
    #rmdir $DIRE
    echo -e "$TIMEANDDATE${DELIM}$SERVER${DELIM}$COMPILE${DELIM}$HARDWARE${DELIM}$TEND${DELIM}$NT${DELIM}$DT${DELIM}$NX${DELIM}$NY${DELIM}$NZ${DELIM}$DX${DELIM}$DY${DELIM}$DZ${DELIM}$RMSU${DELIM}$RMSP${DELIM}$RMST${DELIM}$ERRU${DELIM}$ERRV${DELIM}$ERRW${DELIM}$ERRP${DELIM}$ERRT${DELIM}$TIME${DELIM}$TTIME${DELIM}$CUPS${DELIM}$FILE${DELIM}$NEWNAME${DELIM}${OBSTACLES}${DELIM}${SURFACES}${DELIM}$DYNAMIC${DELIM}$COMMENT${DELIM}$ERRORMES" >> $OUT
    echo -e "\tdt=$DT (t_end=$TEND) size=${NX}x${NY}x${NZ} is finished. Data file: $FILE"
    echo -e "\ttmp data file=$TMP"
    echo "$(date) endtime: $TTIME"
done
echo "output file=$OUT"
