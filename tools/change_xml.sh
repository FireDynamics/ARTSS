#!/bin/bash
HELP="
replace tags in xml. currently available

--da           \t change XML to start ARTSS at a specific time step. variables:\n\t\t\t1. time step\n\t\t\t2. h5 file name
h
 -i | --input  \t xml file to be changed
--logfile      \t specify the new log file name
--o | --output \t file name of where to write the changes
"
unset LOGFILE

while [[ $# -gt 0 ]]
do
  key="$1"
  case $key in 
    -h|--help)
      echo -e "${HELP}"
      exit
      ;;
    --da)
      DA=1
      DA_TIME=$2
      DA_FILE=$3
      DA_PORT=$4
      shift
      shift
      shift
      shift
      ;;
    --HRR)
      HRR=$2
      shift
      shift
      ;;
    -i | --input)
      INPUT=$2
      shift
      shift
      ;;
    --logfile)
      LOGFILE=$2.log
      shift
      shift
      ;;
    --loglevel)
      LOGLEVEL=$2
      shift
      shift
      ;;
    -o | --output)
      OUTPUT=$2
      shift
      shift
      ;;
    --tend)
      t_end=$2
      shift
      shift
      ;;
    --x0)
      x0=$2
      shift
      shift
      ;;
    *)
      echo "unknown parameter $1"
      shift
      ;;
  esac
done

if [ -z $INPUT -o ! -f $INPUT ]
then
  echo "specify a valid input file"
  exit
fi

TMP=.tmp.xml
cp $INPUT $TMP

if [ ! -z $LOGFILE ]
then
  echo "change logging file to $LOGFILE"
  sed 's/<logging file=".*" level="/<logging file="'${LOGFILE}'" level="/g' "${TMP}" > $OUTPUT
  cp $OUTPUT $TMP
fi


if [ ! -z $t_end ]
then
  echo "change t_end to $t_end"
  sed 's/<t_end>\s*[0-9]*\.*[0-9]*\s*<\/t_end>/<t_end> '${t_end}' <\/t_end>/g' "${TMP}" > "${OUTPUT}"
  cp $OUTPUT $TMP
fi

if [ ! -z $HRR ]
then
  echo "change HRR to $HRR"
  sed 's/<HRR>\s*[0-9]*\.*[0-9]*\s*<\/HRR>/<HRR> '${HRR}' <\/HRR>/g' "${TMP}" > "${OUTPUT}"
  cp $OUTPUT $TMP
fi

if [ ! -z $x0 ]
then
  echo "change x0 to $x0"
  sed 's/<x0>\s*[0-9]*\.*[0-9]*\s*<\/x0>/<x0> '${x0}' <\/x0>/g' "${TMP}" > "${OUTPUT}"
  cp $OUTPUT $TMP
fi

if [ ! -z $LOGLEVEL ]
then
  echo "change log level to $LOGLEVEL"
  sed 's/level=".*"/level="'${LOGLEVEL}'"/g' "${TMP}" > "${OUTPUT}"
  cp $OUTPUT $TMP
fi

if [ ! -z $DA ]
then
  echo "set parameter load_file=Yes; file=${DA_FILE}; time=${DA_TIME}; port=${DA_PORT}"
  if ! grep -xq "\s*<time>\s*[0-9]*\.*[0-9]*\s*</time>\s*" ${TMP} -a ! grep -xq "\s*<port>\s*[0-9]\s*</port>\s*" ${TMP}
  then
    sed 's/load_data=.*>/load_data="Yes" file="'${DA_FILE}'">\n    <time> '${DA_TIME}' <\/time>\n    <port> '${DA_PORT}' <\/port>/g' "${TMP}" > "${OUTPUT}"
  else
    if grep -xq "\s*<time>\s*[0-9]*\.*[0-9]*\s*</time>\s*" ${TMP}
    then
      sed 's/<time>\s*[0-9]*\.*[0-9]*\s*<\/time>/<time> '${DA_TIME}' <\/time>/g' ${TMP} > ${OUTPUT}
      cp $OUTPUT $TMP
      sed 's/load_data=.*>/load_data="Yes" file="'${DA_FILE}'">\n    <port> '${DA_PORT}' <\/port>/g' "${TMP}" > "${OUTPUT}"
    else
      sed 's/<port>\s*[0-9]*\s*<\/port>/<port> '${DA_PORT}' <\/port>/g' ${TMP} > ${OUTPUT}
      cp $OUTPUT $TMP
      sed 's/load_data=.*>/load_data="Yes" file="'${DA_FILE}'">\n    <time> '${DA_TIME}' </time>/g' "${TMP}" > "${OUTPUT}"
    fi
  fi
  cp $OUTPUT $TMP
  sed 's/save_vtk="\(Yes\|No\)"/save_vtk="No"/g' "${TMP}" > "${OUTPUT}"
  cp $OUTPUT $TMP
  sed 's/save_csv="\(Yes\|No\)"/save_csv="No"/g' "${TMP}" > "${OUTPUT}"
  cp $OUTPUT $TMP
fi
