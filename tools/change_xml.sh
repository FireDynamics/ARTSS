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
      shift
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
    -o | --output)
      OUTPUT=$2
      shift
      shift
      ;;
    --x0)
      x0=$2
      shift
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


if [ ! -z $x0 ]
then
  echo "change x0 to $x0"
  sed 's/<x0>\s*[0-9]*\.*[0-9]*\s*<\/x0>/<x0> '${x0}' <\/x0>/g' "${TMP}" > "${OUTPUT}"
  cp $OUTPUT $TMP
fi

if [ ! -z $DA ]
then
  echo "set parameter load_file=Yes; file=${DA_FILE}; time=${DA_TIME}"
  if grep -xq "\s*<time>\s*[0-9]*\.*[0-9]*\s*</time>\s*" ${TMP}
  then
    echo "time is there"
    sed 's/<time>\s*[0-9]*\.*[0-9]*\s*<\/time>/<time> '${DA_TIME}' <\/time>/g' ${TMP} > ${OUTPUT}
    cp $OUTPUT $TMP
    sed 's/load_data=.*>/load_data="Yes" file="'${DA_FILE}'">/g' "${TMP}" > "${OUTPUT}"
  else
    echo "time isn't there"
    sed 's/load_data=.*>/load_data="Yes" file="'${DA_FILE}'">\n    <time> '${DA_TIME}' <\/time>/g' "${TMP}" > "${OUTPUT}"
  fi 
  cp $OUTPUT $TMP
fi
