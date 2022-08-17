#!/bin/bash
HELP="
replace tags in xml. currently available

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

if [ -z $OUTPUT -o -f $OUTPUT ]
then
  echo "specify an output file which doesn't exist already"
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
fi
