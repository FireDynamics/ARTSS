#!/bin/bash
FILES=$(ls *.xml)
BUFFER=1
CHECKVALUE=1
EXPANSIONSIZE=1
TIMESTEP=1
TEND=1
EXTRA=""
while [[ $# -gt 0 ]]
do
  key="$1"
  case $key in 
    -f)
      FILES=$2
      shift
      shift
      ;;
    --buffer)
      BUFFER=0
      shift
      ;;
    --checkvalue)
      CHECKVALUE=0
      shift
      ;;
    --expansionsize)
      EXPANSIONSIZE=0
      shift
      ;;
    --name)
      EXTRA="_$2"
      shift
      shift
      ;;
    --tend)
      TEND=0
      shift
      ;;
    --timestep)
      TIMESTEP=0
      shift
      ;;
  esac
done

for var in $FILES
do
  NAME="Test${EXTRA}"
  if [ -f "$var" ]
  then
    if [ $BUFFER -eq 0 ]
    then
      TMP=$(xmllint -xpath "string(//buffer)" $var | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')
      NAME="${NAME}_${TMP}"
    fi
    if [ $CHECKVALUE -eq 0 ]
    then
      TMP=$(xmllint -xpath "string(//check_value)" $var | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')
      NAME="${NAME}_${TMP}"
    fi
    if [ $EXPANSIONSIZE -eq 0 ]
    then
      TMP=$(xmllint -xpath "string(//expansion_size)" $var | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')
      NAME="${NAME}_${TMP}"
    fi
    if [ $TIMESTEP -eq 0 ]
    then
      TMP=$(xmllint -xpath "string(//timestep)" $var | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')
      NAME="${NAME}_${TMP}"
    fi
    if [ $TEND -eq 0 ]
    then
      TMP=$(xmllint -xpath "string(//t_end)" $var | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')
      NAME="${NAME}_${TMP}s"
    fi
    NAME=$NAME.xml
    if [ -f $NAME ]
    then
      echo "$NAME already exists. Skipping $var"
    else
      echo "$var -> $NAME"
      sed "s/$var/$NAME/g" "$var" > $NAME && rm $var
    fi
  else
    echo "File $var does not exist"
  fi
done
