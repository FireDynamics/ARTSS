#!/bin/bash

HELP="
Compare two xml, ignoring comments and trailing whitespaces. If the two files differ, a cleaned version (without comments and same identation) will be created.
"
if [[ $1 == "--help" || $1 == "-h" ]]
then
  echo -e "$HELP"
  exit
fi

YELLOW='\033[1;33m'
GREEN='\e[0;32m'
RED='\e[31m'
NC='\033[0;m'

TMP1=.tmp1.xml
TMP2=.tmp2.xml
TMP1F=.tmp1_format.xml
TMP2F=.tmp2_format.xml
GENERATED=file1.xml
EXAMPLE=file2.xml

xsltproc strip_comments.xsl $1 > ${TMP1}
xsltproc strip_comments.xsl $2 > ${TMP2}
xmllint --format ${TMP1} > ${TMP1F} 
xmllint --format ${TMP2} > ${TMP2F}
tidy -xml -iq ${TMP1F} > "$GENERATED"
tidy -xml -iq ${TMP2F} > "$EXAMPLE"
DIFF=$(diff $GENERATED $EXAMPLE)
COUNT=$(diff $GENERATED $EXAMPLE | wc -l)

if [ $COUNT == 0 ]
then
  echo -e "${GREEN}equal${NC} [$1|$2]"
  rm $EXAMPLE $GENERATED
else
  NAME1=$(xmllint -xpath "string(//xml_filename)" $1 | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')
  NAME2=$(xmllint -xpath "string(//xml_filename)" $2 | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')
  if [[ $COUNT == 4 && $DIFF == *xml_filename* ]]
  then
    echo -e "${YELLOW}difference in xml_name${NC} [$1|$2]"
    echo "filename: $NAME1 [$1]"
    echo "filename: $NAME2 [$2]"
    rm $EXAMPLE $GENERATED
  else
    diff $GENERATED $EXAMPLE
    echo -e "${RED}files differ${NC} [$NAME1|$NAME2]\n"
    mv $GENERATED cleaned_$NAME1
    mv $EXAMPLE cleaned_$NAME2
  fi
fi
rm ${TMP1} ${TMP2} ${TMP1F} ${TMP2F} 
