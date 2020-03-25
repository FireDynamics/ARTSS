#!/bin/bash
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
COUNT=$(echo $diff | wc -l)

if [ $COUNT -eq 0 ]
then
  echo -e "${GREEN}equal${NC} [$1|$2]"
  rm $EXAMPLE $GENERATED
else
  if [[ $COUNT == 1 && $DIFF == *xml_filename* ]]
  then
    echo -e "${YELLOW}difference in xml_name${NC} [$1|$2]"
  else
    echo -e "${RED}files differ${NC} [$1|$2]\n"
  fi
  diff $GENERATED $EXAMPLE
fi
rm ${TMP1} ${TMP2} ${TMP1F} ${TMP2F} 
