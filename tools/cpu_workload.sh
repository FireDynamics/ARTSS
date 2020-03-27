#!/bin/bash
HELP="
Print out the current cpu workload and the current time as seconds since 1970-01-01 (date +%s)
"
COMPILE="serial"
while [[ $# -gt 0 ]]
do
  key="$1"
  case $key in
    -h|--help)
      echo -e "$HELP"
      exit
      ;;
  esac
done

# readout number of cpu cores
CPUCORES=$(grep ^processor /proc/cpuinfo | wc -l)
time=$(date +%s)
# set basic measurement values 
for (( I=0; $I < $CPUCORES; I++ )); do
  TOTAL_LAST[$I]=0
  BUSY_LAST[$I]=0
  CPUDATA=$(grep ^"cpu"$I /proc/stat)
  BUSY_TICKS=$(echo $CPUDATA | awk -F' ' '{printf "%.0f",$2+$3+$4+$7+$8-$BL}')  
  TOTAL_TICKS=$(echo $CPUDATA | awk -F' ' '{printf "%.0f",$2+$3+$4+$5+$6+$7+$8}')

  # save values 
  TOTAL_LAST[$I]=$TOTAL_TICKS
  BUSY_LAST[$I]=$BUSY_TICKS
done

sleep 1

# calculate CPU load for every core 
for (( I=0; $I < $CPUCORES; I++ )); do    
  CPUDATA=$(grep ^"cpu"$I /proc/stat)
  BUSY_TICKS=$(echo $CPUDATA | awk -F' ' '{printf "%.0f",$2+$3+$4+$7+$8-$BL}')  
  TOTAL_TICKS=$(echo $CPUDATA | awk -F' ' '{printf "%.0f",$2+$3+$4+$5+$6+$7+$8}')

  BUSY_1000=$((1000*($BUSY_TICKS-${BUSY_LAST[$I]})/($TOTAL_TICKS-${TOTAL_LAST[$I]})))
  BUSY_INTEGER_FULL=$(($BUSY_1000/10))
  BUSY_INTEGER_FLOAT=$(echo $BUSY_1000 | rev | cut -c1-3 | rev )
  OUTPUT=$(echo -e "$OUTPUT\nCPU_core_"$I"_$BUSY_INTEGER_FULL.$BUSY_INTEGER_FLOAT")
done

# add cpu load of all cores
for i in $(echo "$OUTPUT"); do    
  ((INTEGER_FULL=$INTEGER_FULL+$(echo $i | cut -d "_" -f4 | cut -d . -f1)))
  ((INTEGER_FLOAT=$INTEGER_FLOAT+$(echo $i | cut -d "_" -f4 | cut -d . -f2)))
done

# calculate average cpu load
((INTEGER_FULL=$INTEGER_FULL/$CPUCORES))
((INTEGER_FLOAT=$INTEGER_FLOAT/$CPUCORES+1000))
INTEGER_FLOAT=$(echo $INTEGER_FLOAT | rev | cut -c1-3 | rev)  

# print values
echo "$time"
echo "$CPUCORES" | awk -F " " '{printf "CPU core(s)  : %9s \n", $1} '
echo "$OUTPUT" | grep -v ^$ | awk -F "_" '{printf "%-3s %-4s %-3s : %9s \n", $1, $2, $3, $4 }'
echo "$INTEGER_FULL.$INTEGER_FLOAT" | awk -F " " '{printf "CPU load avg : %9s \n", $1}'

exit 0
