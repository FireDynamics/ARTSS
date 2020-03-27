HELP="
Write out CPU workload for each second

-o|--output\t output file
"

while [[ $# -gt 0 ]]
do
  key="$1"
  case $key in
    -h|--help)
      echo -e "$HELP"
      exit
      ;;
    -o|--output)
      OUTPUT=$2
      shift
      shift
      ;;
  esac
done

while true
do
  echo $(bash cpu_auslastung.sh) >> "${OUTPUT}"
  sleep 1
done

