CASEPATH=$1
CHROMOSOME=$2
#MINSIZE=$3

module load bedtools

echo "##2 indicates the tool being examined for unique entries. 1 indicates the tool was compared against the tool being examined. 0 indicates the tool was not compared"
echo "##Path: $CASEPATH Chrom: $CHROMOSOME"

PINDEL=$(ls $CASEPATH*_somatic_pindel.$CHROMOSOME.vcf.gz)
ABRA=$(ls $CASEPATH*_somatic_abra.$CHROMOSOME.vcf.gz)
RUFUS=$(ls $CASEPATH*_somatic_rufus.$CHROMOSOME.vcf.gz)
PLATYPUS=$(ls $CASEPATH*_somatic_platypus.$CHROMOSOME.vcf.gz)

TOOLS=("PINDEL" "ABRA" "RUFUS" "PLATYPUS")
VCFS=("$PINDEL" "$ABRA" "$RUFUS" "$PLATYPUS")

BASENAME=${PINDEL%pindel*}

HEADER=""
for TOOL in ${TOOLS[@]}
do
  HEADER+="${TOOL}\t"
done
echo -e "${HEADER}UNIQUE_ENTRIES\tSHARED_ENTRIES\tEXCLUSIVE_ENTRIES"
for VCFNUM in ${!TOOLS[@]}
do
  BLENGTH=$((${#TOOLS[@]}-2))
  COMBINATIONS=$((2**$((${#TOOLS[@]}-1))))
  for ((i=0; i<$COMBINATIONS; i++))
  do
    TEMP=$i
    SWITCHES=""
    for ((j=BLENGTH; j>=0; j--))
    do
      if [ $((2**j)) -le $TEMP ]
      then
        SWITCHES+="1"
        TEMP=$(($TEMP-$((2**j))))
      else
        SWITCHES+="0"
      fi
    done
    SWITCHES="${SWITCHES:0:$VCFNUM}2${SWITCHES:$VCFNUM}"
    VCFB=""
    VCFBU=""
    VCFANTI=""
    for VCFTWO in ${!TOOLS[@]}
    do
      if [ ${SWITCHES:$VCFTWO:1} -eq "2" ]
      then
        VCFA=${VCFS[$VCFTWO]}
      elif [ ${SWITCHES:$VCFTWO:1} -eq "1" ]
      then
        VCFB+="${VCFS[$VCFTWO]} "
        VCFBU+="| bedtools intersect -a stdin -b ${VCFS[$VCFTWO]} -header -u "
      else
        VCFANTI+="| bedtools intersect -a stdin -b ${VCFS[$VCFTWO]} -header -v "
      fi
    done
    if [ $i -eq "0" ]
    then
      LINES=$(gunzip -c $VCFA | awk '{if( sqrt((length($4) - length($5))^2) > 5 ) print $0}' | grep "PASS" | wc -l)
      eval "gunzip -c $VCFA $VCFANTI > $BASENAME${TOOLS[$VCFNUM]}.$SWITCHES.EXCLUSIVE.${CHROMOSOME}.vcf"
      echo -e "$(sed "s/./&\t/g" <<<$SWITCHES)$LINES\t$LINES\t$(cat $BASENAME${TOOLS[$VCFNUM]}.$SWITCHES.EXCLUSIVE.${CHROMOSOME}.vcf | awk '{if( sqrt((length($4) - length($5))^2) > 5) print $0}' | grep "PASS" | wc -l)"
    else
      eval "gunzip -c $VCFA $VCFBU > $BASENAME${TOOLS[$VCFNUM]}.$SWITCHES.${CHROMOSOME}.vcf"
      eval "gunzip -c $VCFA $VCFBU $VCFANTI > $BASENAME${TOOLS[$VCFNUM]}.$SWITCHES.EXCLUSIVE.${CHROMOSOME}.vcf"
      echo -e "$(sed "s/./&\t/g" <<<$SWITCHES)$(bedtools intersect -a $VCFA -b ${VCFB::-1} -v | grep "PASS" | awk '{if( sqrt((length($4) - length($5))^2) > 5) print $0}' | wc -l)\t$(cat $BASENAME${TOOLS[$VCFNUM]}.$SWITCHES.${CHROMOSOME}.vcf | awk '{if( sqrt((length($4) - length($5))^2) > 5) print $0}' | grep "PASS" | wc -l)\t$(cat $BASENAME${TOOLS[$VCFNUM]}.$SWITCHES.EXCLUSIVE.${CHROMOSOME}.vcf | awk '{if( sqrt((length($4) - length($5))^2) > 5) print $0}' | grep "PASS" | wc -l)"
    fi
  done
done
  
