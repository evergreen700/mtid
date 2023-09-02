CASEPATH=$1
CHROMOSOME=$2

module load bedtools

echo "Path: $CASEPATH Chrom: $CHROMOSOME"

PINDEL=$(ls $CASEPATH/*_somatic_pindel.$CHROMOSOME.vcf)
ABRA=$(ls $CASEPATH/*_somatic_abra.$CHROMOSOME.vcf)
RUFUS=$(ls $CASEPATH/*_somatic_rufus.$CHROMOSOME.vcf.gz)
PLATYPUS=$(ls $CASEPATH/*_somatic_platypus.$CHROMOSOME.vcf)

TOOLS=("PINDEL" "ABRA" "RUFUS" "PLATYPUS")
VCFS=("$PINDEL" "$ABRA" "$RUFUS" "$PLATYPUS")

bedtools intersect -a "$PINDEL" -b "$ABRA" -v | wc -l

for VCFNUM in ${!TOOLS[@]}
do
  echo "${TOOLS[$VCFNUM]}: ${VCFS[$VCFNUM]}"
  BLENGTH=$((${#TOOLS[@]}-2))
  COMBINATIONS=$((2**$((${#TOOLS[@]}-1))))
  for ((i=1; i<$COMBINATIONS; i++))
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
    echo "$i: ${SWITCHES:0:$VCFNUM}*${SWITCHES:$VCFNUM}"
  done
done
  
