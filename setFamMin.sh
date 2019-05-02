## awk filter bam based on num reads collapsed. (for barcoding)

usage() {
    echo "Usage: $0 [-b <string>] [-a <string>] [-n <int>] "
    echo "        -b bam file"
    echo "        -a Optional filter barcodes"
    echo "        -n Optional, with {a}, number reads min in concensus"
    exit 1;
}
while getopts ":b:a:n:" x; do
    case "${x}" in
        b)
            b=${OPTARG}
            ;;
        a)
            a=${OPTARG}
            ;;
        n)
            n=${OPTARG}
            ;;
    esac
done

if [ -z "${b}" ]; then
    usage
fi
if [ -z "${a}" ]; then
    echo "No Barcode Filter"
fi
if [ -z "${n}" ]; then
    n=2
fi
echo "-b ${b}" > reAlign_log.txt
echo "-a ${a}" >> reAlign_log.txt
echo "-n ${n}" >> reAlign_log.txt

# awk filter
if [ ${a} == "XI" ]; then
    samtools view -h ${b} | awk '{if(/^@/) {print $0} else {if($17 ~ /XI/) split($17,barcodes,":"); else if($16 ~ /XI/) split($16,barcodes,":"); if(barcodes[3] >=2) print }}' | samtools sort -o ${b%.bam}.debar.bam -
fi
if [ ${a} == "XV" ]; then
    samtools view -h ${b} | awk '{if(/^@/) {print $0} else {if($14 ~ /XV/) split($14,barcodes,":"); else if($15 ~ /XV/) split($15,barcodes,":"); if(barcodes[3] >=2) print }}' | samtools sort -o ${b%.bam}.debar.bam -
fi
