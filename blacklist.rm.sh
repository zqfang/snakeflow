

ls macs2/*broadPeak | while read id;
do
    sort -k1,1 -k2,2n ${id} | \
    bedtools subtract -a stdin \
                      -b /home/fangzq/genome/mouse/mm10-blacklist.v2.bed \
                      -A > beds/${id##*/}
done