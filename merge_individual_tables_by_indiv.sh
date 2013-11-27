indir=$1
infiles=$2;
outdir=$3;
indiv=$4;

(
echo -e "indiv_ref\tindiv_test\trank\tchr\tstart\tend\tmu\tp\tadjusted_p\twindow_size"
for f in `echo $infiles | tr ':' ' '` 
do
    fname=`echo $f | awk -F "/" '{print $(NF)}' | sed "s/dCGH_//g"`
    test=`echo $fname | awk -F "." '{print $1}'`
    ref=`echo $fname | awk -F "." '{print $2}'`
    zcat $indir/$f | awk -F '\t' -v OFS="\t" -v test=$test -v ref=$ref '{print ref,test,$0}'
done
) | gzip -c >$outdir/$indiv.segs.gz


