MACS2
```
for i in `ls *sorted`
do
macs2 callpeak -t $i -f BAM -g hs -n $i.macs -B -q 0.05  &
done
```
