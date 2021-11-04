datasets="GSE1730_q"
num_compers="32"

for d in $datasets
do
for c in $num_compers
do
./run.sh ${d} ${c}
done
done
