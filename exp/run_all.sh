datasets="hyves_q patent_q Enron_q youtube_q"
num_compers="1 2 4 8 16 32"
##

#./run.sh Enron_q 32

#datasets="Enron_q hyves_q youtube_q patent_q"
#num_compers="1 2 4 8 16 32"

for d in $datasets
do
for c in $num_compers
do
./run.sh ${d} ${c}
done
done
