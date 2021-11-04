home_path="/home/guimuguo/TthinkerQ"
monitor_path=$home_path/monitor
app_qc_path=$home_path/app_qc_ol
maximal_check_path=$home_path/maximal_check
data_path="/home/guimuguo/quick_data"

# Using best values from table 5
declare -A p1_arr
p1_arr["Enron_q"]="23;0.84"
p1_arr["hyves_q"]="22;0.82"
p1_arr["youtube_q"]="18;0.89"
p1_arr["patent_q"]="20;0.9"
p1_arr["GSE1730_q"]="30;0.85"
p1_arr["GSE10158_q"]="29;0.77"

# clear mem_log and disk_log before start
# Then, mem_sum.sh and disk_sum.sh will remove them after each exp.
#rm -f $monitor_path/mem_log;
#rm -f $monitor_path/disk_log;

f=$1 # dataset
compers=$2 

# make directory for each dataset
mkdir -p ${f}

p1_value=${p1_arr[$f]} # get p1 values by dataset name
arr_split=(${p1_value//;/ }) # split values by semi colon ;
t_size=${arr_split[0]}
gamma=${arr_split[1]}

# Best values from e1
tau_split="200"
tau_time="1"

base=${f}/"${compers}"
log="${base}_log"

(date;time $app_qc_path/run $data_path/${f} ${compers} ${gamma} ${t_size} ${tau_time} ${tau_split};date) 2>&1 | tee ${log};

#sleep 1
#$monitor_path/mem_sum.sh >> ${log};
#$monitor_path/disk_sum.sh >> ${log};

#result="${base}_result"
#cat output_* > ${result}
#wc -l ${result} >> ${log}

#max_result="${base}_max_result"
#($maximal_check_path/quasiCliques ${result} ${max_result}) 2>&1 | tee -a ${log}
#wc -l ${max_result} >> ${log}

# remove output files
#rm output_*

