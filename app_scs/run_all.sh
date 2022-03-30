(./run data/config_amazon.txt 32 0.3) 2>&1 | tee amazon_log
(./run data/config_grqc.txt 32 0.3) 2>&1 | tee grqc_log
(./run data/config_hyves.txt 32 0.3) 2>&1 | tee hyves_log
(./run data/config_youtube.txt 32 0.3) 2>&1 | tee youtube_log