#!/bin/bash
# filepath=/mydata/datasets/
filepath=/mnt/sdd/fkyang/datasets/
dataset=("chromium.zip" "linuxDatasets.zip" "sinadatasets.zip" "imagesubuntu.zip" "sof2019.zip" "synthesis.zip")
# dataset=("linuxDatasets.zip" "sinadatasets.zip" "imagesubuntu.zip" "sof2019.zip" "synthesis.zip")

# deduplication ratio ：
#     chromium.zip :        2.528099017618521
#     linuxDatasets.zip :   2.503153145954249
#     sinadatasets.zip :    6.036346189599833
#     imagesubuntu.zip :    2.209456884885373
#     sof2019.zip :         1.000008217448588
#     synthesis.zip :       11.44383818075717

# filepath=/mnt/sdc/fkyang/
# dataset=("chromium/128.0.6613.164.tar.gz")
current_time=$(date "+%Y-%m-%d")

# ohash_threads_OH_G
# ohash_threads_OH_X
# ohash_threads_OH_GW
# ohash_threads_OD_G
# ohash_threads_OD_X

# ohash_threads_N_G
# ohash_threads_N_X


echo "------------------------start...-------------------"
(date "+%Y-%m-%d_%H:%M:%S")

i=1
for file in ${dataset[@]}   # 1 Feature 1 SF
do
    filebasename=$(basename "$file")
    echo "------------------------starting...-------------------"
    (date "+%Y-%m-%d_%H:%M:%S")
    echo "The file processed is ${filepath}${file}"

    # echo ""
    # echo "------------------------1.run ohash_threads_OH_G using gdelta ----------------------"
    # (date "+%Y-%m-%d_%H:%M:%S")
    # echo ""
    # i=$((i+1))
    # echo "./ohash_threads_OH_G ${filepath}${file} 1 1 _all_datasets | tee -a ./results/all_datasets_ohash_threads_OH_G_${current_time}.log"
    # echo ""
    # # taskset -c $i 
    # ./ohash_threads_OH_G ${filepath}${file} 1 1 _all_datasets | tee -a ./results/all_datasets_ohash_threads_OH_G_${current_time}.log

    echo ""
    echo "------------------------2.run ohash_threads_OH_GW using gdeltaWhash ----------------------"
    (date "+%Y-%m-%d_%H:%M:%S")
    echo ""
    i=$((i+1))
    echo "./ohash_threads_OH_GW ${filepath}${file} 1 1 _all_datasets | tee -a ./results/all_datasets_ohash_threads_OH_GW_${current_time}.log"
    echo ""
    # taskset -c $i 
    ./ohash_threads_OH_GW ${filepath}${file} 1 1 _all_datasets | tee -a ./results/all_datasets_ohash_threads_OH_GW_${current_time}.log

    # echo ""
    # echo "------------------------3.run ohash_threads_OH_X using xdelta ----------------------"
    # (date "+%Y-%m-%d_%H:%M:%S")
    # echo ""
    # i=$((i+1))
    # echo "./ohash_threads_OH_X ${filepath}${file} 1 1 _all_datasets | tee -a ./results/all_datasets_ohash_threads_OH_X_${current_time}.log"
    # echo ""
    # # taskset -c $i 
    # ./ohash_threads_OH_X ${filepath}${file} 1 1 _all_datasets | tee -a ./results/all_datasets_ohash_threads_OH_X_${current_time}.log

    # echo ""
    # echo "------------------------4.run ohash_threads_OD_G using gdelta ----------------------"
    # (date "+%Y-%m-%d_%H:%M:%S")
    # echo ""
    # i=$((i+1))
    # echo "./ohash_threads_OD_G ${filepath}${file} 1 1 _all_datasets | tee -a ./results/all_datasets_ohash_threads_OD_G_${current_time}.log"
    # echo ""
    # # taskset -c $i 
    # ./ohash_threads_OD_G ${filepath}${file} 1 1 _all_datasets | tee -a ./results/all_datasets_ohash_threads_OD_G_${current_time}.log

    # echo ""
    # echo "------------------------5.run ohash_threads_OD_X using xdelta ----------------------"
    # (date "+%Y-%m-%d_%H:%M:%S")
    # echo ""
    # i=$((i+1))
    # echo "./ohash_threads_OD_X ${filepath}${file} 1 1 _all_datasets | tee -a ./results/all_datasets_ohash_threads_OD_X_${current_time}.log"
    # echo ""
    # # taskset -c $i 
    # ./ohash_threads_OD_X ${filepath}${file} 1 1 _all_datasets | tee -a ./results/all_datasets_ohash_threads_OD_X_${current_time}.log

    # echo ""
    # echo "------------------------6.run ohash_ntrans using gdelta ----------------------"
    # (date "+%Y-%m-%d_%H:%M:%S")
    # echo ""
    # i=$((i+1))
    # echo "./ohash_threads_N_G ${filepath}${file} 1 1 _all_datasets | tee -a ./results/all_datasets_ohash_threads_N_G_${current_time}.log"
    # echo ""
    # # taskset -c $i 
    # ./ohash_threads_N_G ${filepath}${file} 1 1 _all_datasets | tee -a ./results/all_datasets_ohash_threads_N_G_${current_time}.log

    # echo ""
    # echo "------------------------7.run ohash_ntrans using xdelta ----------------------"
    # (date "+%Y-%m-%d_%H:%M:%S")
    # echo ""
    # i=$((i+1))
    # echo "./ohash_threads_N_X ${filepath}${file} 1 1 _all_datasets | tee -a ./results/all_datasets_ohash_threads_N_X_${current_time}.log"
    # echo ""
    # # taskset -c $i 
    # ./ohash_threads_N_X ${filepath}${file} 1 1 _all_datasets | tee -a ./results/all_datasets_ohash_threads_N_X_${current_time}.log


    echo ""
    echo "---------------------all_datasets finished----------------"
    (date "+%Y-%m-%d_%H:%M:%S")
    echo "------------------------finished------------------------"
done
echo "------------------------"
echo "using 3 features for 3 sfs........."
echo "------------------------"

# for file in ${dataset[@]}   # 12 Features 3 SFs
# do
#     filebasename=$(basename "$file")
#     echo "------------------------starting...-------------------"
#     (date "+%Y-%m-%d_%H:%M:%S")
#     echo "The file processed is ${filepath}${file}"

#     echo ""
#     echo "------------------------1.run ohash_threads_OH_G using gdelta ----------------------"
#     (date "+%Y-%m-%d_%H:%M:%S")
#     echo ""
#     i=$((i+1))
#     echo "./ohash_threads_OH_G ${filepath}${file}  3 3 _all_datasets | tee -a ./results/all_datasets_ohash_threads_OH_G_${current_time}.log"
#     echo ""
#     # taskset -c $i 
#     ./ohash_threads_OH_G ${filepath}${file}  3 3 _all_datasets | tee -a ./results/all_datasets_ohash_threads_OH_G_${current_time}.log

#     echo ""
#     echo "------------------------2.run ohash_threads_OH_GW using gdeltaWhash ----------------------"
#     (date "+%Y-%m-%d_%H:%M:%S")
#     echo ""
#     i=$((i+1))
#     echo "./ohash_threads_OH_GW ${filepath}${file}  3 3 _all_datasets | tee -a ./results/all_datasets_ohash_threads_OH_GW_${current_time}.log"
#     echo ""
#     # taskset -c $i 
#     ./ohash_threads_OH_GW ${filepath}${file}  3 3 _all_datasets | tee -a ./results/all_datasets_ohash_threads_OH_GW_${current_time}.log

#     echo ""
#     echo "------------------------3.run ohash_threads_OH_X using xdelta ----------------------"
#     (date "+%Y-%m-%d_%H:%M:%S")
#     echo ""
#     i=$((i+1))
#     echo "./ohash_threads_OH_X ${filepath}${file}  3 3 _all_datasets | tee -a ./results/all_datasets_ohash_threads_OH_X_${current_time}.log"
#     echo ""
#     # taskset -c $i 
#     ./ohash_threads_OH_X ${filepath}${file}  3 3 _all_datasets | tee -a ./results/all_datasets_ohash_threads_OH_X_${current_time}.log


#     echo ""
#     echo "------------------------4.run ohash_threads_OD_G using gdelta ----------------------"
#     (date "+%Y-%m-%d_%H:%M:%S")
#     echo ""
#     i=$((i+1))
#     echo "./ohash_threads_OD_G ${filepath}${file}  3 3 _all_datasets | tee -a ./results/all_datasets_ohash_threads_OD_G_${current_time}.log"
#     echo ""
#     # taskset -c $i 
#     ./ohash_threads_OD_G ${filepath}${file}  3 3 _all_datasets | tee -a ./results/all_datasets_ohash_threads_OD_G_${current_time}.log

#     echo ""
#     echo "------------------------5.run ohash_threads_OD_G using xdelta ----------------------"
#     (date "+%Y-%m-%d_%H:%M:%S")
#     echo ""
#     i=$((i+1))
#     echo "./ohash_threads_OD_X ${filepath}${file}  3 3 _all_datasets | tee -a ./results/all_datasets_ohash_threads_OD_X_${current_time}.log"
#     echo ""
#     # taskset -c $i 
#     ./ohash_threads_OD_X ${filepath}${file}  3 3 _all_datasets | tee -a ./results/all_datasets_ohash_threads_OD_X_${current_time}.log


#     # echo ""
#     # echo "------------------------6.run ohash_threads_N_G using gdelta ----------------------"
#     # (date "+%Y-%m-%d_%H:%M:%S")
#     # echo ""
#     # i=$((i+1))
#     # echo "./ohash_threads_N_G ${filepath}${file}  3 3 _all_datasets | tee -a ./results/all_datasets_ohash_threads_N_G_${current_time}.log"
#     # echo ""
#     # # taskset -c $i 
#     # ./ohash_threads_N_G ${filepath}${file}  3 3 _all_datasets | tee -a ./results/all_datasets_ohash_threads_N_G_${current_time}.log


#     # echo ""
#     # echo "------------------------7.run ohash_threads_N_X using xdelta ----------------------"
#     # (date "+%Y-%m-%d_%H:%M:%S")
#     # echo ""
#     # i=$((i+1))
#     # echo "./ohash_threads_N_X ${filepath}${file}  3 3 _all_datasets | tee -a ./results/all_datasets_ohash_threads_N_X_${current_time}.log"
#     # echo ""
#     # # taskset -c $i 
#     # ./ohash_threads_N_X ${filepath}${file}  3 3 _all_datasets | tee -a ./results/all_datasets_ohash_threads_N_X_${current_time}.log


#     echo ""
#     echo "---------------------all_datasets finished----------------"
#     (date "+%Y-%m-%d_%H:%M:%S")
#     echo "------------------------finished------------------------"
# done

# wait
(date "+%Y-%m-%d_%H:%M:%S")
echo "all datasets is finisehd"
