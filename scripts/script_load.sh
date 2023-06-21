if [ $# -lt 2 ]; then
    echo "Provide the board used and the description of the test"
    exit
fi

if [ $1 = "kria" ]; then
    user="ubuntu"
    path="/home/ubuntu/"
    ip="192.168.99.170"
    device="kria170"
elif [ $1 = "ultra96" ]; then
    user="root"
    path="/home/xilinx/"
    ip="192.168.3.1"
    device="${user}@${ip}"
else
    echo "board parameter must be kria or ultra96"
    exit
fi

res_file=results_${1}_$(date +%m%d%H%M).csv
#echo "exec_time_avg,exec_time_std,power_avg,power_std,energy_avg,energy_std" > ${res_file}

#project=subiso_bd

if [ -d overlay ]; then
    rm overlay -r
fi

mkdir overlay
cp ../subiso_bd_kria/subiso_bd_kria.runs/impl_1/design_1_wrapper.bit overlay/design_1.bit
cp ../subiso_bd_kria/subiso_bd_kria.gen/sources_1/bd/design_1/hw_handoff/design_1.hwh overlay/design_1.hwh
cp ./test.txt overlay/test.txt
cp ./host.py overlay/host.py

mkdir overlay/data
cp ../dataset/benchmark/labelled/email-EnronRM.csv overlay/data/
cp ../dataset/benchmark/labelled/musae_githubRM.csv overlay/data/
cp ../dataset/benchmark/labelled/musae_facebookRM.csv overlay/data/
cp ../dataset/benchmark/labelled/twitter_combinedRM.csv overlay/data/
# cp ../dataset/benchmark/labelled/simpleRM.csv overlay/data/
# cp ../dataset/benchmark/labelled/dataEdges2RM.csv overlay/data/
# cp ../dataset/benchmark/labelled/wiki_simpleRM.csv overlay/data/
cp ../dataset/benchmark/queries/query0RM.csv overlay/data/
cp ../dataset/benchmark/queries/query1RM.csv overlay/data/
cp ../dataset/benchmark/queries/query2RM.csv overlay/data/
cp ../dataset/benchmark/queries/query3RM.csv overlay/data/
cp ../dataset/benchmark/queries/query4RM.csv overlay/data/
cp ../dataset/benchmark/queries/query5RM.csv overlay/data/

# upload bitstream to sdcard
# scp -r overlay root@${ip}:/home/xilinx/
scp -r overlay ${device}:${path}

# execute kernel
#cat ./host.py | ssh root@192.168.3.1 'python3 -'
ssh ${device} "source /etc/profile && python3 ${path}overlay/host.py ${path}overlay/"
#>> ${res_file}

# cleanup
scp ${device}:${path}overlay/results.txt ./${res_file}
echo $2 >> ./${res_file}

if [ $1 = "kria" ]; then
    ssh ${device} "rm -r ${path}overlay && xmutil unloadapp k26-starter-kits && xmutil loadapp k26-starter-kits"
else
    ssh ${device} "rm -r ${path}overlay"
fi

rm -r overlay
