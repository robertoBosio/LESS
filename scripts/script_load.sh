if [ $# -eq 0 ]
  then
    echo "No arguments supplied"
    exit
fi

res_file=results_${1}_$(date +%m%d%H%M).csv
#ip="192.168.1.53"
ip="192.168.3.1"
#echo "exec_time_avg,exec_time_std,power_avg,power_std,energy_avg,energy_std" > ${res_file}

#project=subiso_bd

mkdir overlay
#cp ${project}/${project}.runs/impl_1/design_1_wrapper.bit overlay/design_1.bit
#cp ${project}/${project}.gen/sources_1/bd/design_1/hw_handoff/design_1.hwh overlay/design_1.hwh
#cp -R data/ overlay/data
#cp design_1_wrapper.bit overlay/design_1.bit
#cp design_1.hwh overlay/design_1.hwh
#cp bitstream_dynfifov2/design_1_wrapper.bit overlay/design_1.bit
#cp bitstream_dynfifov2/design_1.hwh overlay/design_1.hwh
#cp ./host.py overlay/host.py
cp bitstream_${1}/design_1_wrapper.bit overlay/design_1.bit
cp bitstream_${1}/design_1.hwh overlay/design_1.hwh
cp ./host_${1}.py overlay/host.py
cp -R data/ overlay/data

# upload bitstream to sdcard
scp -r overlay root@${ip}:/home/xilinx/

# execute kernel
#cat ./host.py | ssh root@192.168.3.1 'python3 -'
ssh root@${ip} 'python3 /home/xilinx/overlay/host.py'
#>> ${res_file}

# cleanup
scp root@${ip}:/home/xilinx/overlay/results.txt ./${res_file}
ssh root@${ip} 'rm -r /home/xilinx/overlay'
rm -r overlay
