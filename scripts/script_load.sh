user="root"
path="/home/ubuntu/"
ip="192.168.99.170"
device="kria170"
res_file=results_${1}_$(date +%m%d%H%M).csv

if [ $# -lt 2 ]; then
    echo "Provide the board used and the description of the test"
    exit
fi

# ssh_string="${user}@${ip}"
ssh_string="$device"
if [ -d overlay ]; then
    rm overlay -r
fi

# Input file path
input_file="run_list_fpga.txt"

# Check if the input file exists
if [ ! -f "$input_file" ]; then
  echo "Input file $input_file does not exist."
  exit 1
fi

# Read paths from input file into an array
readarray -t lines < <(grep -v "^#" "$input_file")

# Create a set to store unique paths
declare -A unique_paths

# Iterate over the lines and extract paths
for line in "${lines[@]}"; do

  # Remove leading/trailing whitespace from line
  line=$(echo "$line" | awk '{$1=$1};1')

  # Split the line into two paths
  IFS=' ' read -ra paths <<< "$line"

  # Iterate over the paths and add them to the set
  for file_path in "${paths[@]}"; do

    if [[ $file_path =~ ^[0-9]+$ ]]; then
      continue
    fi

    # Add path to the set
    unique_paths["$file_path"]=1
  done
done


mkdir overlay
cp ../subiso_bd/subiso_bd.runs/impl_1/design_1_wrapper.bit overlay/design_1.bit
cp ../subiso_bd/subiso_bd.gen/sources_1/bd/design_1/hw_handoff/design_1.hwh overlay/design_1.hwh
# cp bitstream/design_1.bit overlay/design_1.bit
# cp bitstream/design_1.hwh overlay/design_1.hwh
cp ./${input_file} overlay/test.txt
cp ./host.py overlay/host.py

mkdir overlay/data
for file_path in "${!unique_paths[@]}"; do
  cp "$file_path" overlay/data/
done

# copy overlay directory to Kria
scp -r overlay ${ssh_string}:${path}
# scp -r overlay kria170:${path}

# execute kernel
ssh ${ssh_string} "source /etc/profile && python3 ${path}overlay/host.py ${path}overlay/"

# copy result
scp ${ssh_string}:${path}overlay/results.txt ./${res_file}
echo $2 >> ./${res_file}

# cleanup and reloading bitstream to control fan speed
ssh ${ssh_string} "rm -r ${path}overlay && xmutil unloadapp k26-starter-kits && xmutil loadapp k26-starter-kits"

rm -r overlay
