user="root"
path="/home/ubuntu/"
ip="192.168.99.170"

if [ -d overlay ]; then
    rm overlay -r
fi

# Input file path
input_file="run_list.txt"

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
cp ./run_list.txt overlay/run_list.txt
cp ./host.py overlay/host.py

mkdir overlay/data
for file_path in "${!unique_paths[@]}"; do
  cp "$file_path" overlay/data/
done

# copy overlay directory to Kria
# scp -r overlay ${user}@${ip}:${path}
scp -r overlay kria170:${path}

# execute kernel
# ssh ${user}@${ip} "source /etc/profile && python3 ${path}overlay/host.py ${path}overlay/"
ssh kria170 "source /etc/profile && python3 ${path}overlay/host.py ${path}overlay/"

# cleanup and reloading bitstream to control fan speed
# ssh ${user}@${ip} "rm -r ${path}overlay && xmutil unloadapp k26-starter-kits && xmutil loadapp k26-starter-kits"
ssh kria170 "rm -r ${path}overlay && xmutil unloadapp k26-starter-kits && xmutil loadapp k26-starter-kits"

rm -r overlay
