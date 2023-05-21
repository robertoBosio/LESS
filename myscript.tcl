set cosim 1

set proj "subgraphiso-kria"
set dir "source"
if { $cosim == 1 } {
    set proj "subgraphiso-kria-nopreproc"
    set dir "source-nopreproc"
}

open_project -reset $proj

add_files "$dir/subisoWrap.cpp" -cflags "-I ./$dir"
add_files -tb "$dir/subiso-test.cpp data dataset" -cflags "-Wno-attributes -Wno-unknown-pragmas" -csimflags "-Wno-unknown-pragmas"

set_top subgraphIsomorphism

open_solution "solution1" -flow_target vivado
config_array_partition -complete_threshold 1
config_interface -m_axi_latency=128
# set_part {xczu3eg-sbva484-1-e}
set_part {xck26-sfvc784-2LV-c}
# set_part {xczu5eg-sfvc784-1-i}
create_clock -period 300MHz -name default
#source "./solution1/directives.tcl"
#csim_design -clean
#csynth_design
#cosim_design
#cosim_design -user_stall stall_file.json -disable_deadlock_detection
#export_design -format ip_catalog
