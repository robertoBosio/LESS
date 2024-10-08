source settings.tcl

set proj "subgraphiso-kria"
set dir "source"

open_project -reset $proj
#open_project $proj

add_files "$dir/logger.cpp"
add_files "$dir/cmdlineparser.cpp"
add_files "$dir/subgraphIsomorphism.cpp" -cflags "-Wno-attributes -Wno-unknown-pragmas"
add_files -tb "$dir/subiso-test.cpp scripts dataset dataset2" -cflags "-O3 -Wno-attributes -Wno-unknown-pragmas"
set_top subgraphIsomorphism

open_solution -reset "solution1" -flow_target vivado
config_array_partition -complete_threshold 1
config_interface -m_axi_latency=1
# set_part {xczu3eg-sbva484-1-e}
set_part $XPART
create_clock -period 300MHz -name default
# set_clock_uncertainty 40%
# set_param hls.enable_fifo_io_regslice true

if {$CSIM == 1} {
    csim_design
}

if {$CSYNTH == 1} {
    csynth_design
}

if {$COSIM == 1} {
    cosim_design
}

if {$EXPORT == 1} {
    export_design
}

# if {$VIVADO_SYN == 1} {
# export_design -flow syn -rtl verilog
# }

# if {$VIVADO_IMPL == 1} {
# export_design -flow impl -rtl verilog
# }

exit

#source "./solution1/directives.tcl"
#csim_design -clean
#csynth_design
#cosim_design -user_stall stall_file.json -disable_deadlock_detection
#cosim_design -tool xsim -trace_level port -rtl verilog -wave_debug -user_stall my_cosim_stall.json
#export_design -format ip_catalog
