source scripts/settings.tcl
set proj "less_hls"
set dir "source"

open_project $proj
add_files "$dir/logger.cpp"
add_files "$dir/cmdlineparser.cpp"
add_files "$dir/subgraphIsomorphism.cpp" -cflags "-Wno-attributes -Wno-unknown-pragmas"
add_files -tb "$dir/subiso-test.cpp scripts dataset_example" -cflags "-O3 -Wno-attributes -Wno-unknown-pragmas"
set_top subgraphIsomorphism

open_solution -reset "solution1" -flow_target vivado
config_array_partition -complete_threshold 1
config_interface -m_axi_latency=1
set_part $XPART
create_clock -period 300MHz -name default

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

if {$VIVADO_SYN == 1} {
    export_design -flow syn -rtl verilog
}

if {$VIVADO_IMPL == 1} {
    export_design -flow impl -rtl verilog
}

exit
