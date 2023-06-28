open_project /home/roberto/Documents/subgraph-iso/subiso_bd/subiso_bd.xpr

source subiso_bd.tcl
#open_bd_design {/home/roberto/Documents/subgraph-iso/subiso_bd/subiso_bd.srcs/sources_1/bd/design_1/design_1.bd}

#upgrading the ip
update_compile_order -fileset sources_1
report_ip_status -name ip_status 
upgrade_ip -vlnv xilinx.com:hls:subgraphIsomorphism:1.0 [get_ips  design_1_subgraphIsomorphism_0_3] -log ip_upgrade.log

#generate
export_ip_user_files -of_objects [get_ips design_1_subgraphIsomorphism_0_3] -no_script -sync -force -quiet
generate_target all [get_files  /home/roberto/Documents/subgraph-iso/subiso_bd/subiso_bd.srcs/sources_1/bd/design_1/design_1.bd]

#synth
reset_run synth_1
set_property strategy Flow_PerfOptimized_high [get_runs synth_1]
set_property STEPS.SYNTH_DESIGN.ARGS.RETIMING true [get_runs synth_1]
launch_runs synth_1 -jobs 8
wait_on_runs synth_1

#implementation
set_property strategy Performance_ExplorePostRoutePhysOpt [get_runs impl_1]
launch_runs impl_1 -jobs 8
wait_on_runs impl_1

#bitstream
launch_runs impl_1 -to_step write_bitstream -jobs 8
wait_on_runs impl_1
