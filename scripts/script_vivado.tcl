set proj_name "less_bd"
set proj_dir "./less_bd"
set proj_dir_hls "./less_hls"

if {[file exists $proj_dir/$proj_name.xpr]} {
    puts "Project $proj_name already exists. Opening the existing project..."
    open_project "$proj_dir/$proj_name.xpr"
    open_bd_design "$proj_dir/$proj_name.srcs/sources_1/bd/design_1/design_1.bd"
    
    # Upgrading the ip
    update_compile_order -fileset sources_1
    report_ip_status
    upgrade_ip -vlnv xilinx.com:hls:subgraphIsomorphism:1.0 [get_ips design_1_subgraphIsomorphism_0_0] -log ip_upgrade.log

} else {
    puts "Project $proj_name does not exists. Creating a new project..."
    create_project $proj_name $proj_dir -part xck26-sfvc784-2LV-c
    set_property board_part xilinx.com:kv260_som:part0:1.4 [current_project]
    set_property ip_repo_paths "$proj_dir_hls/solution1/impl/ip" [current_project]
   
    # Create block design 
    source scripts/less_bd.tcl

    # Make wrapper
    make_wrapper -files [get_files "$proj_dir/$proj_name.srcs/sources_1/bd/design_1/design_1.bd"] -top
    add_files -norecurse "$proj_dir/$proj_name.gen/sources_1/bd/design_1/hdl/design_1_wrapper.v"
}

#generate
export_ip_user_files -of_objects [get_ips design_1_subgraphIsomorphism_0_3] -no_script -sync -force -quiet
generate_target all [get_files "$proj_dir/$proj_name.srcs/sources_1/bd/design_1/design_1.bd"]

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

exit
