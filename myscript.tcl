

open_project -reset subgraphiso

add_files Parameters.hpp
add_files Trie.hpp
add_files queryVertex.hpp
add_files subgraphIsomorphism.hpp

add_files -tb "subiso-test.cpp data" -cflags "-Wno-unknown-pragmas" -csimflags "-Wno-unknown-pragmas"

set_top subgraphIsomorphism

open_solution "solution1" -flow_target vivado
set_part {xcvu11p-flga2577-1-e}
create_clock -period 10 -name default
#source "./solution1/directives.tcl"
#csim_design -clean
#csynth_design
#cosim_design
#export_design -format ip_catalog
