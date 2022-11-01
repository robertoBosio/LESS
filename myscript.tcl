open_project -reset subgraphiso

add_files Parameters.hpp
add_files Trie.hpp
add_files QueryVertex.hpp
add_files subgraphIsomorphism.hpp
add_files subisoWrap.hpp
add_files subisoWrap.cpp

add_files types.hpp
add_files utils.hpp
add_files hash_lookup3.hpp

add_files -tb "subiso-test.cpp data" -cflags "-Wno-unknown-pragmas" -csimflags "-Wno-unknown-pragmas"

set_top subisoWrap

open_solution "solution1" -flow_target vivado
set_part {xcvu9p-flga2104-1-e}
create_clock -period 10 -name default
#source "./solution1/directives.tcl"
#csim_design -clean
#csynth_design
#cosim_design
#export_design -format ip_catalog
