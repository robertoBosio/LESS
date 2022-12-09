#open_project -reset subgraphiso
open_project subgraphiso

add_files "source/Parameters.hpp"
add_files "source/Trie.hpp"
add_files "source/QueryVertex.hpp"
add_files "source/subgraphIsomorphism.hpp"
add_files "source/subisoWrap.hpp"
add_files "source/subisoWrap.cpp"

add_files "source/types.hpp"
add_files "source/utils.hpp"
add_files "source/hash_lookup3.hpp"

add_files -tb "source/subiso-test.cpp data" -cflags "-Wno-unknown-pragmas" -csimflags "-Wno-unknown-pragmas"

set_top subisoWrap

open_solution "solution1" -flow_target vivado
set_part {xczu3eg-sbva484-1-e}
create_clock -period 5 -name default
#source "./solution1/directives.tcl"
#csim_design -clean
#csynth_design
#cosim_design
#cosim_design -user_stall stall_file.json -disable_deadlock_detection
#export_design -format ip_catalog
