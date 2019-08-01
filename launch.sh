#!/bin/bash
# Launch and oomph-convert

make two_d_nst_with_face_element && ./two_d_nst_with_face_element && rm ./RESLT/soln1.vtu && oomph-convert ./RESLT/soln1.dat >/dev/null && echo -e "\n\033[7;32mOomph-convert: done\033[0m \n"
# paraview --state=velocity_plot.pvsm
exit 0