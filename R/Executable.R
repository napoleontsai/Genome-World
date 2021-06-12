
# Make all R scripts executable
dir <- getwd()
setwd(dir)
system('chmod +x RNaseT1.R RNaseA.R RNaseMC1.R RNaseU2.R CID_aw_ions.R CID_bx_ions.R CID_cy_ions.R CID_dz_ions.R CID_internal_wa.R CID_internal_wb.R CID_a-B.R')
q()