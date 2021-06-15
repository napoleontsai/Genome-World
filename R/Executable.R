# Make all R scripts executable and install required R packages
dir <- getwd()
setwd(dir)
system('chmod +x RNaseT1.R RNaseA.R RNaseMC1.R RNaseU2.R CID_aw_ions.R CID_bx_ions.R CID_cy_ions.R CID_dz_ions.R CID_internal_wa.R CID_internal_wb.R CID_a-B.R RNAmod.R')
install.packages ("https://cran.r-project.org/src/contrib/stringr_1.4.0.tar.gz", repos = NULL)
install.packages ("https://cran.r-project.org/src/contrib/seqinr_4.2-8.tar.gz", repos = NULL)
install.packages ("https://cran.rstudio.com/src/contrib/httr_1.4.2.tar.gz", repos = NULL)
install.packages ("https://cran.rstudio.com/src/contrib/jsonlite_1.7.2.tar.gz", repos = NULL)
q()