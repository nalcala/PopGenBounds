## code to prepare `DATASET` dataset goes here
#usethis::use_data(DATASET, overwrite = TRUE)

library(readxl)
library(tidyverse)

# Yellow-bellied toad example from Prohl et al 2021
# data application
toad_raw = read_xlsx("data-raw/ProhlEtAl2021_Bombina variegata Microsat data.xlsx")
toad_genotypes <- microsat_genotypes(toad_raw)

Toad_microsat_freq = vector(mode="list",length = 10)
for(i in 1:10){ # for each locus
  tmp = toad_raw[,(2+i*2):(3+i*2)]
  tmp[tmp==-9] = NA
  alleles = sort(na.exclude(unique(unlist(tmp))))
  tmp[[1]] = factor(tmp[[1]],levels=alleles)
  tmp[[2]] = factor(tmp[[2]],levels=alleles)
  colnames(tmp)=c("A1","A2")
  tmp2 = bind_cols(Population=toad_raw[[2]],tmp) %>% pivot_longer(A1:A2)
  Freq.Mat = table(tmp2$Population,tmp2$value)
  Freq.Mat = sweep(Freq.Mat,1,rowSums(Freq.Mat),"/")
  Toad_microsat_freq[[i]] = matrix(Freq.Mat,nrow = nrow(Freq.Mat),dimnames = list(rownames(Freq.Mat),colnames(Freq.Mat)))
}
names(Toad_microsat_freq) = paste0("locus_",colnames(toad_raw)[seq(4,23,2)])

use_data(Toad_microsat_freq, toad_genotypes, overwrite = T)

