#genemina
library(tidyverse)
library(data.table)

gen.fun <- fread("genemania-functions.txt")
gen.net <- fread("genemania-networks.txt")
gen.par <- fread("genemania-parameters (1).txt")
