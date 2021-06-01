cd analysis/dynGENIE3

```bash

/projects/software/R-3.6.1/bin/R
```

```R
setwd("/projects/fusarium_venenatum/analysis/dynGENIE3")

library ("reshape2")
library ("doRNG")
library ("doParallel")

source ("dynGENIE3.R")

S1<- read.table("C1.txt")
S2<- read.table("C2.txt")
S3<- read.table("C3.txt")
S4<- read.table("C4.txt")


TS1 <- read.expr.matrix("C1_v2.txt",form="rows.are.samples")
TS2 <- read.expr.matrix("C2_v2.txt",form="rows.are.samples")
TS3 <- read.expr.matrix("C3_v2.txt",form="rows.are.samples")
TS4 <- read.expr.matrix("C4_v2.txt",form="rows.are.samples")

time.points <- list(TS1[1,], TS2[1,])
TS.data <- list(TS1[12000:nrow(TS1),], TS2[12000:nrow(TS2),])


time.points <- list(TS1[1,], TS2[1,], TS3[1,], TS4[1,])
TS.data <- list(TS1[2:nrow(TS1),], TS2[2:nrow(TS2),], TS3[2:nrow(TS3),], TS4[2:nrow(TS4),])



# Run dynGENIE3
res <- dynGENIE3(TS.data,time.points)