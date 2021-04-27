rawdata <- read.table("Plot/allgenes_vst.txt",header=T,sep="\t")


dat11 <-as.matrix(rawdata[-1])
dat11 <-t(dat11)
dat11<- scale(dat11, center = TRUE, scale = TRUE)
dat11<-t(dat11)
row.names(dat11)<-rawdata$ID 
write.csv(dat11, "z-scored.csv")

TRI <- read.table("Plot/TRI5_zscore.txt",header=T,sep="\t")


C1 <- read.table("Plot/reshaped_C1.txt",header=T,sep="\t")

# Group samples by time
grouped <- group_by(C1, ID,Timepoint)
# Calculate statistics for each group
stats <- summarise(grouped, N=length(zscore), Average=mean(zscore), StDev=sd(zscore))
# Create a function that calculates 95% confidence intervals for the given data vector using a t-distribution
conf_int95 <- function(C1) {
    n <- length(C1)
    error <- qt(0.975, df=n-1) * sd(C1)/sqrt(n)
    return(error)
}
# Create summary for each group containing sample size, average OD600, and 95% confidence limits
stats <- summarise(grouped, N=length(zscore), Average=mean(zscore), CI95=conf_int95(zscore))
# Plot results 
ggplot(data=stats, aes(x=Timepoint, y=Average, group=ID, color=ID)) + 
geom_ribbon(aes(ymin=Average-CI95, ymax=Average+CI95, fill=ID),
color=NA, alpha=0.3) +
geom_line() +
facet_grid(~ID) +
#scale_fill_manual(values=c('#228B22','#E69F00','#22228b','#8b8b22','#987095','#959870','#987370','#987095','#709873','#709873','#8e7fa9','#7f9aa9')) +
#scale_color_manual(values=c('#228B22','#E69F00','#22228b','#8b8b22','#987095','#959870','#987370','#987095','#709873','#709873','#8e7fa9','#7f9aa9')) +
labs(x="Time (Hours)", y="zscores") +
theme(panel.background = element_rect(fill = NA),panel.grid.major = element_line(colour = "grey90"),panel.ontop = F)



C2 <- read.table("Plot/reshaped_C2.txt",header=T,sep="\t")

# Group samples by time
grouped <- group_by(C2, ID,Timepoint)
# Calculate statistics for each group
stats <- summarise(grouped, N=length(zscore), Average=mean(zscore), StDev=sd(zscore))
# Create a function that calculates 95% confidence intervals for the given data vector using a t-distribution
conf_int95 <- function(C2) {
    n <- length(C2)
    error <- qt(0.975, df=n-1) * sd(C2)/sqrt(n)
    return(error)
}
# Create summary for each group containing sample size, average OD600, and 95% confidence limits
stats <- summarise(grouped, N=length(zscore), Average=mean(zscore), CI95=conf_int95(zscore))
# Plot results 
ggplot(data=stats, aes(x=Timepoint, y=Average, group=ID, color=ID)) + 
geom_ribbon(aes(ymin=Average-CI95, ymax=Average+CI95, fill=ID),
color=NA, alpha=0.3) +
geom_line() +
facet_grid(~ID) +
#scale_fill_manual(values=c('#228B22','#E69F00','#22228b','#8b8b22','#987095','#959870','#987370','#987095','#709873','#709873','#8e7fa9','#7f9aa9')) +
#scale_color_manual(values=c('#228B22','#E69F00','#22228b','#8b8b22','#987095','#959870','#987370','#987095','#709873','#709873','#8e7fa9','#7f9aa9')) +
labs(x="Time (Hours)", y="zscores") +
theme(panel.background = element_rect(fill = NA),panel.grid.major = element_line(colour = "grey90"),panel.ontop = F)

C3 <- read.table("Plot/reshaped_C3.txt",header=T,sep="\t")

# Group samples by time
grouped <- group_by(C3, ID,Timepoint)
# Calculate statistics for each group
stats <- summarise(grouped, N=length(zscore), Average=mean(zscore), StDev=sd(zscore))
# Create a function that calculates 95% confidence intervals for the given data vector using a t-distribution
conf_int95 <- function(C3) {
    n <- length(C3)
    error <- qt(0.975, df=n-1) * sd(C3)/sqrt(n)
    return(error)
}
# Create summary for each group containing sample size, average OD600, and 95% confidence limits
stats <- summarise(grouped, N=length(zscore), Average=mean(zscore), CI95=conf_int95(zscore))
# Plot results 
ggplot(data=stats, aes(x=Timepoint, y=Average, group=ID, color=ID)) + 
geom_ribbon(aes(ymin=Average-CI95, ymax=Average+CI95, fill=ID),
color=NA, alpha=0.3) +
geom_line() +
facet_grid(~ID) +
#scale_fill_manual(values=c('#228B22','#E69F00','#22228b','#8b8b22','#987095','#959870','#987370','#987095','#709873','#709873','#8e7fa9','#7f9aa9')) +
#scale_color_manual(values=c('#228B22','#E69F00','#22228b','#8b8b22','#987095','#959870','#987370','#987095','#709873','#709873','#8e7fa9','#7f9aa9')) +
labs(x="Time (Hours)", y="zscores") +
theme(panel.background = element_rect(fill = NA),panel.grid.major = element_line(colour = "grey90"),panel.ontop = F)

C4 <- read.table("Plot/reshaped_C4.txt",header=T,sep="\t")

# Group samples by time
grouped <- group_by(C4, ID,Timepoint)
# Calculate statistics for each group
stats <- summarise(grouped, N=length(zscore), Average=mean(zscore), StDev=sd(zscore))
# Create a function that calculates 95% confidence intervals for the given data vector using a t-distribution
conf_int95 <- function(C4) {
    n <- length(C4)
    error <- qt(0.975, df=n-1) * sd(C4)/sqrt(n)
    return(error)
}
# Create summary for each group containing sample size, average OD600, and 95% confidence limits
stats <- summarise(grouped, N=length(zscore), Average=mean(zscore), CI95=conf_int95(zscore))
# Plot results 
ggplot(data=stats, aes(x=Timepoint, y=Average, group=ID, color=ID)) + 
geom_ribbon(aes(ymin=Average-CI95, ymax=Average+CI95, fill=ID),
color=NA, alpha=0.3) +
geom_line() +
facet_grid(~ID) +
#scale_fill_manual(values=c('#228B22','#E69F00','#22228b','#8b8b22','#987095','#959870','#987370','#987095','#709873','#709873','#8e7fa9','#7f9aa9')) +
#scale_color_manual(values=c('#228B22','#E69F00','#22228b','#8b8b22','#987095','#959870','#987370','#987095','#709873','#709873','#8e7fa9','#7f9aa9')) +
labs(x="Time (Hours)", y="zscores") +
theme(panel.background = element_rect(fill = NA),panel.grid.major = element_line(colour = "grey90"),panel.ontop = F)




rawdata2 <- read.table("TRI5_vst.txt",header=T,sep="\t")

reshaped2 <- melt(rawdata2, id=c("ID"), variable.name="Timepoint", value.name="vst")
# Write csv file
#annotated <- inner_join(dd, platemap, by="Well")
write.csv(reshaped2, "TRI5_vst_reshaped.csv")


C1 <- read.table("vst/Z1.txt",header=T,sep="\t")

# Group samples by time
grouped <- group_by(C1, ID,Timepoint)
# Calculate statistics for each group
stats <- summarise(grouped, N=length(vst), Average=mean(vst), StDev=sd(vst))
# Create a function that calculates 95% confidence intervals for the given data vector using a t-distribution
conf_int95 <- function(C1) {
    n <- length(C1)
    error <- qt(0.975, df=n-1) * sd(C1)/sqrt(n)
    return(error)
}
# Create summary for each group containing sample size, average OD600, and 95% confidence limits
stats <- summarise(grouped, N=length(vst), Average=mean(vst), CI95=conf_int95(vst))
# Plot results 
ggplot(data=stats, aes(x=Timepoint, y=Average, group=ID, color=ID)) + 
geom_ribbon(aes(ymin=Average-CI95, ymax=Average+CI95, fill=ID),
color=NA, alpha=0.3) +
geom_line() +
#facet_grid(~ID) +
#scale_fill_manual(values=c('#228B22','#E69F00','#22228b','#8b8b22','#987095','#959870','#987370','#987095','#709873','#709873','#8e7fa9','#7f9aa9')) +
#scale_color_manual(values=c('#228B22','#E69F00','#22228b','#8b8b22','#987095','#959870','#987370','#987095','#709873','#709873','#8e7fa9','#7f9aa9')) +
labs(x="Time (Hours)", y="vst") +
theme(panel.background = element_rect(fill = NA),panel.grid.major = element_line(colour = "grey90"),panel.ontop = F)

ggplot(C1, aes(Timepoint,vst, group=ID, color=ID))+
geom_line(alpha=0.01)+
stat_summary(aes(group=ID),
fun=mean, geom="line", size=1)+
#facet_grid(Condition) +
xlab("Timepoints")+
ylab("vst")+
theme(panel.background = element_rect(fill = NA),panel.grid.major = element_line(colour = "grey90"),panel.ontop = F)

C2 <- read.table("vst/Z2.txt",header=T,sep="\t")

# Group samples by time
grouped <- group_by(C2, ID,Timepoint)
# Calculate statistics for each group
stats <- summarise(grouped, N=length(vst), Average=mean(vst), StDev=sd(vst))
# Create a function that calculates 95% confidence intervals for the given data vector using a t-distribution
conf_int95 <- function(C2) {
    n <- length(C2)
    error <- qt(0.975, df=n-1) * sd(C2)/sqrt(n)
    return(error)
}
# Create summary for each group containing sample size, average OD600, and 95% confidence limits
stats <- summarise(grouped, N=length(vst), Average=mean(vst), CI95=conf_int95(vst))
# Plot results 
ggplot(data=stats, aes(x=Timepoint, y=Average, group=ID, color=ID)) + 
geom_ribbon(aes(ymin=Average-CI95, ymax=Average+CI95, fill=ID),
color=NA, alpha=0.3) +
geom_line() +
#facet_grid(~ID) +
#scale_fill_manual(values=c('#228B22','#E69F00','#22228b','#8b8b22','#987095','#959870','#987370','#987095','#709873','#709873','#8e7fa9','#7f9aa9')) +
#scale_color_manual(values=c('#228B22','#E69F00','#22228b','#8b8b22','#987095','#959870','#987370','#987095','#709873','#709873','#8e7fa9','#7f9aa9')) +
labs(x="Time (Hours)", y="vst") +
theme(panel.background = element_rect(fill = NA),panel.grid.major = element_line(colour = "grey90"),panel.ontop = F)

ggplot(C2, aes(Timepoint,vst, group=ID, color=ID))+
geom_line(alpha=0.01)+
stat_summary(aes(group=ID),
fun=mean, geom="line", size=1)+
#facet_grid(Condition) +
xlab("Timepoints")+
ylab("Z-score")+
theme(panel.background = element_rect(fill = NA),panel.grid.major = element_line(colour = "grey90"),panel.ontop = F)

C3 <- read.table("vst/Z3.txt",header=T,sep="\t")

# Group samples by time
grouped <- group_by(C3, ID,Timepoint)
# Calculate statistics for each group
stats <- summarise(grouped, N=length(vst), Average=mean(vst), StDev=sd(vst))
# Create a function that calculates 95% confidence intervals for the given data vector using a t-distribution
conf_int95 <- function(C3) {
    n <- length(C3)
    error <- qt(0.975, df=n-1) * sd(C3)/sqrt(n)
    return(error)
}
# Create summary for each group containing sample size, average OD600, and 95% confidence limits
stats <- summarise(grouped, N=length(vst), Average=mean(vst), CI95=conf_int95(vst))
# Plot results 
ggplot(data=stats, aes(x=Timepoint, y=Average, group=ID, color=ID)) + 
geom_ribbon(aes(ymin=Average-CI95, ymax=Average+CI95, fill=ID),
color=NA, alpha=0.3) +
geom_line() +
#facet_grid(~ID) +
#scale_fill_manual(values=c('#228B22','#E69F00','#22228b','#8b8b22','#987095','#959870','#987370','#987095','#709873','#709873','#8e7fa9','#7f9aa9')) +
#scale_color_manual(values=c('#228B22','#E69F00','#22228b','#8b8b22','#987095','#959870','#987370','#987095','#709873','#709873','#8e7fa9','#7f9aa9')) +
labs(x="Time (Hours)", y="vst") +
theme(panel.background = element_rect(fill = NA),panel.grid.major = element_line(colour = "grey90"),panel.ontop = F)

ggplot(C3, aes(Timepoint,vst, group=ID, color=ID))+
geom_line(alpha=0.01)+
stat_summary(aes(group=ID),
fun=mean, geom="line", size=1)+
#facet_grid(Condition) +
xlab("Timepoints")+
ylab("Z-score")+
theme(panel.background = element_rect(fill = NA),panel.grid.major = element_line(colour = "grey90"),panel.ontop = F)


C4 <- read.table("vst/Z4.txt",header=T,sep="\t")

# Group samples by time
grouped <- group_by(C4, ID,Timepoint)
# Calculate statistics for each group
stats <- summarise(grouped, N=length(vst), Average=mean(vst), StDev=sd(vst))
# Create a function that calculates 95% confidence intervals for the given data vector using a t-distribution
conf_int95 <- function(C4) {
    n <- length(C4)
    error <- qt(0.975, df=n-1) * sd(C4)/sqrt(n)
    return(error)
}
# Create summary for each group containing sample size, average OD600, and 95% confidence limits
stats <- summarise(grouped, N=length(vst), Average=mean(vst), CI95=conf_int95(vst))
# Plot results 
ggplot(data=stats, aes(x=Timepoint, y=Average, group=ID, color=ID)) + 
geom_ribbon(aes(ymin=Average-CI95, ymax=Average+CI95, fill=ID),
color=NA, alpha=0.3) +
geom_line() +
facet_grid(~ID) +
#scale_fill_manual(values=c('#228B22','#E69F00','#22228b','#8b8b22','#987095','#959870','#987370','#987095','#709873','#709873','#8e7fa9','#7f9aa9')) +
#scale_color_manual(values=c('#228B22','#E69F00','#22228b','#8b8b22','#987095','#959870','#987370','#987095','#709873','#709873','#8e7fa9','#7f9aa9')) +
labs(x="Time (Hours)", y="vst") +
theme(panel.background = element_rect(fill = NA),panel.grid.major = element_line(colour = "grey90"),panel.ontop = F)

ggplot(C4, aes(Timepoint,vst, group=ID, color=ID))+
geom_line(alpha=0.01)+
stat_summary(aes(group=ID),
fun=mean, geom="line", size=1)+
#facet_grid(Condition) +
xlab("Timepoints")+
ylab("Z-score")+
theme(panel.background = element_rect(fill = NA),panel.grid.major = element_line(colour = "grey90"),panel.ontop = F)



Z5 <- read.table("vst/Z5.txt",header=T,sep="\t")

ggplot(Z5, aes(Timepoint,vst, group=ID, color=ID))+
geom_line(alpha=0.01)+
stat_summary(aes(group=ID),
fun=mean, geom="line", size=0.5)+
facet_grid(~Condition) +
xlab("Timepoints")+
ylab("vst")+
theme_bw()+
  theme(axis.text= element_text(colour="black", size=7),
        axis.title = element_text(colour = "black", size=12),
        aspect.ratio = 1, legend.title = element_blank())


Z6 <- read.table("Plot/TRI5_zscore.txt",header=T,sep="\t")

ggplot(TRI2, aes(Timepoint,zscore, group=ID, color=ID))+
geom_line(alpha=0.01)+
stat_summary(aes(group=ID),
fun=mean, geom="line", size=1)+
facet_grid(~Condition) +
xlab("Timepoints")+
ylab("vst")+
theme(panel.background = element_rect(fill = NA),panel.grid.major = element_line(colour = "grey90"),panel.ontop = F)

ggplot(TRI2, aes(Timepoint,zscore, group=ID, color=ID))+
geom_line(alpha=0.01)+
stat_summary(aes(group=ID),
fun=mean, geom="line", size=0.5)+
facet_grid(~Condition) +
xlab("Timepoints")+
ylab("z-scores")+
theme_bw()+
  theme(axis.text= element_text(colour="black", size=7),
        axis.title = element_text(colour = "black", size=12),
        aspect.ratio = 1, legend.title = element_blank())

-------------------
rawTRI5 <- read.table("vst_all/vst_all_TRI5.txt",header=T,sep="\t")
res <- melt(rawTRI5, id=c("ID"), variable.name="Timepoints", value.name="vst")
write.table(res, "vst_all/TRI4plot.txt", sep="\t")
L1 <- read.table("vst_all/TRI4plot.txt",header=T,sep="\t")

ggplot(L1, aes(Timepoint,vst, group=ID, color=ID))+
geom_line(alpha=0.01)+
stat_summary(aes(group=ID),
fun=mean, geom="line", size=0.5)+
facet_grid(~Condition) +
xlab("Timepoints")+
ylab("vst")+
theme_bw()+
  theme(axis.text= element_text(colour="black", size=7),
        axis.title = element_text(colour = "black", size=12),
        aspect.ratio = 1, legend.title = element_blank())




dat11 <-as.matrix(rawdata[-1])
dat11 <-t(dat11)
dat11<- scale(dat11, center = TRUE, scale = TRUE)
dat11<-t(dat11)
row.names(dat11)<-rawdata$ID 
write.csv(dat11, "z-scored.csv")

TRI <- read.table("Plot/TRI5_zscore.txt",header=T,sep="\t")



