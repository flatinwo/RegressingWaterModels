library(ggplot2)
mydata <- read.table('RescaledTSEOS.dat',header=TRUE)
theme_set(theme_gray(base_size = 24))
#tiff("AllDataPhaseDiagram.tiff",width=12,height=8,units='in',res=300)
pdf("RescaledTSEOS.pdf")#,width=12,height=8)#,units='in',res=300)
data2 <- subset(mydata, Source %in% c("TSEOS"))
data2a <- subset(data2, Property %in% c("Tmax","Tmin"))
data2b <-subset(data2, !(Property %in% c("Tmax","Tmin")))
data3 <- subset(data2a, Pressure > 400) # to avoid crossing of lines
data4 <- subset(data2a, Pressure <= 400) # to avoid crossing of lines
p <- ggplot(subset(mydata, Source %in% c("Madrid","Evan")),aes(Temperature,Pressure,color=factor(Property),shape=factor(Source)))
p <- p + geom_point(size=4) + scale_color_discrete(parse(text=paste("Property")))
p <- p + ylab("Rescaled Pressure (MPa)")+xlab("Rescaled Temperature (K)") + scale_shape_discrete("Source") 
p <- p + geom_line(data=data3,aes(x=Temperature,y=Pressure)) + geom_line(data=data4,aes(x=Temperature,y=Pressure))
p <- p + geom_line(data=data2b,aes(x=Temperature, y=Pressure))
p <- p + guides(color = guide_legend(override.aes = list(size=4,linetype=0))) 
p <- p + guides(shape = guide_legend(override.aes=list(linetype=c(NA,NA,1),shape=c(16,17,NA),size=4)))
p
dev.off()
