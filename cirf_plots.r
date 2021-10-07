# Codes to create cummulkative impluse response graphs given in Figures 1 to 3 using the .csv files created using TVPcode.R
# use imp16_14c.csv, imp50_14c.csv and imp84_14c.csv for Figure 1,
# use imp16_24c.csv, imp50_24c.csv and imp84_34c.csv for Figure 2 and
# use imp16_34c.csv, imp50_34c.csv and imp84_34c.csv for Figure 3
# Similarly, Figure 4,5,and 6 can be generated using draws from Model 2


install.packages("multipanelfigure")
library("xts")
library("ggplot2")
library("magrittr")
library("multipanelfigure")
library("scales")

dd16 <- read.csv("imp16_34c.csv",header=TRUE)
dd50 <- read.csv("imp50_34c.csv",header=TRUE)
dd84 <- read.csv("imp84_34c.csv",header=TRUE)

Time <- seq(as.Date("2006/7/1"), as.Date("2020/2/1"), "month")

df0 <- data.frame(Time,dd16[,2], dd50[,2],dd84[,2])

colnames(df0) <- c("time","L","M","U")


df3<- data.frame(Time,dd16[,3], dd50[,3],dd84[,3])

colnames(df3) <- c("time","L","M","U")

df6<- data.frame(Time,dd16[,4], dd50[,4],dd84[,4])

colnames(df6) <- c("time","L","M","U")

q0 <- ggplot(df0, aes(x = time, y = M))+geom_line()+
geom_ribbon(data=df0,aes(ymin=L,ymax=U),alpha=0.3)+# adds riddbon
scale_x_date(labels = date_format("%b-%Y"), breaks = date_breaks("3 year"))+
theme(panel.background = element_rect(fill = NA))+ xlab("") + ylab("")+
scale_y_continuous(limits=c(0,1.25))+
theme(legend.position = 'bottom', axis.line.y = element_line(color="black"))+
geom_hline(yintercept = 0)

q3 <- ggplot(df3, aes(x = time, y = M))+geom_line()+
geom_ribbon(data=df3,aes(ymin=L,ymax=U),alpha=0.3) +# adds riddbon
scale_x_date(labels = date_format("%b-%Y"), breaks = date_breaks("3 year"))+
theme(panel.background = element_rect(fill = NA))+ xlab("") + ylab("")+
scale_y_continuous(limits=c(0, 2.5))+
theme(legend.position = 'bottom', axis.line.y = element_line(color="black"))+
geom_hline(yintercept = 0)


q6 <- ggplot(df6, aes(x = time, y = M))+geom_line()+
geom_ribbon(data=df6,aes(ymin=L,ymax=U),alpha=0.3) +# adds riddbon
scale_x_date(labels = date_format("%b-%Y"), breaks = date_breaks("3 year"))+
theme(panel.background = element_rect(fill = NA))+ xlab("") + ylab("")+
scale_y_continuous(limits=c(0, 3.5))+
theme(legend.position = 'bottom', axis.line.y = element_line(color="black"))+
geom_hline(yintercept = 0)


figure1 <- multi_panel_figure(columns = 1, rows =3,panel_label_type = "lower-alpha")

figure1 %<>%
  fill_panel(q0, column = 1, row = 1) %<>%
  fill_panel(q3, column = 1, row = 2) %<>%
  fill_panel(q6, column = 1, row = 3)
figure1
