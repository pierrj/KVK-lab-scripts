## example profile plot ggplot2 script

library(data.table)
library(ggplot2)

file = file
label = label
output = output

df_read <- t(data.frame(fread(file)))

## need to do a lot of rearranging to get it into a format that works for ggplot
df <- data.frame( bin=as.vector(df_read[3:nrow(df_read),1]),
                  val=c(df_read[3:nrow(df_read),2]),
                  label=c(rep('Label',101)))

df$bin <- as.numeric(as.character(df$bin))
df$val <- as.numeric(as.character(df$val))

# plot gc content across tracks

p <- ggplot(df, aes(x=bin, y=val, group=label)) + geom_line(aes(color=label), lwd=0.5) + theme_classic() +
  scale_x_continuous(labels=c('-0.5kb', 'Breakpoint', '0.5kb'), breaks=c(1,51,101)) + xlab('')+
  theme(axis.text=element_text(size=6), axis.title=element_text(size=8), legend.title=element_blank(), legend.text=element_text(size=8)) + ylab(label)+
  theme(legend.position="none")

ggsave(output, plot=p, width=6.5/3, height=6.5/3)