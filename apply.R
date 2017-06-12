# Apply for AMD
# Assuming
# 150 case, 150 control


effect <- seq(0.1, 2, 0.05)

out <- array(0, c(length(effect), 1))
for(i in 1:length(effect)) {
    tmp2 <- fold_change_power.fun(n1=150, n2=150, sig.level=2.5e-06, d=effect[i])
    out[i,1] <- tmp2$power
}

out <- as.data.frame(out)
out <- cbind(effect, out)
names(out)[1:2] <- c("fold_change", "n1=150_&_n2=150")

mdata <- melt(out, id=c("fold_change"))
names(mdata)[2] <- "Sample_size"
names(mdata)[3] <- "Power"
mdata$fold_change <- mdata$fold_change+1

ggplot(data=mdata, aes(x=fold_change, y=Power, group=Sample_size, colour=Sample_size)) + geom_line() + geom_point()
