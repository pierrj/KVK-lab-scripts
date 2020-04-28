library(venneuler)

just1=1572
just2=293
just3=657
oneand2=693
oneand3=1648
twoand3=372
all=2707

v <- venneuler(c(Rep1=just1, Rep2=just2, Rep3=just3, "Rep1&Rep2"=oneand2, "Rep1&Rep3"=oneand3, "Rep2&Rep3"=twoand3, "Rep1&Rep2&Rep3"=all))
v$labels <- rep("", length(v$labels))

par(bg=NA)
plot(v)

dev.copy(png, "venndiagram_all.png")
dev.off()
plot(v)


###EFCs
just1=34
just2=13
just3=8
oneand2=28
oneand3=48
twoand3=9
all=211

v <- venneuler(c(Rep1=just1, Rep2=just2, Rep3=just3, "Rep1&Rep2"=oneand2, "Rep1&Rep3"=oneand3, "Rep2&Rep3"=twoand3, "Rep1&Rep2&Rep3"=all))
v$labels <- rep("", length(v$labels))

par(bg=NA)
plot(v)

dev.copy(png, "venndiagram_efcs.png")
dev.off()
