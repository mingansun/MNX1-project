## Visualize the phastCons curve for different groups of Mnx1 peaks

library(RColorBrewer)
library(ggplot2)

PR = read.table("Mnx1.PR.5kb.cns")[,-1:-6]
IN = read.table("Mnx1.IN.5kb.cns")[,-1:-6]
IG = read.table("Mnx1.IG.5kb.cns")[,-1:-6]

## get averaged profile
PR.ave = apply(PR,2,mean)
IN.ave = apply(IN,2,mean)
IG.ave = apply(IG,2,mean)

## plot using ggplot2
df.cns = data.frame(
  PhastCons = c(PR.ave, IN.ave, IG.ave),
  Distribution = factor(rep(c("Promoter", "Intron", "Intergenic"), each = 10001), levels = c("Promoter", "Intron", "Intergenic")),
  Distance = rep(-5000:5000, 3)
)

## Generate curve plot
png(file="Mnx1_phastCons_vs_distr.curve.png", width = 2400, height = 2000, res = 600)
ggplot(df.cns, aes(x=Distance, y=PhastCons, color = Distribution)) +
  theme_bw() + theme(legend.position=c(0.85,0.8)) +
  xlab("Distance to Mnx1 peak center (bp)") + ylab("PhastCons Score") +
  geom_line(size = I(0.5))
dev.off()






