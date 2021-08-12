## Script used to analyze and visualize differential binding

library(RColorBrewer)
library(DiffBind)
library(scales)

#################################### WT vs. KO  ###################################
## read in the peak sets
k27.wt2hd <- dba(sampleSheet = "K27ac_WT2KO.csv")

## count reads
k27.wt2hd <- dba.count(k27.wt2hd, minOverlap = 1, summits = 500)

## establish a contrast
k27.wt2hd <- dba.contrast(k27.wt2hd, categories = DBA_CONDITION, minMembers = 2)

## perform differential analysis
k27.wt2hd <- dba.analyze(k27.wt2hd)

## retrieve the differentially bound sites
k27.wt2hd.all <- dba.report(k27.wt2hd, contrast = 1, th=1, bFlip = TRUE)
k27.wt2hd.DB  <- k27.wt2hd.all[abs(db.all$Fold) >=  1 & db.all$FDR < 0.1,]
write.table(
  k27.wt2hd.DB,  file = "K27.wt2hd_DB.txt", 
  sep = "\t", quote = FALSE, row.names = FALSE
  )
write.table(
  k27.wt2hd.all, file = "K27.wt2hd_all.txt", 
  sep = "\t", quote = FALSE, row.names = FALSE
  )

# MA plots
# Mnx1 binding status on pooled K27ac peaks
mnx1.cnt <- read.table("K27.wt2hd_all.Hb9.cnt")

#
k27.wt2hd.up <- k27.wt2hd.all[k27.wt2hd.all$Fold >=  1 & k27.wt2hd.all$FDR < 0.1,]
k27.wt2hd.dn <- k27.wt2hd.all[k27.wt2hd.all$Fold <= -1 & k27.wt2hd.all$FDR < 0.1,]
k27.wt2hd.up.bd <- k27.wt2hd.all[
  k27.wt2hd.all$Fold >=  1 & k27.wt2hd.all$FDR < 0.1 & mnx1.cnt[,4] > 0,
  ]
k27.wt2hd.dn.bd <- k27.wt2hd.all[
  k27.wt2hd.all$Fold <= -1 & k27.wt2hd.all$FDR < 0.1 & mnx1.cnt[,4] > 0,
  ]

png(file="K27ac_HD2WT.MAplot.png", width = 3500, height = 2700, res = 600)
par(mfrow = c(1,1), mar = c(4,4,1,1))
plot(
  k27.wt2hd.all$Conc, k27.wt2hd.all$Fold, cex = .1, pch = 20, 
  xlab = "log2(Averaged H3K27ac signal)", ylab = "log2(Mnx1-Mut/WT)"
  )
abline(h = 0, lwd = 3, col = alpha("red", 0.2))
points(
  k27.wt2hd.up$Conc, k27.wt2hd.up$Fold, 
  col = "orange", cex = 0.25, pch = 20
  )
points(
  k27.wt2hd.dn$Conc, k27.wt2hd.dn$Fold, 
  col = "orange", cex = 0.25, pch = 20
  )
points(
  k27.wt2hd.up.bd$Conc, k27.wt2hd.up.bd$Fold, 
  col = "red", cex = 0.5, pch = 20
  )
points(
  k27.wt2hd.dn.bd$Conc, k27.wt2hd.dn.bd$Fold, 
  col = "red", cex = 0.5, pch = 20
  )
legend(
  "bottomright", bty = "n", cex = .8,
  legend = c("MNX1-unbound diff H3K27ac peak", "MNX1-bound diff H3K27ac peak"),
  fill = c("orange", "red")
  )
dev.off()

### compare DB_K27ac vs. Hb9

## vs. Hb9
# Up K27ac pks:  14 out of 40 (35%) are bound by Hb9
# Dn K27ac pks:  13 out of 110 (11.8%) are bound by Hb9
# All K27ac pks: 4466 out of 22440 (19.9%) are bound by Hb9
# FET: p = 0.003

db2hb9 <- matrix(c(14, 40-14, 13, 110-13), nrow = 2)
fisher.test(db2hb9)

c2 <- brewer.pal(3, "Set1")[c(1,3)]

png(file="dbK27ac_vs_Hb9_barplot.png", width = 2200, height = 3200, res = 600)
par(mar=c(4, 4, 1, 1), mfrow = c(1, 1))
barplot(
  t(t(db2hb9) / apply(db2hb9, 2, sum)) * 100, main = "",
  names = c("Mnx1-Mut > WT", "WT > Mnx1-Mut"), 
  ylab = "% of H3K27ac peaks", 
  col = c2, ylim = c(0, 115)
  )
abline(h = 19.9, lwd = 2, lty = 2, col = "grey")
text(x = 0.7, y = 100, labels = "MNX1-", pos = 1, col = 0)
text(x = 0.7, y = 0, labels = "MNX1+", pos = 3, col = 0)
text(x = 1.9, y = 100, labels = "MNX1-", pos = 1, col = 0)
text(x = 1.9, y = 0, labels = "MNX1+", pos = 3, col = 0)

lines(c(0.7, 0.7, 1.9, 1.9), c(102, 105, 105, 102))
text(x = 1.3, y = 105, labels = "p = 0.003", pos = 3)

dev.off()

