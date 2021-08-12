## Visualize the sequence conservation of Mnx1 orthologs from multiple species

library(pheatmap)

d = read.table("MNX1.subset.pep.JD.txt")

plot(d[,1]+1, d[,2], type = "h", ylim = c(0,1))


m = matrix(0, ncol = dim(d), nrow = 7)
row.names(m) = c("Human", "Mouse", "Dog", "Chicken", "Xenopus", "Zebrafish", "Fruit fly")
for(i in 1:dim(m)[2]){
  x = d[i,3]
  v = strsplit(as.character(x), "")[[1]]
  for(j in 1:length(v)){
    if(v[j] != "-"){
      m[j,i] = 1
    }
  }
}

png(file="Mnx1_ortholog_conservation.png", width = 6000, height = 1600, res = 600)
par(mar=c(0,6,0,0))
plot(d[,1], d[,2], type="h", lwd = 3, xlim=c(0,610), ylim=c(0,9), col = "grey", 
     yaxt = "n", xlab = "", ylab = "", bty = "n", axes = FALSE)
axis(side = 2, at = c(8:2, 0.2), tick = FALSE, 
     labels = c("Human", "Mouse", "Dog", "Chicken", "Xenopus", "Zebrafish", "Fruit fly", "Conservation"), las = 2)
for(i in 1:dim(m)[1]){
  for(j in 1:dim(m)[2]){
    if(m[i,j] == 1){
      w = 0.4
      col = "lightblue"
    }
    else{
      w = 0.03
      col = "grey"
    }
    rect(xleft = j, xright = j+1, ybottom = 9 - i + w, ytop = 9 - i - w, col = col, border = NA)
  }
}

dev.off()
