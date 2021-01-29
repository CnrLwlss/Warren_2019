library(deldir)

cdat = data.frame(x=td$xCoord,y=td$yCoord)
meanrad = mean(sqrt(td$Area/(2*pi)))

del = deldir(cdat)

del$delsgs$dist =  sqrt((del$delsgs$x2-del$delsgs$x1)^2+(del$delsgs$y2-del$delsgs$y1)^2)

adj = matrix(0.0,nrow=nrow(cdat),ncol=nrow(cdat)) 
for(i in 1:length(del$delsgs$dist)){
  distance = del$delsgs$dist[i]
  if(distance>6*meanrad) distance = 0.0
  adj[del$delsgs$ind1[i],del$delsgs$ind2[i]] = distance
  adj[del$delsgs$ind2[i],del$delsgs$ind1[i]] = distance
}

adj = round(adj)

colnames(adj) = td$cell_id
rownames(adj) = td$cell_id

candidates = sample(td$cell_id)
unavailable = c()
sampled = c()

pdf("Samples.pdf",width=16,height=8)
while(length(candidates)>0){
 cellid = candidates[1]
 neighb_1 = td$cell_id[adj[,cellid]>0.0]
 neighb_2 = c()
 for(n in neighb_1){
   neighb_2 = unique(c(neighb_2,td$cell_id[adj[,n]>0.0]))
 }
 neighb_2 = neighb_2[neighb_2!=cellid]
 sampled = c(sampled, cellid)
 unavailable = unique(c(unavailable,neighb_1))
 #unavailable = unique(c(unavailable,neighb_1,neighb_2))
 candidates = candidates[!((candidates%in%unavailable)|(candidates%in%sampled))]
}

sdat=cdat[td$cell_id%in%sampled,]
udat=cdat[td$cell_id%in%unavailable,]
op = par(mfrow=c(1,2))
 prop = length(sampled)/length(cdat$x)
 mlab = paste("Proportion sampled:",formatC(100*prop,3),"(%)")
 rang = range(c(cdat$x,cdat$y),xlim=rang,ylim=rang)
 plot(del)
 points(sdat$x,sdat$y,pch=16,col="red",cex=0.5)
 points(udat$x,udat$y,pch=16,col="blue",cex=0.5)
 plot(cdat$x,cdat$y,pch=16,col="grey",xlab="x (px)",ylab="y (px)",cex.axis=1.45,cex.lab=1.45,xlim=rang,ylim=rang,main=mlab)
 points(sdat$x,sdat$y,pch=16,col="red",cex=0.5)
 points(udat$x,udat$y,pch=16,col="blue",cex=0.5)
par(op)

dev.off()
