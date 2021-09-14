reorder_columns<-function(m){

        n=nrow(m)
        k=ncol(m)
        cov=array(-1,k)
        for(j in 1:k){
                cov[j]=sum(m[,j])
        }

        new_order=c()
        while(length(which(cov[]==-1))!=k){
                max=max(cov)
                ids=which(cov[]==max)
                new_order=c(new_order,ids)
                for(i in ids){
                        cov[i]=-1
                }
        }

        new_m=matrix(-1,nrow=n,ncol=k)
        for(j in 1:length(new_order)){
                new_m[,j]=m[,new_order[j]]
        }
        return(new_m)
}



plot_one_matrix<-function(file_cov,tag_order){

        m=read.table(file_cov,skip=1,row.names=1)
        m=as.matrix(m)

        if(tag_order==TRUE){
                a=reorder_columns(m)
                b=reorder_columns(t(a))
                m=t(b)
        }


        nrow_m=nrow(m)
        ncol_m=ncol(m)
        num=2
        colors=colorRampPalette(c("lightgrey","black"))
        for(i in 1:nrow_m){
                if(sum(m[i,])==1){
                        id=which(m[i,]==1)
                        m[i,]=rep(3,ncol_m)
                        m[i,id]=2
                        num=4
                        colors=colorRampPalette(c("lightgrey","black","firebrick","#cfa9a9"))
                }
        }

        #colors=colorRampPalette(c("lightgrey","black","firebrick","#cfa9a9"))
        image(t(m[nrow(m):1,]), col=colors(num),axes = FALSE,xlab="",ylab="")
}


#============================================================================================
#============================================================================================
#   ARGUMENTS
#============================================================================================
#============================================================================================
args = commandArgs(trailingOnly=TRUE)
file_in=args[1]         # path for matrix file
file_out=args[2]        # path for output file

tag_order=FALSE

pdf(paste(file_out,"-",tag_order,".pdf",sep=""),width=10,height=10)
name=paste(file_in)
if(file.exists(name)){
	plot_one_matrix(name,tag_order)
}else{
	print(paste("File not found:",name))
}
dev.off()

