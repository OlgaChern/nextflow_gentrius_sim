#!/usr/local/bin/R
args = commandArgs(trailingOnly=TRUE)

file_in=args[1]
splits_remained=args[2]
splits_ignored=args[3]

library(ape)
library(igraph)

file_out=paste(file_in,"_plot_consensus_tree.pdf",sep="")
pdf(file_out,5,5)
	par(mar=c(0.5,0.5,0.5,0.5),oma=c(0,0,0,0))

	tree<-read.tree(file_in)
	graph=ape::as.igraph.phylo(tree,directed = FALSE,use.labels = FALSE)

	V(graph)$size=ifelse(strength(graph, vids = V(graph), mode = "all")>3,5,0.5)
	V(graph)$color=ifelse(strength(graph, vids = V(graph), mode = "all")>3,"#D93951","#00B7A1")
	
	ids=which(strength(graph, vids = V(graph), mode = "all")>3)
	multi_v=length(ids)
	int_ids=list()
	ext_ids=list()
	int_NUM=0
	ext_NUM=0
	for(v in ids){
		nei_ids=neighbors(graph, V(graph)[v], mode ="all")
                check=length(which(strength(graph, vids = nei_ids, mode = "all")!=1))
                if(check<2){
			ext_NUM=ext_NUM+1
			ext_ids[[ext_NUM]]=v
                }else{  
                        int_NUM=int_NUM+1
			int_ids[[int_NUM]]=v
                }       
	}      

	ext_ids=array(as.numeric(ext_ids))
	int_ids=array(as.numeric(int_ids)) 

	if(length(ext_ids)>0){
		f_cols=colorRampPalette(c("#b4bfe0","royalblue"))
		ext_degrees=degree(graph,v=ext_ids,mode="all")
		ext_cols=f_cols(max(ext_degrees))
		V(graph)[ext_ids]$color=ext_cols[ext_degrees]
	}

	if(length(int_ids)>0){
        	f_cols=colorRampPalette(c("#bf8f8f","firebrick"))
                int_degrees=degree(graph,v=int_ids,mode="all")
                int_cols=f_cols(max(int_degrees))
                V(graph)[int_ids]$color=int_cols[int_degrees]
	}

        plot.igraph(graph,vertex.label=NA,xlim=c(-1.5,1.5),ylim=c(-1.5,1.5)) 
	x=-1.4
	y=1.4
	e=0.1
	sz=0.8
	text(x,y,paste("No.int.branches: ",splits_remained,sep=""),cex=sz,adj=0)
	text(x,y-e,paste("No.mult.nodes: ",multi_v,sep=""),cex=sz,adj=0)

	x=-0.4
	col_nd=ifelse(ext_NUM>0,"royalblue","lightgrey")
	text(x,y,paste("No.mult.ext.nodes: ",ext_NUM,sep=""),cex=sz,adj=0,col=col_nd)
	text(x,y-e,paste("max degree: ",ext_max_degree,sep=""),cex=sz,adj=0,col=col_nd)
	x=0.7
	col_nd=ifelse(int_NUM>0,"firebrick","lightgrey")
	text(x,y,paste("No.mult.int.nodes: ",int_NUM,sep=""),cex=sz,adj=0,col=col_nd)
        text(x,y-e,paste("max degree: ",int_max_degree,sep=""),cex=sz,adj=0,col=col_nd)

	#text_info=paste(r[i,1],"REMAINED",r[i,2]+n,"REMAINED_INT",r[i,2],"DISCARDED",r[i,3],"T_SIZE",9,"NUM_MULTI_NODES",multi_v,"EXTERNAL_MULTI",ext_NUM,"INTERNAL_MULTI",int_NUM,"MULTI_MAX",max(ext_max_degree,int_max_degree),"EXT_MAX",ext_max_degree,"EXT_AVE",ext_ave_degree,"INT_MAX",int_max_degree,"INT_AVE",int_ave_degree,"| EXT_DISTR",paste(ext_distr,collapse=" "),"INT_DISTR",paste(int_distr,collapse=" "),"|",sep=" ")

	#write(text_info,file="output.txt",append=TRUE)
dev.off()

