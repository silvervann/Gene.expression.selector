#######HMAP CLUSTERING #################################################

# Generate column annotations
		Condition= as.character(cat.design)#resultsGrp#
		annotations= data.frame(Condition=cat.design)

		# Specify colors
		
		# The palette with grey:
		cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


		Conditioncol = cbPalette[1:length(unique(Condition))]
		names(Condition) = unique(Condition)#paste(rep("Cluster",nclust),c(1:nclust),sep="_")
		# Groupcol = cbPalette[c(1,8,3,7)]
		# names(Group) = unique(Group)
		# Toxicity_2col = c("gold", "darkgreen")
		# names(Toxicity_2) = unique(Toxicity_2)

		# ann_colors = list(Cluster=Clustercol, Group=Groupcol,Toxicity_2=Toxicity_2col)
		ann_colors = list(Condition=Conditioncol)

		#heatimage_cor <-paste("Sample_Correlation_BG_corrected",".png", sep="_")
			#pdf(file= paste0("./TEST/",heatimage_2[k]))
		png("qNorm_FDR.heatmap_10_samples.png", units="in", width=11, height=8.5, res=300)
		
		mat2amap = mat2[ ,unlist(lapply(c("pos","neg"), function(x) grep(x, colnames(mat2)))) ]
		amap=aheatmap(mat2, annCol = annotations,
			Rowv = TRUE, Colv =TRUE ,
			revC=FALSE,
			distfun = "maximum",
			hclustfun = "complete",
			color = "-RdBu:200",
			scale = "row", ,
			annColors = ann_colors,
			fontsize = 10,
			labRow = NULL, labCol = NULL,
			cexRow=3,cexCol=1)
		dev.off()




##############################################################################
# Plot differentially expressed genes
TTFch_sel$Gene = factor( TTFch_sel$Gene, levels= TTFch_sel$Gene)

  bplotFCh = ggplot(TTFch_sel, aes(x=Gene, y=DDCT, fill=Significant)) +
  theme_minimal() +
  scale_fill_manual(values= c("blue4", "red4") ) +   
    geom_bar(stat="identity", color="black",size=0.3) +
    ylab("DDCt")+ 
    xlab("Gene")+
    theme(panel.background = element_rect(fill = NA),
      axis.title.x = element_text(color="black", size=15, face="bold"),
      axis.title.y = element_text(color="black", size=15, face="bold"),
      axis.text.x = element_text(face="bold", color="black",size=16,angle=90,hjust =0, vjust =1),
      axis.text.y = element_text(face="bold", color="black", size=15,angle=0))

png(file= paste0("-DDCT_Barchart_10_samples.png"), units="in", width=10, height=8, res=200)
print(bplotFCh)
dev.off()

#######ENRICHER#################################################

library(enrichR)
dbsall <- listEnrichrDbs()
dbs <- c("GO_Biological_Process_2018", "WikiPathways_2019_Human", "KEGG_2019_Human")

# gvector = RankedList$GseaScore
# names(gvector) = RankedList$Gene
list2plot =NULL
enrichIn = TTFch_sel
for (i in 1:length(unique(enrichIn$Significant)))
{
	temptenrichIn = enrichIn[enrichIn$Significant == unique(enrichIn$Significant)[i],]	
	glist = as.character(unique(temptenrichIn$Gene))
	enriched <- enrichr(glist, dbs)

	kegg = enriched$KEGG_2019_Human
	# WikiPathways = enriched$WikiPathways_2019_Human
	# GOBP = enriched$GO_Biological_Process_2018
	# GOBP$Term = gsub(" \\(.*", "",GOBP$Term)

	listInput= kegg

	list_temp <- data.frame(Term = gsub("_.*","",listInput$Term),
		log10_Adj_p.val = -log10(listInput$Adjusted.P.value),
		Count = as.numeric(sub("\\/.*", "\\1", listInput$Overlap, perl=TRUE)),
		Significant = unique(enrichIn$Significant)[i],
		Genes = listInput$Genes)

	list_temp = list_temp[order(list_temp$log10_Adj_p.val, decreasing=T),]
	list2plot = rbind(list2plot,list_temp[ c(1:15), ])
}


S1<- ggplot(list2plot, aes(x=Significant, y=Term, size=Count, color=log10_Adj_p.val)) + 
	geom_point(alpha = 0.8) + 
	# facet_grid(Contrasts ~ ., scales = "free_y") + 
	# theme_classic() +
	scale_color_gradient(low = "mediumblue",  high = "red4", space = "Lab", limit = c(min(list2plot$log10_Adj_p), max(list2plot$log10_Adj_p)))+
	theme(axis.title.x = element_text(color="black", size=20, face="bold"),
	  		axis.title.y = element_text(color="black", size=20, face="bold"),
	  		axis.text.x = element_text(face="bold", color="black",size=22, angle=90, hjust =0.5),
	  		axis.text.y = element_text(face="bold", color="black", size=20)) +
	scale_size(range = c(4, 12))

	  	
	png(file= paste("KEGG_DOTPlot_DEgenes_10 samples.png",sep=""), units="in", width=15, height=20, res=300)	
	print(S1)
	dev.off()
	#######################################################################################
