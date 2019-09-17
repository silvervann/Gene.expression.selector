source('/Users/mac/Google Drive/03_miRNA_DrZaffaroni/0.load_packages.R')
setwd("/Users/mac/Google Drive/03_miRNA_DrZaffaroni/Analisi_Rihan2")
library(SLqPCR)
library(plyr)
library(dplyr)
source("./geomMean.R")
source("./collapse.replicates.genes.R")


### Read annotation table
genefiles= read.delim("Annotation_data.Mpphagos_MyTuberculosis.txt")
###Read raw Ctr table
raw=as.matrix(read.table("data.Mpphagos_MyTuberculosis.txt", row.names=1, check.names=FALSE))
raw[raw=="Undetermined"] = as.double(40)
raw = as.data.frame(raw)

#  No Controls
raw = raw[,grep("CTR", colnames(raw), invert=T)]

## Data Celaning
### Split matrix in miRNAs with Non-replicated probes
rep_probes = unique(genefiles$Simbolo[duplicated(genefiles$Simbolo)])
rawmat= raw[!(genefiles$Simbolo %in% rep_probes),]
rownames(rawmat) = unique(genefiles$Simbolo)[unique(genefiles$Simbolo) != rep_probes]

## Matrix of miRNAs with technical replicates
contr_mirs=raw[(genefiles$Simbolo %in% rep_probes),]
contr_mirs$gene = genefiles$Simbolo[(genefiles$Simbolo %in% rep_probes)]

# Table for Kruskall-wallis test for control probes
contr_melt = melt(contr_mirs, id.vars= "gene" )
contr_melt$Delta_CT = as.numeric(contr_melt$value)
contr_melt$group = gsub(".*_","",contr_melt$variable)

# Compute summary statistics by Samples

contr_statistics =as.data.frame(group_by(contr_melt, gene, variable) %>%
  summarise(
    count = n(),
    mean = mean(Delta_CT, na.rm = TRUE),
    sd = sd(Delta_CT, na.rm = TRUE),
    median = median(Delta_CT, na.rm = TRUE),
    IQR = IQR(Delta_CT, na.rm = TRUE)
  ))

tabin_kw = contr_melt[contr_melt$gene==unique(contr_melt$gene)[2],]
kruskal.test(Delta_CT ~ variable, data = tabin_kw)

pairwise.wilcox.test(tabin_kw$Delta_CT, tabin_kw$variable,
                 p.adjust.method = "BH")

# Plot summary statistics by Samples
library(ggpubr)

bpl = ggplot(contr_statistics, aes(x=variable, y=mean,group=gene, color=gene)) + 
  	theme_minimal()+
  	scale_color_brewer(palette="Paired")+
    geom_line(aes(linetype=gene)) + 
    geom_point()+
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.1, position=position_dodge(0.05)) +
	ylab("raw_CT")+ 
	xlab("Sample")+
	rotate_x_text(angle = 45)+
	theme(axis.text.x = element_text(face="bold", color="black",size=6),
	axis.title.x = element_text(color="black", size=8),
	axis.text.y = element_text(face="bold", color="black", size=7),
	axis.title.y = element_text(color="black", size=8))
	# legend.position="none")
	# scale_fill_manual(values=cbbPalette[2:5])+
	# scale_fill_manual(values=cbbPalette[2:5])+ # USE FOR CLINICAL GROUPS
png(file= paste("Summary.controls.png",sep=""), units="in", width=5, height=5, res=200)
print(bpl)
dev.off()
##################### MATRIX FROM WHICH DiffExp FEATRUES ARE EXTRACTED
#######OPTION 2
###Function to control How many NAs could be accepted
#passNa <- function(DF, n) {DF[rowSums(is.na(DF)) <= n,]}
passCt <- function(DF, ct, n) {DF[apply(DF,1, function (x){sum(x>=ct)}) <= n,]}
### Select row Features (miRNAs) with at least 66% samples with detected values
newmat_filt0=passCt(rawmat, ct=35, n=0)
newmat_filt = apply(newmat_filt0,2, as.numeric)
rownames(newmat_filt) = rownames(newmat_filt0)
#############################################################################

newmat_targs = newmat_filt[!(rownames(newmat_filt) %in% as.character(genefiles$Simbolo[genefiles$Descripcion == "Reference_Genes"])) , ]
newmat_RG = newmat_filt[(rownames(newmat_filt) %in% as.character(genefiles$Simbolo[genefiles$Descripcion == "Reference_Genes"])) , ]


Hk_vect = sort(geneStabM(apply(t(newmat_filt),2, as.numeric)), decreasing =TRUE)

Hk_ranking=data.frame(Simbolo=names(Hk_vect), Rank=Hk_vect)
Hk_ranking = merge(Hk_ranking, genefiles, by = "Simbolo", all.x = T)
Hk_ranking$Rank = as.numeric(Hk_ranking$Rank)
Hk_ranking = Hk_ranking[order(Hk_ranking$Rank, decreasing =TRUE), ]

Hk_ranking$Simbolo = factor( Hk_ranking$Simbolo, levels= Hk_ranking$Simbolo)

##############################################################################
# Plot M stability Values for Reference Genes
  bplotMvals = ggplot(Hk_ranking, aes(x=Simbolo, y=Rank, fill=Descripcion)) +
  theme_minimal() +
  scale_fill_manual(values= c("grey80", "red3") ) +   
    geom_bar(stat="identity", color="black",size=0.3) +
    ylab("M stability")+ 
    xlab("Gene")+
    theme(panel.background = element_rect(fill = NA),
      axis.title.x = element_text(color="black", size=15, face="bold"),
      axis.title.y = element_text(color="black", size=15, face="bold"),
      axis.text.x = element_text(face="bold", color="black",size=13,angle=90,hjust =0, vjust =1),
      axis.text.y = element_text(face="bold", color="black", size=15,angle=0))

png(file= paste0("Mvals_stability_Barchart.png"), units="in", width=17, height=10, res=200)
print(bplotMvals)
dev.off()

##############################################################################

### Normalization of expression values using geometric mean of most stable reference genes
Hk=Hk_ranking$Simbolo[Hk_ranking$Descripcion=="Reference_Genes"][4:5]
# Hk=Hk_ranking$Simbolo[1:3]
Hkcol=which(rownames(newmat_filt) %in% Hk==TRUE )##[c(1,3,4)]

# Normalization for Genes
NF <- apply(t(newmat_filt)[, Hk], 1, geomMean)
qNormexprs= newmat_targs - NF

# Normalization for Reference genes
qNormexprs_RG= newmat_RG - NF

##############################################################################
# TTFch_Prob = NULL


cat.design = gsub(".*_","",colnames(qNormexprs))
mat = -qNormexprs ## This Should be negative for the T.test, so we obtain logFCh values suitable for straigthforwar interpretation abput the Relative expression
###DIff Exprss
design <- model.matrix(~0+ cat.design)
colnames(design) <- sub(".*cat.design", "",colnames(design) )
contrasts <- makeContrasts(MTBpos-MTBneg, levels=design)
# contrasts <- makeContrasts(MTBpos-MTBneg,MTBpos-CTR,MTBneg-CTR, levels=design)
fit <- lmFit(mat, design)
# The actual test
fit2 = contrasts.fit(fit, contrasts)
eb = eBayes(fit2)
TTFch= topTable(eb,number=dim(mat)[1])
TTFch$DDCT= TTFch$logFC
TTFch$'2^-DDCT'= 2^(TTFch$DDCT)
TTFch$Gene = rownames(TTFch)

# TTFch$Contrast = "MTBneg vs CTR"

# TTFch_Prob = rbind(TTFch_Prob,  TTFch)

# TTFch = TTFch_Prob

M= log2(1.5)
P= 0.05

Significant= ifelse(TTFch$logFC>M & TTFch$adj.P.Val<P, "Up-modulated",ifelse(TTFch$logFC< -M & TTFch$adj.P.Val<P, "Down-modulated", "Not Sig"))
TTFch$Significant = ifelse(Significant=="Not Sig", "Not Sig",Significant)


write.table(TTFch,"./Differentially_modulated_Genes.txt", quote=FALSE, sep= "\t")

selectedmirs=rownames(TTFch)[TTFch$adj.P.Val < 0.05]
TTFch_sel = TTFch[ TTFch$adj.P.Val < 0.05, ]
TTFch_sel = TTFch_sel[order(TTFch_sel$DDCT), ]


mat2= mat[rownames(mat) %in% selectedmirs,]

##############################################################################

qNorm_melt = melt(qNormexprs)

qNormbpl = ggplot(qNorm_melt, aes(x=Var2, y=value,fill=Var2)) + 
  geom_boxplot() +
  ylab("Delta_CT")+ 
  xlab("Sample")+
  rotate_x_text(angle = 45)+
  theme(axis.text.x = element_text(face="bold", color="black",size=15),
  # axis.title.x = element_blank(),
  axis.text.y = element_text(face="bold", color="black", size=15),
  axis.title.y = element_text(color="black", size=16),
  legend.position="none")
  # scale_fill_manual(values=cbbPalette[2:5])+
  # scale_fill_manual(values=cbbPalette[2:5])+ # USE FOR CLINICAL GROUPS
png(file= paste("gNomr.HKnorm.deltaCt.png",sep=""), units="in", width=11, height=11, res=200)
print(qNormbpl)
dev.off()
