# options(repos = c(CREN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN"))

library('plotgardener')

region_plot <- function(chr, start, end, ymax = 6) {
  
  pageCreate(width = 7.5, height = 4, default.units = "inches")
  
  region <- pgParams(chrom = chr, chromstart = start, chromend = end, assembly = "hg38", 
                     range = c(0, ymax), default.units = "inches", x = 0.5)
  
  plotSignal(data = '/tank/projects/amalygina/SIRT6_splicing/cutrun/nfcore_output_cpm/03_peak_calling/03_bed_to_bigwig/SIRT6_def_R1.bigWig', 
             params = region, y = 0.25, width = 6.5, height = 0.65, fill = "#7ecdbb", linecolor = "#7ecdbb", just = c("left", "top"))
  
  plotSignal(data = '/tank/projects/amalygina/SIRT6_splicing/cutrun/nfcore_output_cpm/03_peak_calling/03_bed_to_bigwig/SIRT6_def_R2.bigWig', 
             params = region, y = 1, width = 6.5, height = 0.65, fill = "#7ecdbb", linecolor = "#7ecdbb", just = c("left", "top"))
  
  plotSignal(data = '/tank/projects/amalygina/SIRT6_splicing/cutrun/nfcore_output_cpm/03_peak_calling/03_bed_to_bigwig/WT_R1.bigWig', 
             params = region, y = 1.75, width = 6.5, height = 0.65, fill = "#37a7db", linecolor = "#37a7db", just = c("left", "top"))
  
  plotSignal(data = '/tank/projects/amalygina/SIRT6_splicing/cutrun/nfcore_output_cpm/03_peak_calling/03_bed_to_bigwig/WT_R2.bigWig', 
             params = region,  y = 2.5, width = 6.5, height = 0.65, fill = "#37a7db", linecolor = "#37a7db", just = c("left", "top"))
  
  plotGenomeLabel(params = region, y = 3.2, length = 6.5, fontsize = 12)
  
  plotText(label = "SIRT6-KO1", fontsize = 17, y = 0.25, just = c("left", "top"), params = region)
  plotText(label = "SIRT6-KO2", fontsize = 17, y = 1, just = c("left", "top"), params = region)
  plotText(label = "WT1", fontsize = 17, y = 1.75, just = c("left", "top"), params = region)
  plotText(label = "WT2", fontsize = 17, y = 2.5, just = c("left", "top"), params = region)
  
  genesPlot <- plotGenes(params = region, fill = "darkgrey", fontcolor = "darkgrey", y = 3.4, width = 6.5, 
                         height = 0.65, just = c("left", "top"), fontsize = 10)
  
  pageGuideHide()
}

region_plot("chr10", 96113905, 96219375) # ZNF185A
region_plot("chr5", 151170124, 151235124) # CCDC69 151180124, 151225124
region_plot("chrX", 129642759, 129657259, 2) # APLN 129644259, 129655759
region_plot("chr1", 205215912, 205285912) # TMCC2
region_plot("chr1", 156666763, 156679450) # NES
region_plot("chr11", 73765872, 73887872) # MRPL48
region_plot("chr1", 201125172, 201181772) # TMEM9

# TANC2 chr17:62,966,235-63,427,703 (461,469) yes small
# TMEM70 chr8:73,972,437-73,982,783 (10,347) yes (all organisms)
# IGSF8 chr1:160,091,339-160,101,109 (9,771) yes
# CNN2 chr19:1,026,547-1,039,077 (12,531) yes
# IDS chrX:149,476,988-149,521,096 (44,109) small (all organisms)
# TMCC2
# NES
# EPN2 chr17:19,214,837-19,336,715 (99,350) yes
# TMEM9

region_plot("chr10", 96113905, 96219375, 8) # ZNF185A
region_plot("chr20", 41126931, 41197801, 11) # PLCG1
region_plot("chr19", 295573, 354815, 9) # MIER2
region_plot("chr19", 40660905, 40695972, 6.5) # NUMBL
region_plot("chr6", 116271760, 116281930, 10) # TSPYL1
region_plot("chr2", 200808094, 200829838, 15) # BZW1


# PLCG1 chr20:41,136,931-41,196,801 (59,871) yes
# MIER2 chr19:305,573-344,815
# NUMBL chr19:40,665,905-40,690,972 (25,068) 
# TSPYL1 chr6:116,267,760-116,279,930 (12,171)
# BZW1 chr2:200,810,594-200,827,338 (16,745)