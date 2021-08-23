
limma_list_liver <- readRDS("results_RNAseqRDS/limma_Liver.RDS")
limma_list_spleen <- readRDS("results_RNAseqRDS/limma_Spleen.RDS")
wgcna_all_spleen <- readRDS("results_RNAseqRDS/WGCNA_MM_GS_all_Spleen.RDS") %>% 
  left_join(grcm38 %>% select(ensgene,symbol)) %>% select(-ensgene)
wgcna_all_liver <- readRDS("results_RNAseqRDS/WGCNA_MM_GS_all_Liver.RDS") %>% 
  left_join(grcm38 %>% select(ensgene,symbol)) %>% select(-ensgene)

y_liver <- readRDS("results_RNAseqRDS/limma_vdl_Liver")
y_spleen <- readRDS("results_RNAseqRDS/limma_vdl_Spleen")



rnaseq_rawcounts %>% 
  filter(ensgene =='ENSMUSG00000032691')


name_any_gene()

my_gene_list <- c("Tmsb4x","Cd74","Ly6e","Tmsb10","H2-D1","Malat1",
                  "Cd52","Ubb","Cyba","Sh3bgrl3",
                  "Apoa1","Serpina1a","Serpina1c",
                  "Rbp4","Ttr","Fabp1","Serpina3k",
                  "Apoa2","Cyp3a11","Alb")

i = 9
for(i in 1:length(my_gene_list)){
  my_gene_list[i] -> gene
  name_any_gene(gene)
}

for(i in 1:length(my_gene_list)){
  my_gene_list[i] -> gene
  png(filename = glue("results_figures/aging_exp_{gene}.png"),
      width = 8, height = 20, units = "cm", res = 600)
  plot_any_gene(gene)
  invisible(dev.off())
}






plot_any_gene("Ly6e")
name_any_gene("Ly6e")




name_any_gene("Il18")
name_any_gene("Nlrp3")
plot_any_gene("Nod2")
plot_any_gene("Nek7") # interesting 
name_any_gene("Nek7")



name_any_gene("Nrf2")
name_any_gene("Myc")
name_any_gene("Nr4a1")
name_any_gene("Cd36")
name_any_gene("Etfdh")

name_any_gene("Hmox1")

limma_list_liver$status %>% 
  filter(grepl("Heme",description, ignore.case = T)) %>% 
  unique()

limma_list_spleen$status %>% 
  # select(-description) %>%
  filter(grepl("Casp",symbol)) %>% 
  unique()

limma_list_liver$status %>% 
  # select(-description) %>%
  filter(grepl("Il18",symbol)) %>% 
  unique()

limma_list_liver$status %>% 
  # select(-description) %>%
  filter(grepl("nlrp",symbol,ignore.case = T)) %>% 
  unique()

limma_list_liver$status %>% 
  filter(grepl("NLR",description, ignore.case = T)) %>% 
  unique()

limma_list_spleen$status %>% 
  filter(grepl("NLR",description, ignore.case = T)) %>% 
  unique()

limma_list_liver$status %>% 
  # select(-description) %>%
  filter(grepl("tnf",symbol,ignore.case = T)) %>% 
  unique()

limma_list_liver$status %>% 
  # select(-description) %>%
  filter(grepl("casp",symbol,ignore.case = T)) %>% 
  unique()

limma_list_liver$status %>% 
  # select(-description) %>%
  filter(grepl("cxcl",symbol,ignore.case = T)) %>% 
  unique()

limma_list_liver$status %>% 
  # select(-description) %>%
  filter(grepl("aim",symbol,ignore.case = T)) %>% 
  unique()

limma_list_liver$status %>% 
  # select(-description) %>%
  filter(grepl("asc",symbol,ignore.case = T)) %>% 
  unique()

limma_list_liver$status %>% 
  filter(grepl("cancer susceptibility",description, ignore.case = T)) %>% 
  unique()

limma_list_liver$status %>% 
  filter(grepl("interferon",description, ignore.case = T)) %>% 
  unique()


plot_any_gene("Casc4")


limma_list_spleen$status %>% 
  # select(-description) %>%
  filter(grepl("Ly6g",symbol)) %>% 
  unique()


limma_list_liver$status %>% 
  unique() %>% 
  filter(P.Value <0.05) %>% 
  arrange(desc(logFC)) %>% 
  nrow
  head(10)


limma_list_liver$status %>% 
  unique() %>% 
  filter(P.Value <0.05) %>% 
  arrange(logFC) %>% 
  head(10)

limma_list_liver$cort %>% 
  unique() %>% 
  filter(P.Value <0.05) %>% 
  arrange(desc(logFC)) %>% 
  head(10)


limma_list_liver$cort %>% 
  unique() %>% 
  filter(P.Value <0.05) %>% 
  arrange(logFC) %>% 
  head(10)



limma_list_spleen$status %>% 
  unique() %>% 
  filter(P.Value <0.05) %>% 
  arrange(desc(logFC)) %>% 
  head(24)


limma_list_spleen$status %>% 
  unique() %>% 
  filter(P.Value <0.05) %>% 
  arrange(logFC) %>% 
  head(10)

limma_list_spleen$cort %>% 
  unique() %>% 
  filter(P.Value <0.05) %>% 
  arrange(desc(logFC)) %>% 
  head(25)


limma_list_spleen$cort %>% 
  unique() %>% 
  filter(P.Value <0.05) %>% 
  arrange(logFC) %>% 
  head(10)


limma_list_spleen$status %>% 
  rename(spleenlogFC = logFC,
         spleenPval = P.Value) %>% 
  select(symbol, spleenlogFC, spleenPval) %>% 
  left_join(limma_list_liver$status %>% 
              rename(liverlogFC = logFC,
                     liverPval = P.Value) %>% 
              select(symbol, liverlogFC, liverPval)) %>%
  filter(spleenPval <0.05 & liverPval < 0.05) %>% 
  filter(spleenlogFC*liverlogFC <0) %>% 
  arrange(spleenlogFC*liverlogFC) %>% 
  left_join(grcm38 %>% select(entrez, symbol)) %>% 
  .$entrez %>% 
  write.table("results_GeneWalk/status_diff_direction.csv", row.names = F)

# 391 genes total , 119 genes different direction   
47/206


limma_list_spleen$cort %>% 
  rename(spleenlogFC = logFC,
         spleenPval = P.Value) %>% 
  select(symbol, spleenlogFC, spleenPval) %>% 
  left_join(limma_list_liver$cort %>% 
              rename(liverlogFC = logFC,
                     liverPval = P.Value) %>% 
              select(symbol, liverlogFC, liverPval)) %>%
  filter(spleenPval <0.05 & liverPval < 0.05) %>% 
  filter(spleenlogFC*liverlogFC <0) %>% 
  arrange(spleenlogFC*liverlogFC) %>% 
  left_join(grcm38 %>% select(entrez, symbol)) %>% 
  .$entrez %>% 
  write.table("results_GeneWalk/CORT_diff_direction.csv", row.names = F)
  
# 391 genes total , 119 genes different direction   
119/391


limma_list_liver$status %>% 
  # select(-description) %>%
  filter(grepl("p53",symbol)) %>% 
  unique()


limma_list_liver$status %>% 
  select(-description) %>%
  filter(grepl("Mup",symbol)) %>% 
  unique()

limma_list_liver$status %>% 
  select(-description) %>%
  filter(grepl("Stat",symbol)) %>% 
  unique()

limma_list_liver$status %>% 
  select(-description) %>%
  filter(grepl("Il",symbol)) %>% 
  unique()

limma_list_liver$status %>% 
  select(-description) %>%
  filter(grepl("Olf",symbol)) %>% 
  # filter(P.Value < 0.05) %>% 
  unique() %>% nrow
16/326

limma_list_liver$status %>% 
  # select(-description) %>%
  filter(grepl("Vmn",symbol)) %>% 
  filter(P.Value < 0.05) %>%
  unique() # %>% nrow
6/137


limma_list_spleen$status %>% 
  select(-description) %>%
  filter(grepl("Olf",symbol)) %>% 
  filter(P.Value < 0.05) %>%
  arrange(logFC) %>% 
  unique() # %>% nrow
136/677

limma_list_spleen$status %>% 
  select(-description) %>%
  filter(grepl("Vmn",symbol)) %>% 
  filter(P.Value < 0.05) %>%
  unique() # %>% nrow
62/204



name_any_gene("Ube3a")
name_any_gene("Decr2")
name_any_gene("Scd1") # Liver cyan hub gene, significantly increased
name_any_gene("Pde3b")

name_any_gene("Eny2")
name_any_gene("Mup1") # GO term negative regulation of insulin 
name_any_gene("Ucp2")
name_any_gene("Hmgcr")


wgcna_all_liver %>% 
  left_join(grcm38 %>% select(symbol, description)) %>% 
  filter(grepl("coagulation",description,ignore.case = T)) %>% 
  unique()


wgcna_all_liver %>% 
  filter(grepl("Mup",symbol)) %>% 
  unique()# apparently not all MUPs are in the pink module! 

wgcna_all_spleen %>% 
  filter(grepl("Vmn",symbol))

wgcna_all_liver %>% 
  filter(grepl("Cdk",symbol))

wgcna_all_liver %>% 
  filter(grepl("Il",symbol))


wgcna_all_liver %>% 
  filter(grepl("Jak",symbol))

wgcna_all_liver %>% 
  filter(grepl("Stat",symbol))


wgcna_all_spleen %>% 
  filter(grepl("Cdk",symbol))

wgcna_all_liver %>% 
  filter(grepl("Cd",symbol)) # apparently not all MUPs are in the pink module! 

wgcna_all_spleen %>% 
  filter(grepl("Cd",symbol)) %>% 
  arrange(module)



wgcna_all_liver %>% 
  filter(grepl("Cyp",symbol)) %>% 
  arrange(abs(GS.status))

wgcna_all_spleen %>% 
  filter(grepl("Cyp",symbol)) %>% 
  arrange(abs(GS.status))

wgcna_all_liver %>% 
  filter(module == "pink") %>% 
  arrange(desc(abs(moduleMembership)))

plot_any_gene("Gm13712")
plot_any_gene("Foxa1")
plot_any_gene("Mup17")
plot_any_gene("Cyp21a1")
plot_any_gene("Cdhr3")
plot_any_gene("Hormad1")
plot_any_gene("Slc25a29")
plot_any_gene("Papln")
plot_any_gene("Cox15")

my_gene = "Papln"


plot_any_gene("Cd80")
plot_any_gene("Cd")

plot_any_gene("Adora2a")
plot_any_gene("Fpr2")

name_any_gene("Mup20")
name_any_gene("Mup3")
name_any_gene("Mup16")
name_any_gene("Mup9")

name_any_gene("Jak2")
plot_any_gene("Jak2")

name_any_gene("Rsl1")
name_any_gene("Rsl2")
name_any_gene("Cyp2a4")
name_any_gene("Ihh")


name_any_gene("Cyp8b1") #liver specific  # important! 
name_any_gene("Apoa4") #liver specific


name_any_gene("Scamp3")
name_any_gene("Avpi1")

name_any_gene("Prdx4")
plot_any_gene("Prdx4")
name_any_gene("Cd36")
plot_any_gene("Cd36")


name_any_gene("Gcnt4")
plot_any_gene("Gcnt4")

name_any_gene("Kirrel3") #none
name_any_gene("Trpc2") #none 


name_any_gene("Hrg")
name_any_gene("Maf")
plot_any_gene("Cd302")

name_any_gene("Mug2")



name_any_gene("Scd1") # upregulated in alpha in both liver and spleen
name_any_gene("Cd55")

name_any_gene("Sik1")

name_any_gene("Nr3c1")

name_any_gene("Cyp8b1") #something to think about 
plot_any_gene("Cyp8b1")

name_any_gene("Cyp2a22")
plot_any_gene("Cyp2a22")

name_any_gene("Ear2")
plot_any_gene("Myc")




name_any_gene("Cyp1a2")
name_any_gene("Cyp1b1")
plot_any_gene("Cyp1b1") # clearly decreasey with cort increase 
name_any_gene("Cyp4a14")
name_any_gene("Cyp26a1")
name_any_gene("Cyp2c9") # no
name_any_gene("Cyp2d9")
name_any_gene("Cyp2c19") # no
name_any_gene("Cyp2d6") # no
name_any_gene("Cyp3a4") # no
name_any_gene("Cyp3a5") # no 


name_any_gene("Cyp2c11")
name_any_gene("Cyp2d10")
name_any_gene("Cyp2d9")
name_any_gene("Cyp3a11") # interaction effect # ~ status both in liver and spleen
plot_any_gene("Cyp3a11")


name_any_gene("Cdkn2a") 
name_any_gene("Cdk8") 

name_any_gene("Cd36")
name_any_gene("Pgc1a") # no
name_any_gene("Pparg")
name_any_gene("Igfals")


# trafficking molecules - Science 2019 paper stress 
name_any_gene("Ccl2")
name_any_gene("Icam1")
name_any_gene("Cxcl10") # all significant in spleen
plot_any_gene("Cxcl10")


name_any_gene("Ghr")
plot_any_gene("Ghr")


wgcna_all_liver %>% 
  filter(grepl("Mup",symbol))

wgcna_all_spleen %>% 
  filter(grepl("Mup",symbol))




name_any_gene("Fbl")
gene = "Hrg"

gene = "Smox"
gene = "Timm10"
gene = "Il1b"
gene = "Il6"
gene = "Il2"
gene = "Il12b"
gene = "Myc"
gene = "Lif"
gene = "Tnfaip8l1"
gene = "Slc10a2"
gene = "Trim5"
gene = "Hsd17b13"
gene = "Mup17"
gene = "Mup16"
# https://www.pnas.org/content/pnas/98/19/10630.full.pdf
gene = "Apoa4"
gene = "Fasn"
gene = "Acaca"
gene = "Acly"
gene = "Mgat1"
gene = "Alad"
gene = "Taldo1"
gene = "Pdha1"
gene = "Lpin2"
gene = "G0s2"
gene = "Ppara"
gene = "Scd1"
gene = "Ubp1"
gene = "Elovl6"
gene = "Ppargc1a"
gene = "Srebf1"
gene = "Cyp3a11"


gene = "Cpt1a"

gene = "Eno2"

gene = "Maf"
gene = "Tgfb1"
gene = "Il1a"

gene = "Cd36" # "Cd36 encodes a fatty acid transporter that facilitates uptake and intracellular trafficking of fatty acids" 
gene = "Apoa4"

gene = "Scd2"
gene = "Lcn2"
gene = "Apom"

gene = "Pdk4"
gene = "Fkbp5"

gene = "Fgf21"
gene = "Cyp8b1"
gene = "G0s2"
gene = "Acaa1b"


gene = "Crybb1"

gene = "Olfr497"
gene = "Cd302"
gene = "Maf"
gene = "Jak2"
gene = "Cr1l"
gene = "Gstm1"
gene = "Serpina6"
gene = "Cdkn2b"
gene = "Junb"
gene = "Cacna1s"


gene = "Tgfb1"
gene = "Ifng"

gene = "Cpt1a"





name_any_gene(gene)

# png(filename = glue("results_figures/expression_{gene}.png"),
#     width = 10, height = 20, units = "cm", res = 600)

gene = "Hrg"
gene = "Hp"
gene = "Hpn"
gene = "Maf"
gene = "Hpx"
gene = "Ube3a"
gene = "Ube3b"
plot_any_gene(gene)
# dev.off()

plot_any_gene("Cox15")

plot_any_gene("Myc")
name_any_gene("Myc")
plot_any_gene("Pmpca")
plot_any_gene("Eci1")
plot_any_gene("Ndufa8")
name_any_gene("Pnpla2")
plot_any_gene("Pnpla2")

name_any_gene("Maob")
name_any_gene("Orm1")
plot_any_gene("Cox15")

plot_any_gene("Eci1")
plot_any_gene("Ndufa8")


limma_list_liver$status %>% 
  filter(symbol %in% focus) %>% # run Codina first 
  # filter(P.Value <0.05) %>%
  select(symbol, logFC, P.Value, description) %>% 
  arrange(desc(abs(logFC)))


# James's inquiry
gois <-c("Cpt1a", "Pdk4", "Gpd1", "Fgf21", "Auh",  "Scd1")
gois <- c("Cpt1a", "Cpt2", "Acadvl", "Hmgcs2", "Hmgcl", "Acat1", "Hadha")
gois <- c("Cpt1a","Srebp1c", "Fasn", "Acly", "Mlycd", "Mcat", "G0s2", "Slc22a5", "Aacs", "Acaca", "Apoa4", "Abhd2","Ldlr","Stat5b", "Acc")

gois <-c("Moxd1", "Cyp4a12","Slp","Elovl3", "Cyp2d9","Cyp3a10")

gois <-c("Apoc1","Cyp8b1")

limma_list_liver$status %>% 
  filter(symbol %in% gois) %>% 
  # filter(P.Value <0.05) %>%
  select(symbol, logFC, P.Value, description) %>% 
  arrange(logFC)


limma_list_liver$status %>% 
  # select(-description) %>%
  filter(grepl("Acsl",symbol)) %>% 
  unique()

limma_list_liver$status %>% 
  # select(-description) %>%
  filter(grepl("Fabp",symbol)) %>% 
  unique()


limma_list_liver$status %>% 
  # select(-description) %>%
  filter(grepl("Apoa",symbol)) %>% 
  unique()
  
limma_list_liver$status %>% 
  # select(-description) %>%
  filter(grepl("malonyl",description, ignore.case = T)) %>% 
  unique()


limma_list_liver$status %>% 
  # select(-description) %>%
  filter(grepl("acetyl-coa",description, ignore.case = T)) %>% 
  unique()



limma_list_liver$status %>% 
  # select(-description) %>%
  filter(grepl("Slc",symbol)) %>% 
  unique() %>% 
  arrange(abs(logFC)) %>% 
  filter(P.Value <0.05) %>% 
  select(-description)


# Spleen ==========================
name_any_gene("Jak1")
name_any_gene("Plk3")
name_any_gene("Stat3")
name_any_gene("Orm1")
plot_any_gene("Orm1")

gois <-c("Sik2","Dnajb14","Ndufa6","Ensa", "Arhgap31",
         "cyp19a1","Foxp3",
         "Cd207",
         "Bicd1","Fgd4", "Cma1", "Tnf","Tnfaip3", "Pde5a")



cbind(tissue = "Liver",limma_list_liver$status) %>% 
  rbind(cbind(tissue = "Spleen",limma_list_spleen$status)) %>%  
  filter(symbol %in% gois) %>% 
  select(tissue, symbol, logFC, P.Value, description) %>% 
  arrange(tissue,logFC)


limma_list_spleen$status %>% 
  # select(-description) %>%
  filter(grepl("Cox",symbol)) %>% 
  unique() %>% 
  arrange(desc(abs(logFC))) %>% 
  # filter(P.Value <0.05) %>% 
  select(-description)

limma_list_spleen$status %>% 
  # select(-description) %>%
  filter(grepl("Il",symbol)) %>% 
  unique() %>% 
  arrange(desc(abs(logFC))) %>% 
  filter(P.Value <0.05) %>%
  select(-description)


limma_list_spleen$status %>% 
  # select(-description) %>%
  filter(grepl("Ccl",symbol)) %>% 
  unique() %>% 
  arrange(desc(abs(logFC))) %>% 
  # filter(P.Value <0.05) %>% 
  select(-description)

limma_list_spleen$status %>% 
  # select(-description) %>%
  filter(grepl("Npy",symbol)) %>% 
  unique() %>% 
  arrange(desc(abs(logFC))) %>% 
  # filter(P.Value <0.05) %>% 
  select(-description)


limma_list_spleen$status %>% 
  # select(-description) %>%
  filter(grepl("Stat",symbol)) %>% 
  unique() %>% 
  arrange(desc(abs(logFC))) %>% 
  # filter(P.Value <0.05) %>%
  select(-description)

limma_list_spleen$status %>% 
  # select(-description) %>%
  filter(grepl("Jak",symbol)) %>% 
  unique() %>% 
  arrange(desc(abs(logFC))) %>% 
  # filter(P.Value <0.05) %>%
  select(-description)



limma_list_spleen$status %>% 
  # select(-description) %>%
  filter(grepl("Hsp",symbol)) %>% 
  unique() %>% 
  arrange(desc(abs(logFC))) %>% 
  filter(P.Value <0.05) %>%
  select(-description)


limma_list_spleen$status %>% 
  # select(-description) %>%
  filter(grepl("Sod",symbol)) %>% 
  unique() %>% 
  arrange(desc(abs(logFC))) %>% 
  # filter(P.Value <0.05) %>%
  select(-description)

# glutathione S-transferease 
limma_list_liver$status %>% 
  # select(-description) %>%
  filter(grepl("Gst",symbol)) %>% 
  unique() %>% 
  arrange(desc(abs(logFC))) %>% 
  filter(P.Value <0.05) %>%
  select(-description)


# glutathione S-transferease 
cbind(tissue = "Liver",limma_list_liver$status) %>% 
  rbind(cbind(tissue = "Spleen",limma_list_spleen$status)) %>% 
  # select(-description) %>%
  filter(grepl("Gpx",symbol, ignore.case = F)) %>% 
  unique() %>% 
  arrange(tissue,desc(abs(logFC))) %>% 
  filter(P.Value <0.05) %>%
  select(-description)



limma_list_spleen$status %>% 
  # select(-description) %>%
  filter(grepl("Gpx",symbol)) %>% 
  unique() %>% 
  arrange(desc(abs(logFC))) %>% 
  # filter(P.Value <0.05) %>%
  select(-description)


limma_list_spleen$cort %>% 
  # select(-description) %>%
  filter(grepl("Gpx",symbol)) %>% 
  unique() %>% 
  arrange(desc(abs(logFC))) %>% 
  # filter(P.Value <0.05) %>%
  select(-description)

limma_list_spleen$cort %>% 
  # select(-description) %>%
  filter(grepl("Gst",symbol)) %>% 
  unique() %>% 
  arrange(desc(abs(logFC))) %>% 
  # filter(P.Value <0.05) %>%
  select(-description)


plot_any_gene("Gsto2")

plot_any_gene("Cd207")
plot_any_gene("Fkbp1a")
plot_any_gene("Tfap2a")


limma_list_spleen$status %>% 
  # select(-description) %>%
  filter(grepl("pola",symbol, ignore.case = T)) %>% 
  unique() %>% 
  arrange(desc(abs(logFC))) %>% 
  # filter(P.Value <0.05) %>% 
  select(-description)
name_any_gene("Pola1")



# what are DEGs in both liver and spleen per CORT?
limma_list_liver$status %>% 
  filter(abs(logFC) >0.15) %>% 
  filter(P.Value <0.05) %>% 
  .$symbol -> x

limma_list_spleen$status %>% 
  filter(abs(logFC) >0.15) %>% 
  filter(P.Value <0.05) %>% 
  .$symbol -> y


x[x %in% y]
y[y %in% x]
plot_any_gene("Mup20")
plot_any_gene("Mup17")



limma_list_liver$status %>% nrow -> liver_n
limma_list_liver$status %>% 
  filter(P.Value < 0.05 & abs(logFC)>0.15) %>% 
  unique() %>% 
  nrow -> liver_status_DEG
limma_list_liver$cort %>% 
  filter(P.Value < 0.05 & abs(logFC)>0.15) %>% 
  unique() %>% 
  nrow -> liver_cort_DEG

liver_n
liver_status_DEG
liver_cort_DEG
liver_status_DEG/liver_n
liver_cort_DEG/liver_n


limma_list_spleen$status %>% nrow -> spleen_n
limma_list_spleen$status %>% 
  filter(P.Value < 0.05 & abs(logFC)>0.15) %>% 
  unique() %>% 
  nrow -> spleen_status_DEG
limma_list_spleen$cort %>% 
  filter(P.Value < 0.05 & abs(logFC)>0.15) %>%
  unique() %>% 
  nrow -> spleen_cort_DEG

spleen_n
spleen_status_DEG
spleen_cort_DEG
spleen_status_DEG/spleen_n
spleen_cort_DEG/spleen_n

