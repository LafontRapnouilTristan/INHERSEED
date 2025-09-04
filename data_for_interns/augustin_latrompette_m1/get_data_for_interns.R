# Load list of Augustin's species
sp_augustin <- readRDS('data_for_interns/augustin_latrompette_m1/results_all.rds')
sp_augustin$code_EPPO <- as.character(sp_augustin$code_EPPO)
sp_augustin$code_EPPO[which(sp_augustin$code_EPPO=="ESCHA")] <- "ESHCA"
eppo_augustin <- unique(sp_augustin$code_EPPO)
sp_augustin$origin <- as.character(sp_augustin$origin)
sp_augustin$origin[which(sp_augustin$origin=="INRAe, Vilmorin")] <- "Vilmorin"
# Load Traitor results
traitor_augustin <- readRDS("./data/final_outs/G0/macroscopy/traitor_measurements/traitor_out_cleaned.rds")
traitor_augustin <- traitor_augustin |>
  dplyr::filter(species%in%eppo_augustin) |>
  dplyr::select("image_name",
                "species",
                "aspect_ratio",
                "surface_structure",
                "solidity",
                "circularity",
                "area_scaled",
                "perimeter_scaled",
                "length_scaled",
                "width_scaled")

# Get eppo-genus correspondance
species_list <- readxl::read_xlsx("data/species_INHERSEED.xlsx")
species_list <- species_list |>
  dplyr::filter(EPPO%in%eppo_augustin) |>
  dplyr::mutate(GENUS=stringr::str_extract(SPECIES,"^[a-zA-Z]+"))

# Load Phylogenetic tree
tree_augustin <- phytools::readNexus("./data/plant_phylogeny/nexus_our_species.nxs")
tree_augustin <- tree_augustin |>
  ape::drop.tip(tip=tree_augustin$tip.label[!tree_augustin$tip.label%in%species_list$GENUS])
tree_augustin$tip.label <- species_list$SPECIES[match(tree_augustin$tip.label,species_list$GENUS)]

tree_augustin <- phytools::bind.tip(tree_augustin,
                   tip.label = "Impatiens balsamina",
                   edge.length = 0,
                   where = 16) 
# save

save(sp_augustin,species_list,traitor_augustin,tree_augustin,file="./data_for_interns/augustin_latrompette_m1/clean_data_for_augustin.RData")



# ACP

traitor_species <- traitor_augustin |>
  dplyr::summarise(dplyr::across(.cols = c("aspect_ratio","surface_structure","solidity","circularity",
                                           "area_scaled","perimeter_scaled","length_scaled","width_scaled" ),
                                 .fns = list(mean = mean,
                                             sd = sd)),
                   .by = species) |>
  dplyr::rename(code_EPPO=species)


num_col <- sp_augustin |>
  dplyr::select(where(is.numeric)) |>
  colnames()

cfu_species <- sp_augustin |>
  dplyr::select(-c("replicate","code_echt","date_maceration","date_etalement","date_lecture_bact")) |>
  dplyr::group_by(code_EPPO) |>
  dplyr::mutate(dplyr::across(.cols = is.numeric,
                              .fns = list(mean = mean,
                                          sd = sd))) |>
  dplyr::select(-all_of(num_col[-c(1,5)])) |>
  dplyr::distinct(.keep_all = TRUE)

data_all <- dplyr::left_join(traitor_species,cfu_species)
data_all  
pca_tmp <- FactoMineR::PCA(data_all[,which(colnames(data_all)%in%c("solidity_mean","circularity_mean","surface_structure_mean","aspect_ratio_mean","imbibition_par_graine_mL_mean","masse_1_graine_mean","area_scaled_mean","width_scaled_mean","length_scaled_mean","CFU_par_graine_mean","dilution_lecture_champi_mean","code_EPPO"))],
                       quali.sup = 1,
                       graph = F)




eig <- pca_tmp$eig
coord_pca <- as.data.frame(pca_tmp$ind$coord)

arrowMul <- function(arrows, data, at = c(0, 0), fill = 0.75) {
  u <- c(range(data[,1], range(data[,2])))
  u <- u - rep(at, each = 2)
  r <- c(range(arrows[, 1], na.rm = TRUE), range(arrows[, 2], na.rm = TRUE))
  rev <- sign(diff(u))[-2]
  if (rev[1] < 0)
    u[1:2] <- u[2:1]
  if (rev[2] < 0)
    u[3:4] <- u[4:3]
  u <- u/r
  u <- u[is.finite(u) & u > 0]
  fill * min(u)
}

coord_pca$species <- data_all$code_EPPO  
mul <- arrowMul(as.data.frame(pca_tmp$var$coord),
                pca_tmp$ind$coord)

df_text_pca <- coord_pca |> 
  dplyr::group_by(species) |> 
  dplyr::summarise(Dim.1=mean(Dim.1),Dim.2=mean(Dim.2))

pca_traitor <- ggplot()+
  geom_point(coord_pca,mapping=aes(x=Dim.1,y=Dim.2,color=species))+
  ggpubr::theme_classic2()+
  geom_vline(xintercept = 0,lty=2)+
  geom_hline(yintercept = 0,lty=2)+
  xlab(paste0("PC1 (",round(eig[1,2],2),"%)"))+
  ylab(paste0("PC2 (",round(eig[2,2],2),"%)"))+
  ggrepel::geom_text_repel(df_text_pca,mapping=aes(x=Dim.1,y=Dim.2,label=species), 
                           na.rm = TRUE,
                           show.legend = F) +
  geom_segment(data= as.data.frame(pca_tmp$var$coord),aes(x = 0, y = 0, xend=mul*Dim.1, yend=mul*Dim.2),
               lineend = "round", 
               linejoin = "round",
               linewidth = .75, 
               arrow = arrow(length = unit(0.2, "inches")),
               colour = "black" )+
  ggrepel::geom_text_repel(data = data.frame(var_names=rownames(pca_tmp$var$coord),pca_tmp$var$coord*mul), # add variable names at the end of arrows
                           aes(x = Dim.1*1.2, # nudge a bit the coordinates so that they're not on the arrows
                               y = Dim.2*1.2,
                               label = var_names),
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 15))+
  # scale_color_manual(values = climate_palette,name="Climatic zone")+
  theme(text=element_text(face="bold"))+
  guides(color = guide_legend(override.aes = list(size = 3) ) )

scree_plot_var_list <- factoextra::fviz_screeplot(pca_tmp,barfill = "lightgrey", barcolor="lightgrey")+
  geom_hline(yintercept = 100/(length(5:11)+1),color="red",lty=2)+
  ggpubr::theme_classic2()+
  theme(text=element_text(face="bold"),
        plot.title = element_blank())+
  scale_y_continuous(expand=c(0,0))

library(patchwork)
pca_traitor+theme(legend.position="none")+scree_plot_var_list
data_all_augustin <- data_all
saveRDS(data_all_augustin,"./data_for_interns/augustin_latrompette_m1/data_species.RDS")

########### Grosse fig ####
load("./data_for_interns/augustin_latrompette_m1/clean_data_for_augustin.RData")
library(patchwork)
library(ggplot2)
ptree <- ggtree::ggtree(tree_augustin)
pdil <- sp_augustin |>
  dplyr::mutate(SPECIES=forcats::fct_relevel(SPECIES, dplyr::arrange(ptree$data[1:32,],y)$label)) |>
  tidyr::pivot_longer(cols=c(bact_dilu_detect_min,dilution_lecture_champi),names_to = "var",values_to = "value") |>
  ggplot(aes(y=SPECIES,x=value,shape = var))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(color=origin))+
  theme_minimal()+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(hjust = .5) ,
        text=element_text(face="bold"))+
  paletteer::scale_color_paletteer_d("tvthemes::AirNomads")


hmaptraitor <- traitor_augustin |>
  dplyr::mutate(dplyr::across(dplyr::all_of(c("aspect_ratio","surface_structure","solidity","circularity",
                                              "area_scaled","perimeter_scaled","length_scaled","width_scaled")),
                              ~scale(.x,center = T,scale = T))) |>
  tidyr::pivot_longer(cols=c("aspect_ratio","surface_structure","solidity","circularity",
                             "area_scaled","perimeter_scaled","length_scaled","width_scaled"),
                      names_to = "Trait",
                      values_to = "Values") |>
  dplyr::summarise(means=mean(Values,na.rm=T),
                   sd=sd(Values,na.rm = T),
                   .by = c(species,Trait)) |>
  dplyr::rename(code_EPPO=species) |>
  dplyr::right_join(dplyr::select(sp_augustin,code_EPPO,SPECIES)) |> 
  dplyr::distinct() |>
  dplyr::mutate(SPECIES=forcats::fct_relevel(SPECIES, dplyr::arrange(ptree$data[1:32,],y)$label)) |>
  ggplot(aes(x=Trait,y=SPECIES,color=means,size=sd))+
  geom_point()+
  theme_minimal()+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle=45,hjust=1,vjust=1))+
  scale_color_gradientn(colors = c("darkblue","lightblue","red"))

ptree+pdil+hmaptraitor+plot_layout(guides="collect")
  

# final #####
cfu_finals <- readRDS("./data_for_interns/augustin_latrompette_m1/finals/results_final.rds")
traitor_finals <- readRDS("./data_for_interns/augustin_latrompette_m1/finals/traitor_final.rds")

sp_augustin <- readRDS('data_for_interns/augustin_latrompette_m1/results_all.rds')
sp_augustin$code_EPPO <- as.character(sp_augustin$code_EPPO)
sp_augustin$code_EPPO[which(sp_augustin$code_EPPO=="ESCHA")] <- "ESHCA"
eppo_augustin <- unique(sp_augustin$code_EPPO)
sp_augustin$origin <- as.character(sp_augustin$origin)
sp_augustin$origin[which(sp_augustin$origin=="INRAe, Vilmorin")] <- "Vilmorin"

traitor_finals <- traitor_finals |>
  dplyr::left_join(sp_augustin[,c(2,4,5,6,7,8,14)]) |>
  dplyr::rename("Longueur"=length_scaled_mean,
                "Largeur"=width_scaled_mean,
                "Surface"=area_scaled_mean,
                "Périmètre"=perimeter_scaled_mean,
                "Poids 1 graine"=masse_1_graine_mean,
                "Aspect ratio"=aspect_ratio_mean,
                "Circularité"=circularity_mean,
                "Convexité"=solidity_mean,
                "Rugosité"=surface_structure_mean,
                "Capacité d'imbibition"=imbibition_par_graine_mL_mean)


pca_tmp <- FactoMineR::PCA(traitor_finals[,which(colnames(traitor_finals)%in%c("Longueur","Largeur","Surface","Périmètre","Poids 1 graine",
                                                                               "Aspect ratio","Circularité","Convexité","Rugosité",
                                                                               "Capacité d'imbibition","code_EPPO"))],
                           quali.sup = 1,
                           graph = F)




eig <- pca_tmp$eig
coord_pca <- as.data.frame(pca_tmp$ind$coord)

coord_pca$species <- traitor_finals$code_EPPO  
coord_pca$Family <- traitor_finals$FAMILY[1:32*4]
mul <- arrowMul(as.data.frame(pca_tmp$var$coord),
                pca_tmp$ind$coord)

df_text_pca <- coord_pca |> 
  dplyr::group_by(species) |> 
  dplyr::summarise(Dim.1=mean(Dim.1),Dim.2=mean(Dim.2))

pca_traitor <- ggplot()+
  geom_point(coord_pca,mapping=aes(x=Dim.1,y=Dim.2,color=Family))+
  ggpubr::theme_classic2()+
  geom_vline(xintercept = 0,lty=2)+
  geom_hline(yintercept = 0,lty=2)+
  xlab(paste0("PC1 (",round(eig[1,2],2),"%)"))+
  ylab(paste0("PC2 (",round(eig[2,2],2),"%)"))+
  ggrepel::geom_text_repel(df_text_pca,mapping=aes(x=Dim.1,y=Dim.2,label=species), 
                           na.rm = TRUE,
                           show.legend = F,
                           size=2) +
  coord_fixed()+
  theme(text=element_text(face="bold"))+
  scale_color_manual(values=colorRampPalette( paletteer::paletteer_d("rcartocolor::Safe"))(20))


ggplot()+
  ggpubr::theme_classic2()+
  geom_segment(data= as.data.frame(pca_tmp$var$coord),aes(x = 0, y = 0, xend=Dim.1, yend=Dim.2),
               lineend = "round", 
               linejoin = "round",
               linewidth = .75, 
               arrow = arrow(length = unit(0.2, "inches")),
               colour = "black" )+
  ggrepel::geom_text_repel(data = data.frame(var_names=rownames(pca_tmp$var$coord),pca_tmp$var$coord), # add variable names at the end of arrows
                           aes(x = Dim.1, # nudge a bit the coordinates so that they're not on the arrows
                               y = Dim.2,
                               label = var_names),
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 15))+
  # scale_color_manual(values = climate_palette,name="Climatic zone")+
  theme(text=element_text(face="bold"))+
  guides(color = guide_legend(override.aes = list(size = 3) ) )+
  xlim(c(-1,1))+
  ylim(c(-1,1))+
  coord_fixed()+
  geom_vline(xintercept = 0,lty=2)+
  geom_hline(yintercept = 0,lty=2)+
  xlab(paste0("PC1 (",round(eig[1,2],2),"%)"))+
  ylab(paste0("PC2 (",round(eig[2,2],2),"%)"))+
  ggforce::geom_circle(aes(x0=0,y0=0,r=1))

scree_plot_var_list <- factoextra::fviz_screeplot(pca_tmp,barfill = "lightgrey", barcolor="lightgrey")+
  geom_hline(yintercept = 100/(length(5:11)+1),color="red",lty=2)+
  ggpubr::theme_classic2()+
  theme(text=element_text(face="bold"),
        plot.title = element_blank())+
  scale_y_continuous(expand=c(0,0))

library(patchwork)
pca_traitor+theme(legend.position="none")+scree_plot_var_list

pca_traitor+theme(legend.position="none")+theme(text=element_text(face="bold"))
