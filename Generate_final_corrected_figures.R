# FINAL CORRECTED PUBLICATION FIGURES
# Preserves the previous layouts and treatment colours while using the amended
# statistical workflow. This script reads raw data and model outputs produced by
# Final_statistical_check_models.R. It never modifies the original figure files.

required <- c("tidyverse", "patchwork")
miss <- setdiff(required, rownames(installed.packages()))
if (length(miss)) stop("Install required packages: ", paste(miss, collapse = ", "))
suppressPackageStartupMessages({library(tidyverse); library(patchwork)})

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
stats_dir <- file.path(project_dir, "Final_statistical_check_outputs")
fig_dir <- file.path(project_dir, "Final_corrected_figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

treat_order <- c("Healthy control", "PGPR A only", "PGPR B only", "PGPR A + B only",
                 "Pathogen only", "PGPR A + Pathogen", "PGPR B + Pathogen",
                 "PGPR A + B + Pathogen")
with_path <- treat_order[5:8]
treat_cols <- c(
  "Healthy control"="black", "PGPR A only"="#1F78B4",
  "PGPR B only"="#33A02C", "PGPR A + B only"="#A6CEE3",
  "Pathogen only"="#E31A1C", "PGPR A + Pathogen"="#FDBF6F",
  "PGPR B + Pathogen"="#B2DF8A", "PGPR A + B + Pathogen"="#6A3D9A")

theme_pub <- theme_bw(base_size=18) + theme(
  plot.title=element_text(size=22,face="bold"),
  axis.title=element_text(size=19,face="bold"),
  axis.text.x=element_text(size=13,angle=45,hjust=1,face="bold"),
  axis.text.y=element_text(size=16), strip.text=element_text(size=17,face="bold"),
  legend.title=element_text(size=18,face="bold"), legend.text=element_text(size=15),
  panel.grid.minor=element_blank(), panel.grid.major=element_line(linewidth=.25),
  plot.margin=margin(8,8,8,8))

design <- function(d) {
  names(d) <- trimws(names(d))
  if (!"Block" %in% names(d) && "ID" %in% names(d)) {
    d$Block <- str_extract(d$ID,"(?i)(?<=R)\\s*[0-9]+")
  }
  d %>% mutate(
  Treatment=factor(str_squish(Treatment),levels=treat_order), Block=factor(Block),
  PGPR=factor(case_when(
    Treatment %in% c("Healthy control","Pathogen only") ~ "None",
    Treatment %in% c("PGPR A only","PGPR A + Pathogen") ~ "A",
    Treatment %in% c("PGPR B only","PGPR B + Pathogen") ~ "B",
    TRUE ~ "A+B"),levels=c("None","A","B","A+B")),
  Pathogen=factor(if_else(Treatment %in% with_path,"With Pathogen","No Pathogen"),
                  levels=c("No Pathogen","With Pathogen")))
}
read_tab <- function(f) design(read.delim(file.path(project_dir,f),check.names=FALSE,
                                          stringsAsFactors=FALSE,na.strings=c("","NA","N/A",".")))
read_letters <- function(f) read_csv(file.path(stats_dir,f),show_col_types=FALSE) %>%
  mutate(Letters=str_remove_all(.group,"\\s+")) %>% select(-.group)

save_all <- function(p,name,w,h) {
  ggsave(file.path(fig_dir,paste0(name,".tiff")),p,width=w,height=h,dpi=600,
         compression="lzw",bg="white",limitsize=FALSE)
  ggsave(file.path(fig_dir,paste0(name,".png")),p,width=w,height=h,dpi=600,
         bg="white",limitsize=FALSE)
  ggsave(file.path(fig_dir,paste0(name,".pdf")),p,width=w,height=h,
         bg="white",limitsize=FALSE)
}

raw_summary <- function(d,response,groups=c("Treatment","PGPR","Pathogen")) d %>%
  filter(!is.na(.data[[response]])) %>% group_by(across(all_of(groups))) %>%
  summarise(Mean=mean(.data[[response]]),SD=sd(.data[[response]]),n=n(),
            SE=SD/sqrt(n),.groups="drop")

bar_plot <- function(s,title,ylabel,letters=NULL,percent=FALSE,ci=FALSE) {
  if (!is.null(letters)) s <- left_join(s,letters,by=intersect(names(s),names(letters)))
  if (!"Letters" %in% names(s)) s$Letters <- ""
  s <- s %>% mutate(Letters=replace_na(Letters,""),
                    Treatment=factor(Treatment,levels=treat_order),
                    low=if(ci) Lower else pmax(0,Mean-SE),
                    high=if(ci) Upper else Mean+SE)
  yr <- max(s$high,na.rm=TRUE); s <- s %>% mutate(label_y=high+.08*max(yr,1))
  ggplot(s,aes(Treatment,Mean,fill=Treatment))+
    geom_col(width=.72,color="black",linewidth=.35)+
    geom_errorbar(aes(ymin=low,ymax=high),width=.18,linewidth=.65)+
    geom_text(aes(y=label_y,label=Letters),size=5,fontface="bold")+
    facet_grid(.~Pathogen,scales="free_x",space="free_x")+
    scale_fill_manual(values=treat_cols,breaks=treat_order,drop=FALSE)+
    labs(title=title,x=NULL,y=ylabel,fill="Treatment")+theme_pub+
    theme(legend.position="none")+
    coord_cartesian(ylim=c(0,if(percent) 115 else yr*1.27),clip="off")
}

line_plot <- function(s,title,ylabel) ggplot(s,aes(Date,Mean,group=Treatment,color=Treatment))+
  geom_line(linewidth=1.15)+geom_point(size=3)+
  geom_errorbar(aes(ymin=Mean-SE,ymax=Mean+SE),width=.10,linewidth=.65)+
  facet_wrap(~Pathogen,scales="free_y")+
  scale_color_manual(values=treat_cols,breaks=treat_order,drop=FALSE)+
  labs(title=title,x="Date",y=ylabel,color="Treatment")+theme_pub

# Figure 1: establishment, growth, leaf number, SPAD --------------------------
shoot <- read_tab("shoot.count.txt") %>% mutate(
  Emerged=if_else(is.na(`20-Nov`),NA_real_,as.numeric(`20-Nov`>0)),
  Shoots=as.numeric(`27-Nov`))
em_w <- read_csv(file.path(stats_dir,"emergence_observed_Wilson95.csv"),show_col_types=FALSE) %>%
  mutate(Mean=100*proportion,Lower=100*lower_95,Upper=100*upper_95)
em_l <- read_letters("emergence_binomial_letters.csv") %>%
  select(PGPR,Pathogen,Letters)
pA1 <- bar_plot(em_w,"A. Plant establishment response","Emergence (%)",em_l,TRUE,TRUE)
sh_s <- raw_summary(shoot,"Shoots")
sh_l <- read_letters("shoot_count_letters.csv") %>% select(PGPR,Pathogen,Letters)
pA2 <- bar_plot(sh_s,NULL,"Shoot no./tuber",sh_l)
pA <- pA1/pA2+plot_layout(heights=c(1,1))

growth_dates <- c("20-Nov","27-Nov","03-Dec","10-Dec","17-Dec","23-Dec")
height <- read_tab("Shoot_length.txt") %>% pivot_longer(all_of(growth_dates),names_to="Date",values_to="Value") %>%
  mutate(Date=factor(Date,levels=growth_dates,ordered=TRUE))
height_s <- raw_summary(height,"Value",c("Treatment","PGPR","Pathogen","Date"))
pB <- line_plot(height_s,"B. Shoot height dynamics","Shoot height (cm)")+theme(legend.position="none")
leaves <- read_tab("Leaves_number.txt") %>% pivot_longer(all_of(growth_dates),names_to="Date",values_to="Value") %>%
  mutate(Date=factor(Date,levels=growth_dates,ordered=TRUE))
leaves_s <- raw_summary(leaves,"Value",c("Treatment","PGPR","Pathogen","Date"))
pC <- line_plot(leaves_s,"C. Leaf number dynamics","Leaf number")+theme(legend.position="none")
spad_dates <- c("03-Dec","10-Dec","23-Dec")
spad <- read_tab("chlorophyll.txt") %>% pivot_longer(all_of(spad_dates),names_to="Date",values_to="Value") %>%
  mutate(Date=factor(Date,levels=spad_dates,ordered=TRUE))
spad_s <- raw_summary(spad,"Value",c("Treatment","PGPR","Pathogen","Date"))
pD <- line_plot(spad_s,"D. Chlorophyll index (SPAD) dynamics","Chlorophyll index (SPAD)")
fig1 <- (pA|pB)/(pC|pD)+plot_layout(guides="collect",widths=c(1.25,1),heights=c(1.15,1)) &
  theme(legend.position="right")
save_all(fig1,"Figure_1_Establishment_Growth_Physiology_CORRECTED",20,16)

# Figure 2: flowering (descriptive only) and tuber number ---------------------
flower_w <- read_csv(file.path(stats_dir,"flowering_observed_Wilson95.csv"),show_col_types=FALSE) %>%
  mutate(Treatment=factor(Treatment,levels=treat_order),Mean=100*proportion,
         Lower=100*lower_95,Upper=100*upper_95)
pFlower <- ggplot(flower_w,aes(Treatment,Mean,fill=Treatment))+
  geom_col(width=.72,color="black",linewidth=.35)+
  geom_errorbar(aes(ymin=Lower,ymax=Upper),width=.18,linewidth=.65)+
  facet_grid(Date~Pathogen,scales="free_x",space="free_x")+
  scale_fill_manual(values=treat_cols,drop=FALSE)+
  labs(title="A. Flowering incidence",x=NULL,y="Pots flowered (%)")+theme_pub+
  theme(legend.position="none")+coord_cartesian(ylim=c(0,105))
tuber <- read_tab("Tuber.txt") %>% mutate(Tubers=as.numeric(Tubers))
tub_s <- raw_summary(tuber,"Tubers")
tub_l <- read_letters("tuber_count_letters.csv") %>% select(PGPR,Pathogen,Letters)
pTuber <- bar_plot(tub_s,"B. Tuber number","Tubers per plant",tub_l)
fig2 <- pFlower/pTuber+plot_layout(heights=c(1.15,1))
save_all(fig2,"Figure_2_Flowering_Tuber_CORRECTED",16,14)

# Figure 3: biomass and root architecture, same 2-column x 3-row layout -------
fresh <- read_tab("fresh.biomass.txt")
dry <- read_tab("dry_biomass_updated.txt")
root <- read_tab("winrhizodata.txt"); names(root) <- make.names(names(root))
trait_plot <- function(d,response,tag,title,ylabel) {
  s <- raw_summary(d,response)
  l <- read_letters(paste0(tag,"_letters.csv")) %>% select(PGPR,Pathogen,Letters)
  bar_plot(s,title,ylabel,l)
}
bp <- list(
  trait_plot(fresh,"total_fresh_weight_gm","fresh_biomass","A. Total fresh biomass","Fresh biomass (g/plant)"),
  trait_plot(dry,"dry_whole_biomass_gm","dry_biomass","B. Total dry biomass","Dry biomass (g/plant)"),
  trait_plot(dry,"root_dry_weight_gm","root_dry_weight","C. Root dry weight","Root dry weight (g/plant)"),
  trait_plot(fresh,"leaf_area_cm2","leaf_area","D. Leaf area","Leaf area (cm2/plant)"),
  trait_plot(root,"TotalLength.cm.","root_total_length","E. Total root length","Root length (cm/plant)"),
  trait_plot(root,"SurfArea.cm2.","root_surface_area","F. Root surface area","Root surface area (cm2/plant)"))
fig3 <- (bp[[1]]|bp[[2]])/(bp[[3]]|bp[[4]])/(bp[[5]]|bp[[6]])+plot_layout(guides="collect") &
  theme(legend.position="right")
save_all(fig3,"Figure_3_Biomass_RootArchitecture_CORRECTED",20,21)

# Figure 4: AUDPC --------------------------------------------------------------
aud <- read_csv(file.path(stats_dir,"AUDPC_per_pot_complete_trajectories.csv"),show_col_types=FALSE) %>% design()
aud_s <- raw_summary(aud,"AUDPC")
aud_l <- read_letters("AUDPC_infected_stratum_letters.csv") %>%
  mutate(Pathogen=factor("With Pathogen",levels=c("No Pathogen","With Pathogen"))) %>%
  select(PGPR,Pathogen,Letters)
# No-pathogen values have zero variance; display descriptively without letters.
pAUD <- bar_plot(aud_s,"AUDPC response to PGPR and pathogen treatments",
                 "AUDPC (disease severity x days)",aud_l)
save_all(pAUD,"Figure_4_AUDPC_CORRECTED",12,7)

writeLines(c(
  "Corrected figure set generated from amended models.",
  "Flowering is descriptive and has no significance letters because the repeated binomial model did not converge.",
  "Letters compare PGPR treatments only within the same pathogen condition (and date where applicable).",
  "Root dry-weight figure uses dry_biomass_updated.txt."
),file.path(fig_dir,"FIGURE_README.txt"))
message("Corrected figures saved to: ",fig_dir)
