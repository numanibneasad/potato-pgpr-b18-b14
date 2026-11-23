# ROBUST REANALYSIS: B14 + B18 seed germination
# Same publication figure layout as the original script:
#   A. Germination percentage (line plot, strain facets)
#   B. Maguire germination index (bar plot, strain facets)
# Statistical changes only: repeated beta-binomial mixed model for cumulative
# germination and factorial HC3-robust inference for the Maguire index.

required <- c("readr","dplyr","tidyr","ggplot2","stringr","purrr",
              "janitor","patchwork","glmmTMB","emmeans","multcomp","sandwich","car")
missing_pkgs <- setdiff(required,rownames(installed.packages()))
if(length(missing_pkgs)) stop("Install required packages: ",paste(missing_pkgs,collapse=", "))
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2)
  library(stringr); library(purrr); library(janitor); library(patchwork)
  library(glmmTMB); library(emmeans); library(multcomp); library(sandwich)
})

N_SEEDS <- 25
out_dir <- "seed_germination_robust_output"
dir.create(out_dir,showWarnings=FALSE,recursive=TRUE)
files <- tibble(file=c("seed_B14.txt","seed_ B18.txt"),
                strain=c("P. sp. strain_B14","P. sp. strain_B18"))
strain_levels <- c("P. sp. strain_B18","P. sp. strain_B14")
cfu_levels <- c("Control","10^5","10^6","10^7")
cfu_labels <- c("Control"="Control","10^5"=expression(10^5),
                "10^6"=expression(10^6),"10^7"=expression(10^7))

read_seed <- function(file,strain_name) {
  read_tsv(file,show_col_types=FALSE,name_repair="unique_quiet") %>% clean_names() %>%
    dplyr::select(any_of(c("treatment","concentration","day_1","day_2","day_3"))) %>%
    mutate(across(c(day_1,day_2,day_3),as.numeric),
           strain=strain_name,
           replicate=row_number(),
           cfu_ml=case_when(
             str_squish(as.character(concentration)) %in% c("0","00") ~ "Control",
             str_squish(as.character(concentration))=="OD1" ~ "10^5",
             str_squish(as.character(concentration))=="OD2" ~ "10^6",
             str_squish(as.character(concentration))=="OD3" ~ "10^7",
             TRUE ~ str_squish(as.character(concentration)))) %>%
    group_by(cfu_ml) %>% mutate(replicate_within_cfu=row_number()) %>% ungroup() %>%
    mutate(plate_id=factor(paste(strain,cfu_ml,replicate_within_cfu,sep="__")),
           strain=factor(strain,levels=strain_levels),
           cfu_ml=factor(cfu_ml,levels=cfu_levels))
}

raw <- map2_dfr(files$file,files$strain,read_seed)
if(any(raw[c("day_1","day_2","day_3")]<0,na.rm=TRUE) ||
   any(raw[c("day_1","day_2","day_3")]>N_SEEDS,na.rm=TRUE))
  stop("Germination counts outside 0..N_SEEDS")

invalid <- raw %>% filter(is.na(day_1)|is.na(day_2)|is.na(day_3)|
                            day_2<day_1|day_3<day_2)
write_csv(invalid,file.path(out_dir,"excluded_invalid_seed_trajectories.csv"))
valid <- anti_join(raw,invalid %>% dplyr::select(plate_id),by="plate_id")
if(!nrow(valid)) stop("No valid plates remain")

calc <- valid %>% mutate(
  germ_pct_d1=100*day_1/N_SEEDS,germ_pct_d2=100*day_2/N_SEEDS,
  germ_pct_d3=100*day_3/N_SEEDS,
  n1=day_1,n2=day_2-day_1,n3=day_3-day_2,
  maguire_gri=n1+n2/2+n3/3)

germ_long <- calc %>%
  dplyr::select(plate_id,strain,cfu_ml,day_1,day_2,day_3) %>%
  pivot_longer(c(day_1,day_2,day_3),names_to="day",values_to="germinated") %>%
  mutate(day=factor(day,levels=c("day_1","day_2","day_3"),
                    labels=c("Day 1","Day 2","Day 3"),ordered=TRUE),
         not_germinated=N_SEEDS-germinated)

# Repeated cumulative binomial response with plate random intercept. A
# beta-binomial conditional distribution accounts for extra-binomial variation.
germ_model_beta <- suppressWarnings(glmmTMB(
  cbind(germinated,not_germinated) ~ strain*cfu_ml + strain*day + cfu_ml*day + (1|plate_id),
  family=betabinomial(link="logit"),data=germ_long,
  control=glmmTMBControl(optCtrl=list(iter.max=100000,eval.max=100000))))
beta_ok <- germ_model_beta$fit$convergence==0 && isTRUE(germ_model_beta$sdr$pdHess)
if(beta_ok) {
  germ_model <- germ_model_beta
  germ_family <- "beta-binomial mixed model"
} else {
  germ_model <- glmmTMB(
    cbind(germinated,not_germinated) ~ strain*cfu_ml + strain*day + cfu_ml*day + (1|plate_id),
    family=binomial(link="logit"),data=germ_long,
    control=glmmTMBControl(optCtrl=list(iter.max=100000,eval.max=100000)))
  germ_family <- "binomial mixed model (beta-binomial candidate rejected because its Hessian was not positive definite)"
}
germ_ok <- germ_model$fit$convergence==0 && isTRUE(germ_model$sdr$pdHess)
if(!germ_ok) stop("Final germination mixed model did not converge")

capture.output(summary(germ_model),file=file.path(out_dir,"germination_model_summary.txt"))
capture.output(car::Anova(germ_model,type=3),file=file.path(out_dir,"germination_TypeIII_tests.txt"))
saveRDS(germ_model,file.path(out_dir,"germination_model.rds"))
writeLines(germ_family,file.path(out_dir,"germination_selected_model.txt"))

germ_emm <- emmeans(germ_model,~cfu_ml|strain*day,type="response")
germ_pairs <- pairs(germ_emm,adjust="tukey")
germ_letters <- as.data.frame(multcomp::cld(germ_emm,adjust="tukey",Letters=letters)) %>%
  mutate(letters=str_remove_all(.group,"\\s+"))
write_csv(as.data.frame(germ_emm),file.path(out_dir,"germination_EMMeans.csv"))
write_csv(as.data.frame(germ_pairs),file.path(out_dir,"germination_pairwise_Tukey.csv"))
write_csv(germ_letters,file.path(out_dir,"germination_letters.csv"))

sum_germ <- germ_letters %>% transmute(
  strain,cfu_ml,day,n=NA_integer_,mean=100*prob,
  lower=100*asymp.LCL,upper=100*asymp.UCL,letters) %>%
  group_by(strain,day) %>%
  mutate(max_upper=max(upper,na.rm=TRUE),
         letter_index=as.integer(cfu_ml)-1,
         y_letters=if_else(max_upper>90,103-3*letter_index,
                           pmin(103,upper+4+2*letter_index))) %>%
  ungroup() %>% dplyr::select(-max_upper,-letter_index)

# Maguire index: strain x concentration factorial model. HC3 covariance makes
# inference less sensitive to unequal variances and the slightly unbalanced cell.
mgi_model <- lm(maguire_gri~strain*cfu_ml,data=calc)
mgi_vcov <- sandwich::vcovHC(mgi_model,type="HC3")
capture.output(summary(mgi_model),file=file.path(out_dir,"MGI_model_summary.txt"))
capture.output(car::Anova(mgi_model,type=3,vcov.=mgi_vcov,white.adjust=FALSE),
               file=file.path(out_dir,"MGI_TypeIII_HC3_tests.txt"))
saveRDS(mgi_model,file.path(out_dir,"MGI_model.rds"))
mgi_emm <- emmeans(mgi_model,~cfu_ml|strain,vcov.=mgi_vcov)
mgi_pairs <- pairs(mgi_emm,adjust="tukey")
mgi_letters <- as.data.frame(multcomp::cld(mgi_emm,adjust="tukey",Letters=letters)) %>%
  mutate(letters=str_remove_all(.group,"\\s+"))
write_csv(as.data.frame(mgi_emm),file.path(out_dir,"MGI_EMMeans_HC3.csv"))
write_csv(as.data.frame(mgi_pairs),file.path(out_dir,"MGI_pairwise_Tukey_HC3.csv"))
write_csv(mgi_letters,file.path(out_dir,"MGI_letters_HC3.csv"))
sum_gri <- mgi_letters %>% transmute(
  strain,cfu_ml,n=NA_integer_,mean=emmean,lower=lower.CL,upper=upper.CL,letters) %>%
  group_by(strain) %>% mutate(y_letters=upper+.08*max(upper,na.rm=TRUE)) %>% ungroup()

# Exact original theme and panel layout ---------------------------------------
theme_pub <- theme_bw(base_size=18)+theme(
  strip.text=element_text(size=20,face="bold"),
  axis.title=element_text(size=22,face="bold"),
  axis.text.x=element_text(size=20,face="bold"),axis.text.y=element_text(size=18),
  legend.title=element_text(size=22,face="bold"),legend.text=element_text(size=19),
  plot.title=element_text(size=22,face="bold"),panel.grid.minor=element_blank(),
  panel.grid.major=element_line(linewidth=.3),plot.margin=margin(8,8,8,8))

p_germ <- ggplot(sum_germ,aes(day,mean,group=cfu_ml,color=cfu_ml,shape=cfu_ml))+
  geom_line(linewidth=1.4)+geom_point(size=4.4)+
  geom_errorbar(aes(ymin=lower,ymax=upper),width=.12,linewidth=.8)+
  geom_text(aes(label=letters,y=y_letters),size=3.5,fontface="bold",show.legend=FALSE)+
  facet_wrap(~strain,nrow=1)+scale_y_continuous(limits=c(0,106),breaks=seq(0,100,20))+
  scale_color_discrete(name="CFU/mL",labels=cfu_labels)+
  scale_shape_discrete(name="CFU/mL",labels=cfu_labels)+
  labs(title="A. Germination percentage",x="Germination time",y="Germination (%)")+theme_pub

p_gri <- ggplot(sum_gri,aes(cfu_ml,mean,fill=cfu_ml))+
  geom_col(width=.72,color="black",linewidth=.5)+
  geom_errorbar(aes(ymin=lower,ymax=upper),width=.20,linewidth=.8)+
  geom_text(aes(label=letters,y=y_letters),size=5.5,fontface="bold")+
  facet_wrap(~strain,nrow=1)+scale_fill_discrete(name="CFU/mL",labels=cfu_labels)+
  scale_x_discrete(labels=cfu_labels)+
  labs(title="B. Maguire germination index",x="CFU/mL",y="MGI")+theme_pub

combined_plot <- p_germ/p_gri+plot_layout(heights=c(1.15,1),guides="collect") &
  theme(legend.position="right")
print(combined_plot)
ggsave(file.path(out_dir,"combined_germination_MGI_B14_B18_ROBUST.png"),combined_plot,
       width=15,height=11,dpi=300,bg="white")
ggsave(file.path(out_dir,"combined_germination_MGI_B14_B18_ROBUST.tiff"),combined_plot,
       width=15,height=11,dpi=300,compression="lzw",bg="white")
ggsave(file.path(out_dir,"combined_germination_MGI_B14_B18_ROBUST.pdf"),combined_plot,
       width=15,height=11,bg="white")
write_csv(sum_germ,file.path(out_dir,"summary_germination_model_based.csv"))
write_csv(sum_gri,file.path(out_dir,"summary_MGI_model_based_HC3.csv"))
writeLines(capture.output(sessionInfo()),file.path(out_dir,"sessionInfo.txt"))
message("Robust seed analysis completed: ",normalizePath(out_dir,winslash="/"))
