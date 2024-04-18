#!/usr/bin/env Rscript

#Rationale: locus prioritisation with Open Target Platform for Asthma - MONDO0004979

library(tidyverse)
library(data.table)
library(readxl)


#input each analysis table:
file_xlsx <- "input/edited_OT-MONDO_0004979-associated-targets-14_04_2024-v24_03.xlsx"
open_target <- read_excel(file_xlsx, sheet = "EFO_MONDO0004979")
locus_gene <- read_excel(file_xlsx, sheet = "Locus2gene_mapping")
locus_gene_in_open_target <- locus_gene %>%
                             left_join(open_target, by = "symbol") %>%
                             relocate(chembl, .after = globalScore) %>%
                             relocate(maxClinicalTrialPhase, .after = chembl)
#save full result:
write_tsv(locus_gene_in_open_target,"output/Locus_prioritisation_full.tsv")

#digest result:
#Summary statistics:
##How many genes in OpenTarget for each locus:
##grep by locus, find the count of gene with globalscore NA and not NA.
locus_gene_in_open_target <- locus_gene_in_open_target %>%
                             group_by(Locus) %>%
                             mutate(Ngene_TOT = n()) %>%
                             ungroup()

locus_Ngene_NOT_OTP  <- locus_gene_in_open_target %>%
    dplyr::group_by(Locus) %>%
    dplyr::summarise_at(vars(globalScore), ~ sum(is.na(.))) %>% rename(Ngene_NOT_OPT = globalScore)

locus_gene_in_open_target <- locus_gene_in_open_target %>% left_join(locus_Ngene_NOT_OTP, by = "Locus")
locus_gene_in_open_target$Ngene_IN_OPT <- locus_gene_in_open_target$Ngene_TOT - locus_gene_in_open_target$Ngene_NOT_OPT
locus_gene_in_open_target$Ngene_IN_OPT_perc <- locus_gene_in_open_target$Ngene_IN_OPT / locus_gene_in_open_target$Ngene_TOT

locus_gene_in_open_target_sumstats <- locus_gene_in_open_target %>%
                                      select(Locus,Ngene_TOT, Ngene_IN_OPT, Ngene_IN_OPT_perc) %>%
                                      unique()

#Locus NOT TO prioritise: locus with gene with precedence for clinical trial in asthma:
locus_with_asthma_clinicaltrial <- locus_gene_in_open_target %>%
                                   filter(chembl != "No data") %>%
                                   select(Locus, symbol, globalScore, chembl, maxClinicalTrialPhase)
locus_with_asthma_clinicaltrial$prioritise <- as.factor("NO")

#Locus TO BE PRIORITISED: locus with NO gene with precedence for clinical trial in asthma:
locus_NOasthma_clinicaltrial <- locus_gene_in_open_target %>%
                                filter(! Locus %in% locus_with_asthma_clinicaltrial$Locus) %>%
                                filter(chembl == "No data")
locus_NOasthma_clinicaltrial <- locus_NOasthma_clinicaltrial %>%
                                     mutate(prioritise = if_else(maxClinicalTrialPhase != "No data",
                                     as.factor("POTENTIAL_IN_CLINICS"), as.factor("POTENTIAL_NO_CLINICS"))) %>%
                                     select(Locus, prioritise, symbol, globalScore, chembl, maxClinicalTrialPhase)

#combine for overall table with summary results on prioritisation:
locus_prioritisation <- rbind(locus_with_asthma_clinicaltrial, locus_NOasthma_clinicaltrial) %>%
                        relocate(prioritise, .after = Locus) %>%
                        arrange(Locus, prioritise, symbol)

##add the stats about gene presence for each locus:
locus_prioritisation <- locus_prioritisation %>%
                        full_join(locus_gene_in_open_target_sumstats, by = "Locus") %>%
                        mutate(prioritise = if_else(Ngene_IN_OPT_perc == 0, "NA_NO_GENE_IN_OPT", prioritise)) %>%
                        arrange(Locus, prioritise, symbol)

#save prioritisation result:
write_tsv(locus_prioritisation,"output/Locus_prioritisation_decision.tsv")