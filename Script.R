rm(list=ls())

#set language
Sys.setenv(LANG = "en")

#Setting working directory
setwd("/R working folder - PAD and metabolomics")

# Loading packages --------------------------------------------------------
library(Hmisc)
library(sjlabelled)
library(matrixStats)

library(pastecs)
library(xtable)
library(dplyr)
library(ggplot2)
library(ggpubr)

#Survival analysis

library(survival)
library(cmprsk)
library(survminer)
library(CMAverse)

#Plotting
library(eoffice)
library(stringr)

#Using ggplot
options(grType='RDevice')
library(broom)
library(flextable)
library(webshot)
library(magick)
library(grid)
library(cowplot)
library(CMAverse)

library(gridExtra)
library(gridGraphics)
library(grid)
library(haven)

library(rms)


#Function for P values in string format
signif.p<-function(x, exponent=F){
  
  result<-c()
  
  for(i in 1:length(x)){
    
    if(exponent==T){
      if (x[i]==0){
        result[i]<-as.character("-16")
        
      } else {
        if(log10(x[i])>(-4)){ 
          result[i]<-NA
        }  else {
          
          result[i]<-as.character(floor(log10(signif(x[i],1))))
        } 
      } 
    }else if (exponent==F) {
      
      
      if(log10(x[i])>(-4)){
        result[i]<-as.character(format(signif(x[i],2),scientific = F)
        )
      } else {
        
        if (x[i]==0){
          q<-"<2"
          
        } else {
          
          q<-signif(x[i],1)*10^-round(log10(x[i]))
          
          if (q<1) {
            
            q<-q*10
          }
        }
        
        result[i]<-paste(as.character(q),"x 10")
      }
    } 
  }
  return(result)
}


# Cleaning up data from PAD and metabolomics ------------------------------

cgps<-cgps[!is.na(cgps$sex),]



#saveRDS(cgps, file="L:/LovbeskyttetMapper/Forskning KBA/Benjamin Nilsson Wadström/cgps_imputed_metabolomics_small_large_VLDL_PAD.Rda")



# Setting up data ------------------------------------------------------------
library(haven)
cgps<-readRDS("cgps_imputed_metabolomics_small_large_VLDL_PAD.Rda")

cgps_TRLC<-read_dta("ObusH_TRLC_CVD_020119.dta")
cgps_afli<-read_dta("AFLI_kort0_OBUSH19012023.dta")
cgps_hopkins<-read.csv("CGPS_ldl_hopkins_17022021.csv",header=TRUE)
cgps_whr_timemeal_leuk<-read.csv("OBUSH_waist_hip_time_since_meal_leu.csv",header=TRUE)[,1:5]

cgps<- merge(cgps,cgps_TRLC[,c("obushnr","TRLC")],by="obushnr", all=TRUE)
cgps<- merge(cgps,cgps_afli[,c("obushnr","af")],by="obushnr", all=TRUE)
cgps<- merge(cgps,cgps_whr_timemeal_leuk[,c("obushnr","time_since_meal")],by="obushnr", all=TRUE)

cgps<-cgps[!is.na(cgps$VLDLC),]

#Creating LDLC variable
cgps$LDLC<-cgps$IDLC+cgps$MLDLC+cgps$LLDLC+cgps$SLDLC

#Recovery correcting cholesterol
cgps$chol_ratio<-cgps$chol/cgps$SerumC
cgps$totalC_recov<-cgps$totalC*cgps$chol_ratio

cgps$VLDLC_recov<-cgps$VLDLC*cgps$chol_ratio
cgps$LDLC_recov<- cgps$LDLC*cgps$chol_ratio

cgps$SLDLC_recov<-cgps$SLDLC*cgps$chol_ratio
cgps$MLDLC_recov<-cgps$MLDLC*cgps$chol_ratio
cgps$LLDLC_recov<-cgps$LLDLC*cgps$chol_ratio
cgps$IDLC_recov<-cgps$IDLC*cgps$chol_ratio

cgps$XSVLDLC_recov<-cgps$XSVLDLC*cgps$chol_ratio
cgps$SVLDLC_recov<-cgps$SVLDLC*cgps$chol_ratio
cgps$MVLDLC_recov<-cgps$MVLDLC*cgps$chol_ratio
cgps$LVLDLC_recov<-cgps$LVLDLC*cgps$chol_ratio
cgps$XLVLDLC_recov<-cgps$XLVLDLC*cgps$chol_ratio
cgps$XXLVLDLC_recov<-cgps$XXLVLDLC*cgps$chol_ratio

#Creating small and large VLDL variables
cgps$VLDLC_S<-cgps$XSVLDLC_recov+cgps$SVLDLC_recov
cgps$VLDLC_L<-cgps$VLDLC_recov-cgps$VLDLC_S



cgps<-cgps[!is.na(cgps$VLDLC) & cgps$m5==0,]
cgps$VLDLC_L_S<-with(cgps,interaction(factor(cut2(VLDLC_S,quantile(VLDLC_S,0.5,na.rm=T)),labels=c("low VLDLC_S", "high VLDLC_S")),factor(cut2(VLDLC_L,quantile(VLDLC_L,0.5,na.rm=T)),labels=c("low VLDLC_L", "high VLDLC_L"))))



#mg per dL conversion
cgps$VLDLC_Smgdl<-cgps$VLDLC_S*38.67
cgps$VLDLC_Lmgdl<-cgps$VLDLC_L*38.67
cgps$VLDLC_recovmgdl<-cgps$VLDLC_recov*38.67
cgps$LDLC_recovmgdl<-cgps$LDLC_recov*38.67

cgps$ldlmgdl<-cgps$ldl*38.67
cgps$remncholcalcmgdl<-cgps$remncholcalc*38.67
cgps$trigmgdl<-cgps$trig*88.50



# Converting data- obush2 for RDR -----------------------------------------

#For regression dilution bias
cgps2_metabolomics<-read.csv("obush2_reg_dil.csv")


for(q in c("","OBUSH1_")){
  cgps2_metabolomics[,paste0(q,"LDLC")]<-(cgps2_metabolomics[paste0(q,"SLDLC")]+cgps2_metabolomics[paste0(q,"MLDLC")]+cgps2_metabolomics[paste0(q,"LLDLC")]+cgps2_metabolomics[paste0(q,"IDLC")])
  cgps2_metabolomics[,paste0(q,"VLDLC")]<-(cgps2_metabolomics[paste0(q,"XXLVLDLC")]+cgps2_metabolomics[paste0(q,"XLVLDLC")]+cgps2_metabolomics[paste0(q,"LVLDLC")]+cgps2_metabolomics[paste0(q,"MVLDLC")]+cgps2_metabolomics[paste0(q,"SVLDLC")]+cgps2_metabolomics[paste0(q,"XSVLDLC")])
  
  cgps2_metabolomics[,paste0(q,"VLDLC_S")]<-(cgps2_metabolomics[paste0(q,"SVLDLC")]+cgps2_metabolomics[paste0(q,"XSVLDLC")])
  cgps2_metabolomics[,paste0(q,"VLDLC_L")]<-(cgps2_metabolomics[paste0(q,"XXLVLDLC")]+cgps2_metabolomics[paste0(q,"XLVLDLC")]+cgps2_metabolomics[paste0(q,"LVLDLC")]+cgps2_metabolomics[paste0(q,"MVLDLC")])
  
  
  }


cgps2_metabolomics<- merge(cgps[c("obushnr","m5","remncholcalc","ldl")],cgps2_metabolomics,by="obushnr", all=TRUE)



lambda.VLDLC_S<-lm(log(VLDLC_S)~log(OBUSH1_VLDLC_S),data=cgps2_metabolomics,subset=m5==0 & statin2==0)$coefficients[2]
lambda.VLDLC_L<-lm(VLDLC_L~OBUSH1_VLDLC_L,data=cgps2_metabolomics,subset=m5==0 & statin2==0)$coefficients[2]
lambda.LDLC_recov<-lm(LDLC~OBUSH1_LDLC,data=cgps2_metabolomics,subset=m5==0 & statin2==0)$coefficients[2]


lambda.apob_mgdl<-0.74
lambda.remncholcalc<-0.48
lambda.ldl<-0.60
lambda.hscrp_log<-0.58 # From article "C-reactive protein concentration and risk of coronary heart disease, stroke, and mortality: an individual participant meta-analysis."


# For lipoprotein subfractions
reg_dil_ratio_all<-read_dta("reg_dil_ratio_all.dta")

reg_dil_ratio_all<-as.data.frame(reg_dil_ratio_all)
reg_dil_ratio_all<-rbind(reg_dil_ratio_all,data.frame(parm=c("apob_mgdl","VLDLC_S","VLDLC_L","hscrp_corr","LDLC_recov"),rdb=c(0.74,lambda.VLDLC_S,lambda.VLDLC_L,lambda.hscrp_log,lambda.LDLC_recov)))

rm("cgps2_metabolomics")



# setting data distribution  ------------------------------------------------------------


#Describing variable distributions to rms package
dd<- datadist(cgps)

dd$limits$VLDLC_S[2]<- quantile(subset(cgps,m5==0)$VLDLC_S,0.001,na.rm=T)
dd$limits$VLDLC_L[2]<- min(subset(cgps,m5==0)$VLDLC_L,0.001,na.rm=T)
dd$limits$LDLC_recov[2]<- 0.5

options(datadist = "dd")
options(contrasts=c("contr.treatment", "contr.treatment"))

options(na.action = na.omit)


quantile(subset(cgps,m5==0)$VLDLC_S,0.001,na.rm=T)


# Descriptive tables ------------------------------------------------------

#Preparation
library(tangram)

EHJ<-nejm
EHJ[['Cell']][["n"]]<-function (n, class = NULL, hdr = FALSE, possible = NULL, ...) 
{
  n <- formatC(n, format="d", big.mark=",")
  cell(n, class = c("cell_n", class), ...)
}
EHJ[['Cell']][['fraction']] <- function(numerator, denominator, format=3, ...)
{ paste0(formatC(numerator, format="d", big.mark=",")," (",render_f(round(100*numerator/denominator, 0), 0),'%)') }
EHJ[['Cell']][['iqr']] <-function (x, format = NA, na.rm = TRUE, names = FALSE, type = 8, 
                                   msd = FALSE, quant = c(0.25, 0.5, 0.75), ...) {
  if (length(quant)%%2 == 0) 
    stop("hmisc_iqr quant argument must be an odd length")
  m <- median(1:length(quant))
  y <- quantile(x, quant, na.rm=TRUE, names, type)
  
  if (median(x, na.rm=T)<10) {
    y<-round(y,1)
  } else {
    y<-round(y,0)
  }
  
  if (is.na(format)) 
    format <- format_guess(y)
  ql <- "-"
  if (sum(!is.na(y)) > 0) {
    ql <- sapply(y, function(x) as.character(x))
    ql <- paste0(ql[m], " (", paste0(ql[1:(m - 1)], 
                                     collapse = "-"), "-", paste0(ql[(m + 
                                                                        1):length(quant)], collapse = "-"), ")")
    if (msd) 
      ql <- paste0(ql, "  ", render_f(mean(x, na.rm = TRUE), 
                                      format), "?", render_f(sd(x, na.rm = TRUE), 
                                                             format))
  }
  cell(ql, ...)
}

EHJ[['Numerical']][['Categorical']] <-function (table, row, column, cell_style, pformat = NULL, msd = FALSE, 
                                                quant = c(0.25, 0.5, 0.75), overall = NULL, test = FALSE, 
                                                useNA = "no", ...) {
  overall_label <- if (is.character(overall)) 
    overall
  else "Overall"
  overall <- column$value != "1" && (isTRUE(overall) || 
                                       is.character(overall))
  datar <- row$data
  datac <- as.categorical(column$data)
  categories <- if (overall) 
    c(levels(datac), overall_label)
  else levels(datac)
  categories <- if (length(categories) == 1) 
    overall_label
  else categories
  format <- ifelse(is.na(row$format), format_guess(datar), 
                   row$format)
  useNA <- useNA == "always" || (sum(is.na(datar)) > 
                                   0 && useNA == "ifany")
  subN <- lapply(levels(datac), FUN = function(cat) {
    cell_style[["n"]](length(datac[datac == cat & !is.na(datac)]), 
                      subcol = cat, hdr = TRUE, possible = length(datac), 
                      ...)
  })
  if (overall) 
    subN[[length(subN) + 1]] <- cell_style[["n"]](sum(!is.na(column$data)), 
                                                  hdr = TRUE, subcol = "Overall", possible = length(column$data), 
                                                  ...)
  stat <- if (length(categories) == 1) {
    ""
  }
  else {
    tst <- suppressWarnings(spearman2(datac, datar, na.action = na.retain))
    cell_style[["fstat"]](f = render_f(tst["F"], 
                                       "%.2f"), df1 = tst["df1"], df2 = tst["df2"], 
                          p = cell_style[["p"]](tst["P"], pformat))
  }
  tbl <- table %>% row_header(derive_label(row))
  if (msd) 
    tbl <- tbl %>% row_header("   Mean?SD")
  if (useNA) 
    tbl <- tbl %>% row_header("   Missing (%)")
  tbl <- if (test) {
    col_header(tbl,  categories, "Test Statistic") %>% 
      col_header("", subN, "")
  }
  else {
    col_header(tbl,  categories) %>% col_header("", 
                                                subN)
  }
  
  tbl <- table_apply(tbl, categories, function(tbl, category) {
    x <- if (category == overall_label) 
      datar
    else datar[datac == category]
    tbl <- tbl %>% add_row(cell_style[["iqr"]](x, 
                                               row$format, subcol = category, na.rm = TRUE, msd = FALSE, 
                                               quant = quant))
    if (msd) {
      tbl <- add_row(tbl, cell(paste0(render_f(mean(x, 
                                                    na.rm = TRUE), row$format), "?", render_f(sd(x, 
                                                                                                 na.rm = TRUE), row$format))))
    }
    if (useNA) {
      tbl <- add_row(tbl, cell_style[["fraction"]](sum(is.na(x)), 
                                                   length(x), format = row$format, subcol = category))
    }
    tbl %>% new_col()
  })
  tbl <- home(tbl) %>% cursor_right(length(categories) + 1)
  if (test) 
    tbl <- add_col(tbl, stat)
  tbl
}

EHJ[['Categorical']][['Categorical']]<-function (table, row, column, cell_style, pformat = NULL, collapse_single = TRUE, 
                                                 overall = NULL, test = FALSE, row_percents = FALSE, useNA = "no", ...) 
{
  grid <- table(as.categorical(row$data), as.categorical(column$data), 
                useNA = useNA)
  if (is.na(colnames(grid)[ncol(grid)])) 
    grid <- grid[, 1:(ncol(grid) - 1)]
  validcol <- which(!apply(grid, 2, FUN = function(x) {
    all(x == 0)
  }))
  validrow <- which(!apply(grid, 1, FUN = function(x) {
    all(x == 0)
  }))
  stat <- if (length(validrow) < 2 || length(validcol) < 2) 
    NA
  else suppressWarnings(chisq.test(grid[validrow, validcol], 
                                   correct = FALSE))
  ncol <- dim(grid)[2]
  nrow <- dim(grid)[1]
  denominators <- if (row_percents) 
    matrix(rep(rowSums(grid), ncol), ncol = ncol, byrow = FALSE)
  else matrix(rep(colSums(grid), nrow), ncol = ncol, byrow = TRUE)
  rowlabels <- rownames(grid)
  subN <- lapply(colnames(grid), FUN = function(cat) cell_style[["n"]](sum(column$data == 
                                                                             cat, na.rm = TRUE), subcol = cat, possible = length(column$data), 
                                                                       hdr = TRUE, ...))
  if (!is.null(overall) && (!is.logical(overall) || overall)) {
    denominators <- cbind(denominators, rep(sum(grid), nrow))
    grid <- cbind(grid, rowSums(grid))
    colnames(grid)[ncol + 1] <- if (is.character(overall)) 
      overall
    else "Overall"
    subN[[ncol + 1]] <- cell_style[["n"]](sum(!is.na(column$data)), 
                                          possible = length(column$data), subcol = "Overall", 
                                          hdr = TRUE, ...)
    ncol <- ncol + 1
  }
  if (collapse_single && dim(grid)[1] <= 2) {
    name <- row$name()
    try({
      l2 <- attr(row$data, "label")
      if (!is.null(l2)) {
        name <- l2
      }
    })
    pos <- dim(grid)[1]
    x <- matrix(grid[pos, ], nrow = 1)
    colnames(x) <- colnames(grid)
    rownames(x) <- paste(name, ":", rownames(grid)[pos])
    grid <- x
    denominators <- matrix(denominators[pos, ], nrow = 1)
    nrow <- 1
  }
  else {
    if (is.na(rownames(grid)[nrow(grid)])) {
      rownames(grid)[nrow(grid)] <- "Missing (%)"
    }
    rownames(grid) <- lapply(rownames(grid), FUN = function(x) paste("  ", 
                                                                     x))
  }
  if (inherits(test, "function")) {
    test_result <- test(row, column, cell_style, ...)
    test <- TRUE
  }
  else if (test) {
    test_result <- if (any(is.na(stat))) 
      cell("NA")
    else cell_style[["chi2"]](render_f(stat$statistic, 
                                       2), stat$parameter, cell_style[["p"]](stat$p.value, 
                                                                             pformat))
  }
  if (test) {
    table <- col_header(table,  colnames(grid), 
                        "Test Statistic")
    table <- col_header(table, subN, "")
  }
  else {
    table <- col_header(table,  colnames(grid))
    table <- col_header(table,  subN)
  }
  if (nrow > 1) 
    table <- row_header(table, derive_label(row, ...))
  for (nm in rownames(grid)) table <- row_header(table, nm)
  for (j in 1:ncol) {
    if (nrow > 1) 
      table <- add_row(table, "")
    format <- if (is.na(row$format) || is.null(row$format)) 
      format_guess(as.vector(grid/denominators))
    else row$format
    for (i in 1:nrow) {
      table <- if (denominators[i, j] == 0) 
        add_row(table, "")
      else add_row(table, cell_style[["fraction"]](grid[i, 
                                                        j], denominators[i, j], format = format, subcol = colnames(grid)[i], 
                                                   subrow = rownames(grid)[j]))
    }
    table <- new_col(table)
  }
  if (test) {
    table <- add_row(table, test_result)
    if (nrow > 1) 
      table <- add_row(table, rep("", nrow))
  }
  table
}




# All individuals

all <- tangram("factor(VLDLC_L_S)~ sex + startage+ VLDLC_S+VLDLC_Smgdl+VLDLC_L+VLDLC_Lmgdl+LDLC_recov+LDLC_recovmgdl+factor(current_smoking)+ cum_smoking+ systolic+diastolic+factor(af)+ bmi+ factor(dmall_pre)+factor(ascvd_pre)+VLDLC_recov+VLDLC_recovmgdl+trig+trigmgdl+remncholcalc+ remncholcalcmgdl+ ldl+ ldlmgdl+apob_mgdl",
               data = cgps, 
               test=FALSE, id="override", transforms = EHJ, msd=FALSE, digits = 1)

csv(all, file="Temp.csv")


# All individuals diagnoses

all <- tangram("ASCVD_infl~ factor(AMI)+factor(IS)+factor(PAD)",
               data = cgps, 
               test=FALSE, id="override", transforms = EHJ, msd=FALSE, digits = 1)
all

csv(all, file="Temp.csv")





#Table 1

cgps$apob_mgdl_cut2<-cut2(cgps$apob_mgdl,100)
table(cgps$apob_mgdl_cut2)


remnchol_desc <- tangram("factor(apob_mgdl_cut2) ~ sex + startage+ remncholcalc+ remncholcalcmgdl+ ldl+ ldlmgdl+VLDLP+LDLP+apob_mgdl+ systolic+diastolic+factor(current_smoking)+ cum_smoking+ bmi+ factor(ascvd_pre)",
                         data = cgps, 
                         test=FALSE, id="override", transforms = EHJ, msd=FALSE, digits = 1)

remnchol_desc

csv(remnchol_desc, file="Temp.csv")

#Table S1

aggregate(cum_smoking~m5, data=subset(cgps, dmall_pre==1 & ever_smoking==T), FUN=function(x) c(median=median(x),iqrlow=quantile(x,0.25),iqrhigh= quantile(x,0.75), count=length(x)))


m5_desc <- tangram("factor(m5,levels=c(1,0)) ~ sex + startage+ remncholcalc+ remncholcalcmgdl+ ldl+ ldlmgdl+ systolic+ ppressure+hscrp+factor(current_smoking)+ cum_smoking+ bmi+ factor(ascvd_pre)",
                   data = cgps, 
                   test=FALSE, id="override", transforms = EHJ, msd=FALSE, digits = 1)


csv(m5_desc, file="Temp.csv")

#Table S2

aggregate(cum_smoking~current_smoking, data=subset(cgps, dmall_pre==1 & ever_smoking==T), FUN=function(x) c(median=median(x),iqrlow=quantile(x,0.25),iqrhigh= quantile(x,0.75), count=length(x)))


m5_desc <- tangram("factor(current_smoking,levels=c(1,0)) ~ sex + startage+ remncholcalc+ remncholcalcmgdl+ ldl+ ldlmgdl+ factor(m5)+ systolic+ ppressure+hscrp+ cum_smoking+ bmi+ factor(ascvd_pre)",
                   data = cgps, 
                   test=FALSE, id="override", transforms = EHJ, msd=FALSE, digits = 1)


csv(m5_desc, file="Temp.csv")


#Table S3
dm_desc <- tangram("factor(dm_pre) ~ sex + startage+ remncholcalc+ remncholcalcmgdl+ ldl+ ldlmgdl+ factor(m5)+systolic+ ppressure+hscrp+factor(current_smoking)+ cum_smoking+ bmi+ factor(ascvd_pre)",
                   data = cgps, 
                   test=FALSE, id="override", transforms = EHJ, msd=FALSE, digits = 1)



csv(dm_desc, file="Temp.csv")

#Table S4
endpoint_descb_all <- tangram("1 ~ PAD + c_diag",
                              data = cgps, 
                              test=FALSE, id="override", transforms = EHJ, msd=FALSE, digits = 1)

csv(endpoint_descb_all, file="Endpoint characteristics.csv")


endpoint_descb <- tangram("factor(current_smoking,levels=c(1,0)) ~ PAD + c_diag",
                          data = cgps, 
                          test=FALSE, id="override", transforms = EHJ, msd=FALSE, digits = 1)



csv(endpoint_descb, file="Endpoint characteristics.csv")



# Table S6
aggregate(cum_smoking~PAD, data=subset(cgps, stopagePAD>startage)[cgps$dmall_pre==FALSE & cgps$ever_smoking==TRUE, ], FUN=function(x) c(median=median(x),iqrlow=iqr(x)/2+median(x),iqrhigh= median(x)-iqr(x)/2, count=length(x)))



PAD_desc <- tangram("PAD ~ sex + startage+ remncholcalc+ remncholcalcmgdl+ ldl+ ldlmgdl+ systolic+ ppressure+hscrp+factor(m5)+factor(current_smoking)+ cum_smoking+ bmi+ factor(ascvd_pre)",
                    data = cgps, 
                    PAD_diabetes_desctest=FALSE, id="override", transforms = EHJ, msd=FALSE, digits = 1)

csv(PAD_desc, file="Temp.csv")

# Extra table, diabetes diagnoses  (deleted)

cgps$glucosemgdl<-cgps$glucose*18


cgps$gluc11<-cut2(cgps$glucose,11.1)


diabetes_desc <- tangram("factor(m5,levels=c(1,0)) ~ factor(dmi_pre) + factor(dmii_pre) + factor(v27) + factor(m11) + factor(m12) + factor(gluc11)",
                         data = cgps, 
                         test=FALSE, id="override", transforms = EHJ, msd=FALSE, digits = 1)
diabetes_desc
csv(diabetes_desc, file="diabetes_desc.csv")





# Figure 1: Correlation matrix and scatterplots particle number and cholesterol content -------

#Correlation matrix

library(reshape2)
library(stringr)
corr_list<-list()

corrmatrix_P<-round(cor(subset(cgps,m5==0)[c(paste0(c("SLDL","MLDL","LLDL","IDL","XSVLDL","SVLDL","MVLDL","LVLDL","XLVLDL","XXLVLDL"),"C_recov"),"trig","remncholcalc","ldl","apob_mgdl")], use = "complete",method="spearman")^2,2)*100


get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)}

corrmatrix_P<-get_upper_tri(corrmatrix_P)

corrmatrix_P <- melt(corrmatrix_P, na.rm = TRUE)

Figure_1<-ggplot(corrmatrix_P, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low ="white", high = "red",  
                       limit = c(0,100), space = "Lab", 
                       name="Spearman\nCorrelation^2") +
  theme_nothing()+ 
  coord_fixed() + 
  geom_text(aes(Var2, Var1, label = str_remove(value, "^0+")), color = "black", size = 4) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1),
        axis.text.y = element_text(size=12),
        legend.justification = c(1, 0),
        legend.position = c(0.6, 0.7),
        legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))
Figure_1

toffice(figure = Figure_1, format = "pptx",title= "Copenhagen General Population Study",
        filename = "Temp.pptx", nr = 1, nc = 1, irow = 1, icol = 1,onsame = FALSE,left = 0.15,top = 0.15,
        append = TRUE, width = 8, height = 8, devsize = FALSE, units = "in")


# Figure 2: RCS plots-----------------------------------------------


fig3<-list()
for (k in c("ASCVD","AMI","IS","PAD"))
  for (i in c("VLDLC_S","VLDLC_L")){
    model<-cph(as.formula(paste0("Surv(startage,stopage",k,",",k,")~rcs(",i,",3)+strat(sex)+systolic+diastolic+strat(status_smoking)+rcs(cum_smoking, 3)+birthdate+rcs(bmi,4)+dmall_pre",
                                 if (grepl("VLDL",i)){"+LDLC_recov"} else  {""},if (grepl("IS",k)){"+af"} else  {""})) ,data=cgps,x=TRUE, y=TRUE)  
    

    pred<-Predict(model,ref.zero = TRUE)
    pred[c("yhat","lower","upper")]<-exp(pred[c("yhat","lower","upper")]/reg_dil_ratio_all[reg_dil_ratio_all$parm==i,"rdb"])
    
    
    col<-if (i=="VLDLC_S"){"darkcyan"} else if (i=="VLDLC_L") {"darkorchid4"} else {"darkblue"}
    
    
    gg<-ggplot(pred[1:200, ] ,
               colfill = col)+
      geom_line(size=1, colour=col,linetype=1) +
      geom_hline(aes(yintercept = 1), linetype = 2)+ coord_cartesian(xlim=c(0,1),ylim= c(0.5, 10))+
      labs(x=i , y=paste("Hazard ratio for",k,"(95% CI)"))+scale_y_continuous(expand = c(0,0),trans="log",breaks = c(0.5,1,2,4,8))+scale_x_continuous(expand = c(0,0))+
      theme(axis.line.x.bottom =element_line(),axis.line.y.left =element_line() ,axis.line.x.top=element_blank() ,plot.margin=margin(0,3,0,3),plot.subtitle = element_blank(), plot.title = element_blank(),axis.text.y = element_text(size=10), 
            panel.grid=element_blank(),strip.background = element_rect(color = "white", fill = "white"), panel.background = element_rect(fill = "white"),axis.ticks.length = unit(.2, "cm"))
    
    
    
    gg_dens<-ggplot(data=cgps) + 
      geom_density(aes_string(x=i),color=col, fill=col, alpha=0.4)+
      theme_nothing()+scale_x_continuous(limits=c(0,1))
    
    
    fig3[paste0(i,k)]<- list(plot_grid(if(k=="ASCVD"){gg_dens} else {ggplot()+theme_nothing()},
                                       gg,ncol=1,nrow=2,
                                       rel_heights = if(k=="ASCVD"){c(0.3,1)}else {c(0.1,1)}))   
    
    print(c(paste(i,k),model$stats[1:2]))}


#Changing function plot_grid
fargs  <- formals(plot_grid)
fargs$ncol           <- 2
formals(plot_grid) <- fargs


figure3<-do.call(plot_grid,fig3)




toffice(figure = figure3, format = "pptx",
        filename = "Temp.pptx", nr = 1, nc = 1, irow = 1, icol = 1,onsame = FALSE,
        append = TRUE, width = 6, height = 6, devsize = FALSE, units = "in")



# Figure 3: Correlation plot ----------------------------------------------

i<-"VLDLC_S"
q<-"VLDLC_L"

gg<-ggplot(data=cgps, aes_string(x=i,y=q))+theme_nothing()+
  theme(axis.line = element_line(), axis.title.x = element_text(),axis.title.y = element_text(angle = 90),axis.ticks = element_line(),axis.text = element_text(),axis.ticks.length = unit(.2, "cm"))+
  #geom_point(aes_string(y=q),color="gray40",alpha=0.1)+
  geom_smooth(color="black",fill="black",method="lm",size=1)+
  geom_hline(yintercept=quantile(cgps[,q],0.5, na.rm=T),linetype=2,size=1,color="black")+geom_vline(xintercept=quantile(cgps[,i],0.5, na.rm=T),linetype=2,size=1,color="black")+
  xlab(i)+ylab(q)+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  coord_cartesian(ylim=c(0,2),xlim=c(0,2))+stat_summary_bin()


gg_raster<-ggplot(data=cgps, aes_string(x=i,y=q,color="VLDLC_L_S"))+theme_nothing()+
  theme(axis.line = element_line(), axis.title.x = element_text(),axis.title.y = element_text(angle = 90),axis.ticks = element_line(),axis.text = element_text(),axis.ticks.length = unit(.2, "cm"))+
  geom_point(alpha=0.1)+scale_color_manual(values = c("black","darkcyan","darkorchid4","darkred"))+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  coord_cartesian(ylim=c(0,2),xlim=c(0,2))


### Calculating mean of 5 closest values for DST 
#df$VLDLC_S_bin<-ntile(df$VLDLC_S, n=nrow(df)/5)

#for(i in 1:nrow(df)){
  #binnr<-df$VLDLC_S_bin[i]
  
  #df$VLDLC_S_bin_mean[i]<-mean(df[df$VLDLC_S_bin==binnr,"VLDLC_S"],na.rm=T)
#}



gg_dens_VLDLC_S<-ggplot(data=cgps) + 
  geom_density(aes_string(x=i),color="darkcyan", fill="darkcyan", alpha=0.4)+
  theme_nothing()+ scale_x_continuous(limits=c(0,2))

gg_dens_VLDLC_L<-ggplot(data=cgps) + 
  geom_density(aes_string(x=q),color="darkorchid4", fill="darkorchid4", alpha=0.4)+
  theme_nothing()+ coord_flip()+ scale_x_continuous(limits=c(0,2))


figure2<-plot_grid(plot_grid(gg_dens_VLDLC_S,gg,ncol=1,nrow=2,rel_heights= c(0.3,1)),gg_dens_VLDLC_L,ncol=2, rel_widths = c(1,0.1))

figure2

toffice(figure = figure2, format = "pptx",
        filename = "Temp.pptx", nr = 1, nc = 1, irow = 1, icol = 1,onsame = FALSE,
        append = TRUE, width = 6, height = 6, devsize = FALSE, units = "in")


figure2<-plot_grid(plot_grid(gg_dens_VLDLC_S,gg_raster,ncol=1,nrow=2,rel_heights= c(0.3,1)),gg_dens_VLDLC_L,ncol=2, rel_widths = c(1,0.1))

ggsave("Temp.png", figure2,device = "png",
       width=6,height = 6,units = c("in"),dpi = 600)


# Figure 4: Discordance analysis by median  ------------------------------------------------------


fig4<-data.frame()
for (i in c("VLDLC_L_S"))
  for (k in c("ASCVD","AMI","IS","PAD")){
    model<-cph(as.formula(paste0("Surv(startage,stopage",k,",",k,")~",i,"+strat(sex)+systolic+diastolic+strat(status_smoking)+rcs(cum_smoking, 3)+birthdate+rcs(bmi,4)+dmall_pre",
                                 if (grepl("VLDL",i)){"+LDLC_recov"} else {""},if (grepl("IS",k)){"+af"} else  {""})) ,data=cgps,x=TRUE, y=TRUE)  
    
    pred<-Predict(model,name=i,ref.zero = T)[,c("VLDLC_L_S","yhat","lower","upper")]
    
    
    # Regression dilutiob bias correction difficult due to categories
    #pred[c("yhat","lower","upper")]<-exp(pred[c("yhat","lower","upper")]/0.65)
    pred[c("yhat","lower","upper")]<-exp(pred[c("yhat","lower","upper")])
    
    labels<-pred[,"VLDLC_L_S"]
    
    fig4[paste(labels,k),"var"]<-paste(k,labels)
    fig4[paste(labels,k),c("yhat","lower","upper")]<-pred[c("yhat","lower","upper")] 
    
    #P-values
    pvals<-pnorm(abs(model$coef[1:3])/sqrt(model$var[1:3,1:3]),lower.tail=F)*2
    fig4[paste(labels,k),"pval"]<-c(1,pvals[1,1],pvals[2,2],pvals[3,3])
    
    
    p_years<-pyears(as.formula(paste0("Surv(startage,stopage",k,",",k,")~",i)) ,data=subset(cgps,!is.na(get(i))),scale=1)
    
    fig4[paste(labels,k),"N"]<-p_years$n
    fig4[paste(labels,k),"N_events"]<-p_years$event
    fig4[paste(labels,k),"events_pyears"]<-p_years$event/p_years$pyear*1000
    
  }




#Formatting P-values for table with blank spaces
#Table
fig4_table <- tibble(cbind(fig4$var,fig4$N,fig4$N_events,round(fig4$events_pyears,1),fig4[grepl("yhat",colnames(fig4))|grepl("pval",colnames(fig4))])) %>% flextable()
fig4_table<- compose(fig4_table,j=paste0("yhat"), value=as_paragraph(paste0(round(fig4[,paste0("yhat")], digits = 2), " (", round(fig4[,paste0("lower")], digits = 2),"-", round(fig4[,paste0("upper")], digits = 2), ")")))
fig4_table<- compose(fig4_table,j=paste0("pval"), value=as_paragraph(signif.p(fig4[,paste0("pval")]),as_sup(signif.p(fig4[,paste0("pval")],exp=T))))

fig4_table<-fontsize(fig4_table, size = 12, part = "body")
fig4_table<-align(autofit(fig4_table, part="all"),align="left", part = "all")
fig4_table


#figure

group.colors <- rep(c("black","darkcyan","darkorchid4","darkred"),4)

fig4$var<-factor(fig4$var,levels=fig4$var)


figure4<-ggplot(fig4, aes(y=var,x=yhat, xmin =lower, xmax = upper))+
  labs(x="Hazard ratio (95% CI)", y=NULL)+ 
  theme_nothing()+
  theme(axis.line.x=element_line(size=1),axis.ticks.x=element_line(size=1),axis.ticks.length.x = unit(0.2, "cm"), axis.text.x =element_text(vjust=0.4), axis.title.x = element_text() , plot.margin=margin(0,3,0,3),plot.title =element_text(face="bold",size=10,hjust = 0.5),strip.text.x = element_text(size = 10))+
  geom_vline(aes(xintercept = 1), linetype = 2, size=0.35) + geom_point(size=3,colour=group.colors)+
  geom_errorbar(width=0.5, cex=1,color=group.colors, position=position_dodge(width=0.1))+ 
  scale_colour_grey(start = 0.75, end = 0)+
  scale_y_discrete(limits=rev)  +coord_cartesian( xlim = c(0.5,2))



toffice(figure = figure4, format = "pptx", title = "Copenhagen General Population Study",
        filename = "Temp.pptx", nr = 1, nc = 1, irow = 1, icol = 1,onsame = FALSE,left = 0.15,top = 0.15,
        append = FALSE, width = 1.5, height = 5,devsize = FALSE, units = "in")


# Figure 5: Continuous discordance analysis ------------------------------


cgps$VLDLC_percentiles<-as.numeric(cut2(cgps[,"XSVLDLC_recov"],g=100))-as.numeric(cut2(cgps$VLDLC_recov,g=100))
dd<- datadist(cgps)
dd$limits$VLDLC_percentiles[1:3]<-c(1,5,11) 
options(datadist = "dd")


fig4<-data.frame()
for (k in c("ASCVD","AMI","IS","PAD")){
  dat<-cgps
  
  
  p_years<-pyears(as.formula(paste0("Surv(startage,stopage",k,",",k,")~1")) ,data=dat,subset=m5==0,scale=1)
  
  fig4["N",paste0("yhat",k)]<-p_years$n
  fig4["N_events",paste0("yhat",k)]<-p_years$event
  fig4["events_pyears",paste0("yhat",k)]<-p_years$event/p_years$pyear*1000
  
  
  for (i in paste0(c("XSVLDL","SVLDL","MVLDL","LVLDL","XLVLDL","XXLVLDL"),"C_recov")){
    
    
    dat$VLDLC_percentiles<-as.numeric(cut2(dat[,i],g=100))-as.numeric(cut2(dat$VLDLC_recov,g=100))
    
    
    model<-cph(as.formula(paste0("Surv(startage,stopage",k,",",k,")~VLDLC_percentiles+VLDLC_recov+LDLC_recov+strat(sex)+systolic+diastolic+strat(status_smoking)+rcs(cum_smoking, 3)+birthdate+rcs(bmi,4)+dmall_pre", 
                                 if (k=="IS"){"+af"} else {""})) ,data=dat,subset=m5==0,x=TRUE, y=TRUE)  
    
    pred<-as.data.frame(summary(model))
    
    fig4[paste0(i),paste0(c("yhat","lower","upper"),k)]<-pred[2,c("Effect","Lower 0.95","Upper 0.95")]
    fig4[paste0(i),paste0("pval",k)]<-pnorm(abs(model$coef[1])/sqrt(model$var[1,1]),lower.tail=F)*2
    
    
    print(paste(i,k,"done!"))
  }}



fig4$var<-factor(rownames(fig4),levels=rownames(fig4))

#Table
fig4_table <- tibble(variable=fig4$var,
                     HR1=fig4$yhatASCVD,P1=paste(signif(fig4$pvalASCVD,2)), 
                     HR2=fig4$yhatAMI, P2=paste(signif(fig4$pvalAMI,2)),
                     HR3=fig4$yhatIS, P3=paste(signif(fig4$pvalIS,2)),
                     HR4=fig4$yhatPAD,P4=paste(signif(fig4$pvalPAD,2)) ) %>% flextable()
fig4_table<- compose(fig4_table,j="HR1", value=as_paragraph(paste0(round(fig4$yhatASCVD, digits = 2), " (", round(fig4$lowerASCVD, digits = 2),"-", round(fig4$upperASCVD, digits = 2), ")")))
fig4_table<- compose(fig4_table,j="HR2", value=as_paragraph(paste0(round(fig4$yhatAMI, digits = 2), " (", round(fig4$lowerAMI, digits = 2),"-", round(fig4$upperAMI, digits = 2), ")")))
fig4_table<- compose(fig4_table,j="HR3", value=as_paragraph(paste0(round(fig4$yhatIS, digits = 2), " (", round(fig4$lowerIS, digits = 2),"-", round(fig4$upperIS, digits = 2), ")")))
fig4_table<- compose(fig4_table,j="HR4", value=as_paragraph(paste0(round(fig4$yhatPAD, digits = 2), " (", round(fig4$lowerPAD, digits = 2),"-", round(fig4$upperPAD, digits = 2), ")")))

align(autofit(fig4_table, part="all"),align="left", part = "all")


#figure

group.colors<-colorRampPalette(colors = c("darkcyan","darkgrey", "darkorchid4"))(6)

figure4a<-ggplot(fig4[!(row.names(fig4) %in% c("N","N_events","events_pyears")),], 
                 aes(y=var,x=yhatASCVD, xmin =lowerASCVD, xmax = upperASCVD))+
  labs(x="Hazard ratio (95% CI)", y=NULL)+ 
  theme_nothing()+
  theme( axis.line.x=element_line(size=1),axis.ticks.x=element_line(size=1),axis.ticks.length.x = unit(0.2, "cm"), axis.text.x =element_text(vjust=0.4), axis.title.x = element_text() , plot.margin=margin(0,3,0,3),plot.title =element_text(face="bold",size=10,hjust = 0.5),strip.text.x = element_text(size = 10))+
  geom_vline(aes(xintercept = 1), linetype = 2, size=0.35) + geom_errorbar(width=0.5, cex=0.5,color="black", position=position_dodge(width=0.1))+ 
  geom_point(size=3,colour=group.colors)+
  scale_colour_grey(start = 0.75, end = 0)+
  scale_y_discrete(limits=rev)+scale_x_continuous(lim=c(0.75,1.5),breaks = c(0.75,1,1.25,1.5))


figure4b<-ggplot(fig4[!(row.names(fig4) %in% c("N","N_events","events_pyears")),], 
                 aes(y=var,x=yhatAMI, xmin =lowerAMI, xmax = upperAMI))+
  labs(x="Hazard ratio (95% CI)", y=NULL)+ 
  theme_nothing()+
  theme( axis.line.x=element_line(size=1),axis.ticks.x=element_line(size=1),axis.ticks.length.x = unit(0.2, "cm"), axis.text.x =element_text(vjust=0.4), axis.title.x = element_text() , plot.margin=margin(0,3,0,3),plot.title =element_text(face="bold",size=10,hjust = 0.5),strip.text.x = element_text(size = 10))+
  geom_vline(aes(xintercept = 1), linetype = 2, size=0.35) + geom_errorbar(width=0.5, cex=0.5,color="black", position=position_dodge(width=0.1))+ 
  geom_point(size=3,colour=group.colors)+
  scale_colour_grey(start = 0.75, end = 0)+
  scale_y_discrete(limits=rev)+scale_x_continuous(lim=c(0.75,1.5),breaks = c(0.75,1,1.25,1.5))


figure4c<-ggplot(fig4[!(row.names(fig4) %in% c("N","N_events","events_pyears")),], 
                 aes(y=var,x=yhatIS, xmin =lowerIS, xmax = upperIS))+
  labs(x="Hazard ratio (95% CI)", y=NULL)+ 
  theme_nothing()+
  theme( axis.line.x=element_line(size=1),axis.ticks.x=element_line(size=1),axis.ticks.length.x = unit(0.2, "cm"), axis.text.x =element_text(vjust=0.4), axis.title.x = element_text() , plot.margin=margin(0,3,0,3),plot.title =element_text(face="bold",size=10,hjust = 0.5),strip.text.x = element_text(size = 10))+
  geom_vline(aes(xintercept = 1), linetype = 2, size=0.35) + geom_errorbar(width=0.5, cex=0.5,color="black", position=position_dodge(width=0.1))+ 
  geom_point(size=3,colour=group.colors)+
  scale_colour_grey(start = 0.75, end = 0)+
  scale_y_discrete(limits=rev)+scale_x_continuous(lim=c(0.75,1.5),breaks = c(0.75,1,1.25,1.5))


figure4d<-ggplot(fig4[!(row.names(fig4) %in% c("N","N_events","events_pyears")),], 
                 aes(y=var,x=yhatPAD, xmin =lowerPAD, xmax = upperPAD))+
  labs(x="Hazard ratio (95% CI)", y=NULL)+ 
  theme_nothing()+
  theme( axis.line.x=element_line(size=1),axis.ticks.x=element_line(size=1),axis.ticks.length.x = unit(0.2, "cm"), axis.text.x =element_text(vjust=0.4), axis.title.x = element_text() , plot.margin=margin(0,3,0,3),plot.title =element_text(face="bold",size=10,hjust = 0.5),strip.text.x = element_text(size = 10))+
  geom_vline(aes(xintercept = 1), linetype = 2, size=0.35) + geom_errorbar(width=0.5, cex=0.5,color="black", position=position_dodge(width=0.1))+ 
  geom_point(size=3,colour=group.colors)+
  scale_colour_grey(start = 0.75, end = 0)+
  scale_y_discrete(limits=rev)+scale_x_continuous(lim=c(0.75,1.5),breaks = c(0.75,1,1.25,1.5))



figure4<-plot_grid(figure4a,figure4b,figure4c,figure4d,ncol=4)


toffice(figure = figure4, format = "pptx",
        filename = "Temp.pptx", nr = 1, nc = 1, irow = 1, icol = 1,onsame = FALSE,
        append = TRUE, width = 5, height = 2.2, devsize = FALSE, units = "in")


# Figure 6: density ridge plots for small-large VLDLC groups ------------------------------------



library(ggridges)

group.colors <- rev(c("grey30","darkcyan","darkorchid4","darkred"))
namesfig5<-c("VLDLC_recov","LDLC_recov","trig", "remncholcalc","ldl", "apob_mgdl" )

fig<-list()
for ( i in namesfig5 ){
    
    #Data frame
    dat<-cgps
    
   dat$VLDLC_L_S <- factor(dat$VLDLC_L_S, levels=rev(levels(dat$VLDLC_L_S)))
    
    #Median (IQR) values (also in mg/dL values for cholesterol)
    qdm1<-round(quantile(subset(dat,VLDLC_L_S==levels(cgps$VLDLC_L_S)[1])[i],c(0.5,0.25,0.75),na.rm=T),1)
    qdm2<-round(quantile(subset(dat,VLDLC_L_S==levels(cgps$VLDLC_L_S)[2])[i],c(0.5,0.25,0.75),na.rm=T),1)
    qdm3<-round(quantile(subset(dat,VLDLC_L_S==levels(cgps$VLDLC_L_S)[3])[i],c(0.5,0.25,0.75),na.rm=T),1)
    qdm4<-round(quantile(subset(dat,VLDLC_L_S==levels(cgps$VLDLC_L_S)[4])[i],c(0.5,0.25,0.75),na.rm=T),1)

    #Violin plot
    fig[paste(i)]<-list(ggplot(data=dat, 
                                           aes_string(x=i, fill="VLDLC_L_S", y="VLDLC_L_S")) + 
                                      geom_density_ridges() +scale_fill_manual(values=group.colors, labels=levels(cgps$VLDLC_L_S))+
                                      labs(x=NULL, y=i,title= paste(i,"\n",paste0(qdm1[1]," (",qdm1[2],"-",qdm1[3],")\n"),paste0(qdm2[1]," (",qdm2[2],"-",qdm2[3],")\n"),
                                                                    paste0(qdm3[1]," (",qdm3[2],"-",qdm3[3],")\n"),paste0(qdm4[1]," (",qdm4[2],"-",qdm4[3],")")))+ 
                                      scale_x_continuous(expand=c(0,0))+ scale_y_discrete(expand=c(0,0)) +coord_cartesian(xlim=c(0,if
                                                                                                                                 (i=="VLDLC_recov"|i=="remncholcalc") {2} else if
                                                                                                                                 (i=="LDLC_recov"|i=="ldl"){6} else if
                                                                                                                                 (i=="trig"){5} else if 
                                                                                                                                 (i=="apob_mgdl"){250}))   + theme_nothing()+
                                      theme(axis.text.x = element_text(),plot.title = element_text(size=10),axis.line.x=element_line(size=1),axis.ticks.length = unit(.2, "cm"),axis.ticks.x=element_line(size=1)))
  }


#Changing function plot_grid
fargs  <- formals(plot_grid)
fargs$ncol           <- 6
fargs$rel_heights   <-c(1,1)
formals(plot_grid) <- fargs


figure5<-do.call(plot_grid,fig[namesfig5])
figure5

mgdl_names<-c("VLDLC_recovmgdl","LDLC_recovmgdl","trigmgdl", "remncholcalcmgdl","ldlmgdl", "apob_mgdl" )
mgdlfig<-list()

for ( i in mgdl_names ){
  
  #Data frame
  dat<-cgps
  
  dat$VLDLC_L_S <- factor(dat$VLDLC_L_S, levels=rev(levels(dat$VLDLC_L_S)))
  
  
  #Median (IQR) values (also in mg/dL values for cholesterol)
  qdm1<-round(quantile(subset(dat,VLDLC_L_S==levels(cgps$VLDLC_L_S)[1])[i],c(0.5,0.25,0.75),na.rm=T),0)
  qdm2<-round(quantile(subset(dat,VLDLC_L_S==levels(cgps$VLDLC_L_S)[2])[i],c(0.5,0.25,0.75),na.rm=T),0)
  qdm3<-round(quantile(subset(dat,VLDLC_L_S==levels(cgps$VLDLC_L_S)[3])[i],c(0.5,0.25,0.75),na.rm=T),0)
  qdm4<-round(quantile(subset(dat,VLDLC_L_S==levels(cgps$VLDLC_L_S)[4])[i],c(0.5,0.25,0.75),na.rm=T),0)
  
  #Violin plot
  mgdlfig[paste(i)]<-list(ggplot(data=dat, 
                                         aes_string(x=i, fill="VLDLC_L_S", y="VLDLC_L_S")) + 
                                    geom_density_ridges() +scale_fill_manual(values=group.colors, labels=levels(cgps$VLDLC_L_S))+
                                    labs(x=NULL, y=i,title= paste(i,"\n",paste0(qdm1[1]," (",qdm1[2],"-",qdm1[3],")\n"),paste0(qdm2[1]," (",qdm2[2],"-",qdm2[3],")\n"),
                                                                  paste0(qdm3[1]," (",qdm3[2],"-",qdm3[3],")\n"),paste0(qdm4[1]," (",qdm4[2],"-",qdm4[3],")")))+ 
                                    scale_x_continuous(expand=c(0,0))+ scale_y_discrete(expand=c(0,0))  + theme_nothing()+
                                    theme(axis.text.x = element_text(),plot.title = element_text(size=10),axis.line.x=element_line(size=1),axis.ticks.length = unit(.2, "cm"),axis.ticks.x=element_line(size=1)))
}


fargs  <- formals(plot_grid)
fargs$ncol           <- 6
fargs$rel_heights   <-c(1,1)
formals(plot_grid) <- fargs

fig5_mgdl<-do.call(plot_grid,mgdlfig[mgdl_names])


toffice(figure5,format = "pptx",
        filename = "Temp.pptx", nr = 1, nc = 1, irow = 1, icol = 1,onsame = FALSE, left = 0.15,top = 0.15,
        append = FALSE, width = 8, height = 4,devsize = FALSE, units = "in")

toffice(fig5_mgdl,format = "pptx",
        filename = "Temp.pptx", nr = 1, nc = 1, irow = 1, icol = 1,onsame = FALSE, left = 0.15,top = 0.15,
        append = FALSE, width = 8, height = 4,devsize = FALSE, units = "in")

# Figure SX: Fasting time -------------------------------------------------




cgps$time_since_meal<-factor(cgps$time_since_meal,levels=c(levels(factor(cgps$time_since_meal))[3:10],"> 8 timer", ""))
quant_025<-function(x){quantile(x,0.25)}
quant_075<-function(x){quantile(x,0.75)}


df<-with(cgps,data.frame(time_since_meal=factor(levels(time_since_meal),levels=levels(time_since_meal)),
                         VLDLC_S=tapply(VLDLC_S, time_since_meal,median),
                         VLDLC_L=tapply(VLDLC_L, time_since_meal,median),
                         VLDLC_S_lower=tapply(VLDLC_S, time_since_meal,quant_025),
                         VLDLC_S_upper=tapply(VLDLC_S, time_since_meal,quant_075),
                         VLDLC_L_lower=tapply(VLDLC_L, time_since_meal,quant_025),
                         VLDLC_L_upper=tapply(VLDLC_L, time_since_meal,quant_075)))


fig<-ggplot(data=subset(df,time_since_meal!=""),aes(x=time_since_meal))+
         geom_point(aes(y=VLDLC_S),color="darkcyan")+geom_errorbar(aes(ymin=VLDLC_S_lower,ymax=VLDLC_S_upper),color="darkcyan")+
  geom_point(aes(y=VLDLC_L),color="darkorchid4")+geom_errorbar(aes(ymin=VLDLC_L_lower,ymax=VLDLC_L_upper),color="darkorchid4")+
  coord_cartesian(ylim=c(0,1))+
  labs(x="Time since last meal, hours" , y="Cholesterol")+
  scale_x_continuous(expand = c(0,0))+scale_x_discrete(expand = c(0,0),labels=c(1:8,'>8 hours'))+
  theme(axis.line.x.bottom =element_line(),axis.line.y.left =element_line() ,axis.line.x.top=element_blank() ,plot.margin=margin(0,3,0,3),plot.subtitle = element_blank(), plot.title = element_blank(),axis.text.y = element_text(size=10), 
        panel.grid=element_blank(),strip.background = element_rect(color = "white", fill = "white"), panel.background = element_rect(fill = "white"),axis.ticks.length = unit(.2, "cm"))

toffice(figure = fig, format = "pptx",
        filename = "Temp.pptx", nr = 1, nc = 1, irow = 1, icol = 1,onsame = FALSE,
        append = TRUE, width = 6, height = 6, devsize = FALSE, units = "in")




# Figure SX: RCS plot LDL-----------------------------------------------


fig3<-list()
for (k in c("ASCVD","AMI","IS","PAD"))
  for (i in c("LDLC_recov")){
    model<-cph(as.formula(paste0("Surv(startage,stopage",k,",",k,")~rcs(",i,",3)+strat(sex)+systolic+diastolic+strat(status_smoking)+rcs(cum_smoking, 3)+birthdate+rcs(bmi,4)+dmall_pre",
                                 if (grepl("VLDL",i)){"+LDLC_recov"} else  {""},if (grepl("IS",k)){"+af"} else  {""})) ,data=cgps,x=TRUE, y=TRUE)  
    
    
    pred<-Predict(model,ref.zero = TRUE)
    pred[c("yhat","lower","upper")]<-exp(pred[c("yhat","lower","upper")]/reg_dil_ratio_all[reg_dil_ratio_all$parm==i,"rdb"])
    
    
    col<-if (i=="VLDLC_S"){"darkcyan"} else if (i=="VLDLC_L") {"darkorchid4"} else {"darkblue"}
    
    
    gg<-ggplot(pred[1:200, ] ,
               colfill = col)+
      geom_line(size=1, colour=col,linetype=1) +
      geom_hline(aes(yintercept = 1), linetype = 2)+ coord_cartesian(xlim=c(0,6),ylim= c(0.5, 10))+
      labs(x=i , y=paste("Hazard ratio for",k,"(95% CI)"))+scale_y_continuous(expand = c(0,0),trans="log",breaks = c(0.5,1,2,4,8))+scale_x_continuous(expand = c(0,0))+
      theme(axis.line.x.bottom =element_line(),axis.line.y.left =element_line() ,axis.line.x.top=element_blank() ,plot.margin=margin(0,3,0,3),plot.subtitle = element_blank(), plot.title = element_blank(),axis.text.y = element_text(size=10), 
            panel.grid=element_blank(),strip.background = element_rect(color = "white", fill = "white"), panel.background = element_rect(fill = "white"),axis.ticks.length = unit(.2, "cm"))
    
    
    
    gg_dens<-ggplot(data=cgps) + 
      geom_density(aes_string(x=i),color=col, fill=col, alpha=0.4)+
      theme_nothing()+scale_x_continuous(limits=c(0,6))
    
    
    fig3[paste0(i,k)]<- list(plot_grid(if(k=="ASCVD"){gg_dens} else {ggplot()+theme_nothing()},
                                       gg,ncol=1,nrow=2,
                                       rel_heights = if(k=="ASCVD"){c(0.3,1)}else {c(0.1,1)}))   
    
    print(c(paste(i,k),model$stats[1:2]))}


#Changing function plot_grid
fargs  <- formals(plot_grid)
fargs$ncol           <- 1
formals(plot_grid) <- fargs


figure3<-do.call(plot_grid,fig3)


toffice(figure = figure3, format = "pptx",
        filename = "Temp.pptx", nr = 1, nc = 1, irow = 1, icol = 1,onsame = FALSE,
        append = TRUE, width = 3, height = 6, devsize = FALSE, units = "in")



# Figure SX: Standard sensitivity analyses --------------------------------


library(data.table)

endpoints<-c("ASCVD","AMI","IS","PAD")
models<-c("original","bmi_ascvd_dmall","time_followup","comprisk","ldl_unadj","VLDLC_adj","remnc_adj")
groups<-c("Low in both","High in small remnants only","High in large TGRL only","High in both")

fig<-data.frame()
  for (q in models)
    for (k in endpoints){
    
    dat<-data.table(cgps[c("VLDLC_L_S","VLDLC_S","VLDLC_L","VLDLC_recov","LDLC_recov","remncholcalc","chol_ratio","sex","af","birthdate","obushnr","startage","systolic","diastolic","cum_smoking","status_smoking",endpoints,paste0("stopage",endpoints),paste0(endpoints,"_comprisk"),paste0("followup",endpoints),"bmi","ascvd_pre","dmall_pre","dm_all")])

    
    if (q=="comprisk"){
      comprisk_dat<-finegray(as.formula(paste0("Surv(startage,stopage",k,",",k,"_comprisk)~.")), etype=k, data=subset(dat,!is.na(VLDLC_L_S) & get(paste0("followup",k))>0),id=obushnr, count=1)  
     
    }else {}
    
      
      adjust<-paste0("VLDLC_L_S", if(q!="ldl_unadj") {"+LDLC_recov"} else {""},if (q=="VLDLC_adj")  {"+VLDLC_recov"} else {""},
                     if (q=="remnc_adj"){"+remncholcalc"} else {""},
                     "+strat(sex)+birthdate+systolic+diastolic+strat(status_smoking)+rcs(cum_smoking, 3)+rcs(bmi,4)+dmall_pre",
                     if (k=="IS"){"+af"} else {""},
                     if (q=="bmi_ascvd_dmall"){"+ascvd_pre"} else {""})
      
      
      if (q=="comprisk"){
        model<-cph(as.formula(paste0("Surv(fgstart,fgstop,fgstatus)~",adjust)) ,data=comprisk_dat,x=TRUE, y=TRUE,weights=fgwt)  
        
      } else if (q=="time_followup"){
        model<-cph(as.formula(paste0("Surv(followup",k,",",k,")~",adjust,"+rcs(startage,4)")) ,data=subset(dat,get(paste0("followup",k))>0),x=TRUE, y=TRUE) 
      } else {
        model<-cph(as.formula(paste0("Surv(startage,stopage",k,",",k,")~",adjust)) ,data=dat,x=TRUE, y=TRUE)  
      }
      
      
      pyrs<-pyears(as.formula(paste0("Surv(startage,stopage",k,",",k,")~VLDLC_L_S")),data=dat)
      
      pred<-Predict(model, VLDLC_L_S,ref.zero = TRUE, fun=exp)
      
      
      fig[paste(k,q,groups),"N_indiv"]<-pyrs$n
      fig[paste(k,q,groups),"N_events"]<-pyrs$event
      
      fig[paste(k,q,groups),c("yhat","lower","upper")]<-pred[,c("yhat","lower","upper")] 
      fig[paste(k,q,groups),c("P")]<-c(NA,
        pnorm(abs(model$coef[1])/sqrt(model$var[1,1]),lower.tail=F)*2,
        pnorm(abs(model$coef[2])/sqrt(model$var[2,2]),lower.tail=F)*2,
        pnorm(abs(model$coef[3])/sqrt(model$var[3,3]),lower.tail=F)*2)
      
      print(paste(k,q,"done!"))
    }

fig$var<-factor(rownames(fig),levels=rownames(fig))

#Creating file
toffice(figure = ggplot(), format = "pptx",title="Sensitivity analyses",
        filename = "Temp.pptx", nr = 1, nc = 1, irow = 1, icol = 1,onsame = FALSE,
        append = F,devsize = F, units = "in")

# Running sensitivity analyses
library(forestploter)

for (k in models){
df<-fig[grepl(k,rownames(fig)),]

df<-data.frame(Endpoint=c("Atherosclerotic","cardiovascular","disease","",
                          "","Myocardial","infarction","",
                          "","Ischemic","stroke","",
                          "","Peripheral","artery disease",""),
               Group=groups,
               N_indiv=format(df[,"N_indiv"], big.mark = ",", scientific = F,trim=T) ,
               N_events=format(df[,"N_events"], big.mark = ",", scientific = F,trim=T) ,
               df[,c("yhat","lower","upper")],
               P=sprintf(c("%3.1g", "%3.2g")[cut(df$P, c(-Inf, 0.10, Inf))], df$P))

#Formatting P values
df$P<-gsub("e-0", "e-",df$P)
df$P<-gsub("e", " x 10",df$P)
df$P<-str_pad(df$P, width=4, side="right", pad="0")
df$P<-gsub("NA00", " ",df$P)

#Creating printed text
df$HR <- ifelse(is.na(df$yhat), "",
                                       sprintf("%.2f (%.2f-%.2f)",
                                               df$yhat, df$lower, df$upper))

#Removing confidence interval for reference
df[grep("(1.00-1.00)",df$HR),"HR"]<-"1.00"


#Creating space for CI
df$` ` <- paste(rep(" ", 20), collapse = " ")

#Set column names here
Endpoint<-"  "
Group<-"\nCholesterol content"
N_individuals<-"No. of\nindividuals"
N_events<-"No. of\nevents"
HR<-if(k=="comprisk") {"Subhazard\nratio (95% CI)"} else {"Hazard\nratio (95% CI)"}
P<-"\nP"

colnames(df)[grep("Group",colnames(df))]<-Group
colnames(df)[grep("Endpoint",colnames(df))]<-Endpoint
colnames(df)[grep("N_indiv",colnames(df))]<-N_individuals
colnames(df)[grep("N_events",colnames(df))]<-N_events
colnames(df)[grep("HR",colnames(df))]<-HR
colnames(df)[grep("P",colnames(df))]<-P


#Changing theme
tm <- forest_theme(base_size = 10,
                   ci_pch = 20,ci_lwd = 2, ci_Theight = 0.4,
                   xaxis_lwd=2,
                   title_col = "gray40",title_fontface = "plain",
                   core=list(fg_params=list(hjust = 0, x = 0),
                             bg_params=list(fill = rep(c("white","gray90"),nrow(df)/2)))
                         )


#Creating plot
p<-forest(df[ ,c(Endpoint,Group,N_individuals,N_events,HR," ",P)],
          est = df$yhat,
          lower = df$lower, 
          upper = df$upper,
          sizes = 1,
          ref_line = 1,
          ci_column = 6,
          xlab = c(if(k=="comprisk") {"Subhazard ratio (95% CI)"} else {"Hazard ratio (95% CI)"}),
          xlim = c(0.5, 2),
          ticks_at = c(0.5,1,1.5,2.0),
          title=paste0(k,""),
          theme=tm)

#Changing colors
colors<-c("black","darkcyan","darkorchid4","darkred")

for (i in colors){
p <- edit_plot(p,
               row = seq(grep(i,colors),nrow(df),4),
               col = 6,
               which = "ci",
               gp = gpar(col = i))
}

#adding line under header
p <- add_border(p, part = "header", row = 1, where = "bottom",gp = gpar(lty=1,lwd=2))

#adding line to mark differebt models
p <- add_border(p, row = seq(4,nrow(df)-1,4), where = "bottom",gp = gpar(lty=1,lwd=1))

#Bold font for first column
p <- edit_plot(p, col = 1, which = "background",gp = gpar(fill = "white"))
p <- edit_plot(p, col = 1, gp = gpar(fontface = "bold"))

#Exporting plot
scale<-get_wh(plot = p, unit = "in")

toffice(figure = p, format = "pptx",
        filename = "Temp.pptx", nr = 1, nc = 1, irow = 1, icol = 1,onsame = FALSE,
        append = T,devsize = F,width=scale[1],height=scale[2], units = "in")

print(paste(k,"done!"))
}



# Figure SX: stratified sensitivity analysis--------------------------------


cgps_int<-cgps
cgps_int$`Sex`<-factor(cgps_int$sex ,labels =c("Women","Men") )
cgps_int$`Age`<-factor(cut2(cgps_int$startage,60) ,labels =c("<60 years",">60 years") )
cgps_int$`LDL cholesterol`<-factor(cut2(cgps_int$LDLC_recov,3.3) ,labels =c("<3.3 mmol/L(128 mg/dL)",">3.3 mmol/L(128 mg/dL)") )
cgps_int$`Current smoking`<-factor(cgps_int$current_smoking ,labels =c("No","Yes") )
cgps_int$`Systolic blood pressure`<-factor(cut2(cgps_int$systolic,140) ,labels =c("<140 mmHg",">140 mmHg") )
cgps_int$`Atrial fibrillation`<-factor(cgps_int$af ,labels =c("No","Yes") )

cgps_int$`Body mass index`<-factor(cut2(cgps_int$bmi,25) ,labels =c("<25 kg/m2",">25 kg/m2") )
cgps_int$`Diabetes`<-factor(cgps_int$dmall_pre ,labels =c("No","Yes") )
cgps_int$`ASCVD `<-factor(cgps_int$ascvd_pre ,labels =c("No","Yes") )

cgps_int$`Triglycerides`<-factor(cut2(cgps_int$trig,1.5) ,labels =c("<1.5 mmol/L(133 mg/dL)",">1.5 mmol/L(133 mg/dL)") )
cgps_int$`Clinical remnant cholesterol`<-factor(cut2(cgps_int$remncholcalc,0.7) ,labels =c("<0.7 mmol/L(27 mg/dL)",">0.7 mmol/L(27 mg/dL)") )
cgps_int$`Clinical LDL cholesterol`<-factor(cut2(cgps_int$ldl,3.3) ,labels =c("<3.3 mmol/L(128 mg/dL)",">3.3 mmol/L(128 mg/dL)"))
cgps_int$`Apolipoprotein B`<-factor(cut2(cgps_int$apob_mgdl,130) ,labels =c("<130 mg/dL",">130 mg/dL") )


endpoints<-c("ASCVD","AMI","IS","PAD")
models<-c("All","Sex","Age","LDL cholesterol","Current smoking","Systolic blood pressure","Atrial fibrillation",
          "Body mass index","Diabetes","ASCVD ",
          "Triglycerides","Clinical remnant cholesterol","Clinical LDL cholesterol","Apolipoprotein B")
          
groups<-c("Only high small remnant-L cholesterol","Only high large remnant-L cholesterol")


fig<-data.frame()
for (k in endpoints)
  for (q in models){
    
    adjust<-paste0(if (grepl("Sex",q)==F) {"+strat(sex)"} else {""},
                    if (grepl("Diabetes",q)==F) {"+dmall_pre"} else {""},
                   if (grepl("smoking",q)==F) {"+strat(status_smoking)+rcs(cum_smoking, 3)"} else {""},
                   if (k=="IS"& q!="Atrial fibrillation"){"+af"} else {""},"+LDLC_recov+systolic+diastolic+birthdate+rcs(bmi,4)")
    
    if (q=="All"){ 
     
      model<-cph(as.formula(paste0("Surv(startage,stopage",k,",",k,")~VLDLC_L_S",adjust)) ,data=cgps,x=TRUE, y=TRUE)  
      
      pyrs<-pyears(as.formula(paste0("Surv(startage,stopage",k,",",k,")~VLDLC_L_S")),data=cgps)
      
      pred<-Predict(model, VLDLC_L_S,ref.zero = TRUE, fun=exp)
      
      fig[paste(q,groups),"stratum"]<-q
      fig[paste(q,groups),"variable"]<-q
      
      fig[paste(q,groups),paste0("N_indiv_",k)]<-pyrs$n[2:3]
      fig[paste(q,groups),paste0("N_events_",k)]<-pyrs$event[2:3]
      
      fig[paste(q,groups),paste(c("yhat","lower","upper"),k,sep = "_")]<-pred[2:3,c("yhat","lower","upper")]
      
      
      model_int<-cph(as.formula(paste0("Surv(startage,stopage",k,",",k,")~VLDLC_L_S",adjust)) ,data=cgps,x=TRUE, y=TRUE) 
      
      fig[paste(q,groups),paste0("P_",k)]<-c(NA,NA)
        
      print(paste(k,q,"done!"))
      
      } else {
        for (i in (2:3)){
    
    level<-as.character(unique(cgps_int[,q])[i])
    
    dat<-cgps_int[cgps_int[,q]==level,]
    
    
    model<-cph(as.formula(paste0("Surv(startage,stopage",k,",",k,")~VLDLC_L_S",adjust)) ,data=dat,x=TRUE, y=TRUE)  
    
    pyrs<-pyears(as.formula(paste0("Surv(startage,stopage",k,",",k,")~VLDLC_L_S")),data=dat)
    
    pred<-Predict(model, VLDLC_L_S,ref.zero = TRUE, fun=exp)
    
    fig[paste(q,groups),"stratum"]<-q
    
    fig[paste(q,groups,level),"variable"]<-q
    fig[paste(q,groups,level),"stratum"]<-level
    
    fig[paste(q,groups,level),paste0("N_indiv_",k)]<-pyrs$n[2:3]
    fig[paste(q,groups,level),paste0("N_events_",k)]<-pyrs$event[2:3]
    
    fig[paste(q,groups,level),paste(c("yhat","lower","upper"),k,sep = "_")]<-pred[2:3,c("yhat","lower","upper")]
    
    cgps_int$newname<-cgps_int[,q]
    
    model_int<-cph(as.formula(paste0("Surv(startage,stopage",k,",",k,")~VLDLC_L_S*newname",adjust)) ,data=cgps_int,x=TRUE, y=TRUE) 

    index1<-grep("high VLDLC_S.low VLDLC_L * ",names(model_int$coef))
    index2<-grep("low VLDLC_S.high VLDLC_L * ",names(model_int$coef))
    
    fig[paste(q,groups,level),paste0("P_",k)]<-if (i==2) {c(pnorm(abs(model_int$coef[index1])/sqrt(model_int$var[index1,index1]),lower.tail=F)*2,
                                     pnorm(abs(model_int$coef[index2])/sqrt(model_int$var[index2,index2]),lower.tail=F)*2)} else {c(NA,NA)}
    
    print(paste(k,q,"done!"))
  }}}


fig$var<-factor(rownames(fig),levels=rownames(fig))


#Creating file
toffice(figure = ggplot(), format = "pptx",title="Sensitivity analyses",
        filename = "Temp.pptx", nr = 1, nc = 1, irow = 1, icol = 1,onsame = FALSE,
        append = F,devsize = F, units = "in")


# Running sensitivity analyses
library(forestploter)

for (k in groups){
  
  dat<-fig[grepl(k,rownames(fig)),]

  
  df<-data.frame(Stratum=dat[,"stratum"])
  
  #Adding hazard ratios
  df<-cbind(df, format(dat[,paste0("N_indiv_",endpoints)], big.mark = ",", scientific = F,trim=T) ,
        format(dat[,paste0("N_events_",endpoints)], big.mark = ",", scientific = F,trim=T),
        dat[,paste0("yhat_",endpoints)],dat[,paste0("lower_",endpoints)],dat[,paste0("upper_",endpoints)])
  
  
  # indent the subgroup if there is a number in the placebo column
  df$Stratum <- ifelse(is.na(df$yhat_ASCVD), 
                       df$Stratum,
                       paste0("   ", df$Stratum))
  
  # NA to blank or NA will be transformed to carachter.
  df[is.na(df)] <- ""
  df[df=="NA"] <- ""
  
  #Converting back variables to numeric
  df[,paste0("yhat_",endpoints)]<-lapply(df[,paste0("yhat_",endpoints)],as.numeric)
  df[,paste0("lower_",endpoints)]<-lapply(df[,paste0("lower_",endpoints)],as.numeric)
  df[,paste0("upper_",endpoints)]<-lapply(df[,paste0("upper_",endpoints)],as.numeric)
  
  #Adding and formatting P values and hazard ratios
  for ( i in endpoints){
    #Adding P-values
    df[,paste0("P_",i)]<-sprintf(c("%3.1g", "%3.2g")[cut(dat[,paste0("P_",i)], c(-Inf, 0.10, Inf))], dat[,paste0("P_",i)])
    
    #Formatting P values
    df[,paste0("P_",i)]<-gsub("e-0", "e-",df[,paste0("P_",i)])
    df[,paste0("P_",i)]<-gsub("e", " x 10",df[,paste0("P_",i)])
    df[,paste0("P_",i)]<-str_pad(df[,paste0("P_",i)], width=4, side="right", pad="0")
    df[,paste0("P_",i)]<-gsub("NA00", " ",df[,paste0("P_",i)])
    
    #Formatting hazard ratios
    df[,paste0("HR_",i)] <- ifelse(is.na(df[,paste0("yhat_",i)]), "",
                    sprintf("%.2f (%.2f-%.2f)",
                            df[,paste0("yhat_",i)], df[,paste0("lower_",i)], df[,paste0("upper_",i)]))
    }

  
  #Creating space for CI
  df$` ` <- paste(rep(" ", 10), collapse = " ")
  df$`  ` <- paste(rep(" ", 10), collapse = " ")
  df$`   ` <- paste(rep(" ", 10), collapse = " ")
  df$`    ` <- paste(rep(" ", 10), collapse = " ")
  

  #Set column names here
  Stratum<-"\nSubgroup"
  
  N_individuals_ASCVD<-"No.\nindiv."
  N_individuals_AMI<-"No. \nindiv."
  N_individuals_IS<-"No.  \nindiv."
  N_individuals_PAD<-"No.   \nindiv."
  
  N_events_ASCVD<-"No.\nevt."
  N_events_AMI<-"No. \nevt."
  N_events_IS<-"No.  \nevt."
  N_events_PAD<-"No.  \nevt. "
  
  HR_ASCVD<- "Hazard\nratio (95% CI)"
  HR_AMI<- "Hazard\nratio (95% CI) "
  HR_IS<- "Hazard\nratio (95% CI)  "
  HR_PAD<- "Hazard\nratio (95% CI)   "
  
  P_ASCVD<-"P\nint."
  P_AMI<-"P \nint."
  P_IS<-"P  \nint."
  P_PAD<-"P   \nint."
  
  colnames(df)[grep("Stratum",colnames(df))]<-Stratum
  
  for ( i in endpoints){
  colnames(df)[grep(paste0("N_indiv_",i),colnames(df))]<-get(paste0("N_individuals_",i))
  colnames(df)[grep(paste0("N_events_",i),colnames(df))]<-get(paste0("N_events_",i))
  colnames(df)[grep(paste0("HR_",i),colnames(df))]<-get(paste0("HR_",i))
  colnames(df)[grep(paste0("P_",i),colnames(df))]<-get(paste0("P_",i))
  }
  
  
  #Changing colors
  color<-if (k=="Only high small remnant-L cholesterol"){"darkcyan"} else {"darkorchid4"}
  
  
  #Changing theme
  tm <- forest_theme(base_size = 9,
                     ci_pch = 20,ci_lwd = 1, ci_Theight = 0.4,
                     xaxis_lwd=2,ci_col = color,
                     title_col = "gray40",title_fontface = "plain",
                     core=list(fg_params=list(hjust = 0, x = 0),
                               bg_params=list(fill = c("white",rep(c(rep("gray90",3),rep("white",3)),nrow(df)-1/3)))))
  
  
  #Creating plot
  p<-forest(df[ ,c(Stratum,N_individuals_ASCVD,N_events_ASCVD,HR_ASCVD," ",P_ASCVD,
                   N_individuals_AMI,N_events_AMI,HR_AMI,"  ",P_AMI,
                   N_individuals_IS,N_events_IS,HR_IS,"   ",P_IS,
                   N_individuals_PAD,N_events_PAD,HR_PAD,"    ",P_PAD)],
            est = with(df,list(yhat_ASCVD,yhat_AMI,yhat_IS,yhat_PAD)),
            lower = with(df,list(lower_ASCVD,lower_AMI,lower_IS,lower_PAD)), 
            upper = with(df,list(upper_ASCVD,upper_AMI,upper_IS,upper_PAD)),
            sizes = 1,
            ref_line = 1,
            ci_column = c(5,10,15,20),
            xlab =rep("Hazard ratio (95% CI)",4),
            xlim = c(0.5, 2),
            ticks_at = c(0.5,1,1.5,2.0),
            title=paste0(k,""),
            theme=tm)
  

  #adding line under header
  p <- add_border(p, part = "header", row = 1, where = "bottom",gp = gpar(lty=1,lwd=2))
  
  #Bold grouping text
  p <- edit_plot(p,
                 row = c(1,seq(2,nrow(df)-1,3)),col=1,
                 gp = gpar(fontface = "bold"))
  
  
  #Exporting plot
  scale<-get_wh(plot = p, unit = "in")
  
  toffice(figure = p, format = "pptx",
          filename = "Temp.pptx", nr = 1, nc = 1, irow = 1, icol = 1,onsame = FALSE,
          append = T,devsize = F,width=scale[1],height=scale[2], units = "in")
  
  print(paste(k,"done!"))
}



