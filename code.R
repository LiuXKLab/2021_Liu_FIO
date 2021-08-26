##################### step01.install packages ######################
rm(list = ls())
options(stringsAsFactors = F)
if(!require("dplyr")) install.packages("dplyr",update = F,ask = F)
if(!require("tidyr")) install.packages("tidyr",update = F,ask = F)
if(!require("plyr")) install.packages("plyr",update = F,ask = F)
if(!require("tidyverse")) install.packages("tidyverse",update = F,ask = F)
if(!require("reshape2")) install.packages("reshape2",update = F,ask = F)
if(!require("scales")) install.packages("scales",update = F,ask = F)
if(!require("stringr")) install.packages("stringr",update = F,ask = F)
if(!require("magrittr")) install.packages("magrittr",update = F,ask = F)
if(!require("ggplot2")) install.packages("ggplot2",update = F,ask = F)
if(!require("doParallel")) install.packages("doParallel",update = F,ask = F)
if(!require("parallel")) install.packages("parallel",update = F,ask = F)
if(!require("foreach")) install.packages("foreach",update = F,ask = F)
if(!require("caret")) install.packages("caret",update = F,ask = F)
if(!require("glmnet")) install.packages("glmnet",update = F,ask = F)
if(!require("leaps")) install.packages("leaps",update = F,ask = F)
if(!require("MASS")) install.packages("MASS",update = F,ask = F)
if(!require("survival")) install.packages("survival",update = F,ask = F)
if(!require("survminer")) install.packages("survminer",update = F,ask = F)
if(!require("pheatmap")) install.packages("pheatmap",update = F,ask = F)
if(!require("ggfortify")) install.packages("ggfortify",update = F,ask = F)
if(!require("glmnet")) install.packages("glmnet",update = F,ask = F)
if(!require("ggpubr")) install.packages("ggpubr",update = F,ask = F)
if(!require("tibble")) install.packages("tibble",update = F,ask = F)
if(!require("cowplot")) install.packages("cowplot",update = F,ask = F)
if(!require("timeROC")) install.packages("timeROC",update = F,ask = F)
if(!require("survivalROC")) install.packages("survivalROC",update = F,ask = F)
if(!require("Hmisc")) install.packages("Hmisc",update = F,ask = F)
if(!require("ggsci")) install.packages("ggsci",update = F,ask = F)
if(!require("pec")) install.packages("pec",update = F,ask = F)
if(!require("riskRegression")) install.packages("riskRegression",update = F,ask = F)
if(!require("sva")) BiocManager::install("sva",update = F,ask = F)
if(!require("limma")) BiocManager::install("limma",update = F,ask = F)
if(!require("DESeq2")) BiocManager::install("DESeq2",update = F,ask = F)
if(!require("edgeR")) BiocManager::install("edgeR",update = F,ask = F)


##################### step02.risk plot ######################
load("train.rda")
load("test.rda")
load("model.rda")
features <- rownames(summary(model)$coefficients)

## training set ##
Data <- train
rs <- Data$riskscore
names(rs) <- rownames(Data)
rs <- rs[rs<=50]
rs_data <- data.frame(x=1:length(rs), rs=as.numeric(sort(rs)))
rs_data$Risk <- ifelse(rs_data$rs >= res.cut, "High", "Low")
rs_data$Risk <- factor(rs_data$Risk, levels = c("Low", "High"), ordered = F)
head(rs_data)

surv_data <- data.frame(x=1:length(rs),
                        t=Data[names(sort(rs)),'times']/365,
                        s=Data[names(sort(rs)),'status']) 
surv_data$Event <- as.factor(ifelse(surv_data$s==0,"No","Yes"))
head(surv_data)

exp_data <- Data[names(sort(rs)),which(colnames(Data) %in% features)]

# A-risk score distribution
plot.A <- ggplot(rs_data, aes(x=x,y=rs)) + 
  geom_point(aes(col=Risk),size=0.5) + 
  scale_color_manual(labels=levels(rs_data$Risk), 
                     name="Risk", 
                     values = c("#00A087FF","#DC0000FF")) + 
  geom_segment(aes(x = sum(rs_data$Risk=="Low"),
                   y = 0, 
                   xend = sum(rs_data$Risk=="Low"), 
                   yend = max(rs_data$rs)), linetype="dashed", size = 0.6) + 
  geom_text(aes(x=sum(rs_data$Risk=="Low")/2,
                y=res.cut+2,
                label=paste0("Cutoff: ",round(res.cut,3))),
            col ="black",size = 4,alpha=0.8) +
  theme(axis.title.x=element_blank()) +
  scale_x_continuous(limits = c(0,NA),expand = c(0,0)) + 
  labs(y="Risk score",x="",fill="Risk") +
  theme_classic() +
  theme(axis.ticks.x=element_blank(),
        axis.line = element_blank(),
        axis.text.x=element_blank())
plot.A

# B-event distribution
plot.B <- ggplot(surv_data,aes(x=x,y=t))+
  geom_point(aes(col=Event), size=0.5)+
  geom_vline(aes(xintercept=sum(rs_data$Risk=="Low")),size=0.6,linetype="dashed")+
  scale_x_continuous(limits = c(0,NA),expand = c(0,0))+
  scale_color_manual(labels = c("No","Yes"),
                     values = c("#00A087FF","#DC0000FF"))+
  labs(y="Time (year)",x="")+
  theme_classic()+
  theme(axis.ticks.x=element_blank(),
        axis.line = element_blank(),
        axis.text.x=element_blank())
plot.B

# C-signature
tmp <- t(exp_data)
tmp[tmp > 1] = 1
tmp[tmp < -1] = -1
reorder_cormat <- function(cormat){
  dd <- dist(cormat)
  hc <- hclust(dd,method = "average")
  cormat <-cormat[hc$order,]
}
tmp1 <- reorder_cormat(tmp)
tmp1 <- rbind(tmp1,ifelse(rs_data$Risk=="Low",-1,1))
tmp.m <- melt(tmp1)

p2 <- ggplot(tmp.m, aes(Var2, Var1),size=0.5) + 
  geom_tile(aes(fill = value))

plot.C <- p2 + scale_fill_gradient2(name="", low="#00A087FF", high="#DC0000FF", mid="white") +
  labs(x = "", y = "")+
  theme_classic()+
  theme(legend.title = element_text(size = 12), legend.position = "right",
        axis.line = element_blank(),
        axis.ticks=element_blank(),
        axis.text.x=element_blank())
plot.C

risk_train <- plot_grid(plot.A, plot.B, plot.C,rel_heights = c(1,1,1),
                        align = 'v',ncol = 1, axis="lr", scale = c(1,1,1), greedy = F)
risk_train


## test set ##
Data <- test
rs <- Data$riskscore
names(rs) <- rownames(Data)
rs <- rs[rs<=50]
rs_data <- data.frame(x=1:length(rs), rs=as.numeric(sort(rs)))
rs_data$Risk <- ifelse(rs_data$rs >= res.cut, "High", "Low")
rs_data$Risk <- factor(rs_data$Risk, levels = c("Low", "High"), ordered = F)
head(rs_data)

surv_data <- data.frame(x=1:length(rs),
                        t=Data[names(sort(rs)),'times']/365,
                        s=Data[names(sort(rs)),'status']) 
surv_data$Event <- as.factor(ifelse(surv_data$s==0,"No","Yes"))
head(surv_data)

exp_data <- Data[names(sort(rs)),which(colnames(Data) %in% features)]

# A-risk score distribution
plot.A <- ggplot(rs_data, aes(x=x,y=rs)) + 
  geom_point(aes(col=Risk),size=0.5) + 
  scale_color_manual(labels=levels(rs_data$Risk), 
                     name="Risk", 
                     values = c("#00A087FF","#DC0000FF")) + 
  geom_segment(aes(x = sum(rs_data$Risk=="Low"),
                   y = 0, 
                   xend = sum(rs_data$Risk=="Low"), 
                   yend = max(rs_data$rs)), linetype="dashed", size = 0.6) + 
  geom_text(aes(x=sum(rs_data$Risk=="Low")/2,
                y=res.cut+2,
                label=paste0("Cutoff: ",round(res.cut,3))),
            col ="black",size = 4,alpha=0.8) +
  theme(axis.title.x=element_blank()) +
  scale_x_continuous(limits = c(0,NA),expand = c(0,0)) + 
  labs(y="Risk score",x="",fill="Risk") +
  theme_classic() +
  theme(axis.ticks.x=element_blank(),
        axis.line = element_blank(),
        axis.text.x=element_blank())
plot.A

# B-event distribution
plot.B <- ggplot(surv_data,aes(x=x,y=t))+
  geom_point(aes(col=Event), size=0.5)+
  geom_vline(aes(xintercept=sum(rs_data$Risk=="Low")),size=0.6,linetype="dashed")+
  scale_x_continuous(limits = c(0,NA),expand = c(0,0))+
  scale_color_manual(labels = c("No","Yes"),
                     values = c("#00A087FF","#DC0000FF"))+
  labs(y="Time (year)",x="")+
  theme_classic()+
  theme(axis.ticks.x=element_blank(),
        axis.line = element_blank(),
        axis.text.x=element_blank())
plot.B

# C-signature
tmp <- t(exp_data)
tmp[tmp > 1] = 1
tmp[tmp < -1] = -1
reorder_cormat <- function(cormat){
  dd <- dist(cormat)
  hc <- hclust(dd,method = "average")
  cormat <-cormat[hc$order,]
}
tmp1 <- reorder_cormat(tmp)
tmp1 <- rbind(tmp1,ifelse(rs_data$Risk=="Low",-1,1))
tmp.m <- melt(tmp1)

p2 <- ggplot(tmp.m, aes(Var2, Var1),size=0.5) + 
  geom_tile(aes(fill = value))

plot.C <- p2 + scale_fill_gradient2(name="", low="#00A087FF", high="#DC0000FF", mid="white") +
  labs(x = "", y = "")+
  theme_classic()+
  theme(legend.title = element_text(size = 12), legend.position = "right",
        axis.line = element_blank(),
        axis.ticks=element_blank(),
        axis.text.x=element_blank())
plot.C

risk_test <- plot_grid(plot.A, plot.B, plot.C,rel_heights = c(1,1,1),
                        align = 'v',ncol = 1, axis="lr", scale = c(1,1,1), greedy = F)
risk_test


##################### step03.survival analysis ######################
res.cut <- surv_cutpoint(train, 
                         time = "times",
                         event = "status", 
                         variables = "riskscore",
                         minprop = 0.3)
res.cut <- res.cut$cutpoint[[1]]

## training set ##
Data <- train
Data$Risk <- as.vector(ifelse(Data$riskscore >= res.cut, "High", "Low"))
Data$Risk <- factor(Data$Risk,levels=c("Low", "High"),ordered = F)

# Calculate HR and 95%CI
Sur <- Surv(Data$times/365, Data$status)
fitd <- survdiff(Sur ~ Risk,
                 data      = Data,
                 na.action = na.exclude)
p.val = 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
HR = (fitd$obs[2]/fitd$exp[2])/(fitd$obs[1]/fitd$exp[1])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/fitd$exp[2]+1/fitd$exp[1]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/fitd$exp[2]+1/fitd$exp[1]))

HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")

fit <- survfit(Sur ~ Risk,
               data      = Data,
               type      = "kaplan-meier",
               error     = "greenwood",
               conf.type = "plain",
               na.action = na.exclude)

names(fit$strata) <- gsub("Risk=", "", names(fit$strata))

survplot.train <- ggsurvplot(fit               = fit,
                      conf.int          = FALSE,
                      risk.table        = TRUE,
                      ncensor.plot      = TRUE,
                      risk.table.col    = "strata",
                      palette           = c("#00A087FF","#DC0000FF"),
                      data              = Data,
                      size              = 1,
                      xlim              = c(0,15),
                      break.time.by     = 3,
                      legend.title      = "",
                      legend.labs=c(paste0(names(fit$strata)[1]," (",fit$n[1],")"),
                                    paste0(names(fit$strata)[2]," (",fit$n[2],")")),
                      pval              = TRUE,
                      surv.median.line  = "hv",
                      xlab              = "Time (year)",
                      ylab              = "Survival probability",
                      risk.table.y.text = FALSE)
survplot.train


## test set ##
Data <- test
Data$Risk <- as.vector(ifelse(Data$riskscore >= res.cut, "High", "Low"))
Data$Risk <- factor(Data$Risk,levels=c("Low", "High"),ordered = F)

# Calculate HR and 95%CI
Sur <- Surv(Data$times/365, Data$status)
fitd <- survdiff(Sur ~ Risk,
                 data      = Data,
                 na.action = na.exclude)
p.val = 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
HR = (fitd$obs[2]/fitd$exp[2])/(fitd$obs[1]/fitd$exp[1])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/fitd$exp[2]+1/fitd$exp[1]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/fitd$exp[2]+1/fitd$exp[1]))

HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")

fit <- survfit(Sur ~ Risk,
               data      = Data,
               type      = "kaplan-meier",
               error     = "greenwood",
               conf.type = "plain",
               na.action = na.exclude)

names(fit$strata) <- gsub("Risk=", "", names(fit$strata))

survplot.test <- ggsurvplot(fit               = fit,
                             conf.int          = FALSE,
                             risk.table        = TRUE,
                             ncensor.plot      = TRUE,
                             risk.table.col    = "strata",
                             palette           = c("#00A087FF","#DC0000FF"),
                             data              = Data,
                             size              = 1,
                             xlim              = c(0,15),
                             break.time.by     = 3,
                             legend.title      = "",
                             legend.labs=c(paste0(names(fit$strata)[1]," (",fit$n[1],")"),
                                           paste0(names(fit$strata)[2]," (",fit$n[2],")")),
                             pval              = TRUE,
                             surv.median.line  = "hv",
                             xlab              = "Time (year)",
                             ylab              = "Survival probability",
                             risk.table.y.text = FALSE)
survplot.test


##################### step04.ROC ######################
## training set ##
data <- train
survivalROC_helper <- function(t) {
  survivalROC(Stime=data$times/365, status=data$status, marker = data$riskscore, 
              predict.time =t, method="KM")
}

library(tidyverse)
library(survivalROC)
survivalROC_dat <- data_frame(t = c(1,2,3,5,7,10)) %>%
  mutate(survivalROC = map(t, survivalROC_helper),
         auc = map_dbl(survivalROC, magrittr::extract2, "AUC"),
         df_survivalROC = map(survivalROC, function(obj) {
           as_data_frame(obj[c("cut.values","TP","FP")])
         })) %>%
  dplyr::select(-survivalROC) %>%
  unnest() %>%
  arrange(t, FP, TP) %>% 
  mutate(auc = sprintf("%.2f",auc)) %>% 
  unite(year, t, auc, sep = " year AUC: ")

AUC <- factor(survivalROC_dat$year,levels = c(unique(survivalROC_dat$year)[1],unique(survivalROC_dat$year)[2],unique(survivalROC_dat$year)[3],unique(survivalROC_dat$year)[4],unique(survivalROC_dat$year)[5],unique(survivalROC_dat$year)[6]))
roc_train <- ggplot(survivalROC_dat, mapping = aes(x = FP, y = TP)) +
  geom_path(aes(color= AUC))+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
  theme_bw() +
  xlab('1−Specificity (False positive rate)') + 
  ylab("Sensitivity (True positive rate)")  + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(plot.title = element_text(colour="black", hjust = 0.5)) +
  theme(axis.text.x = element_text(colour="black", size = 12),
        axis.text.y = element_text(colour="black", size = 12),
        axis.ticks = element_line(size=0.4, color="black"),
        axis.ticks.length = unit(0.2, "cm"), 
        axis.title.x = element_text(colour="black", size = 15), 
        axis.title.y = element_text(colour="black", size = 15), 
        legend.title = element_blank(), 
        legend.text = element_text(colour="black", size = 13),
        legend.position = c(0.75,0.27))
roc_train


## test set ##
data <- test
survivalROC_data <- data_frame(t = c(2,3,5,7,10)) %>%
  mutate(survivalROC = map(t, survivalROC_helper),
         auc = map_dbl(survivalROC, magrittr::extract2, "AUC"),
         df_survivalROC = map(survivalROC, function(obj) {
           as_data_frame(obj[c("cut.values","TP","FP")])
         })) %>%
  dplyr::select(-survivalROC) %>%
  unnest() %>%
  arrange(t, FP, TP) %>% 
  mutate(auc = sprintf("%.2f",auc))%>% 
  unite(year, t, auc, sep = " year AUC: ")

AUC2 <- factor(survivalROC_data$year,levels = c(unique(survivalROC_data$year)[1],unique(survivalROC_data$year)[2],unique(survivalROC_data$year)[3],unique(survivalROC_data$year)[4],unique(survivalROC_data$year)[5]))
roc_test <- ggplot(survivalROC_data, mapping = aes(x = FP, y = TP)) +
  geom_path(aes(color= AUC2))+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
  theme_bw() +
  xlab('1−Specificity (False positive rate)') + 
  ylab("Sensitivity (True positive rate)")  + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(plot.title = element_text(colour="black", hjust = 0.5)) +
  theme(axis.text.x = element_text(colour="black", size = 12),
        axis.text.y = element_text(colour="black", size = 12),
        axis.ticks = element_line(size=0.4, color="black"),
        axis.ticks.length = unit(0.2, "cm"), 
        axis.title.x = element_text(colour="black", size = 15), 
        axis.title.y = element_text(colour="black", size = 15), 
        legend.title = element_blank(), 
        legend.text = element_text(colour="black", size = 13),
        legend.position = c(0.75,0.27))
roc_test



##################### step05.time-dependent C-index ######################
clin <- clin_info %>%
  dplyr::select(PATIENT_ID,
                AGE_AT_DIAGNOSIS,
                TUMOR_SIZE,
                LYMPH_NODES_EXAMINED_POSITIVE,
                CLAUDIN_SUBTYPE,
                DSS.status,
                DSS.time) %>%
  na.omit()
clin$CLAUDIN_SUBTYPE <- ifelse(clin$CLAUDIN_SUBTYPE=="Basal","Basal","Nonbasal")
clin$CLAUDIN_SUBTYPE <- factor(clin$CLAUDIN_SUBTYPE,levels = c("Basal","Nonbasal"),ordered = F)

clin$riskscore <- train[clin$PATIENT_ID,]$riskscore

clin$futime <- clin$DSS.time
clin$fustat <- clin$DSS.status
clin <- clin %>% dplyr::select(-DSS.time,-DSS.status)

data <- clin
data$futime <- data$futime/365

mod1 <- coxph(Surv(futime,fustat==1) ~ TUMOR_SIZE+LYMPH_NODES_EXAMINED_POSITIVE,
              data=data,
              x=T)
mod1
mod2 <- coxph(Surv(futime,fustat==1) ~ riskscore+TUMOR_SIZE+LYMPH_NODES_EXAMINED_POSITIVE,
              data=data,
              x=T)
mod2

pk2 <- pec::cindex(list("Tumor size + positive lymph nodes" = mod1, 
                        "Tumor size + positive lymph nodes + risk score" = mod2),
              formula = Surv(futime,fustat==1)~1,
              data = data,
              eval.times = seq(0,10,1/12),
              splitMethod = "Bootcv",
              B = 1000)

cindex <- data.frame(model = c(rep("Tumor size + positive lymph nodes",length(pk2$AppCindex$`Tumor size + positive lymph nodes`)),rep("Tumor size + positive lymph nodes + risk score",length(pk2$AppCindex$`Tumor size + positive lymph nodes + risk score`))),
                     times = rep(pk2$pred.time,2),
                     Cindex = c(pk2$AppCindex$`Tumor size + positive lymph nodes`,pk2$AppCindex$`Tumor size + positive lymph nodes + risk score`))

t_cindex <- ggplot(data = cindex,
                   aes(x = times,
                       y = Cindex, 
                       color = model)) + 
  geom_line(size=1.2) + 
  theme_classic() + 
  ylim(0.4,1) + 
  scale_colour_manual(values = c("#00A087FF","#DC0000FF")) + 
  scale_x_continuous(
    breaks = c(0,1,2,3,4,5,6,7,8,9,10),
    labels = c(0,1,2,3,4,5,6,7,8,9,10)) + 
  labs(title = "", x = "Time (year)", y = "Concordance index (C-index)") + 
  theme(plot.title = element_text(colour="black", hjust = 0.5, size = 15)) + 
  theme(axis.text.x = element_text(colour="black", size = 12),
        axis.text.y = element_text(colour="black", size = 12),
        axis.ticks = element_line(size=0.4, color="black"),
        axis.ticks.length = unit(0.2, "cm"), 
        axis.title.x = element_text(colour="black", size = 15), 
        axis.title.y = element_text(colour="black", size = 15), 
        legend.title = element_blank(), 
        legend.justification=c(1,0.1),
        legend.position=c(1,0.1))
t_cindex



##################### step06.time-dependent ROC ######################
pk3 <- riskRegression::Score(list("Tumor size + positive lymph nodes" = mod1, 
                                  "Tumor size + positive lymph nodes + risk score" = mod2),
             formula = Surv(futime,fustat==1)~1,
             data = data,
             metrics = "auc", 
             null.model = F, 
             times = seq(0,10,1/12))
sumAUC <- pk3$AUC$contrasts
auc <- plotAUC(pk3)


t_auc <- ggplot(data = auc,
                aes(x = times,
                    y = AUC, 
                    color = model)) +
  geom_line(size=1.2) +
  theme_classic() +
  ylim(0.4,1) +
  scale_colour_manual(values = c("#00A087FF","#DC0000FF")) + 
  scale_x_continuous(
    breaks = c(0,1,2,3,4,5,6,7,8,9,10), #刻度线的位置
    labels = c(0,1,2,3,4,5,6,7,8,9,10)) +
  labs(title = "", x = "Time (year)", y = "Area under curve (AUC)") +
  theme(plot.title = element_text(colour="black", hjust = 0.5, size = 15)) + 
  theme(axis.text.x = element_text(colour="black", size = 12), 
        axis.text.y = element_text(colour="black", size = 12), 
        axis.ticks = element_line(size=0.4, color="black"), 
        axis.ticks.length = unit(0.2, "cm"), 
        axis.title.x = element_text(colour="black", size = 15), 
        axis.title.y = element_text(colour="black", size = 15), 
        legend.title = element_blank(), 
        legend.justification=c(1,0.1), 
        legend.position=c(1,0.1))
t_auc


