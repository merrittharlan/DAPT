#Merritt Harlan
#DAPT PAD
library(bamr)
library(ggplot2)
library(viridis)
library(tidyr)
library(ggsci)
library(scales)
library(ggpubr)

#Test PAD dataset------------
Athabasca_rch1_data <- readRDS(file = "Athabasca_rch1_data.rds")
Athabasca_rch2_data <- readRDS(file = "Athabasca_rch2_data.rds")
Coupe_data <- readRDS(file = "Coupe_data.rds")
Embarras_data <- readRDS(file = "Embarras_data.rds")
Fletcher_data <- readRDS(file = "Fletcher_data.rds")
Limon_data <- readRDS(file = "Limon_data.rds")
Mamawi_data <- readRDS(file = "Mamawi_data.rds")
QF_data <- readRDS(file = "QF_data.rds")
Richardson_data <- readRDS(file = "Richardson_data.rds")
Rochers_rch1_data <- readRDS(file = "Rochers_rch1_data.rds")
Rochers_rch2_data <- readRDS(file = "Rochers_rch2_data.rds")

Athabasca_rch1_Q_test <- test_Q_num(bamdata = Athabasca_rch1_data)
Athabasca_rch2_Q_test <- test_Q_num(bamdata = Athabasca_rch2_data)
Coupe_Q_test <- test_Q_num(bamdata = Coupe_data)
Embarras_Q_test <- test_Q_num(bamdata = Embarras_data)
Fletcher_Q_test <- test_Q_num(bamdata = Fletcher_data)
Limon_Q_test <- test_Q_num(bamdata = Limon_data)
Mamawi_Q_test <- test_Q_num(bamdata = Mamawi_data)
QF_Q_test <- test_Q_num(bamdata = QF_data)
Richardson_Q_test <- test_Q_num(bamdata = Richardson_data)
Rochers_rch1_Q_test <- test_Q_num(bamdata = Rochers_rch1_data)
Rochers_rch2_Q_test <- test_Q_num(bamdata = Rochers_rch2_data)

Athabasca_rch1_Q_num = num_Q_df(num_Q_list = Athabasca_rch1_Q_test, rivname = "Athabasca R1")
Athabasca_rch2_Q_num = num_Q_df(num_Q_list = Athabasca_rch2_Q_test, rivname = "Athabasca R2")
Coupe_Q_num = num_Q_df(num_Q_list = Coupe_Q_test, rivname = "Coupe")
Embarras_Q_num = num_Q_df(num_Q_list = Embarras_Q_test, rivname = "Embarras")
Fletcher_Q_num = num_Q_df(num_Q_list = Fletcher_Q_test, rivname = "Fletcher")
Limon_Q_num = num_Q_df(num_Q_list = Limon_Q_test, rivname = "Limon")
Mamawi_Q_num = num_Q_df(num_Q_list = Mamawi_Q_test, rivname = "Mamawi")
QF_Q_num = num_Q_df(num_Q_list = QF_Q_test, rivname = "QF")
Richardson_Q_num = num_Q_df(num_Q_list = Richardson_Q_test, rivname = "Richardson")
Rochers_rch1_Q_num = num_Q_df(num_Q_list = Rochers_rch1_Q_test, rivname = "Rochers R1")
Rochers_rch2_Q_num = num_Q_df(num_Q_list = Rochers_rch2_Q_test, rivname = "Rochers R2")

Athabasca_rch1_Q_id = Q_ids(Athabasca_rch1_Q_test)
Athabasca_rch2_Q_id = Q_ids(Athabasca_rch2_Q_test)
Coupe_Q_id = Q_ids(Coupe_Q_test)
Embarras_Q_id = Q_ids(Embarras_Q_test)
Fletcher_Q_id = Q_ids(Fletcher_Q_test)
Limon_Q_id = Q_ids(Limon_Q_test)
Mamawi_Q_id = Q_ids(Mamawi_Q_test)
QF_Q_id = Q_ids(QF_Q_test)
Richardson_Q_id = Q_ids(Richardson_Q_test)
Rochers_rch1_Q_id = Q_ids(Rochers_rch1_Q_test)
Rochers_rch2_Q_id = Q_ids(Rochers_rch2_Q_test)


#Plot hydrographs-----
Athabasca_R1_hydro = plot_pad_hydro(num_Q_list = Athabasca_rch1_Q_test, bamdata = Athabasca_rch1_data, 
                                    Q_df = Q_hydro(Athabasca_rch1_Q_test, rivname = "Athabasca R1"),
                                    Q_id = Athabasca_rch1_Q_id, rivname = "Athabasca R1")+
  annotate("text", x = 2, y = 835, label = "NSE = 0.99", size = 3)
Athabasca_R2_hydro = plot_pad_hydro(num_Q_list = Athabasca_rch2_Q_test, bamdata = Athabasca_rch2_data, 
                                    Q_df = Q_hydro(Athabasca_rch2_Q_test, rivname = "Athabasca R2"),
                                    Q_id = Athabasca_rch2_Q_id, rivname = "Athabasca R2")+
  annotate("text", x = 30, y = 1000, label = "NSE = 0.95", size = 3)
Coupe_hydro = plot_pad_hydro(num_Q_list = Coupe_Q_test, bamdata = Coupe_data, 
                             Q_df = Q_hydro(Coupe_Q_test, rivname = "Coupe"),
                             Q_id = Coupe_Q_id,rivname = "Coupe")+
  annotate("text", x = 42, y = 117, label = "NSE = 0.71", size = 3)
Embarras_hydro = plot_pad_hydro(num_Q_list = Embarras_Q_test, bamdata = Embarras_data, 
                                Q_df = Q_hydro(Embarras_Q_test, rivname = "Embarras"),
                                Q_id = Embarras_Q_id, rivname = "Embarras")+
  annotate("text", x = 6, y = 133, label = "NSE = 0.98", size = 3)
Fletcher_hydro = plot_pad_hydro(num_Q_list = Fletcher_Q_test, bamdata = Fletcher_data, 
                                Q_df = Q_hydro(Fletcher_Q_test, rivname = "Fletcher"),
                                Q_id = Fletcher_Q_id,rivname = "Fletcher")+
  annotate("text", x = 34, y = 128, label = "NSE = 0.93", size = 3)
Limon_hydro = plot_pad_hydro(num_Q_list = Limon_Q_test, bamdata = Limon_data, 
                             Q_df = Q_hydro(Limon_Q_test, rivname = "Limon"),
                             Q_id = Limon_Q_id,rivname = "Limon")+
  annotate("text", x = 3, y = 36, label = "NSE = 0.98", size = 3)
Mamawi_hydro = plot_pad_hydro(num_Q_list = Mamawi_Q_test, bamdata = Mamawi_data, 
                              Q_df = Q_hydro(Mamawi_Q_test, rivname = "Mamawi"),
                              Q_id = Mamawi_Q_id,rivname = "Mamawi")+
  annotate("text", x = 7, y = 137, label = "NSE = 0.88", size = 3)
QF_hydro = plot_pad_hydro(num_Q_list = QF_Q_test, bamdata = QF_data, 
                          Q_df = Q_hydro(QF_Q_test, rivname = "Quatres Fourches"),
                          Q_id = QF_Q_id,rivname = "Quatres Fourches")+
  annotate("text", x = 75, y = 390, label = "NSE = 0.97", size = 3)

Richardson_hydro = plot_pad_hydro(num_Q_list = Richardson_Q_test, bamdata = Richardson_data, 
                                  Q_df = Q_hydro(Richardson_Q_test, rivname = "Richardson"),
                                  Q_id = Richardson_Q_id, rivname = "Richardson")+
  annotate("text", x = 16, y = 20, label = "NSE = 0.22", size = 3)
Rochers_R1_hydro = plot_pad_hydro(num_Q_list = Rochers_rch1_Q_test, bamdata = Rochers_rch1_data, 
                                  Q_df = Q_hydro(Rochers_rch1_Q_test, rivname = "Rochers R1"),
                                  Q_id = Rochers_rch1_Q_id,rivname = "Rochers R1")+
  annotate("text", x = 12, y = 1700, label = "NSE = 0.95", size = 3)

Rochers_R2_hydro = plot_pad_hydro(num_Q_list = Rochers_rch2_Q_test, bamdata = Rochers_rch2_data, 
                                  Q_df = Q_hydro(Rochers_rch2_Q_test, rivname = "Rochers R2"),
                                  Q_id = Rochers_rch2_Q_id,rivname = "Rochers R2")+
  annotate("text", x = 18, y = 1535, label = "NSE = 0.95", size = 3)

png("PAD_val_hydro.png",
    width = 6.5,
    height = 5,
    units = 'in',
    res = 300
)
PAD_val_hydro = 
  ggarrange(
    ggarrange(
      Athabasca_R1_hydro,
      Athabasca_R2_hydro,
      Coupe_hydro,
      Embarras_hydro,
      Fletcher_hydro,
      Limon_hydro,
      Mamawi_hydro,
      QF_hydro,
      ncol = 4,
      nrow = 2,
      labels = c("A", "B", "C", "D", "E", "F", "G", "H")
    ),
    ggarrange(
      Richardson_hydro,
      Rochers_R1_hydro,
      Rochers_R2_hydro,
      ncol = 4, 
      common.legend = TRUE,
      legend = "right",
      labels = c("I", "J", "K"),
      widths = c(1,1,1,0.3)
    ),
    nrow = 2,
    heights = c(2,1)
  )
annotate_figure(PAD_val_hydro,
                bottom = text_grob("Number of ADCP Measurements", size = 10),
                left = text_grob(expression(paste(Discharge~(m ^ 3 / s))), size = 10, rot = 90))
dev.off()


#Plot Full hydrographs-------
Athabasca_rch1_full_data <- readRDS(file = "Athabasca_rch1_full_data.rds")
Athabasca_rch1_full_data <- readRDS(file = "Athabasca_rch1_full_data.rds")
Coupe_full_data <- readRDS(file = "Coupe_full_data.rds")
Embarras_full_data <- readRDS(file = "Embarras_full_data.rds")
Fletcher_full_data <- readRDS(file = "Fletcher_full_data.rds")
Limon_full_data <- readRDS(file = "Limon_full_data.rds")
Mamawi_full_data <- readRDS(file = "Mamawi_full_data.rds")
QF_full_data <- readRDS(file = "QF_full_data.rds")
Richardson_full_data <- readRDS(file = "Richardson_full_data.rds")
Rochers_rch1_full_data <- readRDS(file = "Rochers_rch1_full_data.rds")
Rochers_rch2_full_data <- readRDS(file = "Rochers_rch2_full_data.rds")

Athabasca_rch1_WSE <- readRDS(file = "Athabasca_rch1_WSE.rds")
Athabasca_rch1_WSE <- readRDS(file = "Athabasca_rch1_WSE.rds")
Coupe_WSE <- readRDS(file = "Coupe_WSE.rds")
Embarras_WSE <- readRDS(file = "Embarras_WSE.rds")
Fletcher_WSE <- readRDS(file = "Fletcher_WSE.rds")
Limon_WSE <- readRDS(file = "Limon_WSE.rds")
Mamawi_WSE <- readRDS(file = "Mamawi_WSE.rds")
QF_WSE <- readRDS(file = "QF_WSE.rds")
Richardson_WSE <- readRDS(file = "Richardson_WSE.rds")
Rochers_rch1_WSE <- readRDS(file = "Rochers_rch1_WSE.rds")
Rochers_rch2_WSE <- readRDS(file = "Rochers_rch2_WSE.rds")

Athabasca_rch1_corr_Q <- readRDS(file = "Athabasca_rch1_corr_Q.rds")
Athabasca_rch1_corr_Q <- readRDS(file = "Athabasca_rch1_corr_Q.rds")
Coupe_corr_Q <- readRDS(file = "Coupe_corr_Q.rds")
Embarras_corr_Q <- readRDS(file = "Embarras_corr_Q.rds")
Fletcher_corr_Q <- readRDS(file = "Fletcher_corr_Q.rds")
Limon_corr_Q <- readRDS(file = "Limon_corr_Q.rds")
Mamawi_corr_Q <- readRDS(file = "Mamawi_corr_Q.rds")
QF_corr_Q <- readRDS(file = "QF_corr_Q.rds")
Richardson_corr_Q <- readRDS(file = "Richardson_corr_Q.rds")
Rochers_rch1_corr_Q <- readRDS(file = "Rochers_rch1_corr_Q.rds")
Rochers_rch2_corr_Q <- readRDS(file = "Rochers_rch2_corr_Q.rds")

#Estimate Q w/ BAM
Athabasca_rch1_Q_all = test_Q_all(bamdata = Athabasca_rch1_full_data, Q_ADCP = exp(Athabasca_rch1_data$logQ_hat),
                                  Q_time = as.numeric(as.POSIXct(Athabasca_rch1_corr_Q$new_date)), 
                                  WSE_time = as.numeric(na.omit(as.POSIXct(unique(as.Date(Athabasca_rch1_WSE$new_time))))))
Athabasca_rch2_Q_all = test_Q_all(bamdata = Athabasca_rch2_full_data, Q_ADCP = exp(Athabasca_rch2_data$logQ_hat),
                                  Q_time = as.numeric(as.POSIXct(Athabasca_rch2_corr_Q$new_date)), 
                                  WSE_time = as.numeric(na.omit(as.POSIXct(unique(as.Date(Athabasca_rch2_WSE$new_time))))))
Coupe_Q_all = test_Q_all(bamdata = Coupe_full_data, Q_ADCP = exp(Coupe_data$logQ_hat),
                         Q_time = as.numeric(as.POSIXct(Coupe_corr_Q$new_date)), 
                         WSE_time = as.numeric(na.omit(as.POSIXct(unique(as.Date(Coupe_WSE$new_time))))))
Embarras_Q_all = test_Q_all(bamdata = Embarras_full_data, Q_ADCP = exp(Embarras_data$logQ_hat),
                                 Q_time = as.numeric(as.POSIXct(Embarras_corr_Q$new_date)), 
                                 WSE_time = as.numeric(na.omit(as.POSIXct(unique(as.Date(Embarras_WSE$new_time))))))
Fletcher_Q_all = test_Q_all(bamdata = Fletcher_full_data, Q_ADCP = exp(Fletcher_data$logQ_hat),
                            Q_time = as.numeric(as.POSIXct(Fletcher_corr_Q$new_date)), 
                            WSE_time = as.numeric(na.omit(as.POSIXct(unique(as.Date(Fletcher_WSE$new_time, format = "%m/%d/%Y"))))))
Limon_Q_all = test_Q_all(bamdata = Limon_full_data, Q_ADCP = exp(Limon_data$logQ_hat),
                         Q_time = as.numeric(as.POSIXct(Limon_corr_Q$new_date)), 
                         WSE_time = as.numeric(na.omit(as.POSIXct(unique(as.Date(Limon_WSE$new_time))))))
Mamawi_Q_all = test_Q_all(bamdata = Mamawi_full_data, Q_ADCP = exp(Mamawi_data$logQ_hat),
                          Q_time = as.numeric(as.POSIXct(Mamawi_corr_Q$new_date)), 
                          WSE_time = as.numeric(na.omit(as.POSIXct(unique(as.Date(Mamawi_WSE$new_time))))))
QF_Q_all = test_Q_all(bamdata = QF_full_data, Q_ADCP = exp(QF_data$logQ_hat),
                      Q_time = as.numeric(as.POSIXct(QF_corr_Q$new_date)), 
                      WSE_time = as.numeric(na.omit(as.POSIXct(unique(as.Date(QF_WSE$new_time))))))
Richardson_Q_all = test_Q_all(bamdata = Richardson_full_data, Q_ADCP = exp(Richardson_data$logQ_hat),
                              Q_time = as.numeric(as.POSIXct(Richardson_corr_Q$new_date)), 
                              WSE_time = as.numeric(na.omit(as.POSIXct(unique(as.Date(Richardson_WSE$new_time))))))
Rochers_rch1_Q_all = test_Q_all(bamdata = Rochers_rch1_full_data, Q_ADCP = exp(Rochers_rch1_data$logQ_hat),
                                Q_time = as.numeric(as.POSIXct(Rochers_rch1_corr_Q$new_date)), 
                                WSE_time = as.numeric(na.omit(as.POSIXct(unique(as.Date(Rochers_rch1_WSE$new_time, format = "%m/%d/%Y"))))))
Rochers_rch2_Q_all = test_Q_all(bamdata = Rochers_rch2_full_data, Q_ADCP = exp(Rochers_rch2_data$logQ_hat),
                                Q_time = as.numeric(as.POSIXct(Rochers_rch2_corr_Q$new_date)), 
                                WSE_time = as.numeric(na.omit(as.POSIXct(unique(as.Date(Rochers_rch2_WSE$new_time, format = "%m/%d/%Y"))))))

Athabasca_rch1_full_plot = plot_full_hydro(full_BAM = Athabasca_rch1_Q_all,
                                           rivname = "Athabasca R1",
                                           WSE_all = Athabasca_rch1_WSE, Q = Athabasca_rch1_corr_Q)
Athabasca_rch2_full_plot = plot_full_hydro(full_BAM = Athabasca_rch2_Q_all, 
                                           rivname = "Athabasca R2", 
                                           WSE_all = Athabasca_rch2_WSE, Q = Athabasca_rch2_corr_Q)
Coupe_full_plot = plot_full_hydro(full_BAM = Coupe_Q_all, 
                                  rivname = "Coupe", 
                                  WSE_all = Coupe_WSE, Q = Coupe_corr_Q)
Embarras_full_plot = plot_full_hydro(full_BAM = Embarras_Q_all, 
                                          rivname = "Embarras", 
                                          WSE_all = Embarras_WSE, Q = Embarras_corr_Q)
Fletcher_full_plot = plot_full_hydro(full_BAM = Fletcher_Q_all, 
                                     rivname = "Fletcher", 
                                     WSE_all = Fletcher_WSE, Q = Fletcher_corr_Q)
Limon_full_plot = plot_full_hydro(full_BAM = Limon_Q_all, 
                                  rivname = "Limon", 
                                  WSE_all = Limon_WSE, Q = Limon_corr_Q)
Mamawi_full_plot = plot_full_hydro(full_BAM = Mamawi_Q_all, 
                                   rivname = "Mamawi", 
                                   WSE_all = Mamawi_WSE, Q = Mamawi_corr_Q)
QF_full_plot = plot_full_hydro(full_BAM = QF_Q_all, 
                               rivname = "QF", 
                               WSE_all = QF_WSE, Q = QF_corr_Q)
Richardson_full_plot = plot_full_hydro(full_BAM = Richardson_Q_all, 
                                       rivname = "Richardson", 
                                       WSE_all = Richardson_WSE, Q = Richardson_corr_Q)

Rochers_rch1_full_plot = plot_full_hydro(full_BAM = Rochers_rch1_Q_all, 
                                         rivname = "Rochers R1", 
                                         WSE_all = Rochers_rch1_WSE, Q = Rochers_rch1_corr_Q)

Rochers_rch2_full_plot = plot_full_hydro(full_BAM = Rochers_rch2_Q_all, 
                                         rivname = "Rochers R2",
                                         WSE_all = Rochers_rch2_WSE, Q = Rochers_rch2_corr_Q)
PAD_hydro = 
  png("PAD_hydro.png",
      width = 6.5,
      height = 6.5,
      units = 'in',
      res = 300
  )
PAD_full_plot = 
  ggarrange(
    Athabasca_rch1_full_plot,
    Athabasca_rch2_full_plot,
    Coupe_full_plot,
    Embarras_full_plot,
    Fletcher_full_plot,
    Limon_full_plot,
    Mamawi_full_plot,
    QF_full_plot,
    Richardson_full_plot,
    Rochers_rch1_full_plot,
    Rochers_rch2_full_plot,
    nrow = 4,
    ncol = 3,
    labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K"),
    common.legend = TRUE,
    legend = "bottom")

annotate_figure(PAD_full_plot,
                left = text_grob(expression(paste(Discharge~(m ^ 3 / s))), size = 10, rot = 90))
dev.off()


