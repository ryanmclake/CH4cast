#################################################################
# CH4cast version 2                                             #
# Ryan McClure                                                  #
# COMPILE/ANALYZE/PLOT FORECASTS                                #
#################################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# DIRECTORY NEEDS TO BE CHANGED to ./forecast_output/
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

trap_all <- list.files(pattern = "ebullition_hidden_markov_forecast_alltrap_20") %>%
  map(readRDS) %>% 
  data.table::rbindlist(fill = T)%>%
  group_by(forecast_date)%>%
  mutate(days = seq_along(time))%>%
  group_by(forecast_date)%>%
  mutate(days = days-1)%>%
  filter(forecast_date>"2019-06-10")

trap_all_per_null <- list.files(pattern = "ebullition_null_persistence_forecast_alltrap_") %>%
  map(readRDS) %>% 
  data.table::rbindlist() %>%
  group_by(forecast_date)%>%
  mutate(days = seq_along(time))%>%
  group_by(forecast_date)%>%
  mutate(days = days-1)%>%
  filter(forecast_date>"2019-06-10")


trap_compare <- left_join(trap_all, trap_all_per_null, by = "time")


### Gelman rubin doagnostics ###
gelman_temp_model <- list.files(pattern = "temp_scale_gelman_diagnostics_") %>%
  map(readRDS) %>% 
  data.table::rbindlist()

gelman_ebu_model <- list.files(pattern = "ebu_model_gelman_diagnostics_") %>%
  map(readRDS) %>% 
  data.table::rbindlist()

### Parameter estimates ###
trap_all_parameters <- list.files(pattern = "ebullition_parameters_") %>%
  map(readRDS) %>% 
  data.table::rbindlist() %>%
  group_by(forecast_date)%>%
  summarize(mean_process = mean(sd.pro),
            sd_process = sd(sd.pro),
            mean_intercept = mean(mu2),
            sd_intercept = sd(mu2),
            mean_observe = mean(phi),
            sd_observe = sd(phi),
            mean_temp = mean(omega),
            sd_temp = sd(omega))


##### UNCERTAINTY PARTITIONING DATA ####

trap_all_IC <- list.files(pattern = "initial_condition_") %>%
  map(readRDS) %>% 
  data.table::rbindlist()%>%
  select(time, forecast_date, var)%>%
  rename(var_ic = var)

trap_all_PRO <- list.files(pattern = "model_process_") %>%
  map(readRDS) %>% 
  data.table::rbindlist()%>%
  select(time, forecast_date, var)%>%
  rename(var_pro = var)

trap_all_DRI <- list.files(pattern = "driver_data_") %>%
  map(readRDS) %>% 
  data.table::rbindlist()%>%
  select(time, forecast_date, var)%>%
  rename(var_dri = var)

trap_all_PAR <- list.files(pattern = "model_parameter_") %>%
  map(readRDS) %>% 
  data.table::rbindlist()%>%
  select(time, forecast_date, var)%>%
  rename(var_par = var)

partition <- cbind(trap_all_IC, trap_all_PRO[,3], trap_all_DRI[,3], trap_all_PAR[,3])

all_partitioned <- partition %>%
    mutate(sum_var = var_ic+var_pro+var_dri+var_par)%>%
    mutate(IC_proportion = var_ic/sum_var)%>%
    mutate(DRIVE_proportion = var_dri/sum_var)%>%
    mutate(PARA_proportion = var_par/sum_var)%>%
    mutate(PROCESS_proportion = var_pro/sum_var)%>%
    mutate(sum_check = PROCESS_proportion+PARA_proportion+DRIVE_proportion+IC_proportion)%>%
    select(time, forecast_date, IC_proportion, DRIVE_proportion, PARA_proportion, PROCESS_proportion, sum_check)%>%
  mutate(days = seq_along(sum_check))

all_partitioned_melt <- all_partitioned%>%
  select(-sum_check)%>%
  melt(., id = c("time","forecast_date","days"))


trap_all_partition <- left_join(trap_all, full_ebullition_model_alltrap, by = c("time"))
trap_null_partition <- left_join(trap_all_per_null, full_ebullition_model_alltrap, by = c("time"))
  

stats_all_bias_base <- trap_all_partition %>%
  select(time, mean, log_ebu_rate, forecast_date)%>%
  na.omit(.)%>%
  group_by(time)%>%
  filter(row_number(time) == 1) %>%
  mutate(Bias = abs(exp(mean)) - abs(exp(log_ebu_rate)))%>%
  mutate(model = "A: HHM forecast model")%>%
  ungroup(.)

base_model_NSE <- stats_all_bias_base %>%
  select(mean,log_ebu_rate)%>%
  summarize(NSE_model = NSE(exp(mean), exp(log_ebu_rate)))

stats_all_bias_null <- trap_null_partition %>%
  select(time, mean, log_ebu_rate, forecast_date)%>%
  na.omit(.)%>%
  group_by(time)%>%
  filter(row_number(time) == 1) %>%
  mutate(Bias = abs(exp(mean)) - abs(exp(log_ebu_rate)))%>%
  mutate(model = "B: Persistence null forecast model")%>%
  ungroup(.)

null_model_NSE <- stats_all_bias_null %>%
  select(mean,log_ebu_rate)%>%
  summarize(NSE_model = NSE(exp(mean), exp(log_ebu_rate)))

stats_all_bias <- bind_rows(stats_all_bias_base,stats_all_bias_null)

### FIGURES ###
### visualizations of the full ebullition forecasts (Planning to be figure 3)
ebullition_forecasts <- trap_all %>%
  #filter(days != 0)%>%
  ggplot(., aes(x = time, y = exp(mean), group = forecast_date)) +
  # geom_ribbon(aes(ymin = lower_90, ymax = upper_90), alpha = 0.2, fill = "midnightblue") +
  #geom_ribbon(aes(ymin = exp(lower_80), ymax = exp(upper_80)), alpha = 0.2, fill = "midnightblue") +
  geom_ribbon(aes(ymin = exp(lower_70), ymax = exp(upper_70)), alpha = 0.2, fill = "midnightblue") +
  geom_ribbon(aes(ymin = exp(lower_60), ymax = exp(upper_60)), alpha = 0.2, fill = "midnightblue") +
  geom_line(color = "purple4", size = 1, alpha = 0.7)+
  geom_pointrange(data = full_ebullition_model_alltrap, aes(x = time, y = exp(log_ebu_rate), ymin = exp(log_ebu_rate)-exp(log_ebu_rate_sd), ymax = exp(log_ebu_rate)+exp(log_ebu_rate_sd)), inherit.aes = FALSE, pch = 21, color = "red", fill = "red", cex = 0.5) +
  theme_bw()+
  labs(title = "A: HHM forecast model")+
  ylab(expression(paste("Ebullition Rate (mg CH "[4]," ",m^-2,"",d^-1,")")))+
  xlab("")+
  coord_cartesian(xlim=c(as.Date("2019-06-17"),as.Date("2019-11-17")), ylim = c(0,200))+
  theme(axis.text=element_text(size=15, color = "black"),
        axis.title=element_text(size=15, color = "black"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.title = element_blank(),
        title = element_text(size = 15),legend.position = "none",
        legend.text = element_text(size = 16, color = "black"))

null_forecasts <- trap_all_per_null %>%
  #filter(days != 0)%>%
  ggplot(aes(x = time, y = exp(mean), group = forecast_date)) +
  # geom_ribbon(aes(ymin = lower_90, ymax = upper_90), alpha = 0.2, fill = "midnightblue") +
  # geom_ribbon(aes(ymin = lower_80, ymax = upper_80), alpha = 0.2, fill = "midnightblue") +
  geom_ribbon(aes(ymin = exp(lower_70), ymax = exp(upper_70)), alpha = 0.2, fill = "midnightblue") +
  geom_ribbon(aes(ymin = exp(lower_60), ymax = exp(upper_60)), alpha = 0.2, fill = "midnightblue") +
  geom_line(color = "purple4", size = 1, alpha = 0.7)+
  geom_pointrange(data = full_ebullition_model_alltrap, aes(x = time, y = exp(log_ebu_rate), ymin = exp(log_ebu_rate)-exp(log_ebu_rate_sd), ymax = exp(log_ebu_rate)+exp(log_ebu_rate_sd)), inherit.aes = FALSE, pch = 21, color = "red", fill = "red", cex = 0.5) +
  theme_bw()+
  labs(title = "B: Persistence null forecast model")+
  ylab(expression(paste("Ebullition Rate (mg CH "[4]," ",m^-2,"",d^-1,")")))+
  xlab("")+
  coord_cartesian(xlim=c(as.Date("2019-06-17"),as.Date("2019-11-17")), ylim = c(0,200))+
  theme(axis.text=element_text(size=15, color = "black"),
        axis.title=element_text(size=15, color = "black"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.title = element_blank(),
        title = element_text(size = 15),legend.position = "none",
        legend.text = element_text(size = 16, color = "black"))


fig3 <- ebullition_forecasts/null_forecasts

jpeg("FIGURE_3.jpg", width = 600, height = 900)
fig3
dev.off()


daily_variance <- trap_all %>%
  filter(days != 0)%>%
  ggplot(., aes(x = days, y = var, group = days)) +
  geom_boxplot()+
  geom_jitter(aes(color = forecast_date), width = 0.1, size = 2)+
  theme_bw()+
  labs(title = "A: HHM forecast model uncertatiny")+
  ylab("Total forecast variance")+
  xlab("Days into future")+
  scale_x_continuous(limits = c(0,11), 
                   breaks = c(1,2,3,4,5,6,7,8,9,10), 
                   labels = c("1","2","3","4", "5", "6", "7", "8","9","10"))+
  theme(axis.text=element_text(size=15, color = "black"),
        axis.title=element_text(size=15, color = "black"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.title = element_blank(),
        title = element_text(size = 15),legend.position = c(0.1, 0.85),
        legend.text = element_text(size = 16, color = "black"))

daily_variance_null <- trap_all_per_null %>%
  filter(days != 0)%>%
  ggplot(., aes(x = days, y = sd, group = days)) +
  geom_boxplot()+
  geom_jitter(aes(color = forecast_date), width = 0.1, size = 2)+
  theme_bw()+
  labs(title = "B: Persistence null forecast model uncertatiny")+
  ylab("Total forecast variance")+
  xlab("Days into future")+
  scale_x_continuous(limits = c(0,11), 
                     breaks = c(1,2,3,4,5,6,7,8,9,10), 
                     labels = c("1","2","3","4", "5", "6", "7", "8","9","10"))+
  theme(axis.text=element_text(size=15, color = "black"),
        axis.title=element_text(size=15, color = "black"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.title = element_blank(),
        title = element_text(size = 15),legend.position = "none",
        legend.text = element_text(size = 16, color = "black"))


fig4 <- daily_variance/daily_variance_null

jpeg("FIGURE_4.jpg", width = 600, height = 900)
fig4
dev.off()


# Forecast bias figures
forecast_bias_box <- ggplot(stats_all_bias, aes(x = model, y = Bias)) +
  geom_boxplot()+
  geom_jitter(aes(color = time), width = 0.1, size = 2)+
  theme_bw()+
  geom_hline(yintercept = 0, lty = "dashed")+
  labs(title = "A: Forecast bias")+
  ylab("Forecast bias")+
  xlab("")+
  coord_cartesian(ylim = c(-30,30))+
  theme(axis.text=element_text(size=15, color = "black"),
        axis.title=element_text(size=15, color = "black"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.title = element_blank(),
        title = element_text(size = 15), legend.position = c(0.9,0.85),
        legend.text = element_text(size = 16, color = "black"))


ts_bias <- ggplot(stats_all_bias, aes(x = time, y = Bias, group = model)) +
  geom_point(aes(color = model), size = 3)+
  geom_smooth(aes(color = model),method = "loess", lwd = 2)+
  theme_bw()+
  geom_hline(yintercept = 0, lty = "dashed")+
  labs(title = "B: Time series of forecast bias")+
  ylab("Forecast bias")+
  xlab("")+
  coord_cartesian(ylim = c(-30,30))+
  theme(axis.text=element_text(size=15, color = "black"),
        axis.title=element_text(size=15, color = "black"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.title = element_blank(),
        title = element_text(size = 15),legend.position = c(0.7, 0.92),
        legend.text = element_text(size = 16, color = "black"))

fig5 <- forecast_bias_box/ts_bias

jpeg("FIGURE_5.jpg", width = 600, height = 900)
fig5
dev.off()




# PARAMETER ESTIAMTES FROM FORECASTS
process <- ggplot(trap_all_parameters, aes(x = forecast_date, y = mean_process)) +
  geom_ribbon(aes(ymin = mean_process-sd_process, ymax = mean_process+sd_process), alpha = 0.2, fill = "midnightblue") +
  geom_line(color = "black")+
  theme_bw()+
  labs(title = "D: Model process")+
  ylab(expression(paste(epsilon[t])))+
  xlab("")+
  coord_cartesian(xlim=c(as.Date("2019-05-27"),as.Date("2019-11-08")))+
  theme(axis.text=element_text(size=15, color = "black"),
        axis.title=element_text(size=15, color = "black"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.title = element_blank(),
        title = element_text(size = 15),legend.position = "none",
        legend.text = element_text(size = 16, color = "black"))

intercept <- ggplot(trap_all_parameters, aes(x = forecast_date, y = mean_intercept)) +
  geom_ribbon(aes(ymin = mean_intercept-sd_intercept, ymax = mean_intercept+sd_intercept), alpha = 0.2, fill = "midnightblue") +
  geom_line(color = "black")+
  theme_bw()+
  labs(title = "C: Intercept Parameter")+
  ylab(expression(paste(beta[0])))+
  xlab("")+
  coord_cartesian(xlim=c(as.Date("2019-05-27"),as.Date("2019-11-08")))+
  theme(axis.text=element_text(size=15, color = "black"),
        axis.title=element_text(size=15, color = "black"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.title = element_blank(),
        title = element_text(size = 15),legend.position = "none",
        legend.text = element_text(size = 16, color = "black"))

AR <- ggplot(trap_all_parameters, aes(x = forecast_date, y = mean_observe)) +
  geom_ribbon(aes(ymin = mean_observe-sd_observe, ymax = mean_observe+sd_observe), alpha = 0.2, fill = "midnightblue") +
  geom_line(color = "black")+
  theme_bw()+
  labs(title = "B: Autoregressive parameter")+
  ylab(expression(paste(beta[1])))+
  xlab("")+
  coord_cartesian(xlim=c(as.Date("2019-05-27"),as.Date("2019-11-08")))+
  theme(axis.text=element_text(size=15, color = "black"),
        axis.title=element_text(size=15, color = "black"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.title = element_blank(),
        title = element_text(size = 15),legend.position = "none",
        legend.text = element_text(size = 16, color = "black"))

temp <- ggplot(trap_all_parameters, aes(x = forecast_date, y = mean_temp)) +
  geom_ribbon(aes(ymin = mean_temp-sd_temp, ymax = mean_temp+sd_temp), alpha = 0.2, fill = "midnightblue") +
  geom_line(color = "black")+
  theme_bw()+
  labs(title = "A: Temperature parameter")+
  ylab(expression(paste(beta[2])))+
  xlab("")+
  coord_cartesian(xlim=c(as.Date("2019-05-27"),as.Date("2019-11-08")))+
  theme(axis.text=element_text(size=15, color = "black"),
        axis.title=element_text(size=15, color = "black"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.title = element_blank(),
        title = element_text(size = 15),legend.position = "none",
        legend.text = element_text(size = 16, color = "black"))

paramter = (temp+AR)/(intercept+process)


jpeg("FIGURE_6.jpg", width = 800, height = 800)
paramter
dev.off()



  #  PARTITION UNCERTATINY
  c <- ggplot(all_partitioned_melt, aes(x = days, y = value, group=interaction(forecast_date,variable), fill = interaction(forecast_date,variable)))+
           geom_area(aes(fill = variable))+
    scale_fill_manual(values = c("#E69F00", "#D55E00", "#CC79A7", "#56B4E9"))+
    theme_bw()+
    labs(title = "A: Daily forecast uncertainty")+
    ylab("Proportion to total variance")+
    xlab("Forecast cycle")+
    theme(axis.text=element_text(size=15, color = "black"),
          axis.title=element_text(size=15, color = "black"),
          axis.text.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          legend.title = element_blank(),
          title = element_text(size = 15), legend.position = "top",
          legend.text = element_text(size = 10, color = "black"))


  d <- ggplot(all_partitioned_melt, aes(x = variable, y = value, group=variable, fill = variable))+
    geom_boxplot()+
    scale_fill_manual(values = c("#E69F00", "#D55E00", "#CC79A7", "#56B4E9"))+
    theme_bw()+
    labs(title = "B: Aggregate forecast cycle uncertainty")+
    ylab("Proportion to total variance")+
    xlab("Source of Uncertainty")+
    theme(axis.text=element_text(size=15, color = "black"),
          axis.title=element_text(size=15, color = "black"),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          legend.title = element_blank(),
          title = element_text(size = 15), legend.position = "none",
          legend.text = element_text(size = 10, color = "black"))
  
  
  partition = c/d
  
  jpeg("FIGURE_7.jpg", width = 800, height = 800)
  partition
  dev.off()
