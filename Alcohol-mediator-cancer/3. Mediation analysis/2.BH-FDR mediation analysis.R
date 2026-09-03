
  load(file = "immune results/Mediation alcohol to cancer.Rdata")   # res_alcohol_cancer_median: 731性状 x 28对
  load(file = "immune results/alcohol to immune to cancer.Rdata")   # res_median_immune: 名义筛选后的候选(含se1/se2)
  
  # p1、p2 分别做BH-FDR
  res_alcohol_cancer_median$p1_FDR <- NA_real_
  res_alcohol_cancer_median$p2_FDR <- NA_real_
  key_xy <- paste(res_alcohol_cancer_median$X, res_alcohol_cancer_median$Y, sep = "|")
  for (pxy in unique(key_xy)) {
    idx <- which(key_xy == pxy)
    res_alcohol_cancer_median$p1_FDR[idx] <- p.adjust(res_alcohol_cancer_median$p1[idx], method = "BH")
    res_alcohol_cancer_median$p2_FDR[idx] <- p.adjust(res_alcohol_cancer_median$p2[idx], method = "BH")
  }
  
  #合并p1/p2的FDR值，并用delta法重算中介效应SE/CI/P
  res_median_immune <- merge(res_median_immune,
                             res_alcohol_cancer_median[, c("X","Y","M","p1_FDR","p2_FDR")],
                             by = c("X","Y","M"), all.x = TRUE)
  res_median_immune$ind_se_delta <- sqrt(res_median_immune$beta1^2 * res_median_immune$se2^2 +
                                           res_median_immune$beta2^2 * res_median_immune$se1^2)
  res_median_immune$ind_lo_delta <- res_median_immune$indirect_effect - 1.96 * res_median_immune$ind_se_delta
  res_median_immune$ind_up_delta <- res_median_immune$indirect_effect + 1.96 * res_median_immune$ind_se_delta
  res_median_immune$p_delta <- 2 * pnorm(-abs(res_median_immune$indirect_effect / res_median_immune$ind_se_delta))
  
  #对中介效应P值再做BH-FDR
  res_median_immune$p_delta_FDR <- NA_real_
  key_xy2 <- paste(res_median_immune$X, res_median_immune$Y, sep = "|")
  for (pxy in unique(key_xy2)) {
    idx <- which(key_xy2 == pxy)
    res_median_immune$p_delta_FDR[idx] <- p.adjust(res_median_immune$p_delta[idx], method = "BH")
  }
  
  #
  write.csv(res_median_immune, file = "immune results/mediation_revision_all.csv", row.names = FALSE)
  write.csv(res_median_immune[res_median_immune$Cancer == "Esophageal cancer", ],
            file = "immune results/mediation_revision_ESCA.csv", row.names = FALSE)
  