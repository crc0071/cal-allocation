# ============================================================
# 02_metrics.R  评价指标函数库（source 使用）
# 指标定义（论文 Methods 2.5）：
#   主指标：
#   1. set_deviation d(set)：加权归一比例偏差（L1/全变差口径）
#      d = Σ_f w_f × (1/L_f) × Σ_l |p̂_{f,l} − p_{f,l}|
#      p̂ 为 set 内水平比例（缺失水平按 0），p 为全队列理论比例
#   2. all_sets_diff_sd：所有 set 的 d 值标准差（主指标）
#   次要指标（敏感性分析用）：
#   3. chi2_deviation：χ² 口径偏差
#      χ² = Σ_f w_f × Σ_l (n_{s,f,l} − n_s·p_{f,l})² / (n_s·p_{f,l})
#      稀有水平抽样方差小、自动获得更大权重
#   4. hellinger_deviation：Hellinger 距离口径
#      H² = Σ_f w_f × ½ Σ_l (√p̂_{f,l} − √p_{f,l})²，有界且对零计数稳健
#   5. max_d：所有 set 中 d 的最大值（minimax 口径，直接控制最差 set）
#   其余：
#   6. level_rsd：每个特征水平在各 set 中频数的 RSD
#      RSD_l = sd(Freq_l)/mean(Freq_l)，仅对包含该水平的 set 计算
#   7. coverage_penalty：每个 set 的水平覆盖惩罚
#      penalty = Σ_f (L_f − 覆盖水平数_f) × L_f；全覆盖=0
#   8. cohort_integrity：cohort 分裂数（同一 cohort 分布到多个 set）
# ============================================================
# 注意：本文件为库文件，由调用脚本负责工作目录

get_levels <- function(x) {
  if (is.factor(x)) return(levels(x))
  if (is.numeric(x)) return(sort(unique(x)))
  sort(unique(as.character(x)))
}

#' 单个 set 的加权比例偏差 d
#' @param set_data 该 set 的样本子表（数据框）
#' @param data     全队列样本表（用于理论比例与水平全集）
#' @param features 目标特征列名
#' @param weight   命名 list 权重（缺失默认1）
set_deviation <- function(set_data, data, features, weight = NULL) {
  wL <- vapply(features, function(f) {
    lev <- get_levels(data[[f]])
    wf <- if (is.null(weight) || is.null(weight[[f]])) 1 else as.numeric(weight[[f]])
    wf / length(lev)
  }, numeric(1))
  n <- nrow(set_data)
  if (n == 0) {
    return(sum(vapply(features, function(f) wL[f], numeric(1))))
  }
  total <- 0
  for (k in seq_along(features)) {
    f <- features[k]
    lev <- get_levels(data[[f]])
    p_all <- tabulate(match(data[[f]], lev), nbins = length(lev)) / nrow(data)
    p_set <- tabulate(match(set_data[[f]], lev), nbins = length(lev)) / n
    total <- total + wL[k] * sum(abs(p_set - p_all))
  }
  total
}

#' 所有 set 的 d 值与 all_sets_diff_sd
#' @param assignment 含 new_set 列的分配结果数据框
set_deviations <- function(assignment, data, features, weight = NULL) {
  sets <- sort(unique(assignment$new_set))
  dvec <- vapply(sets, function(s)
    set_deviation(assignment[assignment$new_set == s, , drop = FALSE],
                  data, features, weight), numeric(1))
  names(dvec) <- sets
  dvec
}

all_sets_diff_sd <- function(dvec) {
  if (length(dvec) < 2) return(NA_real_)
  sd(dvec)
}

#' χ² 口径偏差（次要指标）：每个 set 的卡方统计量
#' χ² = Σ_f w_f × Σ_l (n_l − n·p_l)² / (n·p_l)
chi2_deviation <- function(set_data, data, features, weight = NULL) {
  n <- nrow(set_data)
  if (n == 0) return(Inf)
  total <- 0
  for (f in features) {
    lev <- get_levels(data[[f]])
    p_all <- tabulate(match(data[[f]], lev), nbins = length(lev)) / nrow(data)
    cnt <- tabulate(match(set_data[[f]], lev), nbins = length(lev))
    wf <- if (is.null(weight) || is.null(weight[[f]])) 1 else as.numeric(weight[[f]])
    total <- total + wf * sum((cnt - n * p_all)^2 / (n * p_all))
  }
  total
}

#' Hellinger 口径偏差（次要指标）：H² = Σ_f w_f × ½ Σ_l (√p̂ − √p)²
hellinger_deviation <- function(set_data, data, features, weight = NULL) {
  n <- nrow(set_data)
  if (n == 0) {
    return(sum(vapply(features, function(f) {
      wf <- if (is.null(weight) || is.null(weight[[f]])) 1 else as.numeric(weight[[f]])
      wf * 0.5
    }, numeric(1))))
  }
  total <- 0
  for (f in features) {
    lev <- get_levels(data[[f]])
    p_all <- tabulate(match(data[[f]], lev), nbins = length(lev)) / nrow(data)
    p_set <- tabulate(match(set_data[[f]], lev), nbins = length(lev)) / n
    wf <- if (is.null(weight) || is.null(weight[[f]])) 1 else as.numeric(weight[[f]])
    total <- total + wf * 0.5 * sum((sqrt(p_set) - sqrt(p_all))^2)
  }
  total
}

#' minimax 口径（次要指标）：所有 set 中 d 的最大值
max_d <- function(dvec) {
  if (length(dvec) == 0) return(NA_real_)
  max(dvec)
}

#' 每个特征水平的跨 set 频数 RSD（默认只对包含该水平的 set 计算）
#' @return data.frame(level, feature, rsd, n_sets_with)
level_rsd <- function(assignment, feature, zero_handling = c("drop", "zero")) {
  zero_handling <- match.arg(zero_handling)
  sets <- sort(unique(assignment$new_set))
  lev <- get_levels(assignment[[feature]])
  out <- data.frame(feature = feature, level = lev, rsd = NA_real_,
                    n_sets_with = NA_integer_, stringsAsFactors = FALSE)
  for (k in seq_along(lev)) {
    l <- lev[k]
    freq <- vapply(sets, function(s)
      sum(assignment[[feature]][assignment$new_set == s] == l), numeric(1))
    if (zero_handling == "drop") {
      keep <- freq > 0
    } else {
      keep <- rep(TRUE, length(freq))
    }
    out$n_sets_with[k] <- sum(keep)
    if (sum(keep) >= 2 && mean(freq[keep]) > 0) {
      out$rsd[k] <- sd(freq[keep]) / mean(freq[keep])
    }
  }
  out
}

#' 每个 set 的水平覆盖惩罚；全覆盖=0
#' @return data.frame(new_set, penalty)
coverage_penalty <- function(assignment, features) {
  sets <- sort(unique(assignment$new_set))
  L <- vapply(features, function(f) length(get_levels(assignment[[f]])), integer(1))
  out <- vapply(sets, function(s) {
    sub <- assignment[assignment$new_set == s, , drop = FALSE]
    pen <- 0
    for (k in seq_along(features)) {
      covered <- length(unique(sub[[features[k]]]))
      pen <- pen + (L[k] - covered) * L[k]
    }
    pen
  }, numeric(1))
  data.frame(new_set = sets, penalty = out)
}

#' cohort 完整性：被分裂到多个 set 的 cohort 列表
#' @return data.frame(cohort, n_sets, n_samples)
cohort_integrity <- function(assignment, cohort_col = "cohortID") {
  if (!(cohort_col %in% colnames(assignment))) {
    return(data.frame(cohort = character(0), n_sets = integer(0),
                      n_samples = integer(0)))
  }
  tab <- table(assignment[[cohort_col]], assignment$new_set) > 0
  nsets <- rowSums(tab)
  nsamples <- table(assignment[[cohort_col]])
  split_ids <- names(nsets)[nsets > 1]
  data.frame(cohort = split_ids,
             n_sets = as.integer(nsets[split_ids]),
             n_samples = as.integer(nsamples[split_ids]),
             stringsAsFactors = FALSE)
}

#' 汇总指标（供模拟对比与案例评价统一调用）
#' @return list(n_sets, set_size_range, mean_d, all_sets_diff_sd,
#'              chi2_sd, hellinger_sd, max_d,
#'              level_rsd 明细表, rsd_median, rsd_max, coverage 明细,
#'              coverage_rate, n_cohort_split)
compute_all_metrics <- function(assignment, data, features, weight = NULL,
                                cohort_col = "cohortID") {
  sets <- sort(unique(assignment$new_set))
  dvec <- set_deviations(assignment, data, features, weight)
  chi2_vec <- vapply(sets, function(s)
    chi2_deviation(assignment[assignment$new_set == s, , drop = FALSE],
                   data, features, weight), numeric(1))
  hell_vec <- vapply(sets, function(s)
    hellinger_deviation(assignment[assignment$new_set == s, , drop = FALSE],
                        data, features, weight), numeric(1))
  rsd_all <- do.call(rbind, lapply(features, function(f) level_rsd(assignment, f)))
  rsd_all$rsd[is.nan(rsd_all$rsd)] <- NA
  cov_all <- coverage_penalty(assignment, features)
  spl <- cohort_integrity(assignment, cohort_col)
  list(n_sets = length(sets),
       set_size_range = range(table(assignment$new_set)),
       mean_d = mean(dvec),
       all_sets_diff_sd = all_sets_diff_sd(dvec),
       chi2_sd = sd(chi2_vec),
       hellinger_sd = sd(hell_vec),
       max_d = max_d(dvec),
       d_per_set = dvec,
       level_rsd = rsd_all,
       rsd_median = median(rsd_all$rsd, na.rm = TRUE),
       rsd_max = max(rsd_all$rsd, na.rm = TRUE),
       coverage = cov_all,
       coverage_rate = mean(cov_all$penalty == 0),
       n_cohort_split = nrow(spl),
       split_detail = spl)
}
