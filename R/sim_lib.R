# ============================================================
# 00_sim_lib.R  模拟研究共享库：场景定义 + 数据生成器（source 使用）
# ============================================================
# 注意：本文件为库文件，由调用脚本负责工作目录

# ---------- 场景数据生成器 ----------
gen_scenario <- function(scen_id, N, set_size, probs_list, cohort = NULL) {
  # probs_list: 命名list，每个特征一个概率向量；cohort: list(donor_probs, rep_range)
  if (is.null(cohort)) {
    data <- data.frame(sampleID = sprintf("S%05d", seq_len(N)))
    for (f in names(probs_list)) {
      data[[f]] <- sample(names(probs_list[[f]]), N, replace = TRUE,
                          prob = probs_list[[f]])
    }
    data$cohortID <- data$sampleID
  } else {
    # 先抽供体（供体级特征），再抽每供体的重复样本（样本级特征独立于供体）
    dprobs <- cohort$donor_probs
    donor_n <- 0L; donor_ids <- character(0); donor_feats <- list()
    target <- N
    while (sum(donor_n) < target) {
      # 逐个抽供体直到达到目标样本数
      nd <- 1L
      if (runif(1) < 0.2) nd <- sample(cohort$rep_range, 1)
      donor_n <- c(donor_n, nd)
      donor_ids <- c(donor_ids, sprintf("D%04d", length(donor_n)))
      df <- list()
      for (f in names(dprobs)) df[[f]] <- sample(names(dprobs[[f]]), 1, prob = dprobs[[f]])
      donor_feats[[length(donor_n)]] <- df
    }
    donor_n <- donor_n[1:which.max(cumsum(donor_n) >= target)]
    donor_ids <- donor_ids[seq_along(donor_n)]
    donor_feats <- donor_feats[seq_along(donor_n)]
    data <- data.frame()
    for (k in seq_along(donor_ids)) {
      for (r in seq_len(donor_n[k])) {
        row <- data.frame(sampleID = sprintf("%s_V%02d", donor_ids[k], r),
                          cohortID = donor_ids[k], stringsAsFactors = FALSE)
        for (f in names(dprobs)) row[[f]] <- donor_feats[[k]][[f]]
        for (f in names(probs_list)) row[[f]] <- sample(names(probs_list[[f]]), 1, prob = probs_list[[f]])
        data <- rbind(data, row)
      }
    }
    data <- data[seq_len(target), ]
  }
  data
}

# 场景定义（与方案一致）
scenarios <- list(
  S1 = list(N = 700, set_size = 35,
            probs = list(f1 = c(A = .5, B = .5), f2 = c(X = .8, Y = .2))),
  S2 = list(N = 700, set_size = 35,
            probs = list(f1 = c(A = .5, B = .5),
                         f2 = c(a = .4, b = .25, c = .15, d = .12, e = .08),
                         f3 = c(P = .6, Q = .4))),
  S3 = list(N = 700, set_size = 35,
            probs = list(region = setNames(rep(0.1, 10), paste0("R", 1:10)),
                         sex = c(M = .5, F = .5),
                         age = c(`30-40` = .10, `41-50` = .20, `51-60` = .30,
                                 `61-70` = .25, `71-80` = .15),
                         disease = c(ami = .15, ct = .10, is = .20, lung = .25, sto = .30))),
  S4 = list(N = 5600, set_size = 35,
            probs = list(region = setNames(rep(0.1, 10), paste0("R", 1:10)),
                         sex = c(M = .5, F = .5),
                         age = c(`30-40` = .10, `41-50` = .20, `51-60` = .30,
                                 `61-70` = .25, `71-80` = .15),
                         disease = c(ami = .15, ct = .10, is = .20, lung = .25, sto = .30))),
  S5 = list(N = 700, set_size = 35,
            probs = list(group = c(G1 = .1, G2 = .1, G3 = .1, G4 = .1, G5 = .1,
                                  G6 = .1, G7 = .1, G8 = .1, G9 = .2),
                         age = setNames(rep(0.1, 10), paste0("A", 1:10))),
            cohort = list(donor_probs = list(center = c(C1 = .35, C2 = .25, C3 = .15,
                                                        C4 = .10, C5 = .08, C6 = .07),
                                             sex = c(F = .55, M = .45)),
                          rep_range = 2:10)),
  S6 = list(N = 5600, set_size = 35,
            probs = list(group = c(G1 = .1, G2 = .1, G3 = .1, G4 = .1, G5 = .1,
                                  G6 = .1, G7 = .1, G8 = .1, G9 = .2),
                         age = setNames(rep(0.1, 10), paste0("A", 1:10))),
            cohort = list(donor_probs = list(center = c(C1 = .35, C2 = .25, C3 = .15,
                                                        C4 = .10, C5 = .08, C6 = .07),
                                             sex = c(F = .55, M = .45)),
                          rep_range = 2:10))
)

