# ============================================================
# sim_cohort_b.R  模拟纵向队列生成器（论文 Cohort B，可完全复现）
# 来源：paper/script/11b_run_simcohort.R 的生成部分（逐点一致）
# 用途：再生成 demo/DEMO.DATA.xlsx 中的数据集（seed=20260813 时
#       与仓库 demo 文件逐位一致，运行 cal_allocation 得到 D=1.69e-05）
# 队列结构（与论文一致）：
#   5,500 样本 / 2,827 供体（约 20% 供体带 2-10 个纵向样本）
#   供体级特征：center（C1-C6，倾斜分布）、sex（F/M）
#   样本级特征：disease（D1-D9，访视特异，比例参照典型纵向研究）、
#               age（18-69，同供体统一为最小年龄后分 4 箱 A1-A4）
# ============================================================

generate_cohort_b <- function(seed = 20260813, target = 5500L) {
  set.seed(seed)
  center_p <- c(C1 = 0.35, C2 = 0.25, C3 = 0.15, C4 = 0.10, C5 = 0.08, C6 = 0.07)
  sex_p <- c(F = 0.55, M = 0.45)
  disease_p <- c(D1 = 0.162, D2 = 0.199, D3 = 0.198, D4 = 0.137, D5 = 0.113,
                 D6 = 0.097, D7 = 0.071, D8 = 0.011, D9 = 0.010)

  # 供体级抽取
  donor_n <- integer(0)
  donor_id <- character(0)
  donor_center <- character(0)
  donor_sex <- character(0)
  k <- 0L
  while (sum(donor_n) < target) {
    k <- k + 1L
    nd <- 1L
    if (runif(1) < 0.20) nd <- sample(2:10, 1)
    donor_n <- c(donor_n, nd)
    donor_id <- c(donor_id, sprintf("SUB%04d", k))
    donor_center <- c(donor_center, sample(names(center_p), 1, prob = center_p))
    donor_sex <- c(donor_sex, sample(names(sex_p), 1, prob = sex_p))
  }
  # 截断到恰好覆盖 target
  cut <- which(cumsum(donor_n) >= target)[1]
  donor_n <- donor_n[1:cut]
  donor_id <- donor_id[1:cut]
  donor_center <- donor_center[1:cut]
  donor_sex <- donor_sex[1:cut]

  # 样本级构建
  sampleID <- character(0); cohortID <- character(0)
  center_v <- character(0); sex_v <- character(0)
  disease_v <- character(0); age_v <- numeric(0)
  for (i in seq_along(donor_id)) {
    for (r in seq_len(donor_n[i])) {
      sampleID <- c(sampleID, sprintf("%s_V%02d", donor_id[i], r))
      cohortID <- c(cohortID, donor_id[i])
      center_v <- c(center_v, donor_center[i])
      sex_v <- c(sex_v, donor_sex[i])
      disease_v <- c(disease_v, sample(names(disease_p), 1, prob = disease_p))
      age_v <- c(age_v, sample(18:69, 1))
    }
  }
  data <- data.frame(sampleID = sampleID, cohortID = cohortID,
                     center = center_v, sex = sex_v,
                     disease = disease_v, age = age_v,
                     stringsAsFactors = FALSE)
  data <- data[seq_len(target), ]
  # 同供体年龄统一为最小年龄（纵向收样时间差异的常规处理）
  multi <- unique(data$cohortID[duplicated(data$cohortID)])
  for (m in multi) {
    ix <- which(data$cohortID == m)
    data$age[ix] <- min(data$age[ix])
  }
  # 年龄分箱
  data$age_bin <- NA_character_
  data$age_bin[data$age <= 27] <- "A1"
  data$age_bin[data$age >= 28 & data$age <= 37] <- "A2"
  data$age_bin[data$age >= 38 & data$age <= 47] <- "A3"
  data$age_bin[data$age >= 48] <- "A4"
  data$age <- NULL
  data[, c("sampleID", "cohortID", "center", "sex", "disease", "age_bin")]
}
