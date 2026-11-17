# ============================================================
# 01_cal_allocation.R  Cal-Allocation 算法核心函数库
# 来源：Cal26-27安定/code2.R（唯一算法事实来源）+ cal-ballance-v2.R
# 本文件只定义函数，source 后使用，不直接运行
#
# 五阶段算法：
#   步骤1：combo 生成（特征水平组合，去空，按样本数降序）
#   步骤2：comb unit 构造（贪心全水平覆盖）
#   步骤3：轮抽关键样本（每空 set 从 unit 每 combo 抽一个 cohort，
#          抽光补位；此阶段样本 key=1，交换优化不动）
#   步骤4：填满（按 d(set) 增加最小的 combo 优先）
#   步骤5：交换优化（等行数、key=0 的 cohort 两两交换，
#          贪心执行使 all_sets_diff_sd 下降的最优交换，直至无改进）
#
# 与 code2.R 的差异（论文版修正）：
#   1. 全局 cohort 原子化：同一 cohortID 的所有行（即使跨 combo，
#      如纵向样本 group 不同）作为一个不可分单元，combo 归属按代表行
#   2. comb unit 全覆盖判定只基于 combo 特征（修复 sample.n 列误入）
#   3. 交换阶段多轮重复直至无改进（code2.R 因 return 位置只做了单轮 s1）
#   4. 加 set.seed，结果完全可复现
#   5. 填满阶段 trial 用整个 cohort 的行（code2.R 只用首行）
#   6. 权重按特征名查找（code2.R 按位置）
# legacy=TRUE 时逐点复刻 code2.R 行为（供 03_parity_check 对拍）
# ============================================================
# 注意：本文件为库文件，由调用脚本负责工作目录

get_levels <- function(x) {
  if (is.factor(x)) return(levels(x))
  if (is.numeric(x)) return(sort(unique(x)))
  sort(unique(as.character(x)))
}

#' Cal-Allocation 主函数
#'
#' @param data              数据框；第1列为样本ID，cohort_col 列标识重复样本
#' @param weight            特征权重（命名 list，如 list(disease=1, sex=0.5)）
#' @param set_size          每 set 目标容量（默认35）
#' @param set.max.size      显式 set 容量向量（默认按 set_size 自动生成，
#'                          前 leftover 个 set 容量 +1）
#' @param seed              随机种子
#' @param cohort_col        cohort 列名；NULL 则每样本独立
#' @param combo_features    参与 combo 构造的特征列名（默认除ID/cohort外全部）
#' @param objective_features 进入目标函数 d 的特征列名（默认同上）
#' @param protect_top       交换阶段跳过偏差最大的前 protect_top 个 set 作为 s1
#'                          （默认 0 = 从偏差最大的 set 开始；原版 code2.R 的
#'                          for(s1 in 10:...) 是交互会话续跑的痕迹，原意为 1）
#' @param legacy            TRUE=复刻 code2.R 行为（对拍用）
#' @param legacy_feature_num legacy 模式下覆盖判定使用的列数（复刻 code2.R 的
#'                           feature_num 参数；默认=目标特征数）
#' @param swap_tol          交换阶段收敛容差（每轮相对改进阈值，默认1e-3；
#'                          大 set 数下尾部轮次改进呈几何衰减，1e-6 过严）
#' @param max_swap_passes   交换阶段最大轮数（默认15；实测超15轮后
#'                          sd 改进 <0.1%，对结果无实质影响）
#' @param adaptive_weight   自适应加权：交换阶段每轮按特征跨 set 离散度
#'                          重算权重 w' = w × disp/mean(disp)（默认FALSE）
#' @param objective_metric  交换阶段优化的目标（默认 "sd" = 论文原版
#'                          sd(all_d)）："mean" = mean(all_d)，"meansq" =
#'                          mean(all_d^2)；后两者为 2026-11-16 pilot 新增，
#'                          仅作用于论文版路径（legacy 只支持 "sd"）
#' @param verbose           打印进度
#' @param subject_block     TRUE=受试者整体不可分（纵向/重复样本）；默认 FALSE。
#'                          需要 cohort_col 指定受试者列，且第一列样本ID唯一。
#'                          初始化时以受试者样本数为 combo，随后执行
#'                          “一个受试者 ↔ 若干受试者”的等样本量子集互换。
#' @param subject_block_seed 受试者子集互换阶段的随机种子；NULL 时用 seed。
#' @param subject_block_max_iter 子集互换最大迭代次数（默认 5e6）。
#' @param subject_block_tol  接受互换所需的最小目标改进（默认 1e-12）。
#' @param subject_block_verbose 是否打印子集互换进度（默认跟随 verbose）。
#' @return list(assignment=数据框(原列+key+new_set),
#'              set_n, seed, runtime, all_sets_diff_sd, mean_d)
cal_allocation <- function(data, weight = NULL, set_size = 35, set.max.size = NULL,
                           seed = 20260813, cohort_col = "cohortID",
                           combo_features = NULL, objective_features = NULL,
                           protect_top = 0, legacy = FALSE,
                           legacy_feature_num = NULL,
                           swap_tol = 1e-3, max_swap_passes = 15,
                           adaptive_weight = FALSE,
                           objective_metric = c("sd", "mean", "meansq"),
                           verbose = TRUE,
                           subject_block = FALSE,
                           subject_block_seed = NULL,
                           subject_block_max_iter = 5e6,
                           subject_block_tol = 1e-12,
                           subject_block_verbose = verbose) {
  objective_metric <- match.arg(objective_metric)
  if (legacy && objective_metric != "sd")
    stop("legacy 模式只支持 objective_metric = 'sd'（对拍口径）")
  set.seed(seed)
  t0 <- Sys.time()
  msg <- function(...) if (verbose) cat(...)
  obj_fun <- switch(objective_metric,
                    sd     = function(d) sd(d),
                    mean   = function(d) mean(d),
                    meansq = function(d) mean(d^2))

  # ---------- 0. 输入整理 ----------
  id_col <- colnames(data)[1]
  if (is.null(cohort_col) || !(cohort_col %in% colnames(data))) cohort_col <- NULL
  all_features <- setdiff(colnames(data), c(id_col, cohort_col))
  if (is.null(combo_features)) combo_features <- all_features
  if (is.null(objective_features)) objective_features <- all_features
  if (is.null(weight)) weight <- as.list(rep(1, length(objective_features)))

  # 目标函数：水平、理论比例、权重
  obj_levels <- lapply(objective_features, function(f) get_levels(data[[f]]))
  names(obj_levels) <- objective_features
  L <- vapply(obj_levels, length, integer(1))
  cnt_all <- lapply(objective_features, function(f)
    tabulate(match(data[[f]], obj_levels[[f]]), nbins = L[f]))
  names(cnt_all) <- objective_features
  p_ratio <- lapply(cnt_all, function(z) z / sum(z))
  w <- vapply(seq_along(objective_features), function(i) {
    if (legacy) {
      v <- weight[[i]]
      if (is.null(v)) 1 else as.numeric(v)
    } else {
      v <- weight[[objective_features[i]]]
      if (is.null(v)) 1 else as.numeric(v)
    }
  }, numeric(1))
  wL <- w / L
  w_orig <- w   # 自适应加权只影响交换阶段内部；最终评估指标用原始权重

  N <- nrow(data)
  if (is.null(set.max.size)) {
    nset0 <- N %/% set_size
    leftover <- N %% set_size
    if (leftover <= nset0) {
      # 余数可作 +1 分配：前 leftover 个 set 容量 set_size+1
      nset <- nset0
      set.max.size <- c(rep(set_size + 1, leftover), rep(set_size, nset - leftover))
    } else {
      # 余数超过 set 数：向上取整，容量为 base 或 base+1
      nset <- ceiling(N / set_size)
      base <- N %/% nset
      r <- N - nset * base
      set.max.size <- c(rep(base + 1, r), rep(base, nset - r))
    }
  } else {
    nset <- length(set.max.size)
  }
  if (sum(set.max.size) != N) stop("set.max.size 总和必须等于样本数 N")

  # ---------- 0b. subject-block 模式 ----------
  # 受试者整体不可分：同 cohort_col 的所有行必须进入同一 set。
  # 初始化：以受试者样本数（subject_size）作为唯一 combo 特征，使同一受试者
  # 的所有访视落在同一 combo / unit；随后执行“一个受试者 ↔ 若干受试者”的
  # 等样本量子集互换，保持每个 set 的样本数精确不变。
  if (subject_block) {
    if (legacy) stop("subject_block 不支持 legacy = TRUE")
    if (is.null(cohort_col))
      stop("subject_block = TRUE 需要 cohort_col 指定受试者列")
    if (anyDuplicated(data[[id_col]]))
      stop("subject_block = TRUE 要求第一列样本ID唯一")

    subj_chr <- as.character(data[[cohort_col]])
    ssize <- as.integer(table(subj_chr)[subj_chr])
    data_sb <- data
    data_sb[[".subject_size"]] <- ssize

    # 初始化沿用标准 Cal 填充/交换，只把 combo 特征换成 subject_size；
    # 这样每个受试者跨访视仍是一个不可分 unit。
    init <- cal_allocation(
      data = data_sb, weight = weight, set_size = set_size,
      set.max.size = set.max.size, seed = seed,
      cohort_col = cohort_col, combo_features = ".subject_size",
      objective_features = objective_features,
      protect_top = protect_top, legacy = FALSE,
      legacy_feature_num = legacy_feature_num,
      swap_tol = swap_tol, max_swap_passes = max_swap_passes,
      adaptive_weight = adaptive_weight,
      objective_metric = objective_metric, verbose = verbose,
      subject_block = FALSE
    )
    asg <- integer(N)
    asg[match(init$assignment[[id_col]], data[[id_col]])] <- init$assignment$new_set
    if (anyNA(asg)) stop("subject_block 初始化结果无法映射回输入行")
    if (any(tapply(asg, subj_chr, function(z) length(unique(z))) > 1))
      stop("subject_block 初始化仍拆分了受试者")

    # 受试者级计数素材
    subs <- sort(unique(subj_chr))
    uidx <- match(subj_chr, subs)
    usz <- as.integer(table(subj_chr)[subs])
    nsub <- length(subs)
    U <- lapply(objective_features, function(f) {
      m <- matrix(0, nrow = nsub, ncol = L[f])
      for (k in seq_len(nsub)) {
        ix <- which(uidx == k)
        m[k, ] <- tabulate(match(data[[f]][ix], obj_levels[[f]]), nbins = L[f])
      }
      m
    })
    names(U) <- objective_features
    d_of <- function(cnt_list, n) {
      sum(vapply(seq_along(cnt_list),
                 function(f) sum(abs(cnt_list[[f]] / n - p_ratio[[f]])) *
                   (1 / L[f]) * w[f], numeric(1)))
    }

    ss <- as.integer(tapply(asg, subj_chr, function(z) z[1])[subs])
    set_n_sb <- as.integer(tapply(usz, ss, sum)[as.character(seq_len(nset))])
    if (anyNA(set_n_sb) || any(set_n_sb != set.max.size))
      stop("subject_block 初始化容量不精确")
    members <- split(seq_along(subs), factor(ss, levels = seq_len(nset)))
    cnt <- lapply(seq_len(nset), function(s)
      lapply(seq_along(U), function(f) colSums(U[[f]][ss == s, , drop = FALSE])))
    all_d <- vapply(seq_len(nset), function(s) d_of(cnt[[s]], set_n_sb[s]),
                    numeric(1))
    start_mean_d <- mean(all_d)
    obj_sb <- obj_fun(all_d)

    set.seed(if (is.null(subject_block_seed)) seed else subject_block_seed)
    best <- obj_sb
    best_ss <- ss
    t_sb <- Sys.time()
    n_accept <- 0L

    for (it in seq_len(subject_block_max_iter)) {
      A <- sample.int(nset, 1); B <- sample.int(nset, 1)
      if (A == B) next
      if (length(members[[A]]) == 0 || length(members[[B]]) == 0) next
      i <- members[[A]][sample.int(length(members[[A]]), 1)]
      target <- usz[i]
      cand <- members[[B]]
      sel <- integer(0); rem <- target
      for (j in cand[sample(length(cand))]) {
        if (usz[j] <= rem) { sel <- c(sel, j); rem <- rem - usz[j] }
        if (rem == 0) break
      }
      if (rem != 0) next
      if (length(sel) == 1 && sel == i) next

      sumSel <- lapply(U, function(M) colSums(M[sel, , drop = FALSE]))
      clA <- vector("list", length(U)); clB <- vector("list", length(U))
      for (f in seq_along(U)) {
        clA[[f]] <- cnt[[A]][[f]] - U[[f]][i, ] + sumSel[[f]]
        clB[[f]] <- cnt[[B]][[f]] - sumSel[[f]] + U[[f]][i, ]
      }
      dA <- d_of(clA, set_n_sb[A]); dB <- d_of(clB, set_n_sb[B])

      if (objective_metric == "mean") {
        delta <- (dA + dB - all_d[A] - all_d[B]) / nset
      } else if (objective_metric == "meansq") {
        delta <- (dA^2 + dB^2 - all_d[A]^2 - all_d[B]^2) / nset
      } else {
        S <- sum(all_d); SS <- sum(all_d^2)
        S2 <- S - all_d[A] - all_d[B] + dA + dB
        SS2 <- SS - all_d[A]^2 - all_d[B]^2 + dA^2 + dB^2
        m1 <- S / nset; m2 <- S2 / nset
        delta <- sqrt((SS2 - nset * m2^2) / (nset - 1)) -
                 sqrt((SS - nset * m1^2) / (nset - 1))
      }

      if (delta < -subject_block_tol) {
        ss[i] <- B; ss[sel] <- A
        members[[A]] <- c(setdiff(members[[A]], i), sel)
        members[[B]] <- c(setdiff(members[[B]], sel), i)
        cnt[[A]] <- clA; cnt[[B]] <- clB
        all_d[A] <- dA; all_d[B] <- dB
        obj_sb <- obj_sb + delta
        n_accept <- n_accept + 1L
        if (obj_sb < best) { best <- obj_sb; best_ss <- ss }
      }
      if (subject_block_verbose && it %% 500000 == 0)
        cat(sprintf("subject-block iter %d obj %.8f best %.8f elapsed %.0fs\n",
                    it, obj_sb, best,
                    as.numeric(difftime(Sys.time(), t_sb, units = "secs"))))
    }
    t_sb <- as.numeric(difftime(Sys.time(), t_sb, units = "secs"))

    # 用 best_ss 重算最终计数（循环结束时 cnt 对应的是当前 ss，不一定等于 best）
    ss_final <- best_ss
    cnt_final <- lapply(seq_len(nset), function(s)
      lapply(seq_along(U), function(f)
        colSums(U[[f]][ss_final == s, , drop = FALSE])))
    final_d <- vapply(seq_len(nset), function(s)
      d_of(cnt_final[[s]], set_n_sb[s]), numeric(1))
    asg_final <- as.integer(ss_final[uidx])
    tab <- tabulate(asg_final, nbins = nset)
    if (!all(tab == set.max.size)) stop("subject_block 最终容量不精确")
    n_split <- sum(tapply(asg_final, subj_chr,
                          function(z) length(unique(z))) > 1)

    ord <- order(asg_final, seq_len(N))
    assignment <- data[ord, , drop = FALSE]
    assignment$key <- 0L
    assignment$new_set <- asg_final[ord]
    runtime <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    msg("subject-block 完成：初始 mean_d=", round(start_mean_d, 8),
        " 最终 mean_d=", round(mean(final_d), 8),
        " 受试者拆分=", n_split, " 耗时=", round(runtime, 1), " 秒\n")
    return(list(
      assignment = assignment,
      set_n = as.integer(tab),
      seed = seed,
      runtime = runtime,
      all_sets_diff_sd = sd(final_d),
      mean_d = mean(final_d),
      final_d = final_d,
      subject_block = TRUE,
      subject_block_seed = if (is.null(subject_block_seed)) seed else subject_block_seed,
      subject_block_iter = subject_block_max_iter,
      subject_block_accepted = n_accept,
      subject_block_runtime = t_sb,
      subject_block_start_mean_d = start_mean_d,
      n_split_subjects = n_split,
      n_multi_subjects = sum(table(subj_chr) > 1)
    ))
  }

  # ---------- 1. combo 生成 ----------
  flist <- lapply(combo_features, function(f) levels(as.factor(data[[f]])))
  names(flist) <- combo_features
  grid <- expand.grid(flist, stringsAsFactors = FALSE)
  grid_key <- do.call(paste, c(unname(grid), list(sep = "|")))
  row_key <- do.call(paste, c(lapply(combo_features, function(f) data[[f]]), list(sep = "|")))
  cnt_map <- table(factor(row_key, levels = unique(row_key)))
  cnt_tmp <- as.integer(cnt_map[grid_key])
  cnt_tmp[is.na(cnt_tmp)] <- 0L   # 空组合计数为0
  grid$sample.n <- cnt_tmp
  keep <- grid$sample.n > 0
  grid <- grid[keep, , drop = FALSE]
  grid_key <- grid_key[keep]      # grid 与 grid_key 必须同步过滤
  ord <- order(grid$sample.n, decreasing = TRUE)
  grid <- grid[ord, , drop = FALSE]
  grid_key <- grid_key[ord]
  rownames(grid) <- paste0("comb", seq_len(nrow(grid)))
  ncomb <- nrow(grid)
  msg("step1 combo 生成完成: ", ncomb, " 个非空组合\n")

  # unit 定义（两种模式一致）：cohort = 同一 cohortID 且所有特征水平
  # 完全相同的样本组（即 combo 内的 cohortID 分组）。同一供体不同访视的
  # 样本（特征不同）是不同 unit，允许分到不同 set；特征完全相同的重复
  # 样本（分管/技术重复）必须同 set。
  if (is.null(cohort_col)) {
    cids <- data[[id_col]]
  } else {
    cids <- data[[cohort_col]]
  }
  units_per_combo <- vector("list", ncomb)
  names(units_per_combo) <- rownames(grid)
  unit_rows <- list(); unit_n <- integer(); unit_combo <- integer()
  for (ci in seq_len(ncomb)) {
    ix <- which(row_key == grid_key[ci])
    # 组名排序复刻 table() 语义（保留原类型：数字按数值排序，字符按字符排序）
    gname <- sort(unique(cids[ix]))
    grp <- lapply(gname, function(g) ix[cids[ix] == g])
    names(grp) <- as.character(gname)
    ids <- paste0("L", ci, "_", as.character(gname))
    units_per_combo[[ci]] <- ids
    for (k in seq_along(grp)) {
      unit_rows[[ids[k]]] <- grp[[k]]
      unit_n[ids[k]] <- length(grp[[k]])
      unit_combo[ids[k]] <- ci
    }
  }
  # 每个 unit 在各目标特征上的水平计数（交换优化与 d 计算的素材）
  unit_cnt <- lapply(unit_rows, function(ix)
    lapply(objective_features, function(f)
      tabulate(match(data[[f]][ix], obj_levels[[f]]), nbins = L[f])))

  # 向量化素材：每个特征一个 unit×levels 计数矩阵（交换阶段批量索引用）
  unit_names <- names(unit_rows)
  um_idx <- setNames(seq_along(unit_names), unit_names)
  UM <- lapply(objective_features, function(f) {
    t(vapply(unit_rows, function(ix)
      tabulate(match(data[[f]][ix], obj_levels[[f]]), nbins = L[f]), numeric(L[f])))
  })
  names(UM) <- objective_features

  # unit 与 combo 一一对应（组内行特征全同），行数口径的 sample.n 与
  # unit 列表天然同步，无需重算

  # set 状态
  set_n <- integer(nset)
  set_rows <- replicate(nset, integer(0), simplify = FALSE)
  set_key <- replicate(nset, integer(0), simplify = FALSE)
  set_cnt <- replicate(nset, lapply(cnt_all, function(z) integer(length(z))), simplify = FALSE)
  set_units <- replicate(nset, character(0), simplify = FALSE)
  set_unit_key <- replicate(nset, integer(0), simplify = FALSE)

  # 注意：按原版 get.set.feature.ratio.diff 的运算顺序逐项计算
  # (sum(|diff|) × (1/L_f) × w_f)，保证浮点舍入与对拍一致
  d_counts <- function(cnt_list, n) {
    sum(vapply(seq_along(cnt_list),
               function(f) sum(abs(cnt_list[[f]] / n - p_ratio[[f]])) * (1 / L[f]) * w[f],
               numeric(1)))
  }
  set_d <- function(i) d_counts(set_cnt[[i]], set_n[i])

  add_unit_to_set <- function(uid, si) {
    ix <- unit_rows[[uid]]
    set_rows[[si]] <<- c(set_rows[[si]], ix)
    set_key[[si]] <<- c(set_key[[si]], rep(0L, length(ix)))
    set_n[si] <<- set_n[si] + length(ix)
    set_units[[si]] <<- c(set_units[[si]], uid)
    set_unit_key[[si]] <<- c(set_unit_key[[si]], 0L)
    for (f in seq_along(objective_features))
      set_cnt[[si]][[f]] <<- set_cnt[[si]][[f]] + unit_cnt[[uid]][[f]]
    ci <- unit_combo[uid]
    grid$sample.n[ci] <<- grid$sample.n[ci] - length(ix)
    units_per_combo[[ci]] <<- setdiff(units_per_combo[[ci]], uid)
  }

  # 从 combo ci 随机序尝试取一个 unit 放入 set si，成功返回 TRUE
  try_pick_one <- function(ci, si, maxsize) {
    us <- units_per_combo[[ci]]
    if (length(us) == 0) return(FALSE)
    ordr <- sample(seq_along(us), length(us))   # 随机顺序（RNG 调用点，legacy 对齐关键）
    for (k in ordr) {
      uid <- us[k]
      if (set_n[si] + unit_n[uid] > maxsize) next
      add_unit_to_set(uid, si)
      return(TRUE)
    }
    return(FALSE)
  }

  # ---------- 2. comb unit 构造（贪心全水平覆盖） ----------
  # legacy 复刻 code2.R：覆盖列 = 前 feature_num 列（code2.R 中 feature_num=4
  # 而 combo 仅 3 特征 → 误含 sample.n 列；Cal26-36 中 feature_num=6=特征数 → 正常）
  if (legacy) {
    if (is.null(legacy_feature_num)) legacy_feature_num <- length(objective_features)
    cover_cols <- seq_len(min(legacy_feature_num, ncol(grid)))
  } else {
    cover_cols <- seq_along(combo_features)
  }
  feature.level <- vapply(cover_cols, function(j) length(unique(grid[[j]])), integer(1))
  cover_check <- function(idx) {
    all(vapply(cover_cols, function(j) length(unique(grid[[j]][idx])), integer(1)) == feature.level)
  }
  contribution_of <- function(cand, unit_idx0) {
    sum(vapply(cover_cols, function(j)
      if (!cand[[j]] %in% grid[[j]][unit_idx0]) feature.level[j] else 0, numeric(1)))
  }
  unit_idx <- 1L
  avail <- setdiff(seq_len(ncomb), unit_idx)
  while (length(avail) > 0) {
    if (cover_check(unit_idx)) break
    contrib <- vapply(avail, function(i) contribution_of(grid[i, , drop = FALSE], unit_idx), numeric(1))
    pick <- avail[which.max(contrib)[1]]
    unit_idx <- c(unit_idx, pick)
    avail <- setdiff(avail, pick)
  }
  unit_org_len <- length(unit_idx)
  if (anyNA(grid$sample.n) || anyNA(unit_idx)) {
    msg("[诊断] NA 出现: sample.n NA 数=", sum(is.na(grid$sample.n)),
        " unit_idx NA 数=", sum(is.na(unit_idx)),
        " sample.n 总计=", sum(grid$sample.n, na.rm = TRUE), "\n")
  }
  msg("step2 comb unit 构造完成: ", unit_org_len, " 个 combo\n")

  # ---------- 3. 轮抽关键样本 ----------
  # 每个空 set 依次从 unit 每个 combo 抽一个 cohort；combo 抽光即从
  # 剩余组合按贡献分补位；循环至无空 set 或无剩余组合
  while (length(avail) > 0 && any(set_n == 0)) {
    if (length(unit_idx) < unit_org_len && length(avail) > 0) {
      contrib <- vapply(avail, function(i) contribution_of(grid[i, , drop = FALSE], unit_idx), numeric(1))
      pick <- avail[which.max(contrib)[1]]
      unit_idx <- c(unit_idx, pick)
      avail <- setdiff(avail, pick)
      next
    }
    while (length(unit_idx) == unit_org_len || length(avail) == 0) {
      empty_i <- which(set_n == 0)[1]
      if (!is.na(empty_i)) {
        for (ci in unit_idx) {
          if (grid$sample.n[ci] > 0) try_pick_one(ci, empty_i, set.max.size[empty_i])
        }
        gone <- unit_idx[grid$sample.n[unit_idx] == 0]
        if (length(gone) > 0) unit_idx <- setdiff(unit_idx, gone)
      } else {
        break
      }
    }
  }
  # 此阶段已入 set 的样本标记为关键样本 key=1，交换优化不动
  for (si in seq_len(nset)) {
    if (set_n[si] > 0) {
      set_key[[si]] <- rep(1L, set_n[si])
      set_unit_key[[si]] <- rep(1L, length(set_units[[si]]))
    }
  }
  msg("step3 轮抽完成: ", sum(set_n), " 个样本入 set（关键样本 key=1）\n")

  # 剩余组合并入 unit
  if (length(avail) > 0) {
    unit_idx <- c(unit_idx, avail)
    avail <- integer(0)
  }

  # ---------- 4. 填满 ----------
  # legacy 分支逐点复刻 code2.R（对拍用，含其固有的行为细节）
  # 论文版分支：cohort 按大小降序尝试（随机平局，减少尾部装箱失败）；
  # 无 cohort 可装入剩余容量时该 set 暂缓；最后对暂缓 set 做回溯精确装箱
  if (legacy) {
    for (si in seq_len(nset)) {
      while (set_n[si] < set.max.size[si]) {
        gone <- unit_idx[grid$sample.n[unit_idx] == 0]
        if (length(gone) > 0) unit_idx <- setdiff(unit_idx, gone)
        nc <- length(unit_idx)
        if (nc == 0) break
        first_u <- vapply(unit_idx, function(ci) {
          us <- units_per_combo[[ci]]
          if (length(us) == 0) NA_character_ else us[1]
        }, character(1))
        valid <- !is.na(first_u)
        contrib <- rep(Inf, nc)
        # 复刻 code2.R：trial 用该 combo 剩余行中行号最小的那一行
        for (k in seq_len(nc)) {
          if (!valid[k]) next
          us <- units_per_combo[[unit_idx[k]]]
          ix1 <- min(vapply(us, function(u) unit_rows[[u]][1], numeric(1)))
          extra <- lapply(objective_features, function(f)
            tabulate(match(data[[f]][ix1], obj_levels[[f]]), nbins = L[f]))
          contrib[k] <- d_counts(Map(`+`, set_cnt[[si]], extra), set_n[si] + 1)
        }
        o <- order(contrib)
        finish <- FALSE
        for (k in o) {
          ci <- unit_idx[k]
          if (grid$sample.n[ci] <= 0) next
          us <- units_per_combo[[ci]]
          if (length(us) == 0) next
          ordr <- sample(seq_along(us), length(us))   # RNG 调用点
          for (kk in ordr) {
            uid <- us[kk]
            if (is.na(unit_n[uid]) || set_n[si] + unit_n[uid] > set.max.size[si]) next
            add_unit_to_set(uid, si)
            finish <- TRUE
            break
          }
          if (finish) break
        }
      }
    }
  } else {
    # 回溯精确装箱（尾部小问题求解）：cohort 大小降序 + 容量降序首配
    # 注意：sizes 由调用方按降序传入（内部不再排序，保证 sol 与调用方
    # 的 cohort 顺序一一对齐）；caps 内部排序并映射回原位置
    pack_exact <- function(sizes, caps) {
      ord_c <- order(caps, decreasing = TRUE)
      caps_s <- caps[ord_c]
      n_it <- length(sizes); n_c <- length(caps_s)
      if (n_it == 0) return(integer(0))
      assign <- integer(n_it)
      rem <- caps_s
      nodes <- 0L
      solve <- function(i, start_b) {
        nodes <<- nodes + 1L
        if (nodes > 5000000L) return(FALSE)   # 防组合爆炸
        if (i > n_it) return(TRUE)
        # 相同大小的物品：后一个只允许放在前一个所选 bin 及之后（对称剪枝）
        b0 <- if (i > 1 && sizes[i] == sizes[i - 1]) start_b else 1L
        tried_cap <- integer(0)
        for (b in b0:n_c) {
          if (rem[b] < sizes[i]) next
          if (rem[b] %in% tried_cap) next    # 等容量 bin 剪枝
          tried_cap <- c(tried_cap, rem[b])
          rem[b] <<- rem[b] - sizes[i]
          assign[i] <<- b
          if (solve(i + 1, b)) return(TRUE)
          rem[b] <<- rem[b] + sizes[i]
        }
        FALSE
      }
      if (solve(1, 1L)) ord_c[assign] else NULL
    }
    # 最小超载回退：大件优先装入剩余容量最大的 bin；装不下时放入剩余
    # 容量最大的 bin（唯一允许超载的途径）
    fallback_pack <- function(sizes, caps) {
      ord_i <- order(sizes, decreasing = TRUE)
      assign <- integer(length(sizes))
      rem <- caps
      for (i in ord_i) {
        b <- which(rem >= sizes[i])
        bb <- if (length(b) > 0) b[which.max(rem[b])] else which.max(rem)
        assign[i] <- bb
        rem[bb] <- rem[bb] - sizes[i]
      }
      assign
    }
    deferred <- integer(0)
    for (si in seq_len(nset)) {
      while (set_n[si] < set.max.size[si]) {
        gone <- unit_idx[grid$sample.n[unit_idx] == 0]
        if (length(gone) > 0) unit_idx <- setdiff(unit_idx, gone)
        nc <- length(unit_idx)
        if (nc == 0) break
        first_u <- vapply(unit_idx, function(ci) {
          us <- units_per_combo[[ci]]
          if (length(us) == 0) NA_character_ else us[1]
        }, character(1))
        valid <- !is.na(first_u)
        contrib <- rep(Inf, nc)
        # 向量化 trial：用每个 combo 第一个 unit 的全部行
        idx_v <- which(valid)
        n_add <- unit_n[first_u[idx_v]]
        contrib_v <- numeric(length(idx_v))
        for (f in seq_along(objective_features)) {
          Mf <- UM[[f]][um_idx[first_u[idx_v]], , drop = FALSE]
          C <- sweep(Mf, 2, set_cnt[[si]][[f]], "+")
          contrib_v <- contrib_v + rowSums(abs(C / (set_n[si] + n_add) -
            rep(p_ratio[[f]], each = length(idx_v)))) * (1 / L[f]) * w[f]
        }
        contrib[idx_v] <- contrib_v
        o <- order(contrib)
        finish <- FALSE
        for (k in o) {
          ci <- unit_idx[k]
          if (grid$sample.n[ci] <= 0) next
          us <- units_per_combo[[ci]]
          if (length(us) == 0) next
          # cohort 大小降序尝试（随机平局）：大 cohort 优先离场，
          # 避免尾部"剩余容量装不下任何 cohort"的死锁
          ordr <- order(-unit_n[us], sample(seq_along(us)))
          for (kk in ordr) {
            uid <- us[kk]
            if (is.na(unit_n[uid]) || set_n[si] + unit_n[uid] > set.max.size[si]) next
            add_unit_to_set(uid, si)
            finish <- TRUE
            break
          }
          if (finish) break
        }
        if (!finish) {
          # 无 cohort 可装入剩余容量：该 set 暂缓，先填其他 set
          deferred <- c(deferred, si)
          break
        }
      }
    }
    # 尾部收口：剩余 cohort 精确装入暂缓 set 的剩余容量
    remaining_uids <- unique(unlist(units_per_combo, use.names = FALSE))
    if (length(remaining_uids) > 0 && length(deferred) > 0) {
      caps <- set.max.size[deferred] - set_n[deferred]
      # cohort 按大小降序排列后再送入 pack_exact（保证 sol 对齐）
      ord_u <- order(unit_n[remaining_uids], decreasing = TRUE)
      uids_sorted <- remaining_uids[ord_u]
      sol <- pack_exact(unit_n[uids_sorted], caps)
      if (is.null(sol)) {
        # 节点上限内未找到精确装箱：改用最小超载分配（记录警告）
        msg("[警告] 精确装箱未在节点上限内找到，改用最小超载分配\n")
        sol <- fallback_pack(unit_n[uids_sorted], caps)
      }
      for (j in seq_along(uids_sorted)) {
        add_unit_to_set(uids_sorted[j], deferred[sol[j]])
      }
    }
  }
  msg("step4 填满完成\n")

  # ---------- 5. 交换优化 ----------
  # 交换执行函数：s1 的第 i 个 unit 与 s2 的第 j 个 unit 互换
  swap_units <- function(s1, i, s2, j) {
    ui <- set_units[[s1]][i]; uj <- set_units[[s2]][j]
    ix1 <- unit_rows[[ui]]; ix2 <- unit_rows[[uj]]
    set_units[[s1]][i] <<- uj; set_units[[s2]][j] <<- ui
    set_unit_key[[s1]][i] <<- 0L; set_unit_key[[s2]][j] <<- 0L
    idx1 <- which(set_rows[[s1]] %in% ix1)
    idx2 <- which(set_rows[[s2]] %in% ix2)
    set_rows[[s1]][idx1] <<- ix2
    set_rows[[s2]][idx2] <<- ix1
    for (f in seq_along(objective_features)) {
      set_cnt[[s1]][[f]] <<- set_cnt[[s1]][[f]] - unit_cnt[[ui]][[f]] + unit_cnt[[uj]][[f]]
      set_cnt[[s2]][[f]] <<- set_cnt[[s2]][[f]] - unit_cnt[[uj]][[f]] + unit_cnt[[ui]][[f]]
    }
  }

  all_d <- vapply(seq_len(nset), set_d, numeric(1))
  pri_sd <- sd(all_d)

  if (legacy) {
    # 复刻 code2.R：仅 s1=偏差第(protect_top+1)名的 set，与其余所有 set
    # 逐一执行最优一次交换；s2 从末位向前
    set_ord <- order(all_d, decreasing = TRUE)
    s1_range <- if (protect_top + 1 <= nset) protect_top + 1 else integer(0)
    for (s1 in s1_range) {
      s1_set <- set_ord[s1]   # 原版: set1 <- set.order[s1]（位置索引 → set 编号）
      for (s2 in rev(set_ord)) {
        if (s1_set == s2) next
        all_d <- vapply(seq_len(nset), set_d, numeric(1))
        pri_sd <- sd(all_d)
        rest_d <- all_d[-c(s1_set, s2)]
        u1 <- set_units[[s1_set]]; u2 <- set_units[[s2]]
        k1 <- set_unit_key[[s1_set]]; k2 <- set_unit_key[[s2]]
        contribution <- c(); contrib_names <- c()
        for (i in seq_along(u1)) {
          for (j in seq_along(u2)) {
            if (unit_n[u1[i]] != unit_n[u2[j]]) next
            if (k1[i] == 1 || k2[j] == 1) next
            cnt1 <- Map(function(a, b) a - unit_cnt[[u1[i]]][[b]] + unit_cnt[[u2[j]]][[b]],
                        set_cnt[[s1_set]], seq_along(objective_features))
            cnt2 <- Map(function(a, b) a - unit_cnt[[u2[j]]][[b]] + unit_cnt[[u1[i]]][[b]],
                        set_cnt[[s2]], seq_along(objective_features))
            tmp_sd <- sd(c(rest_d, d_counts(cnt1, set_n[s1_set]), d_counts(cnt2, set_n[s2])))
            if (tmp_sd < pri_sd) {
              contribution <- c(contribution, tmp_sd)
              contrib_names <- c(contrib_names, paste0(i, "-", j))
            }
          }
        }
        if (length(contribution) == 0) next
        best <- strsplit(contrib_names[which.min(contribution)][1], "-")[[1]]
        i <- as.integer(best[1]); j <- as.integer(best[2])
        # 执行交换：unit 整体互换（行号、计数、unit 列表同步更新）
        swap_units(s1_set, i, s2, j)
      }
    }
  } else {
    # 论文版：多轮重复，每对 (s1,s2) 执行最优一次交换，直至整轮无交换
    s1_idx <- intersect((protect_top + 1):nset, seq_len(nset))
    pass <- 0L
    repeat {
      pass <- pass + 1L
      if (adaptive_weight && pass > 1) {
        # 自适应加权：每轮按特征跨 set 离散度重算权重（boosting 式）
        disp <- vapply(seq_along(objective_features), function(f) {
          d_f <- vapply(seq_len(nset), function(s) {
            if (set_n[s] == 0) return(0)
            sum(abs(set_cnt[[s]][[f]] / set_n[s] - p_ratio[[f]])) * (1 / L[f])
          }, numeric(1))
          sd(d_f)
        }, numeric(1))
        if (all(is.na(disp)) || all(disp <= 0)) {
          disp <- rep(1, length(disp))
        } else {
          disp[is.na(disp) | disp <= 0] <- min(disp[disp > 0], na.rm = TRUE)
        }
        w <- w_orig * disp / mean(disp)
        wL <- w / L
        all_d <- vapply(seq_len(nset), set_d, numeric(1))
      }
      obj_before_pass <- obj_fun(all_d)
      msg("step5 交换轮次 ", pass, " 当前", objective_metric, "=",
          round(obj_before_pass, 8), "\n")
      set_ord <- order(all_d, decreasing = TRUE)
      n_swaps <- 0L
      for (s1 in set_ord[s1_idx]) {
        for (s2 in rev(set_ord)) {
          if (s1 == s2) next
          pri_obj <- obj_fun(all_d)   # 与 code2.R 一致：每对开始时按当前状态重算
          u1 <- set_units[[s1]]; u2 <- set_units[[s2]]
          k1 <- set_unit_key[[s1]]; k2 <- set_unit_key[[s2]]
          # 候选对：等行数且双方 key=0，按 (i,j) 嵌套序
          pi <- rep(seq_along(u1), each = length(u2))
          pj <- rep(seq_along(u2), times = length(u1))
          keep <- (unit_n[u1[pi]] == unit_n[u2[pj]]) & (k1[pi] == 0) & (k2[pj] == 0)
          pi <- pi[keep]; pj <- pj[keep]
          if (length(pi) == 0) next
          rest_d <- all_d[-c(s1, s2)]
          S_rest <- sum(rest_d); SS_rest <- sum(rest_d^2)
          d1p <- numeric(length(pi)); d2p <- numeric(length(pi))
          for (f in seq_along(objective_features)) {
            U1f <- UM[[f]][um_idx[u1[pi]], , drop = FALSE]
            U2f <- UM[[f]][um_idx[u2[pj]], , drop = FALSE]
            # 交换后 s1 计数 = set_cnt[s1] - U1f + U2f（注意 sweep 的减法方向）
            A <- sweep(U1f, 2, set_cnt[[s1]][[f]], FUN = function(x, y) y - x) + U2f
            d1p <- d1p + rowSums(abs(A / set_n[s1] - rep(p_ratio[[f]], each = length(pi)))) * (1 / L[f]) * w[f]
            B <- sweep(U2f, 2, set_cnt[[s2]][[f]], FUN = function(x, y) y - x) + U1f
            d2p <- d2p + rowSums(abs(B / set_n[s2] - rep(p_ratio[[f]], each = length(pi)))) * (1 / L[f]) * w[f]
          }
          tmp_obj <- if (objective_metric == "sd") {
            m <- (S_rest + d1p + d2p) / nset
            sqrt((SS_rest + d1p^2 + d2p^2 - nset * m^2) / (nset - 1))
          } else if (objective_metric == "mean") {
            (S_rest + d1p + d2p) / nset
          } else {                      # meansq
            (SS_rest + d1p^2 + d2p^2) / nset
          }
          better <- which(tmp_obj < pri_obj)
          if (length(better) == 0) next
          b <- better[which.min(tmp_obj[better])[1]]
          swap_units(s1, pi[b], s2, pj[b])
          n_swaps <- n_swaps + 1L
          all_d[s1] <- d1p[b]; all_d[s2] <- d2p[b]
        }
      }
      if (n_swaps == 0) break
      if (pass >= max_swap_passes) break
      if ((obj_before_pass - obj_fun(all_d)) < swap_tol * abs(obj_before_pass)) break
    }
  }
  msg("step5 交换优化完成\n")

  # ---------- 输出 ----------
  # 评估指标一律用原始权重（自适应权重只作用于交换阶段内部）
  w <- w_orig
  wL <- w / L
  runtime <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  final_d <- vapply(seq_len(nset), set_d, numeric(1))
  assignment <- data[unlist(set_rows), , drop = FALSE]
  assignment$key <- unlist(set_key)
  assignment$new_set <- rep(seq_len(nset), times = vapply(set_rows, length, integer(1)))
  msg("完成，耗时 ", round(runtime, 1), " 秒\n")
  list(assignment = assignment,
       set_n = set_n,
       seed = seed,
       runtime = runtime,
       all_sets_diff_sd = sd(final_d),
       mean_d = mean(final_d),
       final_d = final_d)
}
