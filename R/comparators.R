# ============================================================
# 04_comparators.R  对比方法函数库（模拟研究用，source 使用）
# 方法：
#   1. complete_random    完全随机化（全局随机置换后顺序切片）
#   2. stratified_random  分层随机化（2主特征组合为层，层内洗牌轮抽）
#   3. block_random       区组随机化（全特征组合为区组，区内洗牌轮抽）
#   4. minimization_run   最小化（Pocock–Simon 序贯分配，平局随机）
#   5. omixer_run         Omixer 包（omixerRand，wells=set容量，plateNum=set数）
#   6. osat_run           OSAT 包（96孔芯片布局，仅 set 容量匹配时可用，
#                         否则返回 NULL——OSAT 固定芯片布局限制，见 Methods）
# 所有方法返回：data + new_set 列（按输入行序），失败返回 NULL
# 注意：方法 1-4 不感知 cohort（会分裂重复样本），这正是
#       cohort 分裂率指标的对比意义所在
# ============================================================
# 注意：本文件为库文件，由调用脚本负责工作目录

get_levels <- function(x) {
  if (is.factor(x)) return(levels(x))
  if (is.numeric(x)) return(sort(unique(x)))
  sort(unique(as.character(x)))
}

# 轮抽：把洗牌后的行顺序轮流放入容量未满的 set；
# 指针跨层连续（start_s 续传），避免每个层都从 set1 开始造成容量不均
deal_round_robin <- function(n_ordered, cap, start_s = 1L) {
  assign <- integer(n_ordered)
  s <- start_s
  for (i in seq_len(n_ordered)) {
    while (cap[s] <= 0) s <- if (s == length(cap)) 1L else s + 1L
    assign[i] <- s
    cap[s] <- cap[s] - 1L
    s <- if (s == length(cap)) 1L else s + 1L
  }
  list(assign = assign, cap = cap, next_s = s)
}

#' 完全随机化：全局随机置换后顺序切成 set
complete_random <- function(data, set.max.size, seed) {
  set.seed(seed)
  n <- nrow(data)
  ord <- sample(n)
  sizes_ord <- rep(seq_along(set.max.size), times = set.max.size)
  assign <- integer(n)
  assign[ord] <- sizes_ord
  data$new_set <- assign
  data
}

#' 分层随机化：以 strata_cols（默认前2特征）的水平组合为层，
#' 层内洗牌后轮抽入各 set；轮抽指针跨层连续
stratified_random <- function(data, set.max.size, seed, strata_cols) {
  set.seed(seed)
  strata <- do.call(paste, c(lapply(strata_cols, function(f) data[[f]]), list(sep = "|")))
  lev <- unique(strata)
  assign <- integer(nrow(data))
  cap <- set.max.size
  s <- 1L
  for (lv in lev) {
    ix <- which(strata == lv)
    ord <- ix[sample(length(ix))]
    res <- deal_round_robin(length(ord), cap, start_s = s)
    assign[ord] <- res$assign
    cap <- res$cap
    s <- res$next_s
  }
  data$new_set <- assign
  data
}

#' 区组随机化：以全部特征的组合为区组，区内洗牌后轮抽入各 set；
#' 轮抽指针跨区组连续
block_random <- function(data, set.max.size, seed, features) {
  set.seed(seed)
  blocks <- do.call(paste, c(lapply(features, function(f) data[[f]]), list(sep = "|")))
  lev <- unique(blocks)
  assign <- integer(nrow(data))
  cap <- set.max.size
  s <- 1L
  for (lv in lev) {
    ix <- which(blocks == lv)
    ord <- ix[sample(length(ix))]
    res <- deal_round_robin(length(ord), cap, start_s = s)
    assign[ord] <- res$assign
    cap <- res$cap
    s <- res$next_s
  }
  data$new_set <- assign
  data
}

#' 最小化（Pocock–Simon）：样本随机序逐个进入使当前 d(set) 最小的 set，
#' 平局随机；d 与 Cal-Allocation 同口径（加权归一比例偏差）
minimization_run <- function(data, set.max.size, seed, features, weight = NULL) {
  set.seed(seed)
  n <- nrow(data); nset <- length(set.max.size)
  lev <- lapply(features, function(f) get_levels(data[[f]]))
  names(lev) <- features
  L <- vapply(lev, length, integer(1))
  p <- lapply(features, function(f)
    tabulate(match(data[[f]], lev[[f]]), nbins = L[f]) / n)
  wL <- vapply(seq_along(features), function(k) {
    wf <- if (is.null(weight) || is.null(weight[[features[k]]])) 1 else as.numeric(weight[[features[k]]])
    wf / L[k]
  }, numeric(1))
  onehot <- lapply(seq_len(n), function(i)
    lapply(features, function(f) tabulate(match(data[[f]][i], lev[[f]]), nbins = L[f])))
  cnt <- replicate(nset, lapply(features, function(f) integer(L[f])), simplify = FALSE)
  cap <- set.max.size; size <- integer(nset)
  assign <- integer(n)
  ord <- sample(n)
  for (i in ord) {
    avail <- which(cap > 0)
    d <- vapply(avail, function(s) {
      sum(vapply(seq_along(features), function(f)
        wL[f] * sum(abs((cnt[[s]][[f]] + onehot[[i]][[f]]) / (size[s] + 1) - p[[f]])),
        numeric(1)))
    }, numeric(1))
    best <- avail[d == min(d)]
    s <- if (length(best) == 1) best else best[sample.int(length(best), 1)]
    assign[i] <- s; cap[s] <- cap[s] - 1; size[s] <- size[s] + 1
    for (f in seq_along(features)) cnt[[s]][[f]] <- cnt[[s]][[f]] + onehot[[i]][[f]]
  }
  data$new_set <- assign
  data
}

#' Omixer：每个 set 视为一块 35 孔板（wells=set容量，plateNum=set数），
#' omixerRand 跑 iterNum 次随机布局取最优（最小化协变量-位置相关）
omixer_run <- function(data, set.max.size, seed, features, iterNum = 1000) {
  if (!requireNamespace("Omixer", quietly = TRUE)) return(NULL)
  tryCatch({
    nset <- length(set.max.size)
    wells <- unique(set.max.size)
    if (length(wells) != 1) return(NULL)   # 容量不一时不适用
    df <- data.frame(sampleId = data[[1]], block = data[[1]], data[features])
    out <- Omixer::omixerRand(df, sampleId = "sampleId", block = "block",
                              iterNum = iterNum, wells = wells,
                              div = "col", plateNum = nset,
                              randVars = features)
    data$new_set <- as.integer(out$plate[match(data[[1]], out$sampleId)])
    data
  }, error = function(e) NULL)
}

#' OSAT：固定 96 孔芯片布局（IlluminaBeadChip96Plate），仅当 set 容量
#' 与芯片布局匹配时可用；35 样本 set 不适用，返回 NULL（Methods 中说明）
osat_run <- function(data, set.max.size, seed, features, nSim = 1000) {
  if (!requireNamespace("OSAT", quietly = TRUE)) return(NULL)
  tryCatch({
    if (length(unique(set.max.size)) != 1 || unique(set.max.size) != 96) return(NULL)
    nset <- length(set.max.size)
    pheno <- data.frame(SampleID = data[[1]], data[features], stringsAsFactors = FALSE)
    for (f in features) pheno[[f]] <- factor(pheno[[f]])
    gs <- OSAT::setup.sample(pheno, optimal = features)
    gc <- OSAT::setup.container(OSAT::IlluminaBeadChip96Plate, nset, batch = "plates")
    set.seed(seed)
    gSetup <- OSAT::create.optimized.setup(sample = gs, container = gc, nSim = nSim)
    out <- OSAT::get.experiment.setup(gSetup)
    data$new_set <- as.integer(out$plates[match(data[[1]], out$SampleID)])
    data
  }, error = function(e) NULL)
}
