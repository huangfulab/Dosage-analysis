#!/usr/bin/env Rscript

# Guide-averaged bulk RNA-seq gene expression vs post-treatment OCT4 qPCR.
#
# Usage:
#   Rscript 05_gene_vs_qpcr_scatter.R [output_dir]
#       Full DEG batch: union of pairwise OCT4-gRNA vs pooled-control
#       strict DEGs (up/down), one row per gene symbol.
#   Rscript 05_gene_vs_qpcr_scatter.R [output_dir] --batch-deg [max_genes_per_direction]
#       Same union DEG batch, optional per-direction cap (preview).
#   Rscript 05_gene_vs_qpcr_scatter.R [output_dir] --batch-deg-pooled [max_genes_per_direction]
#       Previous pooled 8-vs-3 strict heatmap DEG list.
#   Rscript 05_gene_vs_qpcr_scatter.R [output_dir] --batch-all [max_genes]
#       All genes in the fitted DESeq2 object (alphabetical), optional cap.
#   Rscript 05_gene_vs_qpcr_scatter.R [output_dir] [gene_name_or_id]
#       Single-gene PNG/CSV diagnostics.
#   Add --csv-only anywhere after output_dir to skip PDF pages.
#
# AIC tables include the linear-model slope p-value (Figure 4A Linear vs NA)
# and quadratic coefficients a, b, c of y = a x^2 + b x + c.
# Models are not auto-called as Linear/NA; those columns are for manual calling.
#
# LL.4 matches Dosage-analysis shark: c/d fixed to observed y min/max;
# free parameters are b and e only.

options(stringsAsFactors=F)

script_dir = function() {
  file_arg = grep("^--file=", commandArgs(trailingOnly=F), value=T)
  if (length(file_arg) == 0L) return(normalizePath(getwd(), mustWork=T))
  script_file = sub("^--file=", "", file_arg[[1L]])
  return(dirname(normalizePath(script_file, mustWork=T)))
}
stop_with = function(...) stop(sprintf(...), call.=F)

args = commandArgs(trailingOnly=T)
project_dir = script_dir()
csv_only = any(args == "--csv-only")
args = args[args != "--csv-only"]
output_dir = if (length(args) >= 1L) args[[1L]] else file.path(project_dir, "outputs")

# Default with only output_dir: full pairwise-union DEG batch.
# Optional 2nd arg:
#   --batch-deg [N]         -> union DEG batch, optional per-direction cap N
#   --batch-deg-pooled [N]  -> pooled strict DEG batch, optional cap N
#   --batch-all [N]         -> all genes in DESeq2 object, optional cap N
#   <gene>                  -> single-gene mode
batch_mode = F
batch_scope = NA_character_   # "deg", "deg_pooled", or "all"
gene_query = NA_character_
batch_max_genes = Inf
if (length(args) < 2L) {
  batch_mode = T
  batch_scope = "deg"
  batch_max_genes = Inf
} else if (identical(args[[2L]], "--batch-deg")) {
  batch_mode = T
  batch_scope = "deg"
  if (length(args) >= 3L) {
    batch_max_genes = suppressWarnings(as.integer(args[[3L]]))
    if (!is.finite(batch_max_genes) || batch_max_genes < 1L) {
      stop_with("--batch-deg max_genes_per_direction must be a positive integer")
    }
  }
} else if (identical(args[[2L]], "--batch-deg-pooled")) {
  batch_mode = T
  batch_scope = "deg_pooled"
  if (length(args) >= 3L) {
    batch_max_genes = suppressWarnings(as.integer(args[[3L]]))
    if (!is.finite(batch_max_genes) || batch_max_genes < 1L) {
      stop_with("--batch-deg-pooled max_genes_per_direction must be a positive integer")
    }
  }
} else if (identical(args[[2L]], "--batch-all")) {
  batch_mode = T
  batch_scope = "all"
  if (length(args) >= 3L) {
    batch_max_genes = suppressWarnings(as.integer(args[[3L]]))
    if (!is.finite(batch_max_genes) || batch_max_genes < 1L) {
      stop_with("--batch-all max_genes must be a positive integer")
    }
  }
} else {
  batch_mode = F
  gene_query = args[[2L]]
}

if (!requireNamespace("DESeq2", quietly=T)) stop_with("DESeq2 is required but is not installed")
if (!requireNamespace("ggplot2", quietly=T)) stop_with("ggplot2 is required but is not installed")
if (!requireNamespace("ggrepel", quietly=T)) stop_with("ggrepel is required but is not installed")
if (!requireNamespace("drc", quietly=T)) stop_with("drc is required for LL.4 curve fitting but is not installed")
if (!requireNamespace("matrixStats", quietly=T)) stop_with("matrixStats is required but is not installed")
if (!requireNamespace("devtools", quietly=T)) stop_with("devtools is required for session information")
suppressPackageStartupMessages(library(DESeq2))

alpha = suppressWarnings(as.numeric(Sys.getenv("K2LLPS_ALPHA", unset="0.05")))
if (!is.finite(alpha) || alpha <= 0 || alpha >= 1) {
  stop_with("K2LLPS_ALPHA must be a number strictly between 0 and 1")
}
deg_lfc_threshold = suppressWarnings(as.numeric(Sys.getenv("K2LLPS_DEG_LFC", unset="1")))
if (!is.finite(deg_lfc_threshold) || deg_lfc_threshold < 0) {
  stop_with("K2LLPS_DEG_LFC must be a non-negative number")
}

output_dir = normalizePath(output_dir, mustWork=T)
dds_path = file.path(output_dir, "DESeq2Obj-after-wald-test.rds")
meta_path = file.path(output_dir, "analysis-sample-metadata.csv")
deg_path = file.path(output_dir, "degs-oct4-targeting-vs-pooled-controls-lfc1-padj0.05.csv")
guide_output_dir = file.path(output_dir, "guide-specific")
if (!file.exists(dds_path)) stop_with("Missing %s; run 01_prepare_deseq2.R first", dds_path)
if (!file.exists(meta_path)) stop_with("Missing %s; run 01_prepare_deseq2.R first", meta_path)
if (batch_mode && identical(batch_scope, "deg_pooled") && !file.exists(deg_path)) {
  stop_with("Missing %s; run 01_prepare_deseq2.R first", deg_path)
}
if (batch_mode && identical(batch_scope, "deg") && !dir.exists(guide_output_dir)) {
  stop_with("Missing %s; run 03_run_guide_specific_deseq2.R first", guide_output_dir)
}

guide_levels = c(
  "guide01", "guide02", "guide03", "guide04", "guide05", "guide06",
  "guide14", "guide15", "guide18", "guide19", "guide20"
)
control_guides = c("guide18", "guide19", "guide20")
model_levels = c("linear", "quadratic", "ll4")
model_colors = c(linear="#0072B2", quadratic="#009E73", ll4="#D55E00")
plot_surface = "#fcfcfb"
plot_text = "#0b0b0b"
plot_muted = "#898781"
plot_grid = "#e1e0d9"
plot_axis = "#c3c2b7"
class_colors = c(
  OCT4_targeting="#1baf7a",
  nontargeting_control="#2a78d6",
  safe_harbor_control="#6b5b95"
)

load_analysis_inputs = function() {
  DESeq2Obj = readRDS(dds_path)
  metadata = read.csv(meta_path, check.names=F, na.strings="NA")
  if (NROW(metadata) != NCOL(DESeq2Obj)) {
    stop_with("Metadata rows (%d) do not match DESeq2 samples (%d)", NROW(metadata), NCOL(DESeq2Obj))
  }
  if (!all(metadata$sample_id %in% colnames(DESeq2Obj))) {
    stop_with("Metadata sample_id values are not all present in the DESeq2 object")
  }
  metadata = metadata[match(colnames(DESeq2Obj), metadata$sample_id), , drop=F]
  rownames(metadata) = metadata$sample_id
  if (!setequal(unique(as.character(metadata$guide)), guide_levels)) {
    stop_with("Unexpected guide set in analysis metadata")
  }

  qpcr_by_guide = sapply(guide_levels, function(g) {
    mean(metadata$qpcr_oct4_post_treatment[metadata$guide == g], na.rm=T)
  })
  if (any(!is.finite(qpcr_by_guide))) stop_with("Non-finite guide-mean qPCR values")
  x_scale_ref = mean(qpcr_by_guide[control_guides])
  if (!(is.finite(x_scale_ref) && x_scale_ref > 0)) stop_with("Invalid qPCR control scaling reference")
  x_scaled = qpcr_by_guide / x_scale_ref * 100
  if (any(!(x_scaled > 0))) stop_with("All scaled qPCR values must be > 0 for LL.4")

  guide_meta = do.call(rbind, lapply(guide_levels, function(g) {
    rows = metadata[metadata$guide == g, , drop=F]
    data.frame(
      guide=g,
      guide_class=as.character(rows$guide_class[[1L]]),
      comparison_group=as.character(rows$comparison_group[[1L]]),
      jy_guide_label=as.character(rows$jy_guide_label[[1L]]),
      n_rna_replicates=NROW(rows),
      n_qpcr_replicates=sum(!is.na(rows$qpcr_oct4_post_treatment)),
      mean_qpcr_oct4=unname(qpcr_by_guide[[g]]),
      x_scaled=unname(x_scaled[[g]]),
      stringsAsFactors=F, check.names=F
    )
  }))
  guide_meta$guide_class = factor(
    guide_meta$guide_class,
    levels=c("OCT4_targeting", "nontargeting_control", "safe_harbor_control")
  )

  list(
    DESeq2Obj=DESeq2Obj,
    metadata=metadata,
    guide_meta=guide_meta,
    x_scale_ref=x_scale_ref
  )
}

fit_ll4 = function(fit_data, response_min, response_max) {
  fit_ll4_with_captured_io = function() {
    sink_file = NULL
    sink_connection = NULL
    output_sink_depth = sink.number(type="output")
    message_sink_connection = getConnection(sink.number(type="message"))
    graphics_devices_before = grDevices::dev.list()
    on.exit({
      graphics_devices_after = grDevices::dev.list()
      opened = setdiff(as.integer(graphics_devices_after), as.integer(graphics_devices_before))
      for (device in opened) try(grDevices::dev.off(which=device), silent=T)
      try(sink(message_sink_connection, type="message"), silent=T)
      try({ while (sink.number(type="output") > output_sink_depth) sink(type="output") }, silent=T)
      if (inherits(sink_connection, "connection")) try(close(sink_connection), silent=T)
      if (!is.null(sink_file)) try(unlink(sink_file), silent=T)
    }, add=T)
    sink_file = tempfile(fileext=".ll4-output")
    sink_connection = file(sink_file, open="wt")
    sink(sink_connection, type="output")
    sink(sink_connection, type="message")
    # Match Dosage-analysis shark LL.4: fix c/d to observed y range; free b,e only.
    drc::drm(
      y ~ x,
      data=fit_data,
      fct=drc::LL.4(fixed=c(NA, response_min, response_max, NA)),
      control=drc::drmc(noMessage=T)
    )
  }
  tryCatch(fit_ll4_with_captured_io(), error=function(e) e)
}

direction_text = function(direction_label) {
  if (is.null(direction_label) || is.na(direction_label) || !nzchar(direction_label)) return("single-gene")
  if (identical(direction_label, "all")) return("all genes")
  if (identical(direction_label, "not_deg")) return("not DEG")
  paste0(toupper(direction_label), " DEG")
}

make_zero_input_fallback_plot = function(gene_id, display_name, direction_label, which_axis, detail_note=NULL) {
  which_axis = as.character(which_axis)[[1L]]
  msg = sprintf("The input data of %s axis was all zero.", which_axis)
  detail = if (is.null(detail_note) || is.na(detail_note) || !nzchar(detail_note)) {
    ""
  } else {
    as.character(detail_note)[[1L]]
  }
  dir_txt = direction_text(direction_label)
  label_df = data.frame(
    x=c(0.5, 0.5),
    y=c(0.60, 0.42),
    lab=c(msg, detail),
    face=c("bold", "plain"),
    size=c(5.2, 3.6),
    color=c(plot_text, plot_muted),
    stringsAsFactors=F
  )
  label_df = label_df[nzchar(label_df$lab), , drop=F]
  ggplot2::ggplot(label_df, ggplot2::aes(x=x, y=y, label=lab)) +
    ggplot2::geom_text(
      ggplot2::aes(fontface=face, size=size, color=color),
      show.legend=F
    ) +
    ggplot2::scale_size_identity() +
    ggplot2::scale_color_identity() +
    ggplot2::coord_cartesian(xlim=c(0, 1), ylim=c(0, 1), expand=F) +
    ggplot2::labs(
      title=sprintf("%s (%s)", display_name, gene_id),
      subtitle=sprintf("%s; relative scatter skipped (all-zero input)", dir_txt),
      x=NULL, y=NULL
    ) +
    ggplot2::theme_void(base_family="sans") +
    ggplot2::theme(
      plot.background=ggplot2::element_rect(fill=plot_surface, color=plot_axis, linewidth=0.6),
      plot.title=ggplot2::element_text(color=plot_text, size=13, hjust=0.5, margin=ggplot2::margin(t=14, b=6)),
      plot.subtitle=ggplot2::element_text(color=plot_muted, size=9, hjust=0.5, margin=ggplot2::margin(b=10)),
      plot.margin=ggplot2::margin(16, 16, 16, 16)
    )
}

# Shared schema also supplies zero-row batch tables when one direction is empty.
make_failed_aic_row = function(gene_id, display_name, direction_label, x_scale_ref,
    note, y_scale_ref=NA_real_, y_scale_mode=NA_character_, ll4_c=NA_real_, ll4_d=NA_real_) {
  data.frame(
    gene_id=gene_id, gene_name=display_name, direction=direction_label,
    pvalue_linear=NA_real_,
    linear_intercept=NA_real_, linear_slope=NA_real_,
    hill_n=NA_real_, ll4_ed50=NA_real_,
    quadratic_a=NA_real_, quadratic_b=NA_real_, quadratic_c=NA_real_,
    quadratic_vertex_x=NA_real_, quadratic_vertex_in_range=NA,
    x_min=NA_real_, x_max=NA_real_,
    aic_linear=NA_real_, aic_quadratic=NA_real_, aic_ll4=NA_real_,
    delta_aic_linear=NA_real_, delta_aic_quadratic=NA_real_, delta_aic_ll4=NA_real_,
    weight_linear=NA_real_, weight_quadratic=NA_real_, weight_ll4=NA_real_,
    best_model=NA_character_, ll4_converged=F,
    rmse_linear=NA_real_, rmse_quadratic=NA_real_, rmse_ll4=NA_real_,
    pearson_r=NA_real_, spearman_rho=NA_real_,
    x_scale_ref=x_scale_ref, y_scale_ref=y_scale_ref, y_scale_mode=y_scale_mode,
    ll4_c_fixed=ll4_c, ll4_d_fixed=ll4_d, note=note,
    stringsAsFactors=F, check.names=F
  )
}

# Core per-gene analysis on already guide-averaged log2(norm+1) values.
analyze_gene_points = function(gene_id, display_name, direction_label, y_log2_by_guide, guide_meta, x_scale_ref) {
  y_raw = as.numeric(y_log2_by_guide[guide_levels])
  names(y_raw) = guide_levels
  fail_row = function(...) {
    make_failed_aic_row(gene_id, display_name, direction_label, x_scale_ref, ...)
  }
  fail_zero_axis = function(which_axis, note, y_scale_ref=NA_real_, detail=NULL) {
    list(
      ok=F,
      aic_row=fail_row(note, y_scale_ref=y_scale_ref, y_scale_mode="unavailable"),
      guide_df=NULL,
      pred_df=NULL,
      plot=make_zero_input_fallback_plot(
        gene_id=gene_id,
        display_name=display_name,
        direction_label=direction_label,
        which_axis=which_axis,
        detail_note=detail
      )
    )
  }

  # x-axis all-zero guard (shared qPCR scale reference).
  if (!(is.finite(x_scale_ref) && x_scale_ref > 0) || all(guide_meta$x_scaled == 0)) {
    return(fail_zero_axis(
      "x",
      "all_zero_x_axis",
      detail="Control-mean qPCR scaling reference was zero."
    ))
  }

  if (any(!is.finite(y_raw))) {
    return(list(ok=F, aic_row=fail_row("non_finite_guide_means"), guide_df=NULL, pred_df=NULL, plot=NULL))
  }
  if (all(y_raw == 0)) {
    # Truly no signal in any guide: keep the all-zero fallback page.
    return(fail_zero_axis(
      "y",
      "all_zero_y_axis",
      y_scale_ref=0,
      detail="All 11 guide-mean log2(normalized count + 1) values were zero."
    ))
  }

  # Default: relative y so mean(control guides)=100.
  # Edge case: controls are all zero counts => log2(0+1)=0, so relative-to-ctrl
  # is undefined. Keep the gene by plotting absolute log2(n+1) instead of dropping it.
  y_scale_ref = mean(y_raw[control_guides])
  if (is.finite(y_scale_ref) && y_scale_ref > 0) {
    y_scale_mode = "relative_ctrl100"
    y_plot = y_raw / y_scale_ref * 100
    y_note = "y_scale=relative_ctrl100"
  } else {
    y_scale_mode = "absolute_log2_norm"
    y_scale_ref = 0
    y_plot = y_raw
    y_note = "y_scale=absolute_log2_norm (control guides all zero; cannot set ctrl=100)"
  }

  guide_df = guide_meta
  guide_df$mean_log2_norm_count = unname(y_raw[guide_df$guide])
  guide_df$y_scaled = unname(y_plot[guide_df$guide])
  guide_df$y_scale_ref = y_scale_ref
  guide_df$y_scale_mode = y_scale_mode
  guide_df$x_scale_ref = x_scale_ref

  fit_data = data.frame(
    x=guide_df$x_scaled,
    y=guide_df$y_scaled,
    guide=guide_df$guide,
    guide_class=guide_df$guide_class,
    stringsAsFactors=F, check.names=F
  )
  pearson_r = suppressWarnings(cor(fit_data$x, fit_data$y, method="pearson"))
  spearman_rho = suppressWarnings(cor(fit_data$x, fit_data$y, method="spearman"))

  fit_linear = stats::lm(y ~ x, data=fit_data)
  # Standard (non-orthogonal) quadratic: y = a x^2 + b x + c, matching manuscript methods.
  fit_quadratic = stats::lm(y ~ x + I(x^2), data=fit_data)
  lin_coef = stats::coef(fit_linear)
  linear_intercept = unname(as.numeric(lin_coef["(Intercept)"]))
  linear_slope = unname(as.numeric(lin_coef["x"]))
  lin_coef_mat = summary(fit_linear)$coefficients
  pvalue_linear = if ("x" %in% rownames(lin_coef_mat)) {
    as.numeric(lin_coef_mat["x", "Pr(>|t|)"])
  } else {
    NA_real_
  }
  qcoef = stats::coef(fit_quadratic)
  quadratic_c = unname(as.numeric(qcoef["(Intercept)"]))
  quadratic_b = unname(as.numeric(qcoef["x"]))
  quadratic_a = unname(as.numeric(qcoef["I(x^2)"]))
  x_min = min(fit_data$x)
  x_max = max(fit_data$x)
  quadratic_vertex_x = if (is.finite(quadratic_a) && quadratic_a != 0) {
    -quadratic_b / (2 * quadratic_a)
  } else {
    NA_real_
  }
  quadratic_vertex_in_range = is.finite(quadratic_vertex_x) &&
    quadratic_vertex_x >= x_min && quadratic_vertex_x <= x_max
  response_min = min(fit_data$y)
  response_max = max(fit_data$y)
  ll4_ok_range = is.finite(response_min) && is.finite(response_max) && (response_max > response_min)
  fit_ll4_obj = if (ll4_ok_range) {
    fit_ll4(fit_data, response_min, response_max)
  } else {
    simpleError("max(y) <= min(y); cannot fix LL.4 asymptotes")
  }
  ll4_converged = !inherits(fit_ll4_obj, "error")
  hill_n = NA_real_
  ll4_ed50 = NA_real_
  if (ll4_converged) {
    ll4_coef = stats::coef(fit_ll4_obj)
    b_name = grep("^b", names(ll4_coef), value=T)
    e_name = grep("^e", names(ll4_coef), value=T)
    if (length(b_name) > 0L) hill_n = unname(as.numeric(ll4_coef[[b_name[[1L]]]]))
    if (length(e_name) > 0L) ll4_ed50 = unname(as.numeric(ll4_coef[[e_name[[1L]]]]))
  }

  aic_linear = as.numeric(stats::AIC(fit_linear))
  aic_quadratic = as.numeric(stats::AIC(fit_quadratic))
  aic_ll4 = if (ll4_converged) as.numeric(stats::AIC(fit_ll4_obj)) else NA_real_
  rmse_linear = sqrt(mean(stats::residuals(fit_linear)^2))
  rmse_quadratic = sqrt(mean(stats::residuals(fit_quadratic)^2))
  rmse_ll4 = if (ll4_converged) sqrt(mean(stats::residuals(fit_ll4_obj)^2)) else NA_real_

  aics = c(linear=aic_linear, quadratic=aic_quadratic, ll4=aic_ll4)
  finite = is.finite(aics)
  if (any(finite)) {
    best = min(aics[finite])
    delta = aics - best
    rel = ifelse(finite, exp(-0.5 * delta), 0)
    weights = rel / sum(rel)
    best_model = names(aics)[finite][which.min(aics[finite])]
  } else {
    delta = c(linear=NA_real_, quadratic=NA_real_, ll4=NA_real_)
    weights = delta
    best_model = NA_character_
  }

  ll4_note = if (ll4_converged) {
    sprintf("ll4 fixed c=response_min=%.6g; d=response_max=%.6g", response_min, response_max)
  } else if (inherits(fit_ll4_obj, "error")) {
    paste0("ll4_failed: ", conditionMessage(fit_ll4_obj))
  } else {
    "ll4_not_fit"
  }
  note = paste(y_note, ll4_note, sep=" | ")

  aic_row = data.frame(
    gene_id=gene_id,
    gene_name=display_name,
    direction=direction_label,
    pvalue_linear=pvalue_linear,
    linear_intercept=linear_intercept,
    linear_slope=linear_slope,
    hill_n=hill_n,
    ll4_ed50=ll4_ed50,
    quadratic_a=quadratic_a,
    quadratic_b=quadratic_b,
    quadratic_c=quadratic_c,
    quadratic_vertex_x=quadratic_vertex_x,
    quadratic_vertex_in_range=quadratic_vertex_in_range,
    x_min=x_min,
    x_max=x_max,
    aic_linear=aic_linear,
    aic_quadratic=aic_quadratic,
    aic_ll4=aic_ll4,
    delta_aic_linear=unname(delta[["linear"]]),
    delta_aic_quadratic=unname(delta[["quadratic"]]),
    delta_aic_ll4=unname(delta[["ll4"]]),
    weight_linear=unname(weights[["linear"]]),
    weight_quadratic=unname(weights[["quadratic"]]),
    weight_ll4=unname(weights[["ll4"]]),
    best_model=best_model,
    ll4_converged=ll4_converged,
    rmse_linear=rmse_linear,
    rmse_quadratic=rmse_quadratic,
    rmse_ll4=rmse_ll4,
    pearson_r=pearson_r,
    spearman_rho=spearman_rho,
    x_scale_ref=x_scale_ref,
    y_scale_ref=y_scale_ref,
    y_scale_mode=y_scale_mode,
    ll4_c_fixed=if (ll4_ok_range) response_min else NA_real_,
    ll4_d_fixed=if (ll4_ok_range) response_max else NA_real_,
    note=note,
    stringsAsFactors=F, check.names=F
  )

  x_grid = data.frame(x=seq(min(fit_data$x), max(fit_data$x), length.out=201L))
  pred_list = list(
    data.frame(model="linear", x=x_grid$x, y_hat=as.numeric(stats::predict(fit_linear, x_grid)), stringsAsFactors=F),
    data.frame(model="quadratic", x=x_grid$x, y_hat=as.numeric(stats::predict(fit_quadratic, x_grid)), stringsAsFactors=F)
  )
  if (ll4_converged) {
    pred_list[[3L]] = data.frame(
      model="ll4",
      x=x_grid$x,
      y_hat=as.numeric(stats::predict(fit_ll4_obj, x_grid)),
      stringsAsFactors=F
    )
  }
  pred_df = do.call(rbind, pred_list)
  present_models = intersect(model_levels, unique(pred_df$model))
  label_map = c(
    linear=sprintf("Linear (AIC=%.2f)", aic_linear),
    quadratic=sprintf("Quadratic (AIC=%.2f)", aic_quadratic)
  )
  if (ll4_converged) label_map = c(label_map, ll4=sprintf("LL.4 (c,d fixed; AIC=%.2f)", aic_ll4))
  pred_df$model = factor(pred_df$model, levels=present_models)
  pred_df$model_label = factor(
    unname(label_map[as.character(pred_df$model)]),
    levels=unname(label_map[present_models])
  )
  curve_colors = setNames(unname(model_colors[present_models]), levels(pred_df$model_label))

  # Coefficient table for single-gene mode extras.
  coef_table = rbind(
    data.frame(model="linear", term=names(stats::coef(fit_linear)), estimate=as.numeric(stats::coef(fit_linear)), stringsAsFactors=F),
    data.frame(model="quadratic", term=names(stats::coef(fit_quadratic)), estimate=as.numeric(stats::coef(fit_quadratic)), stringsAsFactors=F)
  )
  if (ll4_converged) {
    ll4_coef = stats::coef(fit_ll4_obj)
    coef_table = rbind(
      coef_table,
      data.frame(
        model="ll4",
        term=c(names(ll4_coef), "c:(fixed_lower)", "d:(fixed_upper)"),
        estimate=c(as.numeric(ll4_coef), response_min, response_max),
        stringsAsFactors=F
      )
    )
  }

  long_metrics = data.frame(
    model=c("linear", "quadratic", "ll4"),
    method=c("stats::lm", "stats::lm", "drc::LL.4(fixed=c(NA, response_min, response_max, NA))"),
    converged=c(T, T, ll4_converged),
    n_points=NROW(fit_data),
    n_parameters=c(length(stats::coef(fit_linear)), length(stats::coef(fit_quadratic)), if (ll4_converged) length(stats::coef(fit_ll4_obj)) else NA_integer_),
    aic=c(aic_linear, aic_quadratic, aic_ll4),
    delta_aic=c(unname(delta[["linear"]]), unname(delta[["quadratic"]]), unname(delta[["ll4"]])),
    akaike_weight=c(unname(weights[["linear"]]), unname(weights[["quadratic"]]), unname(weights[["ll4"]])),
    rmse=c(rmse_linear, rmse_quadratic, rmse_ll4),
    is_best_aic=c("linear", "quadratic", "ll4") == best_model,
    note=c(NA_character_, NA_character_, note),
    stringsAsFactors=F, check.names=F
  )
  long_metrics = long_metrics[order(is.na(long_metrics$aic), long_metrics$aic, long_metrics$model), , drop=F]
  rownames(long_metrics) = NULL

  if (identical(y_scale_mode, "relative_ctrl100")) {
    y_lab = paste0(
      display_name, " log-expression relative to ctrl\n",
      "y = 100 × mean(log2(n+1)) / mean(log2(n+1))_ctrl\n",
      "n = DESeq2 normalized counts"
    )
    scale_txt = "y=100*mean(log2(n+1))/mean(log2(n+1))_ctrl"
  } else {
    y_lab = paste0(
      display_name, " absolute log-expression\n",
      "y = mean(log2(n+1))   [control all zero; NOT relative to ctrl=100]\n",
      "n = DESeq2 normalized counts"
    )
    scale_txt = "y=absolute mean(log2(n+1)); control all zero"
  }
  dir_txt = direction_text(direction_label)

  p = ggplot2::ggplot() +
    ggplot2::geom_line(
      data=pred_df,
      ggplot2::aes(x=x, y=y_hat, color=model_label),
      linewidth=0.95
    ) +
    ggplot2::geom_point(
      data=fit_data,
      ggplot2::aes(x=x, y=y, shape=guide_class),
      color=plot_text, size=2.6, stroke=0.75
    ) +
    ggrepel::geom_text_repel(
      data=fit_data,
      ggplot2::aes(x=x, y=y, label=guide),
      size=2.8,
      color=plot_text,
      box.padding=0.28,
      point.padding=0.12,
      min.segment.length=0,
      segment.size=0.22,
      segment.color=plot_muted,
      max.overlaps=Inf,
      seed=1
    ) +
    ggplot2::scale_shape_manual(
      values=c(OCT4_targeting=17, nontargeting_control=16, safe_harbor_control=15),
      labels=c(
        OCT4_targeting="OCT4 targeting",
        nontargeting_control="Nontargeting control",
        safe_harbor_control="Safe-harbor control"
      ),
      name="Guides"
    ) +
    ggplot2::scale_color_manual(values=curve_colors, name="Model fit", drop=F) +
    ggplot2::guides(
      shape=ggplot2::guide_legend(order=1, direction="horizontal", title.position="left"),
      color=ggplot2::guide_legend(order=2, direction="horizontal", title.position="left")
    ) +
    ggplot2::labs(
      title=sprintf("%s (%s)", display_name, gene_id),
      subtitle=sprintf(
        "%s; best AIC=%s; %s; r=%.2f",
        dir_txt,
        ifelse(is.na(best_model), "none", best_model),
        scale_txt,
        pearson_r
      ),
      x="OCT4 qPCR relative to ctrl\nx = 100 × mean(qPCR) / mean(qPCR)_ctrl",
      y=y_lab
    ) +
    ggplot2::theme_bw(base_family="sans") +
    ggplot2::theme(
      plot.background=ggplot2::element_rect(fill=plot_surface, color=NA),
      panel.background=ggplot2::element_rect(fill=plot_surface, color=NA),
      panel.border=ggplot2::element_rect(color=plot_axis, fill=NA, linewidth=0.6),
      panel.grid.major=ggplot2::element_line(color=plot_grid, linewidth=0.3),
      panel.grid.minor=ggplot2::element_blank(),
      axis.title=ggplot2::element_text(color=plot_text, size=9.5, lineheight=1.05),
      axis.text=ggplot2::element_text(color=plot_text, size=8.5),
      plot.title=ggplot2::element_text(color=plot_text, size=13, hjust=0.5),
      plot.subtitle=ggplot2::element_text(color=plot_muted, size=8.2, hjust=0.5),
      legend.position="top",
      legend.box="vertical",
      legend.box.just="left",
      legend.spacing.y=grid::unit(2, "pt"),
      legend.title=ggplot2::element_text(size=8),
      legend.text=ggplot2::element_text(size=7.5),
      legend.background=ggplot2::element_blank(),
      legend.key=ggplot2::element_blank(),
      plot.margin=ggplot2::margin(10, 12, 10, 12)
    )

  list(
    ok=T,
    aic_row=aic_row,
    guide_df=guide_df,
    pred_df=pred_df,
    coef_table=coef_table,
    long_metrics=long_metrics,
    plot=p
  )
}

resolve_gene = function(DESeq2Obj, gene_query) {
  gene_names = as.character(mcols(DESeq2Obj)$gene_name)
  gene_ids = rownames(DESeq2Obj)
  hit = unique(c(
    which(toupper(gene_names) == toupper(gene_query)),
    which(gene_ids == gene_query)
  ))
  if (length(hit) == 0L) stop_with("Gene '%s' was not found by gene_name or gene_id", gene_query)
  if (length(hit) > 1L) {
    stop_with(
      "Gene query '%s' matched multiple rows: %s",
      gene_query,
      paste(sprintf("%s (%s)", gene_ids[hit], gene_names[hit]), collapse="; ")
    )
  }
  gene_id = gene_ids[[hit]]
  gene_name = gene_names[[hit]]
  if (is.na(gene_name) || gene_name == "") gene_name = gene_id
  list(gene_id=gene_id, gene_name=gene_name)
}

guide_log2_means_for_genes = function(DESeq2Obj, metadata, gene_ids) {
  gene_ids = unique(as.character(gene_ids))
  missing_ids = setdiff(gene_ids, rownames(DESeq2Obj))
  if (length(missing_ids) > 0L) {
    stop_with("%d gene_id values are absent from the DESeq2 object", length(missing_ids))
  }
  norm_mat = counts(DESeq2Obj, normalized=T)[gene_ids, metadata$sample_id, drop=F]
  log2_mat = log2(norm_mat + 1)
  out = sapply(guide_levels, function(g) {
    samples = metadata$sample_id[metadata$guide == g]
    matrixStats::rowMeans2(log2_mat[, samples, drop=F])
  })
  if (is.null(dim(out))) {
    out = matrix(out, nrow=1L, dimnames=list(gene_ids, guide_levels))
  } else {
    rownames(out) = gene_ids
  }
  out
}

# ------------------------------- run modes ------------------------------------
inputs = load_analysis_inputs()
DESeq2Obj = inputs$DESeq2Obj
metadata = inputs$metadata
guide_meta = inputs$guide_meta
x_scale_ref = inputs$x_scale_ref

if (!batch_mode) {
  resolved = resolve_gene(DESeq2Obj, gene_query)
  gene_id = resolved$gene_id
  gene_name = resolved$gene_name
  y_mat = guide_log2_means_for_genes(DESeq2Obj, metadata, gene_id)
  y_by_guide = setNames(as.numeric(y_mat[gene_id, guide_levels]), guide_levels)
  result = analyze_gene_points(
    gene_id=gene_id,
    display_name=gene_name,
    direction_label="",
    y_log2_by_guide=y_by_guide,
    guide_meta=guide_meta,
    x_scale_ref=x_scale_ref
  )
  safe_gene = gsub("[^A-Za-z0-9._-]+", "_", gene_name)
  if (!result$ok) {
    # Still emit a fallback figure for edge cases (e.g. all-zero control y).
    if (is.null(result$plot)) {
      stop_with("Single-gene analysis failed: %s", result$aic_row$note)
    }
    fit_plot_path = file.path(output_dir, sprintf("scatter-gene-vs-qpcr-%s-model-fits.png", safe_gene))
    ggplot2::ggsave(fit_plot_path, result$plot, width=9, height=6.8, units="in", dpi=300, bg=plot_surface)
    message("Gene: ", gene_name, " (", gene_id, ")")
    message("Single-gene analysis skipped relative fit: ", result$aic_row$note)
    message("Fallback plot: ", fit_plot_path)
    devtools::session_info(pkgs="attached")
    quit(save="no", status=0)
  }
  table_path = file.path(output_dir, sprintf("scatter-gene-vs-qpcr-%s-by-guide.csv", safe_gene))
  metrics_path = file.path(output_dir, sprintf("scatter-gene-vs-qpcr-%s-model-metrics.csv", safe_gene))
  coef_path = file.path(output_dir, sprintf("scatter-gene-vs-qpcr-%s-model-coefficients.csv", safe_gene))
  pred_path = file.path(output_dir, sprintf("scatter-gene-vs-qpcr-%s-model-predictions.csv", safe_gene))
  plot_path = file.path(output_dir, sprintf("scatter-gene-vs-qpcr-%s-by-guide.png", safe_gene))
  fit_plot_path = file.path(output_dir, sprintf("scatter-gene-vs-qpcr-%s-model-fits.png", safe_gene))

  write.csv(result$guide_df, table_path, row.names=F, quote=T, na="NA")
  write.csv(result$long_metrics, metrics_path, row.names=F, quote=T, na="NA")
  write.csv(result$coef_table, coef_path, row.names=F, quote=T, na="NA")
  write.csv(result$pred_df, pred_path, row.names=F, quote=T, na="NA")

  # Points-only companion plot (no curves).
  points_only = ggplot2::ggplot(
    result$guide_df,
    ggplot2::aes(x=x_scaled, y=y_scaled, color=guide_class, shape=guide_class)
  ) +
    ggplot2::geom_point(size=3.0, stroke=0.8) +
    ggrepel::geom_text_repel(
      ggplot2::aes(label=guide), size=3.0, color=plot_text,
      box.padding=0.3, point.padding=0.15, min.segment.length=0,
      segment.size=0.25, segment.color=plot_muted, max.overlaps=Inf, seed=1, show.legend=F
    ) +
    ggplot2::scale_color_manual(
      values=class_colors,
      labels=c(OCT4_targeting="OCT4 targeting", nontargeting_control="Nontargeting control", safe_harbor_control="Safe-harbor control"),
      name=NULL
    ) +
    ggplot2::scale_shape_manual(
      values=c(OCT4_targeting=17, nontargeting_control=16, safe_harbor_control=15),
      labels=c(OCT4_targeting="OCT4 targeting", nontargeting_control="Nontargeting control", safe_harbor_control="Safe-harbor control"),
      name=NULL
    ) +
    ggplot2::labs(
      title=sprintf("%s expression vs OCT4 qPCR", gene_name),
      subtitle=sprintf(
        "11 guide means; y_scale=%s; Pearson r=%.2f; Spearman rho=%.2f",
        result$aic_row$y_scale_mode, result$aic_row$pearson_r, result$aic_row$spearman_rho
      ),
      x="OCT4 qPCR relative to ctrl\nx = 100 × mean(qPCR) / mean(qPCR)_ctrl",
      y=if (identical(result$aic_row$y_scale_mode, "relative_ctrl100")) {
        paste0(
          gene_name, " log-expression relative to ctrl\n",
          "y = 100 × mean(log2(n+1)) / mean(log2(n+1))_ctrl\n",
          "n = DESeq2 normalized counts"
        )
      } else {
        paste0(
          gene_name, " absolute log-expression\n",
          "y = mean(log2(n+1))   [control all zero; NOT relative to ctrl=100]\n",
          "n = DESeq2 normalized counts"
        )
      }
    ) +
    ggplot2::theme_bw(base_family="sans") +
    ggplot2::theme(
      plot.background=ggplot2::element_rect(fill=plot_surface, color=NA),
      panel.background=ggplot2::element_rect(fill=plot_surface, color=NA),
      panel.border=ggplot2::element_rect(color=plot_axis, fill=NA, linewidth=0.6),
      panel.grid.major=ggplot2::element_line(color=plot_grid, linewidth=0.35),
      panel.grid.minor=ggplot2::element_blank(),
      axis.title=ggplot2::element_text(color=plot_text, size=10, lineheight=1.1),
      axis.text=ggplot2::element_text(color=plot_text, size=10),
      plot.title=ggplot2::element_text(color=plot_text, size=14, hjust=0.5),
      plot.subtitle=ggplot2::element_text(color=plot_muted, size=8.5, hjust=0.5),
      legend.position="top",
      legend.background=ggplot2::element_blank(),
      legend.key=ggplot2::element_blank(),
      plot.margin=ggplot2::margin(12, 16, 14, 16)
    )

  ggplot2::ggsave(plot_path, points_only, width=8.5, height=6.5, units="in", dpi=300, bg=plot_surface)
  ggplot2::ggsave(fit_plot_path, result$plot, width=9, height=6.8, units="in", dpi=300, bg=plot_surface)

  message("Gene: ", gene_name, " (", gene_id, ")")
  message("Points: 11 guide means")
  message(sprintf("Scaling refs: qPCR=%.4f; log2(norm+1)=%.4f", x_scale_ref, result$aic_row$y_scale_ref))
  message(sprintf("Pearson r=%.3f; Spearman rho=%.3f", result$aic_row$pearson_r, result$aic_row$spearman_rho))
  message("AIC comparison:")
  for (i in seq_len(NROW(result$long_metrics))) {
    row = result$long_metrics[i, , drop=F]
    message(sprintf(
      "  %-12s  AIC=%s  delta=%s  weight=%s  converged=%s",
      row$model,
      ifelse(is.finite(row$aic), sprintf("%.3f", row$aic), "NA"),
      ifelse(is.finite(row$delta_aic), sprintf("%.3f", row$delta_aic), "NA"),
      ifelse(is.finite(row$akaike_weight), sprintf("%.3f", row$akaike_weight), "NA"),
      row$converged
    ))
  }
  message("Table: ", table_path)
  message("Plot: ", plot_path)
  message("Model-fit plot: ", fit_plot_path)
  message("Model metrics: ", metrics_path)
  devtools::session_info(pkgs="attached")
  quit(save="no", status=0)
}

# ------------------------------- batch modes -----------------------------------
annotation_columns = c(
  "gene_biotype", "direction_status",
  "n_guides_deg", "n_guides_up", "n_guides_down", "deg_guides",
  "median_lfc_shrunken_used_for_direction", "median_lfc_raw_pairwise_deg",
  "min_padj_pairwise", "baseMean",
  "n_ids_collapsed", "dropped_gene_ids"
)

merge_annotation_columns = function(aic_table, gene_table) {
  extra = intersect(annotation_columns, names(gene_table))
  if (length(extra) == 0L) return(aic_table)
  lookup = gene_table[, c("gene_id", extra), drop=F]
  lookup = lookup[!duplicated(lookup$gene_id), , drop=F]
  merged = merge(aic_table, lookup, by="gene_id", all.x=T, sort=F)
  merged = merged[match(aic_table$gene_id, merged$gene_id), , drop=F]
  rownames(merged) = NULL
  preferred = c(
    "gene_id", "gene_name", "gene_biotype", "direction", "direction_status",
    "n_guides_deg", "n_guides_up", "n_guides_down", "deg_guides",
    "median_lfc_shrunken_used_for_direction", "median_lfc_raw_pairwise_deg",
    "min_padj_pairwise", "baseMean", "n_ids_collapsed", "dropped_gene_ids",
    "pvalue_linear", "linear_intercept", "linear_slope", "hill_n", "ll4_ed50",
    "quadratic_a", "quadratic_b", "quadratic_c",
    "quadratic_vertex_x", "quadratic_vertex_in_range", "x_min", "x_max"
  )
  front = intersect(preferred, names(merged))
  rest = setdiff(names(merged), front)
  merged[, c(front, rest), drop=F]
}

assert_unique_named_genes = function(aic_table, tag) {
  names_used = as.character(aic_table$gene_name)
  named = !is.na(names_used) & names_used != "" & names_used != as.character(aic_table$gene_id)
  dups = unique(names_used[named][duplicated(names_used[named])])
  if (length(dups) > 0L) {
    stop_with(
      "Duplicate gene_name values remain in %s after collapse: %s",
      tag,
      paste(dups, collapse=", ")
    )
  }
  invisible(TRUE)
}

biotype_rank = function(biotype) {
  biotype = ifelse(is.na(biotype), "", as.character(biotype))
  ifelse(biotype == "protein_coding", 1L,
    ifelse(biotype %in% c("lncRNA", "lincRNA"), 2L, 3L))
}

collapse_duplicate_gene_names = function(gene_table) {
  empty_resolution = data.frame(
    gene_name=character(),
    selected_gene_id=character(),
    selected_gene_biotype=character(),
    selected_baseMean=numeric(),
    n_ids=integer(),
    dropped_gene_ids=character(),
    dropped_gene_biotypes=character(),
    rule=character(),
    stringsAsFactors=F,
    check.names=F
  )
  if (NROW(gene_table) == 0L) {
    gene_table$n_ids_collapsed = integer()
    gene_table$dropped_gene_ids = character()
    return(list(table=gene_table, resolution=empty_resolution))
  }
  gene_table$n_ids_collapsed = 1L
  gene_table$dropped_gene_ids = NA_character_
  key = as.character(gene_table$gene_name)
  empty = is.na(key) | key == ""
  named_idx = which(!empty)
  if (length(named_idx) == 0L) return(list(table=gene_table, resolution=empty_resolution))

  groups = split(named_idx, key[named_idx])
  keep = rep(TRUE, NROW(gene_table))
  resolution_rows = list()
  rule = paste(
    "prefer protein_coding, then higher baseMean,",
    "then smaller min pairwise padj, then gene_id"
  )
  for (nm in names(groups)) {
    idx = groups[[nm]]
    if (length(idx) == 1L) next
    sub = gene_table[idx, , drop=F]
    rank = biotype_rank(sub$gene_biotype)
    bmean = as.numeric(sub$baseMean)
    bmean[!is.finite(bmean)] = -Inf
    padj = as.numeric(sub$min_padj_pairwise)
    padj[!is.finite(padj)] = Inf
    ord = order(rank, -bmean, padj, sub$gene_id, method="radix")
    winner = idx[[ord[[1L]]]]
    losers = idx[ord[-1L]]
    keep[losers] = FALSE
    gene_table$n_ids_collapsed[[winner]] = length(idx)
    gene_table$dropped_gene_ids[[winner]] = paste(gene_table$gene_id[losers], collapse=";")
    resolution_rows[[length(resolution_rows) + 1L]] = data.frame(
      gene_name=nm,
      selected_gene_id=gene_table$gene_id[[winner]],
      selected_gene_biotype=gene_table$gene_biotype[[winner]],
      selected_baseMean=as.numeric(gene_table$baseMean[[winner]]),
      n_ids=length(idx),
      dropped_gene_ids=paste(gene_table$gene_id[losers], collapse=";"),
      dropped_gene_biotypes=paste(gene_table$gene_biotype[losers], collapse=";"),
      rule=rule,
      stringsAsFactors=F,
      check.names=F
    )
  }
  list(
    table=gene_table[keep, , drop=F],
    resolution=if (length(resolution_rows) == 0L) empty_resolution else do.call(rbind, resolution_rows)
  )
}

load_union_pairwise_degs = function(guide_dir, alpha_cutoff, lfc_cutoff) {
  # A rerun of stage 01 may leave old guide tables beside a new pooled object.
  # Check the exact pooled inputs used by this downstream analysis before
  # mixing current expression with guide-specific DEG membership/direction.
  manifest_path = file.path(guide_dir, "pooled-source-md5.csv")
  if (!file.exists(manifest_path)) {
    stop_with("Missing guide source manifest; run 03_run_guide_specific_deseq2.R first")
  }
  source_manifest = read.csv(manifest_path, check.names=F, stringsAsFactors=F)
  required_sources = c(dds_path, meta_path)
  if (!all(c("file", "size_bytes", "md5") %in% names(source_manifest)) ||
      anyDuplicated(source_manifest$file) ||
      !all(basename(required_sources) %in% source_manifest$file)) {
    stop_with("Invalid guide source manifest; rerun 03_run_guide_specific_deseq2.R")
  }
  recorded = source_manifest[match(basename(required_sources), source_manifest$file), , drop=F]
  if (anyNA(recorded$md5) || anyNA(recorded$size_bytes) ||
      !identical(as.character(recorded$md5), unname(tools::md5sum(required_sources))) ||
      !identical(as.numeric(recorded$size_bytes), as.numeric(file.info(required_sources)$size))) {
    stop_with("Stale guide-specific outputs: pooled sources changed; rerun 03_run_guide_specific_deseq2.R")
  }
  target_guides = c(
    "guide01", "guide02", "guide03", "guide04",
    "guide05", "guide06", "guide14", "guide15"
  )
  pieces = lapply(target_guides, function(g) {
    # Complete tables permit a different downstream alpha and retain padj == alpha.
    path = file.path(guide_dir, sprintf("de-results-%s-vs-pooled-controls-complete.csv", g))
    if (!file.exists(path)) stop_with("Missing %s; run 03_run_guide_specific_deseq2.R first", path)
    d = read.csv(path, check.names=F, na.strings="NA")
    required = c(
      "gene_id", "gene_name", "gene_biotype", "baseMean", "padj",
      "log2FoldChange_raw", "log2FoldChange_shrunken"
    )
    missing = setdiff(required, names(d))
    if (length(missing) > 0L) {
      stop_with("%s is missing required columns: %s", path, paste(missing, collapse=", "))
    }
    padj = as.numeric(d$padj)
    lfc_raw = as.numeric(d$log2FoldChange_raw)
    lfc_shr = as.numeric(d$log2FoldChange_shrunken)
    keep = !is.na(padj) & padj <= alpha_cutoff & is.finite(lfc_raw) & abs(lfc_raw) >= lfc_cutoff
    d = d[keep, , drop=F]
    if (NROW(d) == 0L) {
      return(data.frame(
        target_guide=character(), gene_id=character(), gene_name=character(),
        gene_biotype=character(), baseMean=numeric(), padj=numeric(),
        log2FoldChange_raw=numeric(), log2FoldChange_shrunken=numeric(),
        pairwise_direction=character(), stringsAsFactors=F, check.names=F
      ))
    }
    data.frame(
      target_guide=g,
      gene_id=as.character(d$gene_id),
      gene_name=as.character(d$gene_name),
      gene_biotype=as.character(d$gene_biotype),
      baseMean=as.numeric(d$baseMean),
      padj=as.numeric(d$padj),
      log2FoldChange_raw=lfc_raw[keep],
      log2FoldChange_shrunken=lfc_shr[keep],
      pairwise_direction=ifelse(
        lfc_shr[keep] > 0, "up",
        ifelse(lfc_shr[keep] < 0, "down", "zero")
      ),
      stringsAsFactors=F,
      check.names=F
    )
  })
  long = do.call(rbind, pieces)
  rownames(long) = NULL
  if (NROW(long) == 0L) {
    stop_with("Pairwise union DEG list is empty at padj<=%s and |raw LFC|>=%s", alpha_cutoff, lfc_cutoff)
  }
  if (any(long$pairwise_direction == "zero")) {
    stop_with("Pairwise strict DEGs include genes with shrunken LFC == 0")
  }
  split_ids = split(seq_len(NROW(long)), long$gene_id)
  rows = lapply(split_ids, function(idx) {
    sub = long[idx, , drop=F]
    sub = sub[order(sub$target_guide, method="radix"), , drop=F]
    n_up = sum(sub$pairwise_direction == "up")
    n_down = sum(sub$pairwise_direction == "down")
    if (n_up > 0L && n_down == 0L) {
      direction = "up"
      direction_status = "consistent_up"
    } else if (n_down > 0L && n_up == 0L) {
      direction = "down"
      direction_status = "consistent_down"
    } else {
      med = stats::median(sub$log2FoldChange_shrunken)
      if (!is.finite(med) || med == 0) med = stats::median(sub$log2FoldChange_raw)
      direction = if (is.finite(med) && med < 0) "down" else "up"
      direction_status = "mixed"
    }
    data.frame(
      gene_id=sub$gene_id[[1L]],
      gene_name=sub$gene_name[[1L]],
      gene_biotype=sub$gene_biotype[[1L]],
      baseMean=sub$baseMean[[1L]],
      min_padj_pairwise=min(sub$padj, na.rm=T),
      median_lfc_shrunken_used_for_direction=stats::median(sub$log2FoldChange_shrunken),
      median_lfc_raw_pairwise_deg=stats::median(sub$log2FoldChange_raw),
      n_guides_deg=NROW(sub),
      n_guides_up=n_up,
      n_guides_down=n_down,
      deg_guides=paste(sub$target_guide, collapse=";"),
      direction=direction,
      direction_status=direction_status,
      stringsAsFactors=F,
      check.names=F
    )
  })
  gene_table = do.call(rbind, rows)
  rownames(gene_table) = NULL
  gene_table$display_name = ifelse(
    is.na(gene_table$gene_name) | gene_table$gene_name == "",
    gene_table$gene_id,
    gene_table$gene_name
  )
  collapsed = collapse_duplicate_gene_names(gene_table)
  list(table=collapsed$table, resolution=collapsed$resolution, long=long)
}

process_gene_table = function(gene_table, direction_label, pdf_path, aic_path, y_mat, max_genes=Inf, progress_tag=NULL) {
  if (!all(c("gene_id", "display_name") %in% names(gene_table))) {
    stop_with("gene_table must contain gene_id and display_name")
  }
  sub = gene_table
  if ("direction" %in% names(sub) && !is.null(direction_label) && nzchar(direction_label)) {
    sub = sub[sub$direction == direction_label, , drop=F]
  }
  ord = order(sub$display_name, sub$gene_id, method="radix")
  sub = sub[ord, , drop=F]
  if (is.finite(max_genes) && NROW(sub) > max_genes) {
    sub = sub[seq_len(max_genes), , drop=F]
  }
  n = NROW(sub)
  tag = if (is.null(progress_tag) || !nzchar(progress_tag)) {
    if (is.null(direction_label) || !nzchar(direction_label)) "all" else direction_label
  } else {
    progress_tag
  }
  message(sprintf(
    "Batch set %s: %d genes%s -> %s",
    tag,
    n,
    if (is.finite(max_genes)) sprintf(" (preview max=%d)", max_genes) else "",
    if (csv_only) aic_path else pdf_path
  ))
  if (n == 0L) {
    if (!csv_only && file.exists(pdf_path)) {
      unlink(pdf_path)
      if (file.exists(pdf_path)) stop_with("Could not remove stale empty-direction PDF: %s", pdf_path)
    }
    empty = merge_annotation_columns(
      make_failed_aic_row(NA_character_, NA_character_, NA_character_, NA_real_, "")[FALSE, , drop=F],
      sub
    )
    write.csv(empty, aic_path, row.names=F, quote=T, na="NA")
    message(sprintf("Wrote empty AIC table: %s", aic_path))
    return(empty)
  }

  aic_rows = vector("list", n)
  if (!csv_only) {
    grDevices::pdf(pdf_path, width=9, height=6.8, onefile=T, useDingbats=F)
    on.exit(try(grDevices::dev.off(), silent=T), add=T)
  }

  for (i in seq_len(n)) {
    gene_id = sub$gene_id[[i]]
    display_name = sub$display_name[[i]]
    dir_label = if ("direction" %in% names(sub)) {
      as.character(sub$direction[[i]])
    } else if (!is.null(direction_label) && nzchar(direction_label)) {
      direction_label
    } else {
      "all"
    }
    if (i == 1L || i == n || (i %% 100L) == 0L) {
      message(sprintf("  [%s] %d / %d  %s", tag, i, n, display_name))
    }
    y_by_guide = setNames(as.numeric(y_mat[gene_id, guide_levels]), guide_levels)
    result = analyze_gene_points(
      gene_id=gene_id,
      display_name=display_name,
      direction_label=dir_label,
      y_log2_by_guide=y_by_guide,
      guide_meta=guide_meta,
      x_scale_ref=x_scale_ref
    )
    aic_rows[[i]] = result$aic_row
    if (!csv_only) {
      if (!is.null(result$plot)) {
        print(result$plot)
      } else {
        print(make_zero_input_fallback_plot(
          gene_id=gene_id,
          display_name=display_name,
          direction_label=dir_label,
          which_axis="y",
          detail_note=result$aic_row$note[[1L]]
        ))
      }
    }
  }

  aic_table = do.call(rbind, aic_rows)
  rownames(aic_table) = NULL
  if (!identical(aic_table$gene_id, sub$gene_id)) {
    stop_with("AIC table order diverged from gene order for %s", tag)
  }
  aic_table = merge_annotation_columns(aic_table, sub)
  write.csv(aic_table, aic_path, row.names=F, quote=T, na="NA")
  message(sprintf("Wrote AIC table (%d rows): %s", NROW(aic_table), aic_path))
  invisible(aic_table)
}

preview_suffix = paste0(
  if (identical(batch_scope, "deg_pooled")) "-pooled" else "",
  if (is.finite(batch_max_genes)) sprintf("-preview%d", batch_max_genes) else ""
)

run_deg_batch = function(deg_table, label) {
  message(sprintf("Computing guide-mean log2(normalized count + 1) for %d %s DEGs...", NROW(deg_table), label))
  y_mat = guide_log2_means_for_genes(DESeq2Obj, metadata, deg_table$gene_id)

  up_pdf = file.path(output_dir, sprintf("scatter-deg-up-vs-qpcr-model-fits%s.pdf", preview_suffix))
  down_pdf = file.path(output_dir, sprintf("scatter-deg-down-vs-qpcr-model-fits%s.pdf", preview_suffix))
  up_aic = file.path(output_dir, sprintf("scatter-deg-up-vs-qpcr-model-aic%s.csv", preview_suffix))
  down_aic = file.path(output_dir, sprintf("scatter-deg-down-vs-qpcr-model-aic%s.csv", preview_suffix))
  combined_aic = file.path(output_dir, sprintf("scatter-deg-vs-qpcr-model-aic%s.csv", preview_suffix))

  up_table = process_gene_table(
    deg_table, "up", up_pdf, up_aic, y_mat,
    max_genes=batch_max_genes, progress_tag="up"
  )
  down_table = process_gene_table(
    deg_table, "down", down_pdf, down_aic, y_mat,
    max_genes=batch_max_genes, progress_tag="down"
  )
  combined = rbind(up_table, down_table)
  rownames(combined) = NULL
  if (NROW(combined) > 0L) {
    combined = combined[order(combined$gene_name, combined$gene_id, method="radix"), , drop=F]
    rownames(combined) = NULL
    assert_unique_named_genes(combined, label)
  }
  write.csv(combined, combined_aic, row.names=F, quote=T, na="NA")

  message(sprintf("Done %s.", label))
  if (!csv_only) {
    message("UP PDF: ", if (file.exists(up_pdf)) up_pdf else "[no genes; no PDF]")
    message("DOWN PDF: ", if (file.exists(down_pdf)) down_pdf else "[no genes; no PDF]")
  }
  message("UP AIC: ", up_aic)
  message("DOWN AIC: ", down_aic)
  message("COMBINED AIC: ", combined_aic)
  message(sprintf(
    "UP best-model counts: %s",
    paste(sprintf("%s=%d", names(table(up_table$best_model, useNA="ifany")), as.integer(table(up_table$best_model, useNA="ifany"))), collapse=", ")
  ))
  message(sprintf(
    "DOWN best-model counts: %s",
    paste(sprintf("%s=%d", names(table(down_table$best_model, useNA="ifany")), as.integer(table(down_table$best_model, useNA="ifany"))), collapse=", ")
  ))
  if ("pvalue_linear" %in% names(combined) && NROW(combined) > 0L) {
    message(sprintf(
      "Linear p < 0.05: %d / %d",
      sum(is.finite(combined$pvalue_linear) & combined$pvalue_linear < 0.05),
      NROW(combined)
    ))
  }
  invisible(combined)
}

if (identical(batch_scope, "deg")) {
  message(sprintf(
    "Batch mode: union of pairwise OCT4-gRNA vs pooled-control DEGs (padj<=%s and |raw LFC|>=%s)",
    format(alpha, trim=T, scientific=F),
    format(deg_lfc_threshold, trim=T, scientific=F)
  ))
  union_obj = load_union_pairwise_degs(guide_output_dir, alpha, deg_lfc_threshold)
  deg_table = union_obj$table
  membership_path = file.path(
    output_dir,
    sprintf("degs-union-oct4-guides-vs-pooled-controls-lfc%s-padj%s%s.csv", deg_lfc_threshold, alpha, preview_suffix)
  )
  resolution_path = file.path(
    output_dir,
    sprintf("scatter-deg-duplicate-resolution%s.csv", preview_suffix)
  )
  write.csv(deg_table, membership_path, row.names=F, quote=T, na="NA")
  write.csv(union_obj$resolution, resolution_path, row.names=F, quote=T, na="NA")
  message(sprintf(
    "Union DEG list: %d genes after collapsing %d duplicate symbols (%d long-form pairwise hits)",
    NROW(deg_table),
    NROW(union_obj$resolution),
    NROW(union_obj$long)
  ))
  message("Union membership: ", membership_path)
  message("Duplicate resolution: ", resolution_path)
  if (NROW(union_obj$resolution) > 0L) {
    message(sprintf(
      "Collapsed symbols include: %s",
      paste(union_obj$resolution$gene_name, collapse=", ")
    ))
  }
  run_deg_batch(deg_table, "batch-deg-union")
} else if (identical(batch_scope, "deg_pooled")) {
  message("Batch mode: pooled 8-vs-3 strict heatmap DEGs (up/down), alphabetical gene names")
  deg_table = read.csv(deg_path, check.names=F, na.strings="NA")
  if (!all(c("gene_id", "gene_name", "log2FoldChange_raw") %in% names(deg_table))) {
    stop_with("DEG table is missing required columns")
  }
  deg_table$gene_name = as.character(deg_table$gene_name)
  deg_table$display_name = ifelse(
    is.na(deg_table$gene_name) | deg_table$gene_name == "",
    deg_table$gene_id,
    deg_table$gene_name
  )
  deg_table$direction = ifelse(
    deg_table$log2FoldChange_raw > 0, "up",
    ifelse(deg_table$log2FoldChange_raw < 0, "down", "zero")
  )
  if (any(deg_table$direction == "zero")) stop_with("Strict DEG table contains genes with raw LFC == 0")
  if ("padj" %in% names(deg_table)) deg_table$min_padj_pairwise = as.numeric(deg_table$padj)
  collapsed = collapse_duplicate_gene_names(deg_table)
  deg_table = collapsed$table
  resolution_path = file.path(
    output_dir,
    sprintf("scatter-deg-pooled-duplicate-resolution%s.csv", preview_suffix)
  )
  write.csv(collapsed$resolution, resolution_path, row.names=F, quote=T, na="NA")
  message(sprintf(
    "Pooled DEG list: %d genes after collapsing %d duplicate symbols",
    NROW(deg_table),
    NROW(collapsed$resolution)
  ))
  run_deg_batch(deg_table, "batch-deg-pooled")
} else if (identical(batch_scope, "all")) {
  message("Batch mode: ALL genes in DESeq2 object, alphabetical gene names")
  gene_ids = rownames(DESeq2Obj)
  gene_names = as.character(mcols(DESeq2Obj)$gene_name)
  display_name = ifelse(is.na(gene_names) | gene_names == "", gene_ids, gene_names)
  all_table = data.frame(
    gene_id = gene_ids,
    gene_name = gene_names,
    display_name = display_name,
    direction = "all",
    stringsAsFactors=F,
    check.names=F
  )

  # Optional annotation: mark strict DEGs when the DEG table is available.
  if (file.exists(deg_path)) {
    deg_table = read.csv(deg_path, check.names=F, na.strings="NA")
    if (all(c("gene_id", "log2FoldChange_raw") %in% names(deg_table))) {
      deg_dir = ifelse(
        deg_table$log2FoldChange_raw > 0, "up",
        ifelse(deg_table$log2FoldChange_raw < 0, "down", "zero")
      )
      names(deg_dir) = as.character(deg_table$gene_id)
      hit = deg_dir[all_table$gene_id]
      all_table$direction = ifelse(is.na(hit), "not_deg", unname(hit))
    }
  }

  message(sprintf(
    "Computing guide-mean log2(normalized count + 1) for %d genes...",
    NROW(all_table)
  ))
  y_mat = guide_log2_means_for_genes(DESeq2Obj, metadata, all_table$gene_id)

  all_pdf = file.path(output_dir, sprintf("scatter-all-genes-vs-qpcr-model-fits%s.pdf", preview_suffix))
  all_aic = file.path(output_dir, sprintf("scatter-all-genes-vs-qpcr-model-aic%s.csv", preview_suffix))

  all_result = process_gene_table(
    all_table,
    direction_label=NULL,
    pdf_path=all_pdf,
    aic_path=all_aic,
    y_mat=y_mat,
    max_genes=batch_max_genes,
    progress_tag="all"
  )

  message("Done batch-all.")
  message("ALL PDF: ", all_pdf)
  message("ALL AIC: ", all_aic)
  message(sprintf(
    "ALL best-model counts: %s",
    paste(sprintf("%s=%d", names(table(all_result$best_model, useNA="ifany")), as.integer(table(all_result$best_model, useNA="ifany"))), collapse=", ")
  ))
  if ("direction" %in% names(all_result)) {
    message(sprintf(
      "ALL direction counts: %s",
      paste(sprintf("%s=%d", names(table(all_result$direction, useNA="ifany")), as.integer(table(all_result$direction, useNA="ifany"))), collapse=", ")
    ))
  }
} else {
  stop_with("Unknown batch_scope: %s", batch_scope)
}

devtools::session_info(pkgs="attached")
