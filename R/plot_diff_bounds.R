#' Plot differentiation statistics against M with their upper bounds
#'
#' Plots per-locus values of \eqn{F_{ST}}, \eqn{G'_{ST}} and \eqn{D} against
#' \eqn{M}, the mean frequency of the most frequent allele, together with the
#' upper bound of each statistic for `K` populations. Dashed red lines and a
#' red point mark the means over loci. Optionally, confidence intervals of the
#' means are drawn, either as shaded bands or as crossed error bars.
#'
#' @param M Numeric vector of \eqn{M} values, one per locus.
#' @param FST Numeric vector of \eqn{F_{ST}} values, one per locus.
#' @param GpST Optional numeric vector of \eqn{G'_{ST}} values.
#' @param D Optional numeric vector of Jost's \eqn{D} values.
#' @param K Number of populations.
#' @param show_colorbar Logical; if `TRUE`, show the colour bar for the
#'   relative point density (only used when `density_colors = TRUE`).
#' @param density_colors Logical; if `TRUE` (default), colour points by
#'   relative point density, otherwise draw all points in `point_color`.
#' @param relative_density Logical; if `TRUE` (default), colour points by
#'   density relative to the maximum (0-1, with a viridis scale). If `FALSE`,
#'   colour by the raw density and set no colour scale, so that a common scale
#'   can be added to several plots.
#' @param point_color Colour of the points when `density_colors = FALSE`.
#' @param ci Optional data frame with columns `statistic`, `value_lower` and
#'   `value_upper`, e.g. the `mean` element of [boot_stat()]. If `NULL`
#'   (default), no intervals are drawn.
#' @param mean_style How to mark the mean: `"crosshair"` (default; dashed red
#'   lines through the mean, across the whole plot) or `"bound"` (solid red
#'   segment from 0 to the upper bound at the mean of \eqn{M}).
#' @param mean_linetype Line type of the red segment when `mean_style = "bound"`,
#'   e.g. `"solid"` (default), `"dashed"` or `"dotted"`.
#' @param ci_style How to draw the intervals: `"bands"` (shaded bands across
#'   the plot) or `"cross"` (crossed error bars centred on the intervals).
#' @param ci_color Colour of the error bars when `ci_style = "cross"`.
#' @param ci_linewidth Line width of the error bars when `ci_style = "cross"`.
#' @param ci_cap Width of the caps at the ends of the error bars when
#'   `ci_style = "cross"`, in axis units (axes range from 0 to 1, so `0.05`
#'   gives caps 5\% of the axis wide). For the horizontal bar (interval of
#'   \eqn{M}), it sets the cap height. Use `0` for bars without caps.
#' @param ci_alpha Transparency of the confidence bands when `ci_style = "bands"`,
#'   between 0 (invisible) and 1 (opaque).
#' @param ci_names Named character vector giving, for `M`, `FST`, `GpST` and
#'   `D`, the label used for that quantity in the `statistic` column of `ci`.
#' @return A list of three ggplot objects (\eqn{F_{ST}}, \eqn{G'_{ST}}, \eqn{D});
#'   an element is `NULL` if the corresponding statistic was not supplied.
#' @import ggplot2
#' @export
ggbounds_raw <- function(M, FST, GpST = NULL, D = NULL, K = 2, show_colorbar = FALSE,
                         density_colors = TRUE, relative_density = TRUE, density_limits = NULL,
                         point_color = "grey30", ci = NULL, mean_style = c("crosshair", "bound"),
                         mean_linetype = "solid", ci_style = c("bands", "cross"),
                         ci_color = "black", ci_linewidth = 1.5, ci_cap = 0.05, ci_alpha = 0.15,
                         ci_names = c(M = "M", FST = "FST", GpST = "G'ST", D = "D")) {
  ci_style <- match.arg(ci_style)
  mean_style <- match.arg(mean_style)

  M_grid <- seq(0.001, 0.999, 0.001)
  mean_M <- mean(M, na.rm = TRUE)

  # ---- helpers ----------------------------------------------------

  # one row of ci for a given quantity ("M", "FST", "GpST" or "D")
  ci_row <- function(stat) {
    r <- ci[ci$statistic == ci_names[[stat]], ]
    if (nrow(r) != 1) warning("Expected one row in 'ci' for '", ci_names[[stat]], "'; check ci_names.")
    r
  }

  # per-locus points: coloured by density, or in a single colour
  point_layers <- function() {
    if (!density_colors) return(list(geom_point(size = 1.4, col = point_color)))
    colour_aes <- if (relative_density) aes(color = after_stat(density / max(density)))
                  else aes(color = after_stat(density))
    list(ggpointdensity::geom_pointdensity(colour_aes, size = 1.4),
         labs(color = if (relative_density) "relative\ndensity" else "density"),
         scale_color_viridis_c(limits = density_limits, oob = scales::squish))
  }

  # lines marking the mean: dashed crosshair, or solid segment from 0 to the bound at mean M
  mean_lines <- function(mean_y, bound_fun) {
    if (mean_style == "crosshair") {
      list(annotate("segment", x = mean_M, xend = mean_M, y = 0, yend = 1,
                    col = "red", linewidth = 0.4, linetype = "dashed"),
           annotate("segment", x = 0, xend = 1, y = mean_y, yend = mean_y,
                    col = "red", linewidth = 0.4, linetype = "dashed"))
    } else {
      list(annotate("segment", x = mean_M, xend = mean_M, y = 0, yend = bound_fun(K, mean_M),
                    col = "red", linewidth = 0.6, linetype = mean_linetype))
    }
  }

  # CI as shaded bands (drawn first, behind everything)
  ci_bands <- function(stat) {
    if (is.null(ci) || ci_style != "bands") return(NULL)
    m <- ci_row("M"); s <- ci_row(stat)
    list(annotate("rect", xmin = m$value_lower, xmax = m$value_upper, ymin = 0, ymax = 1,
                  fill = "red", alpha = ci_alpha),
         annotate("rect", xmin = 0, xmax = 1, ymin = s$value_lower, ymax = s$value_upper,
                  fill = "red", alpha = ci_alpha))
  }

  # CI as crossed error bars, centred on the middle of the intervals (drawn on top)
  ci_cross <- function(stat) {
    if (is.null(ci) || ci_style != "cross") return(NULL)
    m <- ci_row("M"); s <- ci_row(stat)
    mid_x <- (m$value_lower + m$value_upper) / 2
    mid_y <- (s$value_lower + s$value_upper) / 2
    bars <- list(geom_errorbar(data = data.frame(x = mid_x, ymin = s$value_lower, ymax = s$value_upper),
                               aes(x = x, ymin = ymin, ymax = ymax), inherit.aes = FALSE,
                               width = ci_cap, linewidth = ci_linewidth, col = ci_color))
    if (m$value_upper > m$value_lower)                      # horizontal bar only if M has an interval
      bars <- c(bars, list(geom_errorbar(data = data.frame(y = mid_y, xmin = m$value_lower, xmax = m$value_upper),
                                         aes(y = y, xmin = xmin, xmax = xmax), inherit.aes = FALSE,
                                         width = ci_cap/3, linewidth = ci_linewidth, col = ci_color,
                                         orientation = "y")))
    bars
  }

  # one plot: statistic y against M, with its bound function and axis label
  make_plot <- function(y, stat, bound_fun, y_label) {
    #mean_y <- mean(y, na.rm = TRUE)
    mean_y <- if (is.null(ci)) mean(y, na.rm = TRUE) else ci_row(stat)$value
    p <- ggplot(data.frame(M = M, y = y), aes(x = M, y = y)) +
      point_layers() +
      geom_line(data = data.frame(M = M_grid, y = bound_fun(K, M_grid)),
                aes(x = M, y = y), inherit.aes = FALSE) +
      ci_bands(stat) +
      mean_lines(mean_y, bound_fun) +
      ci_cross(stat) +
      annotate("point", x = mean_M, y = mean_y, col = "red", size = 2) +
      coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
      xlab(expression(italic(M))) +
      ylab(y_label) +
      theme_bw()
    if (!show_colorbar) p <- p + guides(color = "none")
    p
  }

  # ---- plots ------------------------------------------------------
  plotFST  <- make_plot(FST, "FST", Fup, expression(bolditalic(F[ST])))
  plotGpST <- if (!is.null(GpST)) make_plot(GpST, "GpST", Gpup, expression(bolditalic(G*minute[ST]))) else NULL
  plotD    <- if (!is.null(D))    make_plot(D, "D", Dup, expression(bolditalic(D))) else NULL

  list(plotFST, plotGpST, plotD)
}

#' Plots differentiation statistics normalized by their maximal values at corresponding M
#'
#' @inheritParams ggbounds_raw
#' @returns A list of ggplot objects
#' @export
#'
#' @examples
#' library(tidyverse)
#' freqs_locus1 <- matrix(c(1,0.5,0,0.5),nrow=2)
#' freqs_locus2 <- matrix(c(1,0.8,0,0.2),nrow=2)
#' freqs_locus3 <- matrix(c(1,0.2,0,0.8),nrow=2)
#' data  <- rbind(Diff(freqs_locus1),Diff(freqs_locus2),Diff(freqs_locus2))
#' ggbounds_norm(M=data %>% filter(statistic=="M") %>% pull(value),
#'               FST=data %>% filter(statistic=="FST") %>% pull(value))
ggbounds_norm = function(M,FST,GpST=NULL,D=NULL,K=2, show_colorbar=FALSE){
  MFtmp = dplyr::tibble(M= seq(0.001,1-0.001,0.001),
                        FST= Fup(K,seq(0.001,1-0.001,0.001)) )
  MGptmp = dplyr::tibble(M= seq(0.001,1-0.001,0.001),
                         GpST= Gpup(K,seq(0.001,1-0.001,0.001)) )
  MDtmp = dplyr::tibble(M= seq(0.001,1-0.001,0.001),
                        D= Dup(K,seq(0.001,1-0.001,0.001)) )
  # compute normalised FST
  FST_norm <- FST / sapply(M, function(m) Fup(K, m))
  MF_ST_tib <- dplyr::tibble(M=M,FST_n=FST_norm)

  mean_M <- mean(M, na.rm = TRUE)
  mean_FST_n <- mean(FST_norm, na.rm=T)

  nudge = (mean_FST_n<0.65)*0.16-(mean_FST_n>=0.65)*0.16


  plotFST_norm <-
     ggplot(MF_ST_tib,  aes(x = M, y = FST_n, color = after_stat(density / max(density)))) +
     ggpointdensity::geom_pointdensity(size = .8) +
     labs(color = "relative\ndensity") +
     scale_color_viridis() +
     geom_segment(data = dplyr::tibble(M = mean_M, FST = mean_FST_n),
                  aes(x = M, xend = M, y = 0, yend = 1),
                  inherit.aes = FALSE,
                  col = "red", size = 0.8, linetype = "dashed") +
     geom_point(data=dplyr::tibble(M=mean_M,FST_n=mean_FST_n),
                aes(x = M, y = FST_n),
                inherit.aes = FALSE,
                col="red", pch=16,size=3,stroke=2) +
     geom_label(data=dplyr::tibble(M=mean_M, FST_n=mean_FST_n),
                aes(x=M,y=FST_n,
                label = deparse(bquote(bar(italic(F))[ST] == .(format(mean_FST_n, digits = 2))))
                ),
                parse=TRUE,
                inherit.aes = FALSE,
                nudge_x = 0,
                nudge_y = nudge,
                col="red",
                size=3.5) +
     coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = F) +
     xlab(expression(italic(M))) +
     ylab(expression(italic(F[ST]))) +
     theme_bw()

  if(!is.null(GpST)){
    GST_norm <- GpST / sapply(M, function(m) Gpup(K, m))
    MG_ST_tib <- dplyr::tibble(M=M,GST_n=GST_norm)
    mean_GST_n <- mean(GST_norm, na.rm=T)

    plotGpST_norm <-
       ggplot(MG_ST_tib,  aes(x = M, y = GST_n, color = after_stat(density / max(density)))) +
       ggpointdensity::geom_pointdensity(size = .8) +
       labs(color = "relative\ndensity") +
       scale_color_viridis() +
       geom_segment(data = dplyr::tibble(M = mean_M, GST_n = mean_GST_n),
                    aes(x = M, xend = M, y = 0, yend = 1),
                    inherit.aes = FALSE,
                    col = "red", size = 0.8, linetype = "dashed") +
       geom_point(data=dplyr::tibble(M=mean_M,GST_n=mean_GST_n),
                  aes(x = M, y = GST_n),
                  inherit.aes = FALSE,
                  col="red", pch=16,size=3,stroke=2) +
       geom_label(data=dplyr::tibble(M=mean_M,GST_n=mean_GST_n),
                  aes(x=M,y=GST_n,
                      label = deparse(bquote(bar(italic(G*"'"))[ST] == .(format(mean_FST_n, digits = 2)))),
                  ),
                  parse = TRUE,
                  inherit.aes = FALSE,
                  nudge_x = 0,
                  nudge_y = nudge,
                  col="red",
                  size=3.5) +
       coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = F) +
       xlab(expression(italic(M))) +
       ylab(expression(italic(G[ST]))) +
       theme_bw()
  }
  else{
    plotGpST_norm = NULL
  }

  if(!is.null(D)){
    D_norm <- D / sapply(M, function(m) Dup(K, m))
    MD_tib <- dplyr::tibble(M=M,D_n=D_norm)
    mean_D_n <- mean(D_norm, na.rm=T)

    plotD_norm <-
       ggplot(MD_tib,  aes(x = M, y = D_n, color = after_stat(density / max(density)))) +
       ggpointdensity::geom_pointdensity(size = .8) +
       labs(color = "relative\ndensity") +
       scale_color_viridis(name = "Relative Density") +
       geom_segment(data = dplyr::tibble(M = mean_M, D_n = mean_D_n),
                    aes(x = M, xend = M, y = 0, yend = 1),
                    inherit.aes = FALSE,
                    col = "red", size = 0.8, linetype = "dashed") +
       geom_point(data=dplyr::tibble(M=mean_M,D_n=mean_D_n),
                  aes(x = M, y = D_n),  # add this line
                  inherit.aes = FALSE,
                  col="red", pch=16,size=3,stroke=2) +
       geom_label(data=dplyr::tibble(M=mean_M,D_n=mean_D_n),
                  aes(x=M,y=D_n,
                      label = deparse(bquote(bar(italic(D)) == .(format(mean_FST_n, digits = 2)))),
                  ),
                  parse = TRUE,
                  inherit.aes = FALSE,
                  nudge_x = 0,
                  nudge_y = nudge,
                  col="red",
                  size=3.5) +
       coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = F) +
       xlab(expression(italic(M))) +
       ylab(expression(italic(D))) +
       theme_bw()
  }
  else{
    plotD_norm=NULL
  }

  if(!show_colorbar){
    plotFST_norm <- plotFST_norm + guides(color = "none")
    plotGpST_norm <- plotGpST_norm + guides(color = "none")
    plotD_norm <- plotD_norm + guides(color = "none")
  }

  #if(length(M)>2) plot <- plot +  scale_color_viridis_b()
  return(list(plotFST_norm,plotGpST_norm,plotD_norm) )
}


#' Plot differentiation statistics values and their bounds (first column)
#' and the normalized values of statistics (second column)
#'
#'@inheritParams ggbounds_raw
#'
#' @returns A list of ggplot objects
#' @export
#'
#' @examples
#' library(tidyverse)
#' freqs_locus1 <- matrix(c(1,0.5,0,0.5),nrow=2)
#' freqs_locus2 <- matrix(c(1,0.8,0,0.2),nrow=2)
#' freqs_locus3 <- matrix(c(1,0.2,0,0.8),nrow=2)
#' data  <- rbind(Diff(freqs_locus1),Diff(freqs_locus2),Diff(freqs_locus2))
#' ggbounds_new(M=data %>% filter(statistic=="M") %>% pull(value),
#'              FST=data %>% filter(statistic=="FST") %>% pull(value))
ggbounds_new <- function(M,FST,GpST=NULL,D=NULL,K=2){
  gg1 <- ggbounds_raw(M, FST, GpST, D, K)
  gg2 <- ggbounds_norm(M, FST, GpST, D, K)
  return(list(gg1, gg2))
}
