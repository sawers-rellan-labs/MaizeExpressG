# plotly figure builders.

base_layout <- function(p, xlab, ylab, ...) {
  layout(p, xaxis = list(title = xlab, zeroline = FALSE),
         yaxis = list(title = ylab, zeroline = FALSE),
         legend = list(orientation = "h", x = 0, y = 1.02, yanchor = "bottom"),
         margin = list(t = 60), hovermode = "closest", ...) |>
    config(displaylogo = FALSE, toImageButtonOptions = list(format = "svg"))
}

ligule_shape <- function(x0 = 0) list(list(
  type = "line", x0 = x0, x1 = x0, y0 = 0, y1 = 1, yref = "paper",
  line = list(dash = "dot", color = "#888", width = 1)))

profile_plot <- function(d, means, xlab, ylab, show_reps = TRUE, err = "se",
                         ligule_x = NULL, title = NULL, namer = short_name) {
  genes <- unique(means$gene); cols <- setNames(gene_cols(length(genes)), genes)
  p <- plot_ly()
  for (g in genes) {
    mg <- means[means$gene == g, ]
    e <- switch(err, se = mg$se, sd = mg$sd, none = NULL)
    p <- add_trace(p, data = mg, x = ~x, y = ~mean, type = "scatter", mode = "lines+markers",
                   name = namer(g), legendgroup = g,
                   line = list(color = cols[[g]], width = 2), marker = list(color = cols[[g]], size = 7),
                   error_y = if (!is.null(e)) list(array = e, color = cols[[g]], thickness = 1, width = 3),
                   text = ~sprintf("%s<br>section %d (x = %s)<br>mean %.3g, n = %d", namer(g), section, x, mean, n),
                   hoverinfo = "text")
    if (show_reps && any(d$gene == g & !is.na(d$rep))) {
      dg <- d[d$gene == g, ]
      p <- add_trace(p, data = dg, x = ~x, y = ~value, type = "scatter", mode = "markers",
                     name = paste(g, "replicates"), legendgroup = g, showlegend = FALSE,
                     marker = list(color = cols[[g]], size = 5, opacity = 0.45),
                     text = ~sprintf("%s<br>%s<br>value %.3g (FPKM %.3g)", namer(g), sample, value, raw),
                     hoverinfo = "text")
    }
  }
  p <- base_layout(p, xlab, ylab, title = title, showlegend = TRUE)    # keep the name with one gene
  if (!is.null(ligule_x)) p <- layout(p, shapes = ligule_shape(ligule_x),
    annotations = list(list(x = ligule_x, y = 0, yref = "paper", yanchor = "bottom",
                            text = " leaf-2 ligule (approx.)", showarrow = FALSE, xanchor = "left",
                            font = list(size = 10, color = "#666"))))
  p
}

fit_plot <- function(d, curves, focus, xlab, ylab, ligule_x = NULL) {
  p <- plot_ly()
  p <- add_trace(p, data = d, x = ~x, y = ~y, type = "scatter", mode = "markers", name = "data",
                 marker = list(color = "#333", size = 6, opacity = 0.6),
                 text = ~sprintf("%s<br>x = %s, y = %.3g", sample, x, y), hoverinfo = "text")
  mods <- unique(curves$model)
  cols <- setNames(gene_cols(length(mods)), mods)
  for (m in mods) {
    cg <- curves[curves$model == m, ]
    if (m == focus && any(!is.na(cg$lwr))) {
      p <- add_ribbons(p, data = cg, x = ~x, ymin = ~lwr, ymax = ~upr, name = paste(m, "CI"),
                       fillcolor = adjustcolor(cols[[m]], 0.18), line = list(width = 0),
                       hoverinfo = "skip", showlegend = FALSE)
    }
    p <- add_lines(p, data = cg, x = ~x, y = ~fit, name = m,
                   line = list(color = cols[[m]], width = if (m == focus) 3.5 else 1.5,
                               dash = if (m == focus) "solid" else "dash"))
  }
  p <- base_layout(p, xlab, ylab)
  if (!is.null(ligule_x)) p <- layout(p, shapes = ligule_shape(ligule_x))
  p
}

heatmap_plot <- function(Z, xlab, row_labels = TRUE, zlab = "z-score") {
  lim <- max(abs(range(Z, na.rm = TRUE)))
  plot_ly(x = colnames(Z), y = rownames(Z), z = Z, type = "heatmap",
          colors = colorRamp(c("#2166AC", "#F7F7F7", "#B2182B")), zmin = -lim, zmax = lim,
          colorbar = list(title = zlab),
          hovertemplate = "%{y}<br>section %{x}<br>%{z:.2f}<extra></extra>") |>
    layout(xaxis = list(title = xlab, type = "category"),
           yaxis = list(title = "", showticklabels = row_labels, autorange = "reversed",
                        tickfont = list(size = 9))) |>
    config(displaylogo = FALSE)
}
