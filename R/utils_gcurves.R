# for plotting growth curves of a 96-well plate

plotplate <- function(df, dfxy, unsmoothed=TRUE, predicted=FALSE, plate, rows, cols, page, scales = "free_y"){
  dffilt <- dplyr::filter(df, plate_name == {{ plate }})
  xyfilt <- if (!is.null(dfxy)){ left_join(dfxy, distinct(dffilt, id, well, plate_name), by = join_by(id)) %>% 
      drop_na()}
  
  ggplot(dffilt, aes(x = hours)) +
    list(
      ggplot2::geom_line(aes(y=OD600_rollmean), color = "blue"), 
      if (unsmoothed) {ggplot2::geom_line(aes(y=OD600), color = "orange", lty = 2)},
      if (predicted) {ggplot2::geom_line(aes(y=predicted), color = "orange")}, 
      if (!is.null(dfxy)) {ggplot2::geom_point(data = xyfilt, aes(x = x, y = y), color = "red", size = 2)},
      ggplot2::labs(x = "Hours", y = "OD600"), 
      ggplot2::scale_x_continuous(breaks = seq(0, 48, 12), labels = seq(0, 48, 12)), 
      ggforce::facet_wrap_paginate(~ well, nrow = rows, ncol = cols, page = page, scales = scales), 
      ggplot2::theme(axis.text = element_text(size = 5))
    )
}
