# Shared manuscript style. Change colours here, then rerun Figure4.R.
# No analysis or numerical data are changed by this file.
paper_style <- list(
  font = 'Arial',
  text = '#20242B',
  secondary = '#5A626D',
  condition_order = c('LN', 'HN', 'LU', 'HU', 'LA', 'HA'),
  condition = c(LN='#8FBFDD', HN='#2F73AE', LU='#79C8B4',
                HU='#1B927F', LA='#D7A1AD', HA='#A84F65'),
  concentration_shape = c(L=21, H=24),
  colony = c('1'='#E3E6EA', '2'='#B7BDC6', '4'='#778392', '8'='#344050'),
  width_mm = 183,
  height_mm = 164,
  font_size = 7.2,
  axis_text_size = 6.8,
  title_size = 8.1,
  panel_label_size = 9,
  dpi = 600
)

theme_paper <- function(scale = 1) {
  ggplot2::theme_classic(base_size=paper_style$font_size*scale,
                       base_family=paper_style$font) +
    ggplot2::theme(
      axis.line=ggplot2::element_line(linewidth=.28, colour=paper_style$text),
      axis.ticks=ggplot2::element_line(linewidth=.28, colour=paper_style$text),
      axis.ticks.length=grid::unit(1.1, 'mm'),
      axis.text=ggplot2::element_text(size=paper_style$axis_text_size*scale,
                                    colour=paper_style$text),
      axis.title=ggplot2::element_text(size=paper_style$font_size*scale,
                                     colour=paper_style$text),
      axis.title.y=ggplot2::element_text(margin=ggplot2::margin(r=4)),
      axis.text.x=ggplot2::element_text(margin=ggplot2::margin(t=3)),
      plot.title=ggplot2::element_text(size=paper_style$title_size*scale,
                    face='bold', hjust=0, margin=ggplot2::margin(b=3)),
      plot.subtitle=ggplot2::element_text(size=6.3*scale,
                    colour=paper_style$secondary, margin=ggplot2::margin(b=5)),
      plot.tag=ggplot2::element_text(size=paper_style$panel_label_size*scale,
                                    face='bold', colour=paper_style$text),
      plot.tag.position=c(0, 1),
      plot.caption=ggplot2::element_text(size=6.2*scale,
                    colour=paper_style$secondary, hjust=0, lineheight=1.15),
      legend.title=ggplot2::element_text(size=6.5*scale, face='plain'),
      legend.text=ggplot2::element_text(size=6.5*scale),
      legend.key.size=grid::unit(3.2, 'mm'),
      legend.spacing.x=grid::unit(1.4, 'mm'),
      legend.margin=ggplot2::margin(1,0,0,0),
      panel.grid=ggplot2::element_blank(),
      panel.background=ggplot2::element_rect(fill='white', colour=NA),
      plot.background=ggplot2::element_rect(fill='white', colour=NA),
      plot.margin=ggplot2::margin(5, 7, 5, 6),
      legend.position='bottom'
    )
}
