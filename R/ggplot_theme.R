

theme_Job <- ggplot2::theme(
    legend.position = 'bottom',
    legend.direction = 'horizontal',
    text = ggplot2::element_text(size=9, family='Helvetica', face = 'bold'),
    axis.text.x = ggtext::element_markdown(),
    axis.title.x = ggtext::element_textbox_simple(width = NULL, halign = .5),
    axis.title.y = ggtext::element_textbox_simple(size = 8, orientation = 'left-rotated', width = NULL, halign = .5),
    strip.text = ggtext::element_textbox_simple(width = NULL, halign = .5),
    panel.grid.major.x = ggplot2::element_line(colour = 'grey90', linetype = 'dotted'),
    panel.grid.major.y = ggplot2::element_line(colour = '#E5E5E5', linetype = 'dotted'),
    panel.grid.minor.y = ggplot2::element_blank(),
    panel.background = ggplot2::element_rect(fill = NA, colour = 'black'),
    panel.border = ggplot2::element_rect(fill = NA, colour = NA),
    strip.background = ggplot2::element_rect(colour = 'black', fill = 'white'),
    legend.text = ggtext::element_markdown()
)
