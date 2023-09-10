
plot_auc_by_penalty <- function(data){
  ggplot(data, aes(x = penalty, y = mean)) + 
    geom_point(color = 'blue') + 
    geom_line(color = 'blue') + 
    ylab("Area under the ROC Curve") +
    scale_x_log10(labels = scales::label_number()) +
    theme_bw()
}

