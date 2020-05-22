generate_sample_covar <- function(max_n, true_covariance) {
    direction <- runif(max_n) < 0.5
    sign <- 2 * rbinom(max_n, 1, 0.5) -1
    x <- cbind(
        ifelse(direction, sqrt(2), 1 / sqrt(2)) * sign,
        ifelse(direction, sqrt(2), -1 / sqrt(2)) * sign
    )

    sample_covar_elements <- apply(x, 1, function(row) row %*% t(row))
    sample_covar_sums <- apply(sample_covar_elements, 1, cumsum)
    sample_covar_sums
}

get_sample_covar <- function(k, sample_covar_sums) {
    matrix(sample_covar_sums[k,], nrow=2) / k
}

make_covar_ellipse_data <- function(covar_matrix) {
    angles <- seq(0, 2*pi, length.out=100)
    circle_path <- rbind(cos(angles), sin(angles))
    eig <- eigen(covar_matrix)
    U <- eig$vectors
    evalues <- round(eig$values, digits=10) # avoid numerical error causing negative eigenvalues
    covar_sqrt <- U %*% diag(sqrt(evalues)) %*% t(U)
    ellipse_path <- covar_sqrt %*% circle_path
    data.frame(x=ellipse_path[1,], y=ellipse_path[2,])
}

make_covar_bound_data <- function(covar_matrix, operator_norm_radius) {
    eigensystem <- eigen(covar_matrix, TRUE)
    U <- eigensystem$vectors
    inner_eigenvalues <- pmax(0, eigensystem$values - operator_norm_radius)
    inner_matrix <- U %*% diag(inner_eigenvalues) %*% t(U)
    outer_matrix <- (
        U %*% diag(eigensystem$values + operator_norm_radius) %*% t(U))
    inner <- make_covar_ellipse_data(inner_matrix)
    outer <- make_covar_ellipse_data(outer_matrix)
    ind <- c(1:nrow(inner), 1)
    polygon_data <- do.call(rbind, lapply(1:nrow(inner), function(i) {
        data.frame(
            x=c(inner$x[ind[i]], outer$x[ind[i]], outer$x[ind[i+1]],
                inner$x[ind[i+1]]),
            y=c(inner$y[ind[i]], outer$y[ind[i]], outer$y[ind[i+1]],
                inner$y[ind[i+1]]),
            group=i
        )
    }))
    polygon_data
}

get_radius <- function(n, Sigma_op, b, alpha=0.05) {
    boundary <- discrete_mixture_bound(
       b * n * Sigma_op, alpha=alpha, mixture_density=lil_density, c=2 * b / 3,
       v_min=b * Sigma_op)
    boundary / n
}

get_step_covar_bound_data <- function(n, sample_covar_sums, true_sigma_op, b) {
    make_covar_bound_data(get_sample_covar(n, sample_covar_sums),
                          get_radius(n, true_sigma_op, b))
}

ellipse_plot <- function(n, max_n=100000, save=FALSE) {
    set.seed(DEFAULT_SEED)

    b <- 4
    true_covariance <- matrix(c(5/4, 3/4, 3/4, 5/4), nrow=2)
    true_sigma_op <- max(eigen(true_covariance)$values)
    sample_covar_sums <- generate_sample_covar(max_n, true_covariance)
    polygon_data <- get_step_covar_bound_data(n, sample_covar_sums,
                                              true_sigma_op, b)
    truth_data <- make_covar_ellipse_data(true_covariance)

    colors <- brewer.pal(3, 'Dark2')
    light_orange <- lighten(colors[2], 1.5)

    time_label <- prettyNum(as.integer(n), big.mark=',')

    plot <- (
        ggplot(polygon_data, aes(x, y))
        + geom_polygon(aes(group=group), fill=light_orange, color=light_orange)
        + geom_path(data=truth_data, color=colors[3], size=1)
        + ggtitle(bquote(paste(italic(t) == .(time_label))))
        + scale_x_continuous('First coordinate', limits=c(-1.7, 1.7),
                             breaks=NULL)
        + scale_y_continuous(
              'Second coordinate',
              limits=c(-1.7, 1.8),
              breaks=NULL
          )
        + report_theme
        + margins_pt(right=20)
    )
    plot <- (
        plot
        + text_(label="paste('True ', Sigma)",
                parse=TRUE, x=0, y=-1.65, color=colors[3])
        + arrow_(color=colors[3],
                 x=0,
                 y=-1.4,
                 size=0.7,
                 xend=truth_data$x[75],
                 yend=truth_data$y[75])
        + text_(label='Confidence set', x=0, y=max(polygon_data$y) + 0.3,
                color=colors[2])
    )
    if (save) {
        save_plot(
            plot,
            sprintf('covariance_%s', n),
            c(1.6, 1.8))
    }
    plot
}

save_ellipse_plots <- function() {
    ellipse_plot(200, save=TRUE)
    ellipse_plot(500, save=TRUE)
    ellipse_plot(2000, save=TRUE)
}
