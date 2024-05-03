sscpp = Module("simple_sim_cpp")

.onUnload <- function(libpath) {
    library.dynam.unload("cspaper", libpath)
}

make_config <- function(num_walks=1000, num_steps=100000, num_threads=4) {
    get_width_times <- unique(round(
        exp(seq(log(10), log(num_steps), length.out=300))
    ))
    new(sscpp$SimulationConfig, num_walks, num_steps, num_threads,
        get_width_times)
}

make_stopping_time_plot <- function(results, alpha, name_map=NULL,
                                    x_label=TRUE, y_label=TRUE) {
    num_steps = results[[1]]$num_steps
    stopping_time_dfs <- lapply(names(results), function(name) {
        times <- results[[name]]$stopping_times
        data.frame(
            name=name,
            stopping_time=ifelse(times > -1, times, num_steps + 1)
        )
    })
    max_false_positive_rate <- max(
        sapply(results, function(r) mean(r$stopping_times > -1))
    )
    stopping_time_data <- do.call(rbind, stopping_time_dfs)
    if (!is.null(name_map)) {
        stopping_time_data$name <- plyr::revalue(stopping_time_data$name,
                                                 name_map)
    }
    label_map <- c(`Hoeffding`='Hoeffding     ',
                   `Beta-Binomial`='Beta\u00adBinomial     ',
                   `Naive SN`='Naive SN     ',
                   `Pointwise Bernoulli`='Pointwise Bernoulli     ')
    stopping_time_data <- mutate(stopping_time_data,
                                 name=plyr::revalue(name, label_map))
    stopping_time_plot <- (
        ggplot(stopping_time_data, aes(stopping_time, color=name))
        + geom_hline(yintercept=alpha, color='gray')
        + stat_ecdf()
        + scale_x_log10(
              if (x_label) 'Number of samples' else NULL,
              breaks=scales::trans_breaks("log10", function(x) 10^x),
              labels=scales::trans_format("log10", scales::math_format(10^.x)))
        + ylab(if (y_label) 'False positive rate' else NULL)
        + scale_color_manual(NULL, values=COLORS)
        + coord_cartesian(
              xlim=c(10, num_steps),
              ylim=c(0, max(1.05 * alpha, max_false_positive_rate)),
              expand=F
          )
        # + right_theme(format='paper')
        + theme(legend.position='top')
    )
    list(data=stopping_time_data, plot=stopping_time_plot)
}

plot_results <- function(results, alpha, title=NULL, use_legend=FALSE,
                         scale_width=FALSE, y_label=TRUE, right_margin_pt=0) {
    stopping_time_results <- make_stopping_time_plot(
        results, alpha, x_label=FALSE, y_label=y_label)
    stopping_time_plot <- (
        stopping_time_results$plot
        + aes(linetype=name)
        + scale_linetype_manual(NULL, values=c(5, 3, 2, 4, 1))
        + guides(color=guide_legend(nrow=2, byrow=TRUE))
        + margins_pt(right=right_margin_pt)
        + theme(legend.key.width=unit(1.7, 'line'),
                legend.key.height=unit(0.2, 'line'))
    )

    num_steps <- max(stopping_time_results$data$stopping_time)
    width_dfs <- lapply(names(results), function(name) {
        avg_widths <- apply(results[[name]]$widths, 2, mean, na.rm=TRUE)
        data.frame(
            name=name,
            t=results[[name]]$evaluate_width_times,
            avg_width=avg_widths
        )
    })
    width_data <- dplyr::filter(do.call(rbind, width_dfs), !is.na(avg_width))

    if (scale_width) {
        max_final_width <- with(dplyr::filter(width_data, t == max(width_data$t)),
                                max(avg_width * sqrt(t)))
        min_width <- with(width_data, min(avg_width * sqrt(t)))
        max_width <- with(width_data, max(avg_width * sqrt(t)))
        y_max <- min(max_width, max(max_final_width, 5 * min_width))
        width_plot <- (
            ggplot(width_data, aes(t, avg_width * sqrt(t), color=name))
            + geom_line()
            + scale_x_log10(
                  bquote(paste('Number of samples, ', italic(t))),
                  breaks=scales::trans_breaks("log10", function(x) 10^x),
                  labels=scales::trans_format("log10", scales::math_format(10^.x)))
            + scale_y_continuous(
                  bquote('Mean CI width' %*% sqrt(italic(t)))
              )
            + scale_color_manual(NULL, values=COLORS)
            # + right_theme(format='paper')
            + coord_cartesian(
                  xlim=c(10, num_steps),
                  ylim=c(0, y_max),
                  expand=F
              )
            + theme(legend.position='none')
        )
    } else {
        max_final_width <- with(dplyr::filter(width_data, t == max(width_data$t)),
                                max(avg_width))
        min_width <- with(width_data, min(avg_width))
        max_width <- with(width_data, max(avg_width))
        y_max <- min(max_width, max(max_final_width, 5 * min_width))
        width_plot <- (
            ggplot(width_data, aes(t, avg_width, color=name, linetype=name))
            + geom_line()
            + scale_x_log10(
                  'Number of samples',
                  breaks=scales::trans_breaks("log10", function(x) 10^x),
                  labels=scales::trans_format("log10", scales::math_format(10^.x)))
            + scale_y_log10(if (y_label) 'CI width' else NULL)
            + annotation_logticks(sides='l')
            + scale_color_manual(NULL, values=COLORS)
            + scale_linetype_manual(NULL, values=c(5, 3, 2, 4, 1))
            # + right_theme(format='paper')
            + coord_cartesian(
                  xlim=c(10, num_steps),
                  ylim=c(0.8 * min_width, 1.1 * max_width),
                  expand=F
              )
            + margins_pt(right=right_margin_pt)
            + theme(legend.position='none')
        )
    }

    legend <- extract_legend(stopping_time_plot)
    stopping_time_plot <- stopping_time_plot + theme(legend.position='none')
    top <- NULL
    if (!is.null(title)) {
        top <- grid::textGrob(title,
                              gp=grid::gpar(fontfamily=GG_FAMILY, fontsize=10))
    }
    plot <- gridExtra::arrangeGrob(stopping_time_plot, width_plot, ncol=1,
                                   right=' ', top=top, heights=c(1, 1.15))
    if (use_legend) {
        plot <- gridExtra::arrangeGrob(plot, legend, ncol=1, heights=c(8, 2))
    }
    gridExtra::grid.arrange(plot)
    return(list(plot=plot, legend=legend))
}

DEFAULT_ALPHA <- 0.05

# all for [-1, 1]-bounded observations

z_test_radius <- function(t, sample_mean, alpha=DEFAULT_ALPHA)
{
    sample_p <- (sample_mean - -1) / 2
    std_dev <- sqrt(t * sample_p * (1 - sample_p)) * 2
    z_factor <- qnorm(1 - alpha / 2)
    z_factor * std_dev
}

pointwise_hoeffding_radius <- function(t, alpha=DEFAULT_ALPHA)
{
    sqrt(2 * t * log(2 / alpha))
}

linear_boundary <- function(t, lambda, alpha=DEFAULT_ALPHA) {
    log(2 / alpha) / lambda + lambda * t / 2
}

two_sided_normal_mixture <- function(t, rho, alpha=DEFAULT_ALPHA) {
    sqrt((t + rho) * log((t + rho) / (alpha * alpha * rho)))
}

make_intro_plot <- function(results, alpha, title=NULL, use_legend=TRUE) {
    set.seed(1984228733)

    stopping_time_results <- make_stopping_time_plot(
        results, alpha,
        name_map=c('Z-test'='Pointwise CLT      ',
                   `Pointwise Hoeffding`='Pointwise Hoeffding      ',
                   `Linear boundary`='Linear bound      ',
                   'Normal mixture'='Curved bound'))
    stopping_time_data <- stopping_time_results$data
    stopping_time_plot <- (
        stopping_time_results$plot
        + aes(linetype=name)
        + scale_linetype_manual(NULL, values=c(3, 4, 2, 1))
        + ylab('Cumulative miscoverage prob.')
        + guides(color=guide_legend(nrow=2, byrow=TRUE),
                 linetype=guide_legend(nrow=2, byrow=TRUE))
        + theme(legend.key.width=unit(1, 'line'),
                legend.key.height=unit(0.2, 'line'))
    )

    num_steps <- max(stopping_time_data$stopping_time)

    random_walk_data <- (
        data.frame(t=1:num_steps,
                   walk=cumsum(2 * rbinom(num_steps, 1, 0.5) - 1))
        %>% filter(t %in% stopping_time_data$stopping_time)
        %>% mutate(
                z_test=z_test_radius(t, walk / t),
                hoeffding=pointwise_hoeffding_radius(t),
                linear=linear_boundary(t, lambda=0.121),
                curved=two_sided_normal_mixture(t, rho=81.535)
            )
        %>% gather('name', 'radius', -t, -walk)
        %>% mutate(mean=walk / t, upper_bound=mean + radius / t,
                   lower_bound=mean - radius / t,
                   name=factor(name, levels=c('z_test', 'hoeffding', 'linear',
                                              'curved')))
    )
    random_walk_plot <- (
        ggplot(random_walk_data, aes(t))
        + geom_line(y=0, color='gray')
        + geom_line(aes(y=mean),
                    filter(random_walk_data, name == random_walk_data$name[1]),
                    linetype='dashed')
        + geom_line(aes(y=upper_bound, color=name, linetype=name))
        + geom_line(aes(y=lower_bound, color=name, linetype=name))
        + text_('Empirical mean', x=1000, y=-0.5, format='paper')
        + arrow_(x=200, y=-0.4, xend=100, yend=0.05, format='paper')
        + scale_x_log10(
              bquote(paste('Number of samples, ', italic(t))),
              breaks=scales::trans_breaks("log10", function(x) 10^x),
              labels=scales::trans_format("log10", scales::math_format(10^.x)))
        + ylab('Confidence bounds')
        + scale_color_manual(NULL, values=COLORS)
        + scale_linetype_manual(NULL, values=c(3, 4, 2, 1))
        # + right_theme(format='paper')
        + coord_cartesian(
              xlim=c(10, num_steps),
              ylim=c(-1, 1),
              expand=F
          )
        + theme(legend.position='none')
    )

    legend <- extract_legend(stopping_time_plot)
    stopping_time_plot <- stopping_time_plot + theme(legend.position='none')
    top <- NULL
    if (!is.null(title)) {
        top <- grid::textGrob(title,
                              gp=grid::gpar(fontfamily=GG_FAMILY, fontsize=10))
    }
    plot <- gridExtra::arrangeGrob(random_walk_plot, stopping_time_plot,
                                   #width_plot,
                                   ncol=2, right=' ', top=top)
    if (use_legend) {
        plot <- gridExtra::arrangeGrob(plot, legend, ncol=1, heights=c(8, 1.1))
    }
    gridExtra::grid.arrange(plot)
    plot
}

save_results <- function(results, alpha, filename) {
    filename <- sprintf('build/%s', filename)
    cat(sprintf('Saving results to %s...\n', filename))
    save(results, alpha, file=filename)
}

load_and_plot <- function(filename, save=FALSE, format='paper', ...) {
    full_filename <- sprintf('build/%s', filename)
    cat(sprintf('Loading results from %s...\n', full_filename))
    envir <- new.env()
    load(full_filename, envir=envir)
    plot <- plot_results(envir$results, envir$alpha, ...)
    if (save) {
        label <- paste0(filename, '.plot')
        height_factor <- if (use_legend) 1.25 else 1
        save_plot(plot$plot, label, format,
                  list(paper=c(5, 1.8 * height_factor),
                       slides=c(5.5, 2.5 * height_factor)))
    }
    return(plot)
}

# which = three_point, ber_50, ber_99
bounded_simulations <- function(which, alpha=0.05, v_opt=500, save=FALSE, ...) {
    config <- make_config(...)

    if (which == 'three_point') {
        x_min <- -1 - 20*.01/.49
        x_max <- 20
        config$discrete_steps(c(x_min, 1, 20), c(.495, .495, .01))
        variance <- .49*x_min^2 + .49*1 + .01*20
    } else if (which == 'ber_50') {
        x_min <- -1
        x_max <- 1
        config$discrete_steps(c(-1, 1), c(.5, .5))
        variance <- 1
    } else if (which == 'ber_01') {
        x_min <- -.01
        x_max <- .99
        config$discrete_steps(c(-.01, .99), c(.99, .01))
        variance <- .01 * .99
    } else stop()

    results = list()

    config$beta_binomial_mixture(abs(x_min) * x_max * v_opt, alpha)
    results[['Beta-Binomial']] <- sscpp$run_simulation(config)

    config$pointwise_bernoulli(alpha)
    results[['Pointwise Bernoulli']] <- sscpp$run_simulation(config)

    config$normal_mixture(v_opt, alpha);
    results[['Hoeffding']] <- sscpp$run_simulation(config)

    config$naive_self_normalized(v_opt * variance, alpha)
    results[['Naive SN']] <- sscpp$run_simulation(config)

    config$empirical_bernstein(v_opt * variance, alpha)
    results[['Empirical Bernstein']] <- sscpp$run_simulation(config)

    if (save) {
        save_results(results, alpha, sprintf('ev_%s_simulations.robj', which))
    } else {
        plot_results(results, alpha)
    }
}

intro_simulations <- function(alpha=0.05, t_opt=500, save=FALSE,
                              num_walks=10000, num_steps=100000, ...) {
    config <- make_config(num_steps=num_steps, num_walks=num_walks, ...)
    config$discrete_steps(c(-1, 1), c(.5, .5))
    results = list()

    config$z_test(alpha);
    results[['Z-test']] <- sscpp$run_simulation(config)

    config$pointwise_hoeffding(alpha)
    results[['Pointwise Hoeffding']] <- sscpp$run_simulation(config)

    config$linear_bound(t_opt, alpha)
    results[['Linear boundary']] <- sscpp$run_simulation(config)

    config$normal_mixture(t_opt, alpha)
    results[['Normal mixture']] <- sscpp$run_simulation(config)

    if (save) {
        save_results(results, alpha, 'intro_simulations.robj')
    } else {
        make_intro_plot(results, alpha)
    }
}

plot_intro <- function(save=FALSE) {
    envir <- new.env()
    load('build/intro_simulations.robj', envir=envir)
    plot <- make_intro_plot(envir$results, envir$alpha)
    if (save) {
        save_plot(plot, 'intro_simulations_plot', 'paper',
                  list(paper=c(5, 2.6)))
    }
}

save_all_plots <- function() {
    p1 <- load_and_plot('ev_ber_50_simulations.robj', title='Bernoulli(0.5)')
    p2 <- load_and_plot('ev_ber_01_simulations.robj', title='Bernoulli(0.01)',
                        y_label=FALSE)
    p3 <- load_and_plot('ev_three_point_simulations.robj', title='Three point',
                        y_label=FALSE)
    final <- gridExtra::arrangeGrob(p1$plot, p2$plot, p3$plot, ncol=3,
                                    widths=c(4.3, 4, 4))
    final <- gridExtra::arrangeGrob(final, p3$legend, ncol=1, heights=c(7, 1))
    save_plot(final, 'final_simulations', 'paper', list(paper=c(5, 3)))
}
