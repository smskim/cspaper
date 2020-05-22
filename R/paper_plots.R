psi_G <- function(lambda, c_=0) lambda^2 / (2 * (1 - c_ * lambda))

psi_G_inverse <- function(u, c_=0) {
    stopifnot(u > 0)
    upper <- 1.01 * sqrt(2 * u)
    if (c_ > 0 && upper > 1 / c_) {
        upper <- 1 / c_ - 1e-9
    }
    uniroot(function(lambda) psi_G(lambda) - u, c(0, upper))$root
}

stitching_linear_bound <- function(v, k, h_fn, alpha, c_=0, eta=2) {
    r_k <- log(h_fn(k) / alpha)
    lambda_k <- psi_G_inverse(r_k / eta^(k + 1/2), c_=c_)
    (r_k + psi_G(lambda_k, c_=c_) * v) / lambda_k
}

stitching_bound <- function(v, h_fn, alpha, c_=0, eta=2) {
    l_fn <- function(v) log(h_fn(log(v, eta))) + log(1 / alpha)
    k1 <- (eta^(1/4) + eta^(-1/4)) / sqrt(2)
    k2 <- sqrt(eta) + 1
    k1 * sqrt(v * l_fn(v)) + c_ * k2 * l_fn(v)
}

poly_stitching_plot <- function(h_fn=function(k) (k + 1)^1.4 * zeta(1.4),
                           eta=2.9, error_rate=0.05, max_k=2,
                           save=FALSE) {
    plot_v <- seq(1, eta^(max_k + 1/2), length.out=1000)
    linear_data <- do.call(
        rbind,
        lapply(0:max_k, function(k) {
            data.frame(
                group=paste('Linear', k),
                which='1 linear',
                v=plot_v,
                y=stitching_linear_bound(plot_v, k, h_fn, error_rate, eta=eta)
            )
        })
    )
    piecewise_data <- do.call(
        rbind,
        lapply(0:max_k, function(k) {
            v_in_epoch <- plot_v[plot_v >= eta^k & plot_v <= eta^(k+1)]
            data.frame(
                group=paste('Piecewise', k),
                which='2 piecewise',
                v=v_in_epoch,
                y=stitching_linear_bound(v_in_epoch, k, h_fn, error_rate,
                                         eta=eta)
            )
        })
    )
    final_data <- data.frame(
        group='Final',
        which=rep('3 final', length(plot_v)),
        v=plot_v,
        y=stitching_bound(plot_v, h_fn, error_rate, eta=eta)
    )

    colors <- c('black', brewer.pal(3, 'Dark2')[2])

    plot_data <- rbind(linear_data, piecewise_data)
    y_max <- max(final_data$y)
    plot <- (
        ggplot(plot_data, aes(v, y, group=group, color=which, alpha=which))
        + geom_vline(xintercept=eta^(0:max_k), color='gray',
                     linetype='dashed')
        + scale_color_manual(NULL, values=c('black', 'black'))
        + scale_alpha_manual(NULL, values=c(1/5, 1, 1))
        + scale_y_continuous(bquote(paste('Boundary for ', italic(S[t]))),
                             breaks=0)
        + coord_cartesian(
              xlim=c(1, eta^(max_k + 1/2)),
              ylim=c(0, 1.1 * y_max),
              expand=FALSE
          )
        + report_theme
        + theme(panel.grid.minor.x=element_blank(), legend.position='none')
        + scale_x_continuous(
              expression(italic(V[t])),
              breaks=eta^(0:max_k),
              labels=sapply(0:max_k, function(k) bquote(eta^.(k)))
          )
    )
    plot <- (
        plot
        + geom_line(data=final_data, color=colors[2])
        + text_(label='Final boundary', color=colors[2], hjust='right',
                x=final_data$v[920] * 0.9, y=final_data$y[920] * 1.01)
    )
    stitching_arrow <- function(xindex, yindex, curve_index) {
        arrow_(x=final_data$v[xindex], y=final_data$y[yindex],
               xend=piecewise_data$v[curve_index],
               yend=piecewise_data$y[curve_index],
               color=colors[1])
    }
    plot <- (
        plot
        + geom_line(size=0.5)
        + text_(label='Linear uniform\nChernoff bounds', color=colors[1],
                vjust='top', x=final_data$v[400], y=final_data$y[50])
        + stitching_arrow(290, 30, 50)
        + stitching_arrow(400, 60, 370)
        + stitching_arrow(490, 60, 800)
    )
    if (save) {
        save_plot(
            plot,
            'poly_stitching_final',
            c(PAPER_WIDTH, 2))
    }
    plot
}

plot_normalized_boundaries <- function(max_v=1e4, alpha=0.025, y_max=7,
                                       border_size=0.2, save=FALSE) {
    v <- exp(seq(0, log(max_v), length.out=500))
    bound_fns <- list(
        `paste('Polynomial stitching, ', italic(c) == 1, ', ', italic(m) == 100)`=
            partial(poly_stitching_bound, eta=2, v_min=100, c=1),
        `paste('Polynomial stitching, ', italic(c) == 0, ', ', italic(m) == 100)`=
            partial(poly_stitching_bound, eta=2, v_min=100),
        `paste('Discrete mixture LIL, ', italic(m) == 50)`=
            partial(discrete_mixture_bound, mixture_density=lil_density,
                    lambda_max=0.38),
        `paste('Gamma mixture, ', italic(c) == 1, ', ', italic(m) == 300)`=
            partial(gamma_exponential_bound, c=1, rho=best.rho(300, alpha)),
        `paste('Normal mixture, ',  italic(m) == 300)`=
            partial(one_sided_normal_mixture, rho=best.rho(300, alpha)),
        `paste('Gamma mixture, ', italic(c) == 1, ', ', italic(m) == '5,000')`=
            partial(gamma_exponential_bound, c=1, rho=best.rho(5000, alpha)),
        `paste('Normal mixture, ', italic(m) == '5,000')`=
            partial(one_sided_normal_mixture, rho=best.rho(5000, alpha))
    )
    plot_datas <- lapply(1:length(bound_fns), function(index) {
        data.frame(
            index=index,
            type=names(bound_fns)[index],
            v=v,
            y=bound_fns[[index]](v, alpha) / sqrt(v)
        )
    })
    plot_data <- (
        do.call(rbind, plot_datas)
        %>% mutate(type=factor(as.character(type),
                               levels=unique(names(bound_fns))))
    )
    spaces <- do.call(paste0, as.list(rep(' ', 70)))
    min_v <- 5
    plot <- (
        ggplot(plot_data, aes(v, y, color=type, linetype=type))
        + geom_vline(xintercept=c(min_v, max_v), size=border_size)
        + annotate('segment', x=min_v, y=0, xend=max_v, yend=0,
                   size=border_size)
        + annotate('segment', x=min_v, y=y_max, xend=max_v, yend=y_max,
                   size=border_size)
        + geom_line()
        + scale_color_manual(values=COLORS)
        + scale_linetype_manual(values=c(2, 2, 2, 1, 1, 1, 1))
        + scale_x_log10(bquote(paste(italic(v), .(spaces))),
                        breaks=10^(1:4),
                        labels=lapply(1:4, function(x) bquote(10^.(x))))
        + scale_y_continuous(bquote(italic(u(v) / sqrt(v))))
        + coord_cartesian(xlim=c(min_v, 15000 * max_v), ylim=c(0, y_max),
                          expand=FALSE)
        + report_theme
        + theme(panel.grid.major=element_blank(),
                panel.grid.minor=element_blank(),
                panel.border=element_blank(),
                legend.position='none')
    )
    for (index in 1:length(bound_fns)) {
        label <- names(bound_fns)[index]
        label_y <- y_max - index * 0.8 - 0.7
        plot <- (
            plot
            + text_(label=label, parse=TRUE, x=2.2 * max_v, y=label_y,
                    hjust='left')
            + annotate('segment',
                       x=2 * max_v, y=label_y,
                       xend=max_v,
                       yend=tail(plot_data[plot_data$type == label,]$y, 1),
                       color=COLORS[index])
        )
    }
    if (save) {
        save_plot(
            plot,
            'normalized_boundaries',
            c(5, 2))
    }
    plot
}

bernoulli_ate_plot <- function(max_n=1e5, treat_p=1/2, alpha=0.05, save=FALSE,
                               label='Corollary 1') {
    stopifnot(treat_p <= 1/2)
    set.seed(DEFAULT_SEED)

    n <- 1:max_n

    y0 <- rbinom(max_n, 1, 0.5)
    y1 <- ifelse(y0 == 1, 1, rbinom(max_n, 1, 0.2))

    treated <- rbinom(max_n, 1, treat_p)
    Y_obs <- ifelse(treated, y1, y0)

    full_data <- (
        data.frame(n=n, Tn=treated, Y_obs=Y_obs) %>%
        mutate(
            ate=cummean(y1 - y0),
            n_treated=cumsum(Tn),
            n_control=n - n_treated,
            naive_estimate=(
                cumsum(Tn * Y_obs) / n_treated
                - cumsum((1 - Tn) * Y_obs) / n_control
            ),
            Vn_treated_naive=(
                cumsum(Tn * Y_obs^2) - (1/n_treated) * cumsum(Tn * Y_obs)^2
            ),
            Vn_control_naive=(
                cumsum((1-Tn)*Y_obs^2) - (1/n_control) * cumsum((1-Tn)*Y_obs)^2
            ),
            diff_var_naive=(
                Vn_treated_naive / (n_treated - 1) / n_treated
                + Vn_control_naive / (n_control - 1) / n_control
            ),
            naive_bound=naive_estimate + qnorm(1-alpha/2) * sqrt(diff_var_naive)
        )
        %>% mutate(naive_bound=ifelse(is.finite(naive_bound), naive_bound, 1))
    )

    plot_n <- unique(round(exp(seq(0, log(max_n), length.out=300))))
    data <- full_data[plot_n,]
    density <- get_half_normal_density(50, alpha)
    generic_bound_data <- generic_ate_bound(Y_obs, treated, treat_p, alpha,
                                            100, plot_n)

    data$generic_bound <- generic_bound_data$upper
    data$point_estimate <- generic_bound_data$point_estimate
    with(data, cat(sprintf('Naive min UCB/truth: %.2f, uniform: %.2f\n',
                           min(naive_bound / ate), min(generic_bound / ate))))
    with(data, cat(sprintf('Naive failed n: %d/%d, uniform: %d\n',
                           sum(naive_bound < ate),
                           nrow(data),
                           sum(generic_bound < ate))))

    plot_data <- (
        data
        %>% gather('type', 'bound', naive_bound, generic_bound, point_estimate)
        %>% select(n, bound, type)
    )
    plot_data <- rbind(
        plot_data,
        data.frame(n=data$n, type='ate', bound=data$ate)
    )
    plot_data <- mutate(
        plot_data,
        type=factor(as.character(type),
                    levels=c('ate', 'point_estimate', 'naive_bound',
                             'generic_bound'))
    )

    colors <- c('gray', 'black', 'black', COLORS[2], COLORS[1])
    x_axis <- scale_x_log10(
        bquote(paste(italic(t), ' (log scale)')),
        breaks=10^(1:6),
        labels=scales::trans_format("log10", scales::math_format(10^.x))
    )
    ci_plot <- (
        ggplot(plot_data, aes(n, bound, color=type, linetype=type))
        + geom_line(size=0.7)
        + scale_color_manual(NULL, values=colors)
        + scale_linetype_manual(NULL, values=c(1, 1, 3, rep(1, 8)))
        + x_axis
        + scale_y_continuous(bquote(paste('UCB for ', ATE[italic(t)])),
                             breaks=c(0, 0.1, 0.2, 0.3))
        + coord_cartesian(xlim=c(85, max_n), ylim=c(-0.05, 0.35), expand=FALSE)
        + report_theme
        + theme(legend.position='none')
        + margins_pt(right=10)
        #+ text_(label='ATE[italic(t)]', parse=T, x=data$n[220],
        #        y=data$ate[220]-.02, color=colors[1], vjust='top',
        #        format=format)
    )
    ratio_data <- (
        data
        %>% mutate(ratio=((generic_bound - point_estimate) /
                          (naive_bound - naive_estimate)))
    )
    ratio_plot <- (
        ggplot(ratio_data, aes(n, ratio))
        + geom_hline(yintercept=1, color=colors[2], linetype='dotted')
        + geom_line(color = colors[4])
        + x_axis
        + scale_y_continuous(
              'Ratio of UCB radius to CLT',
              breaks=0:3
          )
        + coord_cartesian(xlim=c(100, max_n), ylim=c(0,3), expand=FALSE)
        + report_theme
        + theme(legend.position='none')
    )

    ci_plot <- ci_plot + margins_pt(right=10)
    ratio_plot <- ratio_plot + margins_pt(left=10, right=10)

    plot <- arrangeGrob(ci_plot, ratio_plot, ncol=2)
    grid.arrange(plot)

    if (save) {
        save_plot(
            plot,
            'bernoulli_ate_ci_uniform',
            c(PAPER_WIDTH, 2.2))
    }
    plot
}

discrete_mixture_plot <- function(
        s=30, n=400, error_rate=0.05, psi=psi_hoeffding,
        mixture_density=lil_density, gamma=1.3, lambda_max=1/2, save=FALSE) {
    k <- 0:round(log(10 * max(n), gamma))
    lambda_k <- lambda_max / gamma^(k + 1/2)
    weight_k <- (
        mixture_density(lambda_max / gamma^k, lambda_max) * lambda_max * (gamma - 1) / gamma^(k+1)
    )

    x_grid <- seq(0, lambda_max, length.out=400)
    density_data <- data.frame(
        x=x_grid,
        y=mixture_density(x_grid, lambda_max)
    )
    box_data <- data.frame(
        x=lambda_max / 2 * (gamma + 1) / gamma^(k+1),
        y=mixture_density(lambda_k * sqrt(gamma), lambda_max)
    )
    integrand <- function(lambda) exp(lambda * s - n * psi(lambda))
    integrand_data <- data.frame(
        x=x_grid,
        y=integrand(x_grid)
    )
    point_mass_data <- data.frame(
        x=lambda_k,
        y=integrand(lambda_k)
    )
    colors <- brewer.pal(8, 'Dark2')
    plot <- (
        ggplot(data=density_data, aes(x, y))
        + geom_col(
              data=box_data,
              width=lambda_max * (gamma - 1) / gamma^(k+1),
              color=colors[1],
              fill=NA,
          )
        + geom_col(
              data=box_data[6,],
              width=lambda_max * (gamma - 1) / gamma^(5+1),
              color=colors[1],
              fill=colors[1],
              alpha=0.5
          )
        + geom_line(data=density_data)
        + geom_line(data=integrand_data, color=colors[2])
        + geom_segment(aes(x=x, y=0, xend=x, yend=y), data=point_mass_data, color=colors[3])
        + geom_point(data=point_mass_data, color=colors[3])
        + scale_x_continuous(
              bquote(lambda),
              breaks=c(0, lambda_k[1:30], lambda_max),
              labels=c(0, lapply(0:9, function(x) bquote(lambda[.(x)])),
                       bquote(cdots * phantom(0)), rep(NA, 19), bquote(lambda[max]))
          )
        + scale_y_continuous(NULL, breaks=0)
        + coord_cartesian(xlim=c(0, 1/2), ylim=c(0, 4), expand=F)
        + report_theme
        + margins_pt(right=10)
        + text_(label='italic(f(lambda))', parse=TRUE,
                x=density_data$x[250], y=density_data$y[250] + 0.1,
                hjust='left', vjust='bottom')
        + text_(label='exp(italic(lambda * s - psi(lambda) * v))', parse=TRUE,
                x=integrand_data$x[60], y=integrand_data$y[60] + 0.2,
                hjust='left', vjust='bottom', color=colors[2])
        + text_(label='italic(w)[5]', parse=TRUE, x=point_mass_data$x[6] * 1.3, y=1.4,
                hjust='left', vjust='bottom', color=colors[1])
        + arrow_(x=box_data$x[6] * 1.29, y=1.39, xend=box_data$x[6] * 1.05, yend=box_data$y[6],
                 color=colors[1])
        + text_(label='exp(italic(lambda[5] * s - psi(lambda[5]) * v))', parse=TRUE,
                x=point_mass_data$x[6] * 1.1, y=point_mass_data$y[6] + 0.5,
                hjust='left', color=colors[3])
        + arrow_(x=point_mass_data$x[6] * 1.09, y=point_mass_data$y[6] + 0.4,
                 xend=point_mass_data$x[6], yend=point_mass_data$y[6],
                 color=colors[3])
    )
    if (save) {
        save_plot(
            plot,
            'discrete_mixture',
            c(PAPER_WIDTH, 2.5))
    }
    plot
}

make_finite_lil_plot_simplified <- function(
        max_n=1e5, error_rate=0.025, A=1.7, save=FALSE, t_opt=10000) {
    cat(sprintf('Normal mixture rho = %.2f\n', best.rho(t_opt, error_rate)))
    bound_fns <- list(
        `Polynomial stitching (ours)`=partial(poly_stitching_bound, A=A,
                                              v_min=t_opt),
        # .82 = fraction of alpha spent by polynomial stitching in 42 epochs
        # with s = 1.4
        # sum((1:42)^(-1.4)) / zeta(1.4)
        `Inverted stitching (ours)`=partial(inverted_stitching_lil_bound, A=A,
                                            alpha_fraction=0.82, v_min=t_opt),
        `Discrete mixture (ours)`=partial(
            discrete_mixture_bound, mixture_density=lil_density,
            lambda_max=0.04444
        ),
        `Normal mixture`=partial(
            one_sided_normal_mixture,
            rho=best.rho(t_opt, error_rate)
        )
    )
    colors <- brewer.pal(8, 'Dark2')
    colors <- c(colors, colors[1:2])
    linetypes <- rep('solid', 10)

    plot_n <- unique(round(exp(seq(0, log(max_n), length.out=1000))))
    plot_datas <- lapply(1:length(bound_fns), function(index) {
        data.frame(
            index=index,
            type=names(bound_fns)[index],
            n=plot_n,
            y=bound_fns[[index]](plot_n, error_rate)
        )
    })
    plot_data <- (
        do.call(rbind, plot_datas)
        %>% mutate(type=factor(as.character(type),
                               levels=unique(names(bound_fns))))
    )
    hoeffding_data <- data.frame(
        index=-1, n=plot_n, y=sqrt(2 * plot_n * log(1/error_rate))
    )
    clt_data <- data.frame(
        index=-1, n=plot_n, y=qnorm(1-error_rate) * sqrt(plot_n)
    )
    lil_data <- data.frame(
        index=-1, n=plot_n, y=sqrt(2 * plot_n * log(pmax(1, log(plot_n))))
    )

    y_max <- 4 * sqrt(max_n)
    border_size <- 0.2
    spaces <- do.call(paste0, as.list(rep(' ', 60)))
    plot <- (
        ggplot(plot_data, aes(n, y, group=index))
        + geom_vline(xintercept=c(0, max_n), size=border_size)
        + geom_segment(x=0, y=0, xend=max_n, yend=0, size=border_size)
        + geom_segment(x=0, y=y_max, xend=max_n, yend=y_max, size=border_size)
        + geom_line(data=clt_data, linetype='dotted')
        + geom_line(data=hoeffding_data, linetype='dotted')
        + geom_line(aes(color=type, linetype=type), size=0.7)
        + scale_color_manual(NULL, values=colors)
        + scale_linetype_manual(NULL, values=linetypes)
        + scale_x_continuous(bquote(paste(italic(V[t]), .(spaces))),
                             breaks=c(0, 10^5),
                             labels=c(0, bquote(10^5)))
        + scale_y_continuous('Boundary', breaks=c(0, 2000), labels=c(0, '2,000'))
        + coord_cartesian(xlim=c(0, 1.77 * max_n), ylim=c(0, y_max),
                          expand=FALSE)
        + report_theme
        + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
                panel.border=element_blank(), legend.position='none',
                legend.key.width=unit(2.5, 'line'))
        + text_(label='Hoeffding bound', x=max_n * .55, y=290, hjust='left')
        + arrow_(x=max_n * .7, y=360, xend=max_n * .65, yend=660)
        + text_(label='CLT bound', x=max_n * .3, y=190)
    )
    for (index in 1:length(bound_fns)) {
        label <- names(bound_fns)[index]
        label_y <- y_max - index * 150 + 100
        plot <- (
            plot
            + text_(label=label, x=1.1 * max_n, y=label_y, hjust='left')
            + annotate('segment',
                       x=1.09 * max_n, y=label_y,
                       xend=1.0 * max_n, yend=max(plot_data[plot_data$type == label,]$y),
                       color=colors[index])
        )
    }

    if (save) {
        save_plot(
            plot,
            'finite_lil_plot_simplified',
            c(PAPER_WIDTH, 2))
    }
    plot
}

make_finite_lil_plot <- function(
        max_n=1e5, error_rate=0.025, A=1.7, save=FALSE) {
    bound_fns <- list(
        `Jamieson et al. (2013)`=partial(jamieson_bound, A=A),
        `Balsubramani (2014)`=balsubramani_bound,
        `Zhao et al. (2016)`=partial(zhao_bound, A=A),
        `Darling & Robbins (1967b)`=partial(dr_ili_bound, A=A),
        `Kaufmann et al. (2016)`=partial(kaufmann_exact_bound, A=A),
        `Darling & Robbins (1968)`=partial(dr_mixture_lil_bound, A=A),
        `Polynomial stitching (ours)`=partial(poly_stitching_bound, A=A),
        # .82 = fraction of alpha spent by polynomial stitching in 42 epochs
        # with s = 1.4
        # sum((1:42)^(-1.4)) / zeta(1.4)
        `Inverted stitching (ours)`=partial(inverted_stitching_lil_bound, A=A,
                                            alpha_fraction=0.82),
        `Discrete mixture (ours)`=partial(
            discrete_mixture_bound, mixture_density=lil_density, lambda_max=4
        )
    )
    colors <- brewer.pal(8, 'Dark2')
    colors <- c(colors, colors[1:2])
    linetypes <- rep('solid', 10)

    plot_n <- unique(round(exp(seq(0, log(max_n), length.out=1000))))
    plot_datas <- lapply(1:length(bound_fns), function(index) {
        data.frame(
            index=index,
            type=names(bound_fns)[index],
            n=plot_n,
            y=bound_fns[[index]](plot_n, error_rate)
        )
    })
    plot_data <- (
        do.call(rbind, plot_datas)
        %>% mutate(type=factor(as.character(type),
                               levels=unique(names(bound_fns))))
    )
    lil_data <- data.frame(
        index=-1, n=plot_n, y=sqrt(2 * plot_n * log(pmax(1, log(plot_n))))
    )

    y_max <- 6.5 * sqrt(max_n)
    border_size <- 0.2
    spaces <- do.call(paste0, as.list(rep(' ', 60)))
    plot <- (
        ggplot(plot_data, aes(n, y, group=index))
        + geom_vline(xintercept=c(0, max_n), size=border_size)
        + geom_segment(x=0, y=0, xend=max_n, yend=0, size=border_size)
        + geom_segment(x=0, y=y_max, xend=max_n, yend=y_max, size=border_size)
        + geom_line(aes(color=type, linetype=type), size=0.7)
        + scale_color_manual(NULL, values=colors)
        + scale_linetype_manual(NULL, values=linetypes)
        + scale_x_continuous(bquote(paste(italic(V[t]), .(spaces))),
                             breaks=c(0, 10^5),
                             labels=c(0, bquote(10^5)))
        + scale_y_continuous('Boundary', breaks=c(0, 2000), labels=c(0, '2,000'))
        + coord_cartesian(xlim=c(0, 1.55 * max_n), ylim=c(0, y_max),
                          expand=FALSE)
        + report_theme
        + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
                panel.border=element_blank(), legend.position='none',
                legend.key.width=unit(2.5, 'line'))
    )
    for (index in 1:length(bound_fns)) {
        label <- names(bound_fns)[index]
        label_y <- y_max - index * 150 + 110
        plot <- (
            plot
            + text_(label=label, x=1.1 * max_n, y=label_y, hjust='left')
            + annotate('segment',
                       x=1.09 * max_n, y=label_y,
                       xend=1.0 * max_n, yend=max(plot_data[plot_data$type == label,]$y),
                       color=colors[index])
        )
    }

    if (save) {
        save_plot(
            plot,
            'finite_lil_plot',
            c(PAPER_WIDTH, 3))
    }
    plot
}
