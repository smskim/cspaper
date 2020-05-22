report_abc <- function(label, A, B, C) {
    cat(sprintf('%s: A=%.3f, B=%.3f, C=%.3f\n', label, A, B, C))
}

zhao_bound <- function(v, error_rate, A, a=NULL, c=1.1) {
    if (!is.null(a)) {
        A <- 2*sqrt(a)
    } else {
        a <- A^2 / 4
    }
    stopifnot(a > c/2)
    B <- c
    C <- c / (2 * a) * log(zeta(2 * a / c) / (error_rate * log(c)^(2 * a / c)))
    report_abc('Zhao', A, B, C)
    A * sqrt(v * (log(log(B*v)) + C))
}

jamieson_bound <- function(v, error_rate, A, epsilon=NULL) {
    if (!is.null(epsilon)) {
        A <- (1 + sqrt(epsilon)) * sqrt(2 * (1 + epsilon))
    } else {
        epsilon <- uniroot(function(e) (1 + sqrt(e)) * sqrt(2*(1 + e)) - A, c(0, 100))$root
    }
    B <- 1 + epsilon
    C <- (
        1 / (1 + epsilon) *
        log((2 + epsilon) / (error_rate * epsilon * log(1 + epsilon)^(1 + epsilon)))
    )
    report_abc('Jamieson', A, B, C)
    A * sqrt(v * (log(log(B*v)) + C))
}

balsubramani_bound <- function(v, error_rate) {
    v0 <- 173 * log(2/error_rate)
    A <- sqrt(6)
    B <- 173 * 3 / 2 * log(2 / error_rate)
    C <- log(1 / error_rate) / 2
    report_abc('Balsubramani', A, B, C)
    ifelse(v < 1, NA, pmax(1, A * sqrt(v * (log(log(B*v)) + C))))
    #sapply(v, function(vv) uniroot(function(m) m - A * sqrt(vv * (log(log(B*vv / m)) + C)), c(sqrt(vv), 40*sqrt(vv)))$root)
}

kaufmann_exact_bound <- function(v, error_rate, A, eta=NULL) {
    if (is.null(eta)) {
        eta <- A^2 / 2
    } else {
        A <- sqrt(2 * eta)
    }
    error_prob <- function(x) {
        exp(1/2) * zeta(eta * (1 - 1/(2*x))) * (sqrt(x / 8) + 1)^eta * exp(-x)
    }
    optimal_x <- uniroot(function(x) error_prob(x) - error_rate, c(8 / (exp(1) - 1)^2, 20))$root
    A <- sqrt(2 * eta)
    B <- exp(1)
    C <- optimal_x / eta
    report_abc('Kaufmann exact', A, B, C)
    A * sqrt(v * (log(log(B*v)) + C))
}

best_nu <- function(t, alpha) {
    objective <- function(nu) {
        gamma <- 1 + sqrt(1 + t^2 / nu) / (2 * alpha)
        2 * nu * gamma * log(gamma) / (gamma - 1) - t^2
    }
    uniroot(objective, c(t^2 / 100, t^2 * 100))$root
}

exact_normal_mixture_z <- function(t, error_rate, t0) {
    nu = best_nu(t0, error_rate)
    supermg <- function(s, v) {
        sqrt(4*nu / (v + nu)) * exp(s^2 / (2 * (v + nu))) * pnorm(s / sqrt(v + nu))
    }
    sapply(t^2, function(v) {
        uniroot(function(s) supermg(s, v) - 1/error_rate, c(0, 1), extendInt='upX')$root / sqrt(v)
    })
}

# Darling & Robbins (1967) Iterated Logarithm Inequalities
dr_ili_bound <- function(v, error_rate, A, eta=NULL, c=1.4, j=1) {
    stopifnot(c > 1)
    if (is.null(eta)) {
        eta <- uniroot(function(e) (1 + e) * sqrt(c / (2 * e)) - A, c(1, 100))$root
    } else {
        A <- (1 + eta) * sqrt(c / (2 * eta))
    }
    B <- eta^j
    C <- 1 / c * log(1 / (error_rate * (c - 1) * (j - 1/2)^(c - 1) * log(eta)^c))
    report_abc('D&R 1967', A, B, C)
    A * sqrt(v * (log(log(B*v)) + C))
}

# Darling & Robbins (1968) Some Further Remarks on Inequalities for Sample Sums
upper_incomplete_gamma <- function(s, x) {
    pgamma(x, s, lower=FALSE) * gamma(s)
}

dr_mixture_lil_bound <- function(v, error_rate, A, v0=3, t0=sqrt(1/3), optimistic_factor=FALSE)
{
    m <- 3
    B <- 3
    error_prob <- function(C) {
        (
            exp(1/2) / 2 * A * exp(-C) * (2 / (A^2 - 2))^(3/2) *
            upper_incomplete_gamma(3/2, (A^2 - 2) / 2 * (log(log(m)) + C))
        ) / ifelse(optimistic_factor, 2, 1)
    }
    C <- uniroot(function(C) error_prob(C) - error_rate, c(0.1, 100))$root
    report_abc('D&R 1968', A, B, C)
    A * sqrt(v * (log(log(B*v)) + C))
}

# Stitching

poly_stitching_bound <- function(v, error_rate, A=NULL, eta=NULL, s=1.4, c=0,
                                     v_min=1) {
    if (is.null(eta)) {
        eta <- uniroot(function(e) (e^.25 + e^(-.25)) * sqrt(s/2) - A,
                       c(1, 100))$root
    }
    k1 <- (eta^{1/4} + eta^{-1/4}) / sqrt(2)
    k2 <- (sqrt(eta) + 1) / 2
    l_v <- (s * log(log(eta * pmax(v, v_min) / v_min))
            + log(zeta(s) / (error_rate * log(eta)^s)))
    second_term <- k2 * c * l_v
    sqrt(k1^2 * pmax(v, v_min) * l_v + second_term^2) + second_term
}

inverted_stitching_lil_p_bound <- function(A, C, t_horizon) {
    bound_for_eta <- function(eta) {
        max_k <- ceiling(log(t_horizon, eta))
        xk <- eta^(0:max_k)
        yk <- A * sqrt(xk * (log(log(exp(1)*xk)) + C))
        ykp1 <- yk[-1]
        yk <- head(yk, -1)
        xk <- head(xk, -1)
        exponents <- 2*(ykp1 - yk) * (eta*yk - ykp1) / (xk * (eta - 1)^2)
        stopifnot(all(exponents > 0))
        sum(exp(-exponents))
    }
    result <- optimize(bound_for_eta, c(1.001, 10))
    stopifnot(result$minimum < 9.9)
    cat(sprintf('A=%.3f C=%.3f best eta=%.3f\n', A, C, result$minimum))
    result$objective
}

get_lil_beta <- function(A, t_horizon, error_rate=0.05) {
    uniroot(
        function(C) {
            sapply(C, function(c) {
                inverted_stitching_lil_p_bound(A, c, t_horizon) - error_rate
            })
        },
        c(1, 10)
    )$root
}

inverted_stitching_lil_bound <- function(v, error_rate, A, t_horizon=1e20,
                                         alpha_fraction=1, v_min=1) {
    error_rate <- error_rate * alpha_fraction
    B <- exp(1)
    C <- get_lil_beta(A, t_horizon, error_rate=error_rate)
    report_abc('Inverted stitching', A, B, C)
    v_effective <- pmax(1, v / v_min)
    A * sqrt(v_min * v_effective * (log(log(B*v_effective)) + C))
}

best.rho <- function(v, alpha) {
    stopifnot(0 < alpha && alpha < 1)
    v / (log(1 / (4 * alpha^2)) + log(log(exp(1) / (4 * alpha^2))))
}

two_sided_best_rho <- function(v, alpha) {
    stopifnot(0 < alpha && alpha < 1)
    v / (2 * log(1 / alpha) + log(log(exp(1) / alpha^2)))
}

# assumes potential outcomes are [0,1]-valued
generic_ate_bound <- function(y_obs, treated, treat_p, alpha, v_opt, indices) {
    stopifnot(all(0 <= y_obs && y_obs <= 1))
    y1_mean <- cumsum(y_obs * treated) / pmax(1, cumsum(treated))
    y0_mean <- cumsum(y_obs * !treated) / pmax(1, cumsum(!treated))
    y1_prediction <- c(1/2, head(y1_mean, -1))
    y0_prediction <- c(1/2, head(y0_mean, -1))
    prediction <- ifelse(treated, y1_prediction, y0_prediction)
    Xhat <- y1_prediction - y0_prediction
    X <- (
        Xhat
        + (treated - treat_p) / (treat_p * (1 - treat_p)) * (y_obs - prediction)
    )
    p_min <- if (treat_p <= 0.5) treat_p else 1 - treat_p
    St <- cumsum(X)[indices]
    mean <- St / indices
    Vt <- cumsum((X - Xhat)^2)[indices]
    confidence_radius <- (
        gamma_exponential_bound(Vt, alpha/2, two_sided_best_rho(v_opt, alpha),
                                c=2 / min(treat_p, 1 - treat_p))
        / indices)
    data.frame(n=indices, point_estimate=mean, upper=mean + confidence_radius,
               lower=mean - confidence_radius)
}

######################################################################

generic_discrete_mixture_lcb <- function(
    Sn, n, error_rate, phi, mixture_density, eta=1.1, mu_lower_bound=NA,
    max_Vn=NA, v_min=1/10
) {
    stopifnot(length(Sn) == length(n))
    stopifnot(min(n) > 0)
    if (is.na(max_Vn)) {
        max_Vn <- max(n)
    }
    lambda_max = 1 / (1 + sqrt(v_min / (2 * log(1 / error_rate))))
    k_max = ceiling(log(
        lambda_max * (1 + sqrt(5 * max_Vn / log(1 / error_rate))),
        eta
    ))
    k <- 0:k_max
    lambda_k <- lambda_max / eta^(k + 1/2)
    weight_k <- (
        mixture_density(lambda_k * sqrt(eta), lambda_max)
        * lambda_max * (eta - 1) / eta^(k+1)
    )
    stopifnot(sum(weight_k) <= 1)

    mixture <- function(mu, index) {
        min(
            1/error_rate + 1,
            sum(
                weight_k
                * exp(lambda_k * Sn[index]
                      - n[index] * phi(lambda_k, mu, index))
            )
        )
    }
    get_bound <- function(index) {
        mean <- Sn[index] / n[index]
        use_lower_bound <- (
            !is.na(mu_lower_bound)
            && mixture(mu_lower_bound, index) < 1 / error_rate
        )
        if (use_lower_bound) {
            mu_lower_bound
        } else {
            uniroot(
                function(mu) mixture(mu, index) - 1 / error_rate,
                c(mean - 1e-6, mean),
                extendInt='downX',
                tol=1e-6
            )$root
        }
    }
    sapply(1:length(Sn), get_bound)
}

generic_discrete_mixture_ucb <- function(
    Sn, n, error_rate, phi, mixture_density, ...
) {
    1 - generic_discrete_mixture_lcb(n - Sn, n, error_rate, phi,
                                     mixture_density, ...)
}

# radius of bound on Sn
discrete_mixture_radius <- function(Vn, error_rate, psi, mixture_density, ...) {
    phi <- function(lambda, mu, n) psi(lambda) * Vn[n] / n
    n <- 1:length(Vn)
    n * generic_discrete_mixture_ucb(0, n, error_rate, phi, mixture_density,
                                     max_Vn=max(Vn), ...)
}

# legacy interface
discrete_mixture_z <- function(t, error_rate, psi, mixture_density, ...) {
    discrete_mixture_radius(t^2, error_rate, psi, mixture_density, ...) / t
}

psi_hoeffding <- function(lambda) lambda^2/2

get_half_normal_density <- function(v_opt, error_rate) {
    rho <- best_nu(sqrt(v_opt), error_rate)
    function(lambda, lambda_max) {
        total_mass <- 2 * (1 - pnorm(lambda_max, sd=sqrt(1 / rho), lower=FALSE))
        2 * dnorm(lambda, sd=sqrt(1 / rho)) / total_mass
    }
}

lil_density <- function(l, lambda_max, s=1.4) {
    stopifnot(max(l) <= lambda_max)
    (s - 1) / (l * log(lambda_max * exp(1) / l)^s)
}

# compare discrete mixture to normal mixture
normal_mixture_comparison <- function(opt_n=100, max_n=1e6, error_rate=0.05,
                                      check_all_n=FALSE) {
    cat(sprintf('rho = %.3f\n', best_nu(sqrt(opt_n), error_rate)))
    if (check_all_n) {
        n <- 1:max_n
    } else {
        n <- unique(round(exp(seq(log(100), log(1e6), length.out=1e4))))
    }
    exact_z <- exact_normal_mixture_z(sqrt(n), error_rate, sqrt(opt_n))
    density <- get_half_normal_density(opt_n, error_rate)
    discrete_z <- discrete_mixture_z(sqrt(n), error_rate, psi_hoeffding, density,
                                     lambda_max=100, eta=1.05)
    cat('Exact z-values:\n')
    print(summary(exact_z))
    cat('Approx z-values:\n')
    print(summary(discrete_z))
    cat('Approx / exact:\n')
    print(summary(discrete_z / exact_z), digits=5)
    if (!check_all_n) {
        qplot(n, discrete_z / exact_z) + scale_x_log10()
    }
}

lil_density <- function(l, lambda_max, s=1.4) {
    stopifnot(max(l) <= lambda_max)
    (s - 1) / (l * log(lambda_max * exp(1) / l)^s)
}

discrete_mixture_bound_nonvec <- function(
        v, alpha, mixture_density, c=0, v_min=1, eta=1.1, lambda_max=NULL
) {
    psi <- function(lambda) lambda^2 / (2 * (1 - c * lambda))
    if (is.null(lambda_max)) {
        lambda_max <- 1 / (c + sqrt(v_min / (2 * log(1/alpha))))
        cat(sprintf('Discrete mixture using lambda_max = %.3f\n', lambda_max))
    }
    k_max <- max(
        ceiling(log(10, eta)),
        ceiling(log(lambda_max * (c + sqrt(5 * v / log(1/alpha))), eta))
    )
    k <- 0:k_max
    lambda <- lambda_max / eta^(k + 1/2)
    weight <- (
        lambda_max * (eta - 1) * mixture_density(lambda * sqrt(eta), lambda_max)
        / eta^(k + 1)
    )
    root_function <- function(s) {
        sum(weight * exp(lambda * s - psi(lambda) * v)) - 1 / alpha
    }
    uniroot(root_function, c(sqrt(v), 10 * sqrt(v)), extendInt='upX')$root
}
discrete_mixture_bound <- function(vs, ...) {
    sapply(vs, function(v) discrete_mixture_bound_nonvec(v, ...))
}

two_sided_normal_mixture <- function(v, alpha, rho) {
    sqrt((v + rho) * log((v + rho) / (alpha^2 * rho)))
}

one_sided_normal_mixture <- function(v, alpha, rho) {
    sqrt(2 * (v + rho) * log(1 / (2 * alpha) * sqrt((v + rho) / rho) + 1))
}

lower_incomplete_gamma <- function(a, x) {
    pgamma(x, a)
}

gamma_exponential_log_mixture <- function(s, v, rho, c=1) {
    rho_csq = rho / c^2
    K2 = (v + rho) / c^2
    constant_term = (
        rho_csq * log(rho_csq) + lgamma(K2) - lgamma(rho_csq)
        - log(lower_incomplete_gamma(rho_csq, rho_csq))
    )
    K1 = (c * s + v) / c^2
    (K1 -
     K2 * log(K1 + rho_csq) +
     log(lower_incomplete_gamma(K2, K1 + rho_csq)) +
     constant_term)
}

gamma_exponential_bound <- function(v, alpha, rho, c=1) {
    sapply(v, function(vv) {
        uniroot(function(s) {
            gamma_exponential_log_mixture(s, vv, rho, c) - log(1/alpha)
        }, c(0, vv), extendInt='upX', tol=1e-10)$root
    })
}
