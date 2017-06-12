
#############################
#   Fold change function    #
#############################


fold_change_power.fun <- function(n1 = NULL, n2 = NULL, d = NULL, sig.level = 0.05, power = NULL, 
    alternative = c("two.sided", "less", "greater")) 
{
    if (sum(sapply(list(n1, n2, d, power, sig.level), is.null)) != 
        1) 
        stop("exactly one of n1, n2, d, power, and sig.level must be NULL")
    if (!is.null(sig.level) && !is.numeric(sig.level) || any(0 > 
        sig.level | sig.level > 1)) 
        stop(sQuote("sig.level"), " must be numeric in [0, 1]")
    if (!is.null(power) && !is.numeric(power) || any(0 > power | 
        power > 1)) 
        stop(sQuote("power"), " must be numeric in [0, 1]")
    if (!is.null(n1) && n1 < 2) 
        stop("number of observations in the first group must be at least 2")
    if (!is.null(n2) && n2 < 2) 
        stop("number of observations in the second group must be at least 2")
    alternative <- match.arg(alternative)
    tsample <- 2
    ttside <- switch(alternative, less = 1, two.sided = 2, greater = 3)
    tside <- switch(alternative, less = 1, two.sided = 2, greater = 1)
    if (tside == 2 && !is.null(d)) 
        d <- abs(d)
    if (ttside == 1) {
        p.body <- quote({
            nu <- n1 + n2 - 2
            pt(qt(sig.level/tside, nu, lower = TRUE), nu, ncp = d * 
                (1/sqrt(1/n1 + 1/n2)), lower = TRUE)
        })
    }
    if (ttside == 2) {
        p.body <- quote({
            nu <- n1 + n2 - 2
            qu <- qt(sig.level/tside, nu, lower = FALSE)
            pt(qu, nu, ncp = d * (1/sqrt(1/n1 + 1/n2)), lower = FALSE) + 
                pt(-qu, nu, ncp = d * (1/sqrt(1/n1 + 1/n2)), 
                  lower = TRUE)
        })
    }
    if (ttside == 3) {
        p.body <- quote({
            nu <- n1 + n2 - 2
            pt(qt(sig.level/tside, nu, lower = FALSE), nu, ncp = d * 
                (1/sqrt(1/n1 + 1/n2)), lower = FALSE)
        })
    }
    if (is.null(power)) 
        power <- eval(p.body)
    else if (is.null(n1)) 
        n1 <- uniroot(function(n1) eval(p.body) - power, c(2 + 
            1e-10, 1e+07))$root
    else if (is.null(n2)) 
        n2 <- uniroot(function(n2) eval(p.body) - power, c(2 + 
            1e-10, 1e+07))$root
    else if (is.null(d)) {
        if (ttside == 2) {
            d <- uniroot(function(d) eval(p.body) - power, c(1e-07, 
                10))$root
        }
        if (ttside == 1) {
            d <- uniroot(function(d) eval(p.body) - power, c(-10, 
                5))$root
        }
        if (ttside == 3) {
            d <- uniroot(function(d) eval(p.body) - power, c(-5, 
                10))$root
        }
    }
    else if (is.null(sig.level)) 
        sig.level <- uniroot(function(sig.level) eval(p.body) - 
            power, c(1e-10, 1 - 1e-10))$root
    else stop("internal error")
    METHOD <- c("t test power calculation")
    structure(list(n1 = n1, n2 = n2, d = d, sig.level = sig.level, 
        power = power, alternative = alternative, method = METHOD), 
        class = "power.htest")
}


# Power of statistical tests
# joseph powell: joseph.powell@uq.edu.au
# 2/3/12

# 1.Function to calculate the power for unqual case control examples


pow_unequal_CC.fun <- function (
    h = NULL,           # Effect size - listed as odds ratio
    n1 = NULL,          # Number cases 
    n2 = NULL,          # Number of controls    
    sig.level = 0.05,   # Defult sig level - p.value    
    power = NULL,       # Leave blank if you want to calc the   
    alternative = c("two.sided", "less", "greater")) 
{
    if (sum(sapply(list(h, n1, n2, power, sig.level), is.null)) != 
        1) 
        stop("exactly one of h, n1, n2, power, and sig.level must be NULL")
    if (!is.null(n1) && n1 < 2) 
        stop("number of observations in the first group must be at least 2")
    if (!is.null(n2) && n2 < 2) 
        stop("number of observations in the second group must be at least 2")
    if (!is.null(sig.level) && !is.numeric(sig.level) || any(0 > 
        sig.level | sig.level > 1)) 
        stop(sQuote("sig.level"), " must be numeric in [0, 1]")
    if (!is.null(power) && !is.numeric(power) || any(0 > power | 
        power > 1)) 
        stop(sQuote("power"), " must be numeric in [0, 1]")
    alternative <- match.arg(alternative)
    tside <- switch(alternative, less = 1, two.sided = 2, greater = 3)
    if (tside == 2 && !is.null(h)) 
        h <- abs(h)
    if (tside == 3) {
        p.body <- quote({
            pnorm(qnorm(sig.level, lower = FALSE) - h * sqrt((n1 * 
                n2)/(n1 + n2)), lower = FALSE)
        })
    }
    if (tside == 1) {
        p.body <- quote({
            pnorm(qnorm(sig.level, lower = TRUE) - h * sqrt((n1 * 
                n2)/(n1 + n2)), lower = TRUE)
        })
    }
    if (tside == 2) {
        p.body <- quote({
            pnorm(qnorm(sig.level/2, lower = FALSE) - h * sqrt((n1 * 
                n2)/(n1 + n2)), lower = FALSE) + pnorm(qnorm(sig.level/2, 
                lower = TRUE) - h * sqrt((n1 * n2)/(n1 + n2)), 
                lower = TRUE)
        })
    }
    if (is.null(power)) 
        power <- eval(p.body)
    else if (is.null(h)) {
        if (tside == 2) {
            h <- uniroot(function(h) eval(p.body) - power, c(1e-10, 
                10))$root
        }
        if (tside == 1) {
            h <- uniroot(function(h) eval(p.body) - power, c(-10, 
                5))$root
        }
        if (tside == 3) {
            h <- uniroot(function(h) eval(p.body) - power, c(-5, 
                10))$root
        }
    }
    else if (is.null(n1)) 
        n1 <- uniroot(function(n1) eval(p.body) - power, c(2 + 
            1e-10, 1e+05))$root
    else if (is.null(n2)) 
        n2 <- uniroot(function(n2) eval(p.body) - power, c(2 + 
            1e-10, 1e+05))$root
    else if (is.null(sig.level)) 
        sig.level <- uniroot(function(sig.level) eval(p.body) - 
            power, c(1e-10, 1 - 1e-10))$root
    else stop("internal error")
    structure(list(h = h, n1 = n1, n2 = n2, sig.level = sig.level, 
        power = power, alternative = alternative, class = "power.htest"))
}



# Usage
# pow_unequal_CC.fun(0.1, 390, n2=400, sig.level=0.05)





#2. Power to detect if something is expressed significantly

#/- Function to calcualte power of chi-sq test -\#
pow_chisq.fun <- function (
    w = NULL,       # Effect size
    N = NULL,       # Number of samples
    df = NULL,      # degrees of freedom in the test 
    sig.level=0.05, # alpha (pvalue) level
    power = NULL)   # leave blank if wanting to calculate the power) 
{

    if (sum(sapply(list(w, N, df, power, sig.level), is.null)) != 
        1) 
        stop("exactly one of w, N, df, power, and sig.level must be NULL")

    if (!is.null(w) && w < 0) 
        stop("w must be positive")

    if (!is.null(N) && N < 1) 
        stop("number of observations must be at least 1")

    if (!is.null(sig.level) && !is.numeric(sig.level) || any(0 > 
        sig.level | sig.level > 1)) 
        stop(sQuote("sig.level"), " must be numeric in [0, 1]")

    if (!is.null(power) && !is.numeric(power) || any(0 > power | 
        power > 1)) 
        stop(sQuote("power"), " must be numeric in [0, 1]")

    p.body <- quote({
        k <- qchisq(sig.level, df = df, lower = FALSE)
        pchisq(k, df = df, ncp = N * w^2, lower = FALSE)
    })

    if (is.null(power)) 
        power <- eval(p.body)

    else if (is.null(w)) 
        w <- uniroot(function(w) eval(p.body) - power, c(1e-10, 
            1e+05))$root

    else if (is.null(N)) 
        N <- uniroot(function(N) eval(p.body) - power, c(1 + 
            1e-10, 1e+05))$root

    else if (is.null(sig.level)) 
        sig.level <- uniroot(function(sig.level) eval(p.body) - 
            power, c(1e-10, 1 - 1e-10))$root

    else stop("internal error")
    structure(list(w = w, N = N, df = df, sig.level = sig.level, 
        power = power), class = "power.htest")
}






