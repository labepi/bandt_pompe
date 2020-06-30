
# Metrics to be computed

# Shannon entropy 'H'
#
# probs: vector of probabilities from the Bandt-Pompe distribution
shannon_entropy = function(probs, normalized=FALSE)
{
    # consider only the non-zero values
    p = which(probs > 1e-30)
    
    # considering log base n
    entropy = -sum(probs[p]*log(probs[p]))
    
    # if normalized, the log base is irrelevant
    if(normalized)
    {
        entropy = entropy/log(length(probs))
    }

    return(entropy)
}

# Statistical Complexity
#
# probs: vector of probabilities from the Bandt-Pompe distribution
# entropy: the Shannon entropy for probs
complexity = function(probs, entropy = NULL, normalized=TRUE)
{
    if (is.null(entropy))
    {
        entropy = shannon_entropy(probs, normalized=normalized)
    }

    # the length of the probabilities, 
    N = length(probs)

    # the reference distribution (uniform)
    P_u = rep(1/N, N)

    # the Jensen-shannon divergence
    JS = shannon_entropy( (probs + P_u) / 2  ) - 
        shannon_entropy(probs)/2 - shannon_entropy(P_u)/2

    # the statistical complexity
    aux = (   ((N+1)/N) * log(N + 1) - 2*log(2*N) + log(N)    )
    
    # checking the case ofr  asingle probability being computed
    if (aux != 0)
    {
        Q_0 = -2*(1/aux)
    } else {
        Q_0 = -2
    }
    Q = Q_0 * JS
    C = Q*entropy

    # "piggypacking" the JS variable as an 'attr'
    attr(C, "JS") = JS

    return(C)
}

# Fisher Information
# probs: vector of probabilities from the Bandt-Pompe distribution
fisher_information = function(probs)
{
    # the number of probabilities
    N = length(probs)

    # the normalization constant
    if (probs[1] == 1 | probs[N] == 1)
    {
        F0 = 1
    } else {
        F0 = 1/2
    }

    #print((sqrt(probs[2:N]) - sqrt(probs[1:(N-1)]))^2)

    # the fisher information
    f = F0 * sum((sqrt(probs[2:N]) - sqrt(probs[1:(N-1)]))^2)

    return(f)
}

# computes the probability of self-transitions (pst) from a given
# transition graph
pst = function(g)
{
    # self-transitions
    st = E(g)$weight[which_loop(g)]

    # pst
    return(sum(st))
}

# probability of self transitions
prob_st = function(g)
{
    st = E(g)$weight[which_loop(g)]
    nst = E(g)$weight[!which_loop(g)]

    res = list(pst=sum(st), pnst=sum(nst))

    attr(res, 'st') = st
    attr(res, 'nst') = nst

    return(res)
}


# fast computation of PST using matrix trace
pst_trace = function(x, D=3, tau=1)
{
    res = sum(diag(bandt_pompe_transition(x, D=D, tau=tau)))

    return(res)
}



# Ideas


# test of a modification in the C_JS function altering Pe to be another
# distribution
complexity_distance = function(x, y, D=4, tau=1)
{
    # computing the B-P distribution
    xp = bandt_pompe_distribution(x, D, tau)
    yp = bandt_pompe_distribution(y, D, tau)


    # the length of the probabilities, 
    N = length(yp$probabilities)

    # the Jensen-shannon divergence
    JS = shannon_entropy( (xp$probabilities + yp$probabilities) / 2  ) 
            - 
          shannon_entropy(xp$probabilities)/2
            -
          shannon_entropy(yp$probabilities)/2

    # the statistical complexity
    #aux = (   ((N+1)/N) * log(N + 1) - 2*log(2*N) + log(N)   )
    #Q_0 = -2*(1/aux)
    #Q = Q_0 * JS
    #complexity = Q*entropy
    
    print(JS)
    #print(complexity)
}



# returns a probability distribution that extremizes C for a fixed H
extremizes = function(p, N, m, n)
{
    P = rep(0, N)

    #     { 0                1 <= j <= m
    # P = { p                m+1 <= j <= m+n
    #     { (1-p*n)/(N-m-n)  m+n+1 <= j <= N

    # the case 0 was already set
    P[(m+1):(m+n)] = p
    P[(m+n+1):N]   = (1-p*n)/(N-m-n)

    return (P)
}


# computes the Jensen-Shannon divergence between two distributions
# - distributions must follow the BPD format
JSDiv = function(p, q)
{
    # computing the B-P distribution
    #xp = bandt_pompe_distribution(x, D, tau)
    #yp = bandt_pompe_distribution(y, D, tau)

    #p = xp$probabilities
    #q = yp$probabilities
    
    # the Jensen-shannon divergence
    #JS = shannon_entropy( (xp$probabilities + yp$probabilities) / 2  ) 
    #    - shannon_entropy(xp$probabilities)/2
    #    - shannon_entropy(yp$probabilities)/2

    m = 0.5 * (p + q)
    JS = 0.5 * (sum(p * log2(p / m)) + sum(q * log2(q / m)))

    # the distance
    # The square root of the Jensen–Shannon divergence is a metric often
        # referred to as Jensen-Shannon distance (wikipedia).
    attr(JS, 'd') = sqrt(JS)

    return(JS)

}


