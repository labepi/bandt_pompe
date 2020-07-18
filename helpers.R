
# Helper function to call the Bandt-Pompe functions
# - data is a vector
# - D and tau: parameters
# - equal: if TRUE consider equal sequences independently
complexity_entropy = function(data, D=4, tau=1, equal=FALSE)
{   
    # the bandt-pompe distribution
    da = bandt_pompe_distribution(data, D, tau, equal=equal)

    # computing information theory quantifiers

    # shannon entropy
    H = shannon_entropy(da$probabilities, normalized=TRUE)
    
    # statistical complexity
    C = complexity(da$probabilities,H)

    return(
           c(
             H=H,
             C=C,
             JS=attr(C, 'JS')
             )
            )
}

# helper function to compute BP from the distribution
shannon_complexity = function(bpd, D=3, tau=1, normalized=TRUE)
{
    # shannon entropy
    H = shannon_entropy(bpd$probabilities, normalized)
    
    # statistical complexity
    C = complexity(bpd$probabilities, H, normalized)

    return(
           c(
             H=H,
             C=C,
             JS=attr(C, 'JS')
             )
            )
}

# Function to plot the CCEP plane
# - the points may be informed (H,SC), or empty to plot just the limits
plot.ccep = function(H=NULL, SC=NULL, D=4, main='', N=NULL,
                   xlim=c(0,1), ylim=c(0,1), lwd=1, col=1, add=FALSE)
{
    # added an option to first check by the N
    if (is.null(N))
    {
        N = factorial(D)
    }
    
    # limits
    d = read.table(paste(bandt_pompe_path, 'limits/limits_N',N,'.dat', sep=''), header=T)

    if (add==FALSE)
    {
        if ( !is.null(H) & !is.null(SC) )
        {
            plot(H, SC, 
                main=main, xlim=xlim, ylim=ylim, 
                xlab="Normalized Shannon Entropy", 
                ylab="Statistical Complexity"
            )
        } else {
            plot(NA, 
                main=main, xlim=xlim, ylim=ylim, 
                xlab="Normalized Shannon Entropy", 
                ylab="Statistical Complexity"
            )
        }
    }
  
    lines(d$H, d$SC, lwd=lwd, col=col)
}


# Function to plot the CCEP plane with ggplot
# - the points may be informed (H,SC), or empty to plot just the limits
# NOTE: it returns the ggplot object p
gplot.ccep = function(H=NULL, SC=NULL, D=4, main='', N=NULL,
                   xlim=c(0,1), ylim=c(0,1), lwd=1, col=1, shp=1)
{
    # added an option to first check by the N
    if (is.null(N))
    {
        N = factorial(D)
    }
    
    # limits
    d = read.table(paste(bandt_pompe_path, 'limits/limits_N',N,'.dat', sep=''), header=T)
    d.df = as.data.frame(d)

    p = ggplot() +
        #geom_line(aes(H, SC), data=d.df) +
        geom_path(aes(H, SC), data=d.df) +
        labs(x="Normalized Shannon Entropy (H)",
             y="Statistical Complexity (C)",
             title=main) +
        coord_cartesian(xlim=xlim, ylim=ylim) +
        theme_bw()

    if ( !is.null(H) & !is.null(SC) )
    {
        H_SC = data.frame(x=H, y=SC)
        p = p + geom_point(aes(x, y), data=H_SC, color=col, shape=shp)
    }

    return(p)
}

# checking if the parameters are valid
# - m is the length of the series
# - D is the embedding dimension
# - tau is the embedding delay
# - lim is the minimum number of patterns allowed
# NOTE:
# - the maximum possible value for tau is
#   tau < m / (D-1)
checkParameters = function(m, D, tau, lim=1)
{
    if ( ( m - (D-1)*tau ) >= lim )
    {
        return(TRUE)
    }

    return(FALSE)
}



# checking the maximum possible value for tau
# - m is the length of the series
# - D is the embedding dimension
# - lim is the minimum number of patterns allowed
#
# - the maximum possible value for tau, for having at least one pattern,
#   must satisfy the condition:
#       tau < m / (D-1)
checkMaxTau = function(m, D, lim=1)
{
    return(ceiling((m-(lim-1))/(D-1)-1))
}

