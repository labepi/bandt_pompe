

# Implementation of the Bandt-Pompe transformations

source('includes.R')

###################
# regular functions
###################

# return an empty list of symbols
# D - the embedding dimension
# equal - if TRUE, it will compute equal sequences separately
bandt_pompe_empty = function(D = 4, equal=FALSE)
{
    # to get the index of the permutation pi
    perms = sort(apply(permutations(D), 1, paste, collapse=''))

    if (equal==TRUE)
    {
        perms = c(perms, paste(rep(0, D), collapse=''))
    }

    # the list of empty elements
    res = rep(0, length(perms))

    # defining name index
    names(res) = perms

    return(res)
}

# compute a single pattern
# w - the sliding window being permutted
# equal - if TRUE, it will compute equal sequences separately
bandt_pompe_pattern = function(w, equal=FALSE)
{
    if (equal==TRUE)
    {
        if(sd(w) == 0)
        {
            return(paste(rep(0, length(w)), collapse=''))
        }
    }

    # the current permutation pattern
    pattern = paste(order(w), collapse='')

    return(pattern)
}

# Bandt-Pompe transformation
#
# data: the time series (univariate vector)
# D:    the embedding dimension (size of sliding window)
# tau:  the embedding delay ('step' value)
# equal: if TRUE, it will compute equal sequences separately
# amplitude: if TRUE, it computes the amplitude of each sliding window
#           (largest difference between elements)
#
# return: The list of symbols from the transformed dataset,
#         the list of differences (if diff is TRUE)
bandt_pompe = function(data, D=4, tau=1, by=1, equal=FALSE, amplitude=FALSE)
{
    # the list of symbols to be returned
    symbols = c()

    # the amplitude, the largest differences of elements for each sub
    # (sliding window)
    amps = c()
    
    # dicovering the sequences of order n
    for (s in seq(1, length(data)-(D-1)*tau, by=by))
    {
        # the indices for the subsequence
        ind = seq( s, s+(D-1)*tau, by=tau)
        
        # get the sub-timeseries (sliding window)
        sub = data[ind]

        # computing the differences (amplitudes)
        amps = c(amps, (max(sub) - min(sub)))

        # the current permutation pattern
        pattern = bandt_pompe_pattern(sub, equal=equal)

        # adding the current pattern to the list of symbols
        symbols = c(symbols, pattern)
    }

    # adding the amplitude as an attribute
    if (amplitude == TRUE)
    {
        attr(symbols, 'amplitudes') = amps
    }

    return(symbols)
}

# Bandt-Pompe distribution
#
# Parameters:
# data: the time series (univariate vector) or the pre-computed symbols
# D:    the embedding dimension (size of sliding window)
# tau:  the embedding delay ('step' value)
# numred: numerosity reduction, similar to BOSS algorithm, do not count
#         repetitions of the same symbol
# equal - if TRUE, it will compute equal sequences separately
# useSymbols - if TRUE, the symbols were already computed, and passed as 'data'
bandt_pompe_distribution = function(data, D=4, tau=1, numred=FALSE, by=1, 
                                    equal=FALSE, useSymbols=FALSE)
{
    # check if the symbols were already computed
    if (useSymbols == TRUE)
    {
        # symbols were passed as data
        symbols = data
    }
    else
    {
        # the list of symbols from the BP transformation
        symbols = bandt_pompe(data, D=D, tau=tau, by=by, equal=equal)
    }

    # the distribution of permutations (pi)
    dpi = rep(0, factorial(D))
    
    # to get the index of the permutation pi
    perms = bandt_pompe_empty(D, equal=equal)

    # the names for the permutations
    perms_n = names(perms)

    # there is no reduction
    if (numred == FALSE)
    {
        # dicovering the sequences of order n
        for (i in perms_n)
        {
            # counting the pattern
            perms[i] = sum(symbols == i)
        }
    }
    else
    {
        # only counts different subsequence of symbols
        oldsymbol = symbols[1]
        
        for (i in 2:length(symbols))
        {
            if (oldsymbol != symbols[i])
            {
                perms[symbols[i]] = perms[symbols[i]] + 1
            }
            oldsymbol = symbols[i]
        }
    }
    
    # the returning format
    output = data.frame(patterns = names(perms), 
                        frequencies = perms,
                        probabilities = perms/sum(perms))
    return(output)
}



