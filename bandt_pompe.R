

# Implementation of the Bandt-Pompe transformations

# required packages

suppressMessages(library(igraph))

require(e1071, quietly=TRUE)

suppressMessages(library(entropy))

library(Rcpp)


if (!exists('bandt_pompe_path'))
{
    bandt_pompe_path = '.'
}

# Rcpp implementation of bandt-pompe functions
sourceCpp(paste(bandt_pompe_path,'bandt_pompe.cpp', sep='/'))


###################
# regular functions
###################

# return an empty list of symbols
# D - the embedding dimension
# equal - if TRUE, it will compute equal sequences separately
# na_aware - if TRUE, it will compute NA patterns separately
bandt_pompe_empty = function(D = 4, equal=FALSE, na_aware=FALSE)
{
    # to get the index of the permutation pi
    perms = sort(apply(permutations(D), 1, paste, collapse=''))

    if (equal==TRUE)
    {
        perms = c(perms, paste(rep(0, D), collapse=''))
    }
    
    if (na_aware==TRUE)
    {
        perms = c(perms, 'NA')
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
# TODO: make this function calls its Cpp version
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
    
    # discovering the sequences of order n
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
# na_aware - if TRUE, the symbols with only NAs will be counted separated
# na_rm - if TRUE and na_aware=TRUE, the "NA patterns" are not counted
bandt_pompe_distribution = function(data, D=4, tau=1, numred=FALSE, by=1, 
                                    equal=FALSE, useSymbols=FALSE, 
                                    na_aware=FALSE, na_rm=FALSE)
{
    # check if the symbols were already computed
    if (useSymbols == TRUE)
    {
        # symbols were passed as data
        #symbols = data

        # counting the patterns in Rcpp
        perms = bandt_pompe_distribution_symbols_c(data, D, tau)
    }
    else
    {
        # the list of symbols from the BP transformation
        # NOTE: changing the main bandt-pompe function to its cpp version
        #symbols = bandt_pompe_c(data, D=D, tau=tau)
        perms = bandt_pompe_distribution_c(data, D, tau, na_aware)
    }

    # not counting the "NA patterns"
    if (na_aware == TRUE & na_rm == TRUE)
    {
        perms = perms[names(perms) != 'NA']
    }

    # the names for the permutations
    perms_n = names(perms)
    
    # the returning format
    output = data.frame(patterns = perms_n, 
                        frequencies = perms,
                        probabilities = perms/sum(perms))
    return(output)
}




# Bandt-Pompe distribution (old version in R)
#
# Parameters:
# data: the time series (univariate vector) or the pre-computed symbols
# D:    the embedding dimension (size of sliding window)
# tau:  the embedding delay ('step' value)
# numred: numerosity reduction, similar to BOSS algorithm, do not count
#         repetitions of the same symbol
# equal - if TRUE, it will compute equal sequences separately
# useSymbols - if TRUE, the symbols were already computed, and passed as 'data'
bandt_pompe_distribution2 = function(data, D=4, tau=1, numred=FALSE, by=1, 
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
        #symbols = bandt_pompe(data, D=D, tau=tau, by=by, equal=equal)

        # NOTE: changing the main bandt-pompe function to its cpp version
        symbols = bandt_pompe_c(data, D=D, tau=tau)
    }

    # the distribution of permutations (pi)
    #dpi = rep(0, factorial(D))
    
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

# Bandt-Pompe transformation
# TODO: make this function calls its Cpp version
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
bandt_pompe_tie = function(data, D=4, tau=1)
{
    # compute the symbols, removing the tied patterns
    symbols = bandt_pompe_tie_c(data, D, tau, NULL)

    # compute the probability distribution of patterns from these symbols
    bpd = bandt_pompe_distribution_symbols_c(symbols, 3, 1)

    # normalize distribution
    probs = bpd/sum(bpd)

    # compute the symbols, imputing tied patterns with the a priori
    # distribution of patterns without the tied patterns
    symbols = bandt_pompe_tie_c(data, D, tau, prob=probs)


    return(symbols)
}

