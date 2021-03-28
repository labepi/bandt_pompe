#library(Rcpp)
#
#if (!exists('bandt_pompe_path'))
#{
#    bandt_pompe_path = '.'
#}
#
## Rcpp implementation of bandt-pompe functions
#sourceCpp(paste(bandt_pompe_path,'bandt_pompe_graph.cpp', sep='/'))

############################
# transition graph functions
############################


# The Band-Pompe Transition (adjacency matrix)
#
# Parameters:
# data: the time series (univariate vector) or the pre-computed symbols
# D:    the embedding dimension (size of sliding window)
# tau:  the embedding delay ('step' value)
# normalized: returns as percentage
# markov: normalization as a stochastic matrix
# loop: enable (TRUE) the creation of loopings in the graph
# vect: returns as a vector instead of matrix
# sthocastic: if TRUE, it will create a sthocastic matrix (cols sum up to 1)
# equal: if TRUE, it will compute equal sequences separately
# useSymbols - if True, the symbols were already computed, and passed as 'data'
# amplitude: if TRUE, it computes the amplitude of each sliding window
#           (largest difference between elements)
# na_aware - if TRUE, the symbols with only NAs will be counted separated
# na_rm - if TRUE and na_aware=TRUE, the "NA patterns" are not counted
# hop - the hop to consider for the consecutive pattern in transition        
bandt_pompe_transition = function(data, D=4, tau=1, 
                                  normalized=TRUE, markov=FALSE,
                                  loop=TRUE, vect=FALSE, by=1,
                                  sthocastic=FALSE, 
                                  equal=FALSE, useSymbols=FALSE,
                                  amplitude=FALSE,
                                  na_aware=FALSE, na_rm=FALSE, hop=1)
{
    # the number of permutations
    dfact = factorial(D)

    # plus 000 pattern
    if(equal==TRUE | na_aware==TRUE)
    {
        dfact = dfact+1
    }
    
    # to get the index of the permutation pi
    perms_n = names(bandt_pompe_empty(D, equal=equal, na=na_aware))

    # the adjacency matrix
    M = matrix(0, ncol=dfact, nrow=dfact)

    # labeling the patterns
    rownames(M) = perms_n
    colnames(M) = perms_n

    # the Rcpp functions return a list of edges (lists)

    #buildTime = Sys.time()
    # check if the symbols were already computed
    if (useSymbols == TRUE)
    {
        # symbols were passed as data
        # counting transitions in Rcpp
        L = bandt_pompe_transition_symbols_c(data, D, tau, hop)
    }
    else
    {
        # computing symbols and counting transitions in Rcpp
        L = bandt_pompe_transition_c(data, D, tau, na_aware)
    }
    #buildTime = difftime(Sys.time(), buildTime, units='sec')
    #cat('\tTIME TG_IN:',buildTime,'\n')

    # filling the matrix with the counting of patterns
    for(pattern in names(L))
    {
        M[names(L[[pattern]]),pattern] = L[[pattern]]
    }

    # not counting the "NA patterns"
    if (na_aware == TRUE & na_rm == TRUE)
    {
        M = M[-dfact,-dfact]
    }
    
    # consider the amplitude as the weights
    if (amplitude == TRUE)
    {
        amps = attr(bp, 'amplitude')
        
        # NOTE: for the case of a constant series, it is not possible to
        # compute amplitudes
        if (sd(amps) == 0)
        {
            amplitude = FALSE
        }
    }

    #buildTime = Sys.time()
    if (normalized == TRUE)
    {
        # similar to a markov chain probabilities
        if (markov == TRUE)
        {
            # normalization by columns: stochastic matrix
            M = scale(M, center=FALSE, scale=colSums(M))
            if (sthocastic==FALSE)
            {
                M = t(M)
            }
            M[is.nan(M)] = 0 # for the 0 columns
        }
        else
        {
            # doing a general normalization
            M = M/sum(M)
        }
    }
    #buildTime = difftime(Sys.time(), buildTime, units='sec')
    #cat('\tTIME NORM:',buildTime,'\n')

    if (vect == TRUE)
    {
        M = c(M)
    }
    
    return(M)
}


# The Band-Pompe Transition (adjacency matrix)
#
# Parameters:
# data: the time series (univariate vector) or the pre-computed symbols
# D:    the embedding dimension (size of sliding window)
# tau:  the embedding delay ('step' value)
# normalized: returns as percentage
# markov: normalization as a stochastic matrix
# loop: enable (TRUE) the creation of loopings in the graph
# vect: returns as a vector instead of matrix
# sthocastic: if TRUE, it will create a sthocastic matrix (cols sum up to 1)
# equal: if TRUE, it will compute equal sequences separately
# useSymbols - if True, the symbols were already computed, and passed as 'data'
# amplitude: if TRUE, it computes the amplitude of each sliding window
#           (largest difference between elements)
bandt_pompe_transition2 = function(data, D=4, tau=1, 
                                  normalized=TRUE, markov=FALSE,
                                  loop=TRUE, vect=FALSE, by=1,
                                  sthocastic=FALSE, 
                                  equal=FALSE, useSymbols=FALSE,
                                  amplitude=FALSE)
{
    # the number of permutations
    dfact = factorial(D)

    # plus 000 pattern
    if(equal==TRUE)
    {
        dfact = dfact+1
    }

    # to get the index of the permutation pi
    perms = names(bandt_pompe_empty(D, equal=equal))
    
    # the transitions matrix
    M = matrix(0, ncol=dfact, nrow=dfact)
    rownames(M) = perms
    colnames(M) = perms

    # check if the symbols were already computed
    if (useSymbols == TRUE)
    {
        # symbols were passed as data
        symbols = data
    }
    else
    {
        # the list of symbols from the BP transformation
        #symbols = bandt_pompe(data, D=D, tau=tau, by=by, equal=equal, amplitude=amplitude)
        # NOTE: changing the main bandt-pompe function to its cpp version
        symbols = bandt_pompe_c(data, D=D, tau=tau)
    }

    # consider the amplitude as the weights
    if (amplitude == TRUE)
    {
        amps = attr(bp, 'amplitude')
        
        # NOTE: for the case of a constant series, it is not possible to
        # compute amplitudes
        if (sd(amps) == 0)
        {
            amplitude = FALSE
        }
    }

    # discovering the sequences of order n
    for (i in 2:length(symbols))
    {
        # the previous permutation pattern
        from = symbols[i-1]

        # the current permutation pattern
        to = symbols[i]
        
        # checking if the creation of loopings are allowed
        if (from != to | loop == TRUE)
        {
            # incrementing the counting for this transition
    
            # considering the amplitudes as weights
            if (amplitude == TRUE)
            {
                M[to, from] = M[to, from] + abs(amps[i-1] - amps[i])
            }
            else
            {
                # considering a sthocastic matrix direction
                M[to, from] = M[to, from] + 1
                #M[from, to] = M[from, to] + 1
            }
        }
    }

    if (normalized == TRUE)
    {
        # similar to a markov chain probabilities
        if (markov == TRUE)
        {
            # normalization by columns: stochastic matrix
            M = scale(M, center=FALSE, scale=colSums(M))
            if (sthocastic==FALSE)
            {
                M = t(M)
            }
            M[is.nan(M)] = 0 # for the 0 columns
        }
        else
        {
            # doing a general normalization
            M = M/sum(M)
        }
    }

    if (vect == TRUE)
    {
        M = c(M)
    }
    
    # labeling the patterns
    rownames(M) = perms
    colnames(M) = perms
    
    return(M)
}

# get the BP transition graph
# - D: embedded dimension
# - tau: embedded delay
# - empty: TRUE to remove the vertices with 0 degree
# - markov: normalization proposed by quiroz, similar to a markov-chain
# - sthocastic: if TRUE, it will create a sthocastic matrix (cols sum up to 1)
# - equal: if TRUE, it will compute equal sequences separately
bandt_pompe_transition_graph = function(x, D=4, tau=1, empty=TRUE, 
                                        markov=FALSE, sthocastic=FALSE,
                                        equal=FALSE)
{
    # adjacency matrix from bandt pompe transition
    A = bandt_pompe_transition(x, D=D, tau=tau, 
                               markov=markov, sthocastic=sthocastic,
                               equal=equal)
    
    # to order the graph edges
    #if (markov==TRUE)
    #{
    #    A = t(A)
    #}
    
    # the graph
    gA = graph_from_adjacency_matrix(A, mode="directed", weighted=TRUE)
    
    # removing vertices without transition
    if (empty == TRUE)
    {
        gA = delete_vertices(gA, which(igraph::degree(gA) == 0))
    }

    return(gA)
}



