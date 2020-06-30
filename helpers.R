
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
    d = read.table(paste(bp_path, 'limits/limits_N',N,'.dat', sep=''), header=T)

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
    d = read.table(paste(bp_path, 'limits/limits_N',N,'.dat', sep=''), header=T)
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

# plot a region withing the min and max limits for the CCEP
# the values (alread equalized) for the limits (continua and trozos)
# must be informed
plot.region = function(alpha=0.5, col=1, lty=1, cont, troz)
{
    midy = alpha*cont[,2] + (1-alpha)*troz[,2]
    midx = alpha*cont[,1] + (1-alpha)*troz[,1]

    lines(midx, midy, col=col, lty=lty)
}


plot.tg = function(g, layout=layout.circle, weights=TRUE)
{
    plot(g, 
         layout=layout, 
         edge.label=round(E(g)$weight, 3)
         )
}


# function to plot the BP transition graph (to DEBUG)
# - x: a first time series
# - y: an optional second time series
# - D: embedded dimension
# - tau: embedded delay
# - filename: if specified, saves the plot to the 'filename' path
# - empty: TRUE to remove the vertices with 0 degree
# - xlab, ylab: labels to the plots for the two time series
# - xfactors, yfactors: vectors informinig (edge, vertex) sizes to t.graph
plot.tg2 = function(x, y=NULL, D=4, tau=1, 
                   filename=NULL, empty=TRUE,
                   xlab="x", ylab="y", axislim=NULL,
                   xfactors=c(5,2), yfactors=c(5,2),
                   quiroz=FALSE
                   )
{
    x = as.numeric(x)

    # get the transition graph
    gx =  bandt_pompe_transition_graph(x, D=D, tau=tau, empty=empty)
    # the BP distributions
    bpx = bandt_pompe_distribution(x, D=D, tau=tau)
    # entropies
    Hx = shannon_entropy(bpx$prob, norm=T)

    # different behavior if y is provided
    if (!is.null(y))
    {
        y = as.numeric(y)
        gy = bandt_pompe_transition_graph(y, D=D, tau=tau, empty=empty)
        bpy = bandt_pompe_distribution(y, D=D, tau=tau)
        Hy = shannon_entropy(bpy$prob, norm=T)
    
        # adjusting figure sizes
        par(mfrow=c(2,4))
    }
    else
    {
        # adjusting figure sizes
        par(mfrow=c(1,4))
    }

    if (!is.null(filename))
    {
        pdf(filename, width=11)
    }

    # plotting the time series
    # order of figures:
    # - raw time series
    # - histogram of time series
    # - barplots of BP distribution
    # - BP transition graphs

    # just to easy the visualization
    if (is.null(axislim))
    {
        axislim = length(x)
    }

    # plotting figures for x
    plot(1:axislim, x[1:axislim], type='b', 
         pch=19, main=paste("Time series of",xlab))
    hist(x, main=paste("Histogram of",xlab))
    barplot(bpx$prob, main=paste("Bandt-Pompe dist. of",xlab), 
            names.arg=bpx$patterns, las=2)
    #plot(gx, layout=layout.kamada.kawai, 
    plot(gx, #layout=layout.circle, 
            edge.label=round(E(gx)$weight, 3),
            main=paste("Transition Graph of",xlab), sub=paste('H:',Hx),
            edge.width=xfactors[1]*E(gx)$weight + 1-min(E(gx)$weight), 
            vertex.size=degree(gx) * xfactors[2])

    if (! is.null(y))
    {
        # plotting figures for y
        plot(1:axislim, y[1:axislim], type='b', pch=19, 
             main=paste("Time series of",ylab))
        hist(y, main=paste("Histogram of",ylab))
        barplot(bpy$prob, main=paste("Bandt-Pompe dist. of",ylab), 
                names.arg=bpy$patterns, las=2)
        #plot(gy, layout=layout.kamada.kawai, 
        plot(gy, #layout=layout.circle, 
                edge.label=round(E(gy)$weight, 3),
                main=paste("Transition Graph of",ylab), sub=paste('H:',Hy),
                edge.width=yfactors[1]*E(gy)$weight + 1-min(E(gy)$weight), 
                vertex.size=degree(gy) * yfactors[2])
    }

    if (!is.null(filename))
    {
        suppressMessages(dev.off())
    }

    par(mfrow=c(1,1))

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


######
# function to plot the BP distribution
# bpd = recevies as input the result from bandt_pompe_distribution()
gplot.bpd = function(bpd, title='')
{
    bpd.df = data.frame(x=bpd$patterns, y=bpd$prob)
    p = ggplot(aes(x=x, y=y), data=bpd.df) +
            geom_bar(stat='identity', width=0.7, 
                     fill='darkgray', colour='black') + 
            xlab('Patterns') + ylab('Probabilities') + theme_bw() + 
            theme(axis.text.x = element_text(angle = 90, hjust = 1),
                  plot.title=element_text(hjust=0.5)) +
            ggtitle(title)
    return(p)
}


# function to merge two Bandt-Pompe distributions in one
merge.bpd = function(p1, p2)
{
    # checking identical patterns
    if(sum(p1$patterns == p2$patterns) != length(p1$patterns))
    {
        print('error')
        return(NULL)
    }

    # merging
    p1$frequencies = p1$frequencies + p2$frequencies
    p1$probabilities = p1$frequencies/sum(p1$frequencies)

    return(p1)
}
































# TODO: [deprecated] old functions

# Function to plot the CCEP plane
# - the points may be informed (H,SC), or empty to plot just the limits
plot.ccep.old = function(H=NULL, SC=NULL, D=4, main='', 
                   xlim=c(0,1), ylim=c(0,1), lwd=1, col=1)
{
    # NOTE: var. names was inherited from original code of Rosso et al.

    cont_name = paste(bp_path, '/limits/rosso/continua-N',factorial(D),'.q1', sep='')
    trozos_name = paste(bp_path, '/limits/rosso/trozos-N',factorial(D),'.q1', sep='')

    # upper and lower limits
    continua = read.table(cont_name, skip=7)
    trozos = read.table(trozos_name, skip=7)

    maxY = max(continua$V2)
      
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
  
    lines(trozos, lwd=lwd, col=col)
    lines(continua, lwd=lwd, col=col)
}


# Function to plot the CCEP plane with ggplot
# - the points may be informed (H,SC), or empty to plot just the limits
# NOTE: it returns the ggplot object p
gplot.ccep.old = function(H=NULL, SC=NULL, D=4, main='', 
                   xlim=c(0,1), ylim=c(0,1), lwd=1, col=1, shp=1)
{
    # NOTE: var. names was inherited from original code of Rosso et al.

    cont_name = paste(bp_path, '/limits/rosso/continua-N',factorial(D),'.q1', sep='')
    trozos_name = paste(bp_path, '/limits/rosso/trozos-N',factorial(D),'.q1', sep='')

    # upper and lower limits
    continua = read.table(cont_name, skip=7)
    trozos = read.table(trozos_name, skip=7)

    # renaming the columns
    colnames(continua) = c('x', 'y')
    colnames(trozos) = c('x', 'y')

    maxY = max(continua$y)

    continua.df = as.data.frame(continua)
    trozos.df = as.data.frame(trozos)

    p = ggplot() +
        geom_line(aes(x, y), data=continua.df) +
        geom_line(aes(x, y), data=trozos.df) + 
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




##### deprecated
# get the BP transition graph
# - D: embedded dimension
# - tau: embedded delay
# - empty: TRUE to remove the vertices with 0 degree
get.tg = function(x, D=4, tau=1, empty=TRUE, quiroz=FALSE)
{
    # adjacency matrix from bandt pompe transition
    A = bandt_pompe_transition(x, D=D, tau=tau, quiroz=quiroz)
    
    # the graph
    gA = graph_from_adjacency_matrix(A, mode="directed", weighted=TRUE)
    
    # removing vertices without transition
    if (empty == TRUE)
    {
        gA = delete_vertices(gA, which(degree(gA) == 0))
    }

    return(gA)
}


##### deprecated
# compute and print some features from the transition graph
# - x: time series
# - D: embedded dimension
# - tau: embedded delay
# - xlab: time series label
# - empty: TRUE to remove the vertices with 0 degree
# - debug: TRUE to print detailed features
# - header: TRUE to print the header of the table
feat.tg = function(x, D=4, tau=1, lab='x', empty=FALSE, debug=FALSE, header=FALSE)
{
    # the features
    #feats = data.frame()
    feats = list()
    
    # handling the case of 
    if (length(x) >= D)
    {
        # getting the transition graph
        gx = get.tg(x, D=D, tau=tau, empty=empty)

        # computing the features vector
        feats$weight.avg          = mean(E(gx)$weight)
        feats$degree.avg          = mean(degree(gx))
        #feats$degree.max          = max(degree(gx))
        #feats$degree.min          = min(degree(gx))
        feats$density             = edge_density(gx, loops=TRUE)
        feats$distance.avg        = mean_distance(gx, directed=TRUE)
        feats$diameter            = diameter(gx)
        feats$transitivity.avg    = mean(transitivity(gx, type='weighted'), 
                                    na.rm=TRUE)
        feats$betweeness.edge.avg = mean(edge_betweenness(gx))
        feats$betweeness.vert.avg = mean(betweenness(gx))
        feats$closeness           = mean(closeness(gx, mode='out'))
        feats$assortativity       = assortativity_degree(gx)
        feats$reciprocity         = reciprocity(gx, ignore.loops=FALSE)

        # adding more features

        # degree
        feats$degree.shannon = shannon_entropy(degree_distribution(gx),
                                normalized=T)
        feats$degree.complexity = complexity(degree_distribution(gx))
        feats$degree.fisher = fisher_information(degree_distribution(gx))
        
        # transitivity
        transitivity.dist = transitivity(gx, type='weighted')
        transitivity.dist = transitivity.dist[!is.na(transitivity.dist)]
        if (sum(transitivity.dist) != 0)
        {
            transitivity.dist = transitivity.dist[!is.infinite(
                                 transitivity.dist)]
            transitivity.dist = transitivity.dist/sum(transitivity.dist)
            feats$transitivity.sha = shannon_entropy(transitivity.dist,
                                        normalized=T)
            feats$transitivity.com = complexity(transitivity.dist)
            feats$transitivity.fis = fisher_information(transitivity.dist)
        }
        else
        {
            feats$transitivity.sha = 0
            feats$transitivity.com = 0
            feats$transitivity.fis = 0
        }

        # edge betweenness
        edge_betweenness.dist     = edge_betweenness(gx)
        if (sum(edge_betweenness.dist) != 0)
        {
            edge_betweenness.dist = edge_betweenness.dist[!is.na(
                                        edge_betweenness.dist)]
            edge_betweenness.dist = edge_betweenness.dist[!is.infinite(
                                        edge_betweenness.dist)]
            edge_betweenness.dist = edge_betweenness.dist/sum(
                                        edge_betweenness.dist)
            feats$betweeness.edge.sha = shannon_entropy(
                                            edge_betweenness.dist,
                                            normalized=T)
            feats$betweeness.edge.com = complexity(edge_betweenness.dist)
            feats$betweeness.edge.fis = fisher_information(
                                            edge_betweenness.dist)
        }
        else
        {
            feats$betweeness.edge.sha = 0
            feats$betweeness.edge.com = 0
            feats$betweeness.edge.fis = 0
        }

        # vertex betweenness
        betweenness.dist = betweenness(gx)
        if (sum(betweenness.dist) != 0)
        {
            betweenness.dist = betweenness.dist[!is.na(betweenness.dist)]
            betweenness.dist = betweenness.dist[!is.infinite(
                                betweenness.dist)]
            betweenness.dist = betweenness.dist/sum(betweenness.dist)
            feats$betweeness.vert.sha = shannon_entropy(betweenness.dist,
                                        normalized=T)
            feats$betweeness.vert.com = complexity(betweenness.dist)
            feats$betweeness.vert.fis = fisher_information(betweenness.dist)
        }
        else
        {
            feats$betweeness.vert.sha = 0
            feats$betweeness.vert.com = 0
            feats$betweeness.vert.fis = 0
        }
        
        # closeness
        closeness.dist            = closeness(gx, mode='out')
        if (sum(closeness.dist) != 0)
        {
            closeness.dist = closeness.dist[!is.na(closeness.dist)]
            closeness.dist = closeness.dist[!is.infinite(closeness.dist)]
            closeness.dist = closeness.dist/sum(closeness.dist)
            feats$closeness.vert.sha = shannon_entropy(closeness.dist,
                                        normalized=T)
            feats$closeness.vert.com = complexity(closeness.dist)
            feats$closeness.vert.fis = fisher_information(closeness.dist)
        }
        else
        {
            feats$closeness.vert.sha = 0
            feats$closeness.vert.com = 0
            feats$closeness.vert.fis = 0
        }

    }
    else
    {
        feats$weight.avg          = 0
        feats$degree.avg          = 0
        #feats$degree.max          = 0
        #feats$degree.min          = 0
        feats$density             = 0
        feats$distance.avg        = 0
        feats$diameter            = 0
        feats$transitivity.avg    = 0
        feats$betweeness.edge.avg = 0
        feats$betweeness.vert.avg = 0
        feats$closeness           = 0
        feats$assortativity       = 0
        feats$reciprocity         = 0
        feats$degree.shannon      = 0
        feats$degree.complexity   = 0
        feats$degree.fisher       = 0
        feats$transitivity.sha    = 0
        feats$transitivity.com    = 0
        feats$transitivity.fis    = 0
        feats$betweeness.edge.sha = 0
        feats$betweeness.edge.com = 0
        feats$betweeness.edge.fis = 0
        feats$betweeness.vert.sha = 0
        feats$betweeness.vert.com = 0
        feats$betweeness.vert.fis = 0
        feats$closeness.vert.sha  = 0
        feats$closeness.vert.com  = 0
        feats$closeness.vert.fis  = 0
    
    }
    
    # converting to data.frame (???)
    feats = as.data.frame(feats)

    # printing the header (features names)
    if (header==TRUE)
    {
        cat(names(feats), '\n')
    }

    #dist = distance_table(gx)
    #cat('max distance:', max(dist$res),'\n')
    #cat('min distance:', min(dist$res),'\n')

    # TODO: tvz nao tenha muito sentido essas metricas serem as medias,
    # é provavel que o mais adequado seja a entropia delas, como uma
    # distribuicao de uma dada quantidade/caracteristica pelos nós!!! 

    # clustering coefficient
    #cat('transitivity global:', mean(transitivity(gx, type='global'), na.rm=TRUE), '\n')
    # betweeness

    # ditribuicao dos graus dos nos?
    # distribuicao da betweeness?
    
    feats_res = as.numeric(feats)

    attr(feats_res, 'cols') = names(feats)

    return(feats_res)
}


# Helper function to call the Bandt-Pompe functions
complexity_entropy_old = function(data, D=6, tau=1, dist='bp', 
                              entropy='shannon', breaks=40, method='deg')
{   
    # the distribution method
    # bandt-pompe
    if (dist == 'bp')
    {
        #da = bandt_pompe_distribution(data, D,tau)
        #da = bandt_pompe_gap(data, D,tau)
        da = bandt_pompe_distribution(data, D,tau)

    # visibility graph (nodes degree distribution)
    } 
    else if (dist == 'visibility')
    {
        da = visibility_graph(data, breaks=breaks)
    } 
    else if (dist == 'hvg')
    {
        if (method == 'dist')
        {
            unit=TRUE
        }
        else
        {
            unit = FALSE
        }

        da = hvg(data, method=method, unit=unit)
    }

    if(entropy == 'shannon')
    {
      H = shannon_entropy(da$probabilities, normalized=TRUE)
    }
    C = complexity(da$probabilities,H)
    return(c(H=H,C=C[1],JS=C[2]))
}
