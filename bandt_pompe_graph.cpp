#include <Rcpp.h>
#include <iostream>
#include <algorithm>
#include <map>
using namespace Rcpp;

// Returns an empty matrix (list of lists) of symbols
// 
// The list is a map, to counting the patterns
//
// @param D The embedding dimension
// 
// [[Rcpp::export]]
std::map<std::string, std::map<std::string, int>> bandt_pompe_empty_matrix_c(int D=3)
{
    // the matrix of maps
    std::map<std::string, std::map<std::string, int>> M;

    int i, itens[D];

    // initializing the list of itens
    for(i = 0; i < D; i++)
    {
        itens[i] = i+1;
    }

    String s;
    
    // create each permutation and initialize it with 0
    do {
        s = "";
        for(i = 0; i < D; i++)
        {
            s += std::to_string(itens[i]);
        }
        M[s] = bandt_pompe_empty_c(D);
    } while ( std::next_permutation(itens,itens+D) );

    return M;
}

//' The Band-Pompe Transition (adjacency matrix)
//'
//' Parameters:
//' data: the time series (univariate vector) or the pre-computed symbols
//' D:    the embedding dimension (size of sliding window)
//' tau:  the embedding delay ('step' value)
//
// [[Rcpp::export]]
std::map<std::string, std::map<std::string, int>> 
bandt_pompe_transition_c(NumericVector x, int D=3, int tau=1)
{
    // the transition matrix
    std::map<std::string, std::map<std::string, int>> M = bandt_pompe_empty_matrix_c(D);

    // computing the symbols
    StringVector symbols = bandt_pompe_c(x, D, tau);

    // counting the transitions
    for(int i = 1; i < symbols.size(); i++)
    {
        ++M[(String)symbols[i-1]][(String)symbols[i]];
    }

    return M;
}

//' The Band-Pompe Transition (adjacency matrix)
//'
//' Parameters:
//' @symbols the pre-computed symbols
//' @D    the embedding dimension (size of sliding window)
//' @tau  the embedding delay ('step' value)
//
// [[Rcpp::export]]
std::map<std::string, std::map<std::string, int>> 
bandt_pompe_transition_symbols_c(StringVector symbols, int D=3, int tau=1)
{
    // the transition matrix
    std::map<std::string, std::map<std::string, int>> M = bandt_pompe_empty_matrix_c(D);

    // counting the transitions
    for(int i = 1; i < symbols.size(); i++)
    {
        ++M[(String)symbols[i-1]][(String)symbols[i]];
    }

    return M;
}

