#include <Rcpp.h>
#include <iostream>
#include <algorithm>
#include <map>
using namespace Rcpp;

#include <chrono> 
using namespace std;
using namespace std::chrono; 

//' Ordering Permutation
//' This function is a replacement to the R's order() function, but only
//'  works for the ascending order.
//'
//' @param x a numeric vector to order
//' @return A numeric vector with the order of the vector fields
// 
// NOTE: the built-in function order() beats orderc() in execution time
// above n=100, wich is very impractal for an embedding dimension D.
//
// [[Rcpp::export]]
NumericVector orderc(NumericVector x)
{
    // the initial vector size
    int n = x.size();

    // auxiliar counters
    int i,j;
    
    // the minimun index
    int mi;
    
    // an array to control the already checked values
    bool checked[n] = {false};

    // the return order
    NumericVector out(n);

    // TODO: this is a naive try

    for(i = 0; i < n; ++i) 
    {
        // set mi as the first checked false
        for(j = 0; j < n; j++)
        {
            if (checked[j] == false)
            {
                mi = j;
                break;
            }
        }

        for(j = 0; j < n; j++)
        {
            if (x[j] < x[mi] && checked[j] == false)
            {
                mi = j;
            }
        }

        //Rprintf("mi:%d and val:%f\n", mi+1, x[mi]);
        
        // setting the value as checked and its position
        checked[mi] = true;
        out[i] = mi+1;
    }

    return out;
}

//' Compute a single pattern
//'
//' @param w The sliding window to be permutted
//' @return Returns the permutation pattern as a string
//
// [[Rcpp::export]]
String bandt_pompe_pattern_c(NumericVector x)
{
    // get the order of x
    NumericVector ord = orderc(x);

    // converting to a string
    String res = "";
    for(int i = 0; i < x.size(); i++)
    {
        res += std::to_string((int)ord[i]);
    }
    
    return res;
}

//' Bandt-Pompe transformation
//'
//' @param x The time series (univariate vector)
//' @param D The embedding dimension (size of sliding window)
//' @param tau The embedding delay ('step' value)
//' @param narm If true, do not count symbols if all elements in w are nan
//
// [[Rcpp::export]]
StringVector bandt_pompe_c(NumericVector x, int D=3, int tau=1)
{
    // aux
    int i, j, s;

    // the lenght of the vector
    int n = x.size();

    // the number of symbols
    int sym_num = n-(D-1)*tau;

    //Rprintf("D,tau: %d %d\n",D,tau);
    //Rprintf("symnum: %d\n",sym_num);
    
    // the list of symbols to be returned
    StringVector symbols(sym_num);

    // the sliding window
    NumericVector w(D);

    // discovering the sequences of order n
    for (s = 0; s < (n-(D-1)*tau); s++)
    {
        // getting the subsequence for the respective window
        j = 0;
        for(i = s; i <= (s+(D-1)*tau); i+=tau)
        {
            //Rprintf("%d ",i);
            w[j] = x[i];
            j++;
        }    

        //Rprintf("%d %f,%f,%f\n",s,w[0],w[1],w[2]);

        // copying the pattern to the symbol list
        symbols[s] = bandt_pompe_pattern_c(w);
    }

    return symbols;
}

//' Bandt-Pompe transformation (na-removal)
//'
//' Do not count symbols if all elements in w are NAN.
//'
//' @param x The time series (univariate vector)
//' @param D The embedding dimension (size of sliding window)
//' @param tau The embedding delay ('step' value)
//
// [[Rcpp::export]]
StringVector bandt_pompe_na_c(NumericVector x, int D=3, int tau=1)
{
    // aux
    int i, j, s;

    // the lenght of the vector
    int n = x.size();

    // the number of symbols
    int sym_num = n-(D-1)*tau;

    // to count or not symbols in case of narm=true
    bool all_nan;

    // the list of symbols to be returned
    StringVector symbols(sym_num);

    // the sliding window
    NumericVector w(D);

    // discovering the sequences of order n
    for (s = 0; s < (n-(D-1)*tau); s++)
    {
        // getting the subsequence for the respective window
        j = 0;
        // 
        all_nan = true;
        for(i = s; i <= (s+(D-1)*tau); i+=tau)
        {
            //Rprintf("%d ",i);
            w[j] = x[i];
        
            // any number of non-nan makes this symbols to bee count
            if (isnan(w[j]) == false)
            {
                all_nan = false;
                //Rprintf("%d isnan\n",s);
            } 
            
            j++;
        }

        // do not count this symbol, set as NA
        if (all_nan == true)
        {
            //Rprintf("%d is all nan\n",s);
            symbols[s] = NA_STRING;
        }
        else
        {
            //Rprintf("%d %f,%f,%f\n",s,w[0],w[1],w[2]);

            // copying the pattern to the symbol list
            symbols[s] = bandt_pompe_pattern_c(w);
        }
    }

    return symbols;
}


// Returns an empty list of symbols
// 
// The list is a map, to counting the patterns
//
// @param D The embedding dimension
// 
// [[Rcpp::export]]
std::map<std::string, int> bandt_pompe_empty_c(int D=3)
{
    std::map<std::string, int> perms;

    int i, itens[D];

    // initializing the list of itens
    for(i = 0; i < D; i++)
    {
        itens[i] = i+1;
    }

    String s;
    
    //auto start = high_resolution_clock::now(); 
    
    // create each permutation and initialize it with 0
    do {
        s = "";
        for(i = 0; i < D; i++)
        {
            s += std::to_string(itens[i]);
        }
        perms[s] = 0;
    } while ( std::next_permutation(itens,itens+D) );

    //auto stop = high_resolution_clock::now(); 
    //auto duration = duration_cast<microseconds>(stop - start); 
    //Rcout << "\t\tTIME EMPTY: " << duration.count() << " seconds" << endl;

    return perms;
}

//' Bandt-Pompe distribution
//'
//' @param x The time series (univariate vector)
//' @param D The embedding dimension (size of sliding window)
//' @param tau The embedding delay ('step' value)
// 
// [[Rcpp::export]]
std::map<std::string, int> 
bandt_pompe_distribution_c(NumericVector x, int D=3, int tau=1, bool na_aware=false)
{
    // create the empty permutations
    std::map<std::string, int> perms = bandt_pompe_empty_c(D);
    
    // computing the symbols
    StringVector symbols;
    if (na_aware == true)
    {
        symbols = bandt_pompe_na_c(x, D, tau);
    }
    else
    {
        symbols = bandt_pompe_c(x, D, tau);
    }

    // counting the patterns 
    for(int i = 0; i < symbols.size(); i++)
    {
        ++perms[(String)symbols[i]];
    }

    return perms;
}


//' Bandt-Pompe distribution (symbols)
//'
//' @param symbols The symbols already computed
//' @param D The embedding dimension (size of sliding window)
//' @param tau The embedding delay ('step' value)
// 
// [[Rcpp::export]]
std::map<std::string, int> 
bandt_pompe_distribution_symbols_c(StringVector symbols, int D=3, int tau=1)
{
    // create the empty permutations
    std::map<std::string, int> perms = bandt_pompe_empty_c(D);
    
    // symbols already computed

    // counting the patterns 
    for(int i = 0; i < symbols.size(); i++)
    {
        ++perms[(String)symbols[i]];
    }

    return perms;
}

/////////////////////////////
// transition graph functions
/////////////////////////////

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
    
    //auto start = high_resolution_clock::now(); 
    
    // create each permutation and initialize it with 0
    do {
        s = "";
        for(i = 0; i < D; i++)
        {
            s += std::to_string(itens[i]);
        }
        M[s] = bandt_pompe_empty_c(D);
    } while ( std::next_permutation(itens,itens+D) );
    
    //auto stop = high_resolution_clock::now(); 
    //auto duration = duration_cast<microseconds>(stop - start); 
    //Rcout << "\t\tTIME EMPTY MATRIX: " << duration.count() << " seconds" << endl;

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
bandt_pompe_transition_c(NumericVector x, int D=3, int tau=1, bool na_aware=false)
{
    // the transition matrix
    // NOTE: the matrix initialization as removed, and this adjustment is made in R
    //std::map<std::string, std::map<std::string, int>> M = bandt_pompe_empty_matrix_c(D);
    std::map<std::string, std::map<std::string, int>> M;

    // computing the symbols
    StringVector symbols;
    if (na_aware == true)
    {
        symbols = bandt_pompe_na_c(x, D, tau);
    }
    else
    {
        symbols = bandt_pompe_c(x, D, tau);
    }
    //StringVector symbols = bandt_pompe_c(x, D, tau);

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
    //auto start = high_resolution_clock::now(); 
    
    // the transition matrix
    // NOTE: the matrix initialization as removed, and this adjustment is made in R
    //std::map<std::string, std::map<std::string, int>> M = bandt_pompe_empty_matrix_c(D);
    std::map<std::string, std::map<std::string, int>> M;

    //auto stop = high_resolution_clock::now(); 
    //auto duration = duration_cast<microseconds>(stop - start); 
    //Rcout << "\t\tTIME MAP: " << duration.count() << " seconds" << endl;

    //auto start2 = high_resolution_clock::now(); 
    //
    // counting the transitions
    for(int i = 1; i < symbols.size(); i++)
    {
        ++M[(String)symbols[i-1]][(String)symbols[i]];
    }

    //auto stop2 = high_resolution_clock::now(); 
    //auto duration2 = duration_cast<microseconds>(stop2 - start2); 
    //Rcout << "\t\tTIME LOOP: " << duration2.count() << " seconds" << endl;

    return M;
}

