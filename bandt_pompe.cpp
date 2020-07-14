#include <Rcpp.h>
using namespace Rcpp;

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

// TODO: needs to return a string vector of symbols
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

