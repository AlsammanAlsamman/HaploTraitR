#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double calculate_ld_between_two_snps(IntegerVector snp1, IntegerVector snp2) {
    if (snp1.size() != snp2.size()) {
        stop("SNP vectors must be of the same length.");
    }

    // Variables to hold counts
    long n = 0, naa = 0, naA = 0, nAA = 0;
    long nbb = 0, nbB = 0, nBB = 0;
    long nAABB = 0, naabb = 0, naaBB = 0, nAAbb = 0;

    // Iterate over SNP pairs, skipping NA values (3)
    for (int i = 0; i < snp1.size(); i++) {
        int val1 = snp1[i];
        int val2 = snp2[i];

        // Skip if either SNP value is NA (3)
        if (val1 == 3 || val2 == 3) continue;

        n++;  // Count valid pairs

        // Update counts based on SNP values
        if (val1 == 0) naa++;
        else if (val1 == 1) naA++;
        else if (val1 == 2) nAA++;

        if (val2 == 0) nbb++;
        else if (val2 == 1) nbB++;
        else if (val2 == 2) nBB++;

        // Update joint counts
        if (val1 == 2 && val2 == 2) nAABB++;
        if (val1 == 0 && val2 == 0) naabb++;
        if (val1 == 0 && val2 == 2) naaBB++;
        if (val1 == 2 && val2 == 0) nAAbb++;
    }

    // Check if there are enough data points
    if (n == 0) return NA_REAL;

    // Calculate allele frequencies
    double pa = (2.0 * naa + naA) / (2.0 * n);
    double pb = (2.0 * nbb + nbB) / (2.0 * n);
    double pA = 1.0 - pa;
    double pB = 1.0 - pb;

    // Calculate delta
    double delta = (nAABB + naabb - naaBB - nAAbb) / (2.0 * n) -
                   (double)(naa - nAA) * (nbb - nBB) / (2.0 * n * n);

    // Calculate variances for each allele
    double DA = ((double)nAA / n) - pA * pA;
    double DB = ((double)nBB / n) - pB * pB;

    // Calculate r^2 (LD measure)
    double denom = (pA * pa + DA) * (pB * pb + DB);
    if (denom > 0) {
        return delta * delta / denom;
    }

    return NA_REAL;  // return NA if calculation fails
}

// [[Rcpp::export]]
NumericMatrix calculate_ld_matrix(IntegerMatrix snp_matrix) {
    int num_snps = snp_matrix.nrow();
    NumericMatrix ld_matrix(num_snps, num_snps);

    // Calculate LD for each unique SNP pair
    for (int i = 0; i < num_snps; i++) {
        for (int j = i + 1; j < num_snps; j++) {
            IntegerVector snp1 = snp_matrix(i, _);  // Row i
            IntegerVector snp2 = snp_matrix(j, _);  // Row j

            // Calculate LD between SNPs i and j
            double ld_value = calculate_ld_between_two_snps(snp1, snp2);
            ld_matrix(i, j) = ld_value;
            ld_matrix(j, i) = ld_value;  // Symmetric matrix
        }
    }

    // Diagonal elements should be 1 as they are perfectly correlated with themselves
    for (int i = 0; i < num_snps; i++) {
        ld_matrix(i, i) = 1.0;
    }

    return ld_matrix;
}

