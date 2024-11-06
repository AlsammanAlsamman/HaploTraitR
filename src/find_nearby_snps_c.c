#include <R.h>
#include <Rinternals.h>

// Function to get SNPs within threshold distance
SEXP find_nearby_snps_c(SEXP list1, SEXP list2, SEXP threshold) {
  // Coerce threshold to an integer
  int thresh = INTEGER(threshold)[0];

  // Get the length of the lists
  int len1 = LENGTH(list1);
  int len2 = LENGTH(list2);

  // Initialize an R list to store the results
  SEXP result = PROTECT(allocVector(VECSXP, len1));

  // Iterate over each SNP in list1
  for (int i = 0; i < len1; i++) {
    int pos1 = INTEGER(list1)[i];

    // Initialize a temporary list to collect nearby SNPs from list2
    SEXP nearby_snps = PROTECT(allocVector(INTSXP, len2));
    int count = 0;

    // Iterate over each SNP in list2
    for (int j = 0; j < len2; j++) {
      int pos2 = INTEGER(list2)[j];

      // Check if the absolute distance is within the threshold
      if (abs(pos1 - pos2) <= thresh) {
        INTEGER(nearby_snps)[count++] = pos2;
      }
    }

    // Resize nearby_snps to the actual count of nearby SNPs
    SET_VECTOR_ELT(result, i, lengthgets(nearby_snps, count));

    // Release the temporary nearby_snps list
    UNPROTECT(1);
  }

  // Release the main result list
  UNPROTECT(1);
  return result;
}
