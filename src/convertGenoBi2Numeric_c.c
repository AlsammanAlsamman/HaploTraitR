#include "R.h"
#include "Rinternals.h"

SEXP convertGenoBi2Numeric_c (SEXP matrixOfGenotypes){
  int numberOfRows = nrows(matrixOfGenotypes);
  int numberOfColumns = ncols(matrixOfGenotypes);
  int i,j;
  // length of matrix
  int length = length(matrixOfGenotypes);
  // convert matrixOfGenotypes to char**
  char** matrixOfGenotypesChar = (char**)R_alloc(length,sizeof(char*));
  // memset the memory
  memset(matrixOfGenotypesChar,0,length*sizeof(char*));
  // copy the data
  for(int c=0;c<numberOfColumns;c++){
    for(int r=0;r<numberOfRows;r++){
      matrixOfGenotypesChar[r+c*numberOfRows] = (char*)R_alloc(2,sizeof(char));
      // copy the two characters
      for(int k=0;k<2;k++){
        matrixOfGenotypesChar[r+c*numberOfRows][k] = CHAR(STRING_ELT(matrixOfGenotypes,r+c*numberOfRows))[k];
      }
      // add the null character
      matrixOfGenotypesChar[r+c*numberOfRows][2] = '\0';
    }
  }
  // create matrix of 0,1,2 with the same size as matrixOfGenotypes
  SEXP matrixOfSNPsNumeric;
  PROTECT(matrixOfSNPsNumeric = allocMatrix(INTSXP,numberOfRows,numberOfColumns));
  // convert matrixOfGenotypes to numeric 0,1,2
  for(i=0;i<numberOfRows;i++){
    // we need to determine the allele A, which is the first allele in the first SNP in each row
    char alleleA ;
    for(j=0;j<numberOfColumns;j++){
      if(matrixOfGenotypesChar[i+j*numberOfRows][0] != '0'){
        alleleA = matrixOfGenotypesChar[i+j*numberOfRows][0];
        break;
      }
    }
    // convert matrixOfGenotypes to numeric 0,1,2
    // now with alleleA we can convert the genotypes to numeric 0,1,2
    // restart j
    for(j=0;j<numberOfColumns;j++){
      // if both alleles are alleleA
      if(matrixOfGenotypesChar[i+j*numberOfRows][0] == alleleA && matrixOfGenotypesChar[i+j*numberOfRows][1] == alleleA){
        // assign 0 in matrixOfSNPsNumeric,take care of the order of the indexes in matrix is by columns
        INTEGER(matrixOfSNPsNumeric)[i+j*numberOfRows] = 0;
      }
      // if both alleles are not the same
      else if(matrixOfGenotypesChar[i+j*numberOfRows][0] != matrixOfGenotypesChar[i+j*numberOfRows][1]){
        INTEGER(matrixOfSNPsNumeric)[i+j*numberOfRows] = 1;
      }
      // if its is missing data
      else if(matrixOfGenotypesChar[i+j*numberOfRows][0] == 'N' && matrixOfGenotypesChar[i+j*numberOfRows][1] == 'N'){
        INTEGER(matrixOfSNPsNumeric)[i+j*numberOfRows] = -1;
      } else{
        INTEGER(matrixOfSNPsNumeric)[i+j*numberOfRows] = 2;
      }

    }
  }
  UNPROTECT(1);
  return(matrixOfSNPsNumeric);
}

