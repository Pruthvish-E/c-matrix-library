typedef struct matrix{
	int c;
	int r;
	double **m;
}matrix;


//Creating a matrix of order = r*c with r*c elements
matrix * create_matrix(int c,int r,double * elements); // DONE

//Display a Matrix
void displayMat(const matrix *m); //DONE


//Check matrices are compatible for multiplication returns 1 if possible else false
int multiply_compat(const matrix *m1,const matrix* m2); //DONE

//Multiply two compatible matrices returns NULL if the matrices are not compatible
matrix * multiply_mat(const matrix* m1,const matrix* m2); //DONE

//Delete a matrix
matrix* destroy_matrix(matrix *m1); //DONE


//Add two matrices m1 and m2 returns null if incompatible
matrix* add(const matrix *m1,const matrix *m2); //DONE

//Subtract two matrices m1 and m2 returns null if incompatible
matrix* subtract(const matrix *m1,const matrix *m2); //DONE


//Create Identitiy Matrix of order n
matrix *identity(int n); //DONE



//Create Zero Matrix of order r*c
matrix *zeroMatrix(int r, int c); //DONE


//Swap rows a and b of matrix m returns 1 on success
int row_swap(matrix *m, int a, int b); //DONE


//Swap columns a and b of matrix returns 1 on success
int col_swap(matrix *m, int a, int b); //DONE


//Multiply matrix by a Scalar returns 1 on success
int scalar_multiply(matrix *m, double f); //DONE


//Reduce row b by a factor*a times returns 1 on success
int reduceRow(matrix *m, int b, int a, double factor); //DONE


//Reduce column b by a factor*a times returns 1 on success
int reduceCol(matrix *m, int b, int a, double factor); //DONE


//Check for equality Of Two Matrices returns 1 if equal 0 otherwise
int equals(matrix *m1, matrix *m2); //DONE


//clone a matrix m
matrix *clone(const matrix *m);  //DONE


//Transpose a matrix m
matrix *transpose(const matrix *m); //DONE


//Generate a random matrix
matrix *rand_matrix(int r, int c,float maxNo); //DONE


//Augment a matrix m1 with another matrix m2 returns NULL if incompatible
matrix *aug_matrix(const matrix *m1,const matrix * m2); //DONE


//To find Cofactor of a matrix
matrix * cofactor(const matrix * m1); //DONE



//To find adjoint of a matrix
matrix * adjoint(const matrix * m1); //DONE



//Find Determinant of a matrix m
double determinant(const matrix *m,int n); //DONE

//Find Inverse of the matrix returns NULL if matrix is singular

matrix *inverse(const matrix *m1); //DONE


//minimum matrix Chain Multiplication for a system of compatible matrices

int minmatrixOrder(int *orderDim,int n); //DONE
int recursiveMinOrder(int *orderDim,int start,int end); //DONE


//Gaussian elimination to solve linear system equations of a augmented matrix returns NULL if no solution unique exists
double* GaussianELimationSolutions(matrix * m); //DONE



//forward Elimination of a augmented matrix m
int forwardElimation(matrix *m); //DONE


//backward substitution of the matrix to find the solution returns NULL if no solution unique exists
double * backSubstitution(matrix *m); //DONE



//Find eigenvalues of a matrix m
double *eigenvalues(matrix *m); //DONE


//Finding minimum eigen value and corresponding eigen vector using Rayleighs power method
matrix * eigenVectorPower(matrix *m); //DONE



//LU factorization of a matrix m parameters are m and two matrices o same order of m to store L and U
void LUfactorize(matrix *m,matrix *L,matrix * U); //DONE


//Strassen's matrix multiplication 
matrix * strassen_mul(const matrix *m1,const matrix* m2); //DONE
