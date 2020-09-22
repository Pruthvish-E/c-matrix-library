#include "mymatrix.h"
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <math.h>
//Creating a matrix of order = r*c with r*c elements
matrix * create_matrix(int c,int r,double * elements)
{
	int n = c*r;
	matrix * m = (matrix *)malloc(sizeof(matrix));
	m->c = c;
	m->r = r;
	m->m = (double **)malloc(sizeof(double *)*r);
	int k = 0;
	for(int i = 0;i<r;i++)
	{
		m->m[i] = (double*)malloc(sizeof(double)*c);
		for(int j = 0;j<c;j++)
		{
			m->m[i][j] = elements[k++];
		}
	}
	return m;
}

//Display a Matrix
void displayMat(const matrix *m)
{
	int r = m->r;
	int c = m->c;
	for(int i = 0;i<r;i++)
	{
		for(int j = 0;j<c;j++)
		{
			printf("%lf ",m->m[i][j]);
		}
		printf("\n");
	}
}

//Check matrices are compatible for multiplication returns 1 if possible else false
int multiply_compat(const matrix *m1,const matrix *m2)
{
	if(m1->c == m2->r)
	{
		return 1;
	}
	return 0;
}


//Multiply two compatible matrices
matrix * multiply_mat(const matrix* m1,const matrix* m2)
{
	matrix *m = (matrix *)malloc(sizeof(matrix));
	m->r = m1->r;
	m->c = m2->c;
	m->m = (double **)malloc(sizeof(double *)*(m1->r));
	if(!multiply_compat(m1,m2))
	{
		return NULL;
	}
	for(int i = 0;i<m1->r;i++)
	{
		m->m[i] = (double *)calloc(m2->c,sizeof(double));

	}


	for(int i=0; i<m1->r; ++i)
	{
        for(int j=0; j<m2->c; ++j)
        {
            for(int k=0; k<m1->c; ++k)
            {
                m->m[i][j]+=m1->m[i][k]*m2->m[k][j];
            }
        }
	}

	return m;

}


//Delete a matrix
matrix* destroy_matrix(matrix *m1)
{
	if(m1 == NULL)
	{
		return NULL;
	}
	int r = (m1)->r;
	if(m1->m == NULL)
	{
		free(m1);
		return NULL;
	}
	for(int i = 0;i<r;i++)
	{
		if((m1)->m[i] == NULL){
			free(m1);
			return NULL;
				
		} 
		free((m1)->m[i]);
	}
	free(m1);
	return NULL;
}


//Add two matrices m1 and m2 returns null if incompatible
matrix* add(const matrix *m1,const matrix *m2)
{
	
	if(m1->r != m2->r || m1->c != m2->c)
	{
		return NULL;
	}

	matrix * m = (matrix *)malloc(sizeof(matrix));
	m->r = m1->r;
	m->c = m1->c;

	m->m = (double **)malloc(sizeof(double *)*(m1->r));
	for(int i = 0;i<m1->r;i++)
	{
		m->m[i] = (double *)calloc(m2->c,sizeof(double));

	}

	for(int i = 0; i < m1->c; i++){
		for(int j = 0; j < m1->r; j++)
			m->m[i][j] = m1->m[i][j] +  m2->m[i][j];
	}
	return m;
}


//Subtract two matrices m1 and m2 returns null if incompatible
matrix* subtract(const matrix *m1,const matrix *m2)
{
	if(m1->r != m2->r || m1->c != m2->c)
	{
		return NULL;
	}

	matrix * m = (matrix *)malloc(sizeof(matrix));
	m->r = m1->r;
	m->c = m1->c;

	m->m = (double **)malloc(sizeof(double *)*(m1->r));
	for(int i = 0;i<m1->r;i++)
	{
		m->m[i] = (double *)calloc(m2->c,sizeof(double));

	}

	for(int i = 0; i < m1->c; i++){
		for(int j = 0; j < m1->r; j++)
			m->m[i][j] = m1->m[i][j] -  m2->m[i][j];
	}
	return m;
}


//Create Identitiy Matrix of order n
matrix *identity(int n)
{
	matrix * m = (matrix *)malloc(sizeof(matrix));
	m->r = n;
	m->c = n;

	m->m = (double **)malloc(sizeof(double *)*n);
	for(int i = 0;i<n;i++)
	{
		m->m[i] = (double *)calloc(n,sizeof(double));
		m->m[i][i] = 1;

	}
	return m;
}

//Create Zero Matrix of order r*c
matrix *zeroMatrix(int r, int c)
{
	matrix * m = (matrix *)malloc(sizeof(matrix));
	m->r = r;
	m->c = c;

	m->m = (double **)malloc(sizeof(double *)*r);
	for(int i = 0;i<r;i++)
	{
		m->m[i] = (double *)calloc(c,sizeof(double));

	}
	return m;
}


//Swap rows a and b of matrix m returns 1 on success
int row_swap(matrix *m, int a, int b)
{
	if(a>=m->r || a<0 || b>=m->r || b<0)
	{
		return 0;
	}

	double * temp = m->m[a];
	m->m[a] = m->m[b];
	m->m[b] = temp;
	return 1; 
}


//Swap columns a and b of matrix returns 1 on success
int col_swap(matrix *m, int a, int b)
{
	if(a>=m->c || a<0 || b>=m->c || b<0)
	{
		return 0;
	}

	double temp;

	for(int i = 0;i<m->r;i++)
	{
		temp = m->m[i][a];
		m->m[i][a] = m->m[i][b];
		m->m[i][b] = temp;
	}
	return 1;
}


//Multiply matrix by a Scalar returns 1 on success
int scalar_multiply(matrix *m, double f)
{
	if(m == NULL || m->m == NULL)
	{
		return 0;
	}
	for(int i = 0;i<m->r;i++)
	{
		for (int j = 0; j <m->c ; ++j)
		{
			m->m[i][j]*=f;
		}
	}
	return 1;
}


//Reduce row b by a factor*a times returns 1 on success
int reduceRow(matrix *m, int b, int a, double factor)
{
	if(b>=m->r || b<0 || a>=m->r || a<0)
	{
		return 0;
	}
	for(int i = 0;i<m->c;i++)
	{
		m->m[b][i] -= m->m[a][i]*factor;
	}
	return 1;
}


//Reduce column b by a factor*a times returns 1 on success
int reduceCol(matrix *m, int b, int a, double factor)
{
	if(b>=m->c || b<0 || a>=m->c || a<0)
	{
		return 0;
	}
	for (int i = 0; i < m->r; ++i)
	{
		m->m[i][b] -= m->m[i][a]*factor;	
	}
	return 1;
}


//Check for equality Of Two Matrices returns 1 if equal 0 otherwise
int equals(matrix *m1, matrix *m2)
{
	if(m1->r != m2->r || m1->c != m2->c)
	{
		return 0;
	}
	for(int i = 0;i<m1->r;i++)
	{
		for (int j = 0; j < m1->c; ++j)
		{
			if(m1->m[i][j] != m2->m[i][j])
			{
				return 0;
			}
		}
	}
	return 1;
}


//clone a matrix m
matrix *clone(const matrix *m)
{
	matrix * Mclone = (matrix *)malloc(sizeof(matrix));
	Mclone->r = m->r;
	Mclone->c = m->c;

	Mclone->m = (double **)malloc(sizeof(double *)*(m->r));
	for(int i = 0;i<(m->r);i++)
	{
		Mclone->m[i] = (double *)calloc(m->c,sizeof(double));
	}

	for (int i = 0; i < m->r; ++i)
	{
		for (int j = 0; j < m->c; ++j)
		{
			Mclone->m[i][j] = m->m[i][j];
		}
	}

	return Mclone;
}


//Transpose a matrix m
matrix *transpose(const matrix *m)
{
	matrix * Mclone = (matrix *)malloc(sizeof(matrix));
	Mclone->r = m->c;
	Mclone->c = m->r;

	Mclone->m = (double **)malloc(sizeof(double *)*(m->c));
	for(int i = 0;i<(m->c);i++)
	{
		Mclone->m[i] = (double *)calloc(m->r,sizeof(double));
	}

	for (int i = 0; i < m->r; ++i)
	{
		for (int j = 0; j < m->c; ++j)
		{
			Mclone->m[j][i] = m->m[i][j];
		}
	}

	return Mclone;
}


//Generate a random matrix
matrix *rand_matrix(int r, int c,float maxNo)
{
	matrix * m = (matrix *)malloc(sizeof(matrix));
	m->r = r;
	m->c = c;

	m->m = (double **)malloc(sizeof(double *)*r);
	for(int i = 0;i<r;i++)
	{
		m->m[i] = (double *)calloc(c,sizeof(double));
		for (int j = 0; j < c; ++j)
		{
			m->m[i][j] = ((float)rand()/(float)(RAND_MAX)) *(maxNo);
		}
	}
	return m;	
}

//Augment a matrix m1 with another matrix m2 returns NULL if incompatible
matrix *aug_matrix(const matrix *m1,const matrix *m2)
{
	if(m1->r != m2->r)
	{
		return NULL;
	}

	matrix * m = (matrix *)malloc(sizeof(matrix));
	m->r = m1->r;
	m->c = m1->c + m2->c;

	m->m = (double **)malloc(sizeof(double *)*(m1->r));
	for(int i = 0;i<m1->r;i++)
	{
		m->m[i] = (double *)calloc(m->c,sizeof(double));
		for (int j = 0; j < m->c; ++j)
		{
			if(j<m1->c)
			{
				m->m[i][j] = m1->m[i][j];	
			}
			else
			{
				m->m[i][j] = m2->m[i][j-m1->c];
			}
			
		}
	}
	return m;	
}


//To find Cofactor of a matrix returns NULL if not square
matrix * getCofactor(const matrix * m1,int p ,int q,int n)
{
	if(m1->r != m1->c)
	{
		return NULL;
	}
    int i = 0, j = 0; 
  	matrix * m = (matrix *)malloc(sizeof(matrix));
	m->r = n;
	m->c = n;
	
	m->m = (double **)malloc(sizeof(double *)*n);
	for(int i = 0;i<n;i++)
	{
		m->m[i] = (double *)calloc(n,sizeof(double));
	}

    for (int row = 0; row < n; row++) 
    { 
        for (int col = 0; col < n; col++) 
        { 
             
            if (row != p && col != q) 
            { 
                m->m[i][j++] = m1->m[row][col]; 
  
               	if (j == n - 1) 
                { 
                    j = 0; 
                    i++; 
                } 
            } 
        } 
    }
    return m; 

} 


//Find Determinant of a matrix m
double determinant(const matrix *m,int n)
{
	double D = 0; 
    if (n == 1) 
        return m->m[0][0]; 
  
    matrix* temp;  
  
    int sign = 1;   
  
      
    for (int f = 0; f < n; f++) 
    { 
        // Getting Cofactor of A[0][f] 
        temp = getCofactor(m, 0, f,n); 
        D += sign * m->m[0][f] * determinant(temp, n - 1); 
  
        // terms are to be added with alternate sign 
        sign = -sign; 
    } 
  
    return D; 
}


//To find adjoint of a matrix
matrix * adjoint(const matrix * m1)
{
	int n = m1->r;
	matrix * m = (matrix *)malloc(sizeof(matrix));
	m->r = n;
	m->c = n;
	
	m->m = (double **)malloc(sizeof(double *)*n);
	for(int i = 0;i<n;i++)
	{
		m->m[i] = (double *)calloc(n,sizeof(double));
	}
	if(n == 1)
	{
		m->m[0][0] = 1;
		return m;
	}

	int sign = 1;
	matrix * temp; 
  
    for (int i=0; i<n; i++) 
    { 
        for (int j=0; j<n; j++) 
        { 
            temp = getCofactor(m1, i, j, n); 
  
            // sign of adj[j][i] positive if sum of row 
            // and column indexes is even. 
            sign = ((i+j)%2==0)? 1: -1; 
  
            m->m[j][i] = (sign)*(determinant(temp, n-1)); 
        } 
    }
    return m;

}

//To find Cofactor of a matrix
matrix * cofactor(const matrix * m1)
{
	matrix *adj = adjoint(m1);
	return transpose(adj);
}


//Find Inverse of the Matrix returns NULL if matrix is singular
matrix *inverse(const matrix *m1)
{

	int n = m1->r;
	// Find determinant of A[][] 
    double det = determinant(m1, n); 
    if (det == 0) 
    { 
        return NULL; 
    } 
  	matrix * m = (matrix *)malloc(sizeof(matrix));
	m->r = n;
	m->c = n;
	
	m->m = (double **)malloc(sizeof(double *)*n);
	for(int i = 0;i<n;i++)
	{
		m->m[i] = (double *)calloc(n,sizeof(double));
	}
    // Find adjoint 
    matrix * adj = adjoint(m1); 
  
    // Find Inverse using formula "inverse(A) = adj(A)/det(A)" 
    for (int i=0; i<n; i++) 
        for (int j=0; j<n; j++) 
            m->m[i][j] = adj->m[i][j]/det; 
  
    return m; 
}


//minimum matrix Chain Multiplication for a system of compatible matrices
int recursiveMinOrder(int *orderDim,int start,int end)
{
	if(start == end) 
        return 0; 
    int k; 
    int min = INT_MAX; 
    int count; 
  
  	for (k = start; k <end; k++) 
    { 
        count = recursiveMinOrder(orderDim, start, k) + 
                recursiveMinOrder(orderDim, k+1, end) + 
                orderDim[start-1]*orderDim[k]*orderDim[end];

  		
        if (count < min){ 
            min = count; 

    } 
  
    // Return minimum count 
    return min; 
}

int minmatrixOrder(int *orderDim,int n)
{
	int * order = (int *)malloc(sizeof(int)*n);

	int ans =  recursiveMinOrder(orderDim,1,n-1,order);

}



//Gaussian elimination to solve linear system equations of a augmented matrix returns NULL if no solution unique exists
double* GaussianELimationSolutions(matrix * aug)
{

    
    int singular_flag = forwardElimation(aug); 
  
    
    if (singular_flag != -1) 
    { 
        
        return NULL; 
    } 
  
    
    return backSubstitution(aug); 
}

//forward Elimination of a augmented matrix m
int forwardElimation(matrix *aug)
{
	int N = aug->r;
	for (int k=0; k<N; k++) 
    { 
        int i_max = k; 
        int v_max = aug->m[i_max][k]; 
  
        for (int i = k+1; i < N; i++) 
            if (abs(aug->m[i][k]) > v_max) 
                v_max = aug->m[i][k], i_max = i; 
  
        if (!aug->m[k][i_max]) 
            return k; // Matrix is singular 
  
        if (i_max != k) 
            row_swap(aug, k, i_max); 
  
  
        for (int i=k+1; i<N; i++) 
        { 
            double f = aug->m[i][k]/aug->m[k][k]; 
  
            for (int j=k+1; j<=N; j++) 
                aug->m[i][j] -= aug->m[k][j]*f; 
  
            aug->m[i][k] = 0; 
        } 
  
        //print(mat);        //for matrix state 
    } 
    //print(mat);            //for matrix state 
    return -1; 
}

//backward substitution of the matrix to find the solution returns NULL if no solution unique exists
double * backSubstitution(matrix *m)
{
	int N = m->r;
	 double* x = (double *)malloc(sizeof(double)*N);  // An array to store solution 
  
    for (int i = N-1; i >= 0; i--) 
    { 
        x[i] = m->m[i][N]; 
  
        for (int j=i+1; j<N; j++) 
        { 
            x[i] -= m->m[i][j]*x[j]; 
        } 
  
        x[i] = x[i]/m->m[i][i]; 
    } 
  	double* result = (double *)malloc(sizeof(double)*N);  // An array to store solution 
  	
    return x;
}



//Find eigenvalues of a matrix m 
double *eigenvalues(matrix *m)
{
	double *values, factor;
	matrix *red;
	unsigned int i, j, l;
	if(m == NULL)
		return NULL;
	if(m->r != m->c)
		return NULL;
	values = (double *) malloc(sizeof(double)*m->r);
	red = clone(m);
	/* reduce each of the rows to get a lower triangle */	
	for(i = 0; i < red->c; i++){
		for(j = i + 1; j < red->r; j++){
			if(red->m[i][i] == 0){
				for(l = i+1; l < red->r; l++){
					if(red->m[l][l] != 0){
						row_swap(red, i, l);
						break;
					}
				}
				continue;
			}
			factor = red->m[i][j]/(red->m[i][i]);
			reduceRow(red, j, i, factor);
		}
	}
	for(i = 0; i < red->c; i++)
		values[i] = red->m[i][i];
	return values;
}


//Finding minimum eigen value and corresponding eigen vector using Rayleighs power method
matrix * eigenVectorPower(matrix *m)
{
	int i,j,k;
	int n = m->r;
	int len =m->r;
	float x[len],z[len],e[len],zmax,emax;
	do
    {
        for(i=0; i<n; i++)
        {
            z[i]=0;
            for(j=0; j<n; j++)
            {
                z[i]=z[i]+m->m[i][j]*x[j];
            }
        }
        zmax=fabs(z[1]);
        for(i=1; i<n; i++)
        {
            if((fabs(z[i]))>zmax)
                zmax=fabs(z[i]);
        }
        for(i=0; i<n; i++)
        {
            z[i]=z[i]/zmax;
        }
        for(i=0; i<n; i++)
        {
            e[i]=0;
            e[i]=fabs((fabs(z[i]))-(fabs(x[i])));
        }
        emax=e[1];
        for(i=1; i<n; i++)
        {
            if(e[i]>emax)
                emax=e[i];
        }
        for(i=0; i<n; i++)
        {
            x[i]=z[i];
        }
    }while(emax>0.001);

    matrix * result = (matrix *)malloc(sizeof(matrix));
    result->r = n;
    result->c = 1;
    result->m = (double **)malloc(sizeof(double *)*n);
    
    for(i=0; i<n; i++)
    {
    	result->m[i] = (double *)malloc(sizeof(double *)*1);
        result->m[i][0] = z[i];
    }
    return result;
}



//LU factorization of a matrix m parameters are m and two matrices o same order of m to store L and U
void LUfactorize(matrix *a,matrix *l,matrix * u)
{
	int n = a->r;
	int i = 0, j = 0, k = 0;
   for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
         if (j < i)
         l->m[j][i] = 0;
         else {
            l->m[j][i] = a->m[j][i];
            for (k = 0; k < i; k++) {
               l->m[j][i] = l->m[j][i] - l->m[j][k] * u->m[k][i];
            }
         }
      }
      for (j = 0; j < n; j++) {
         if (j < i)
            u->m[i][j] = 0;
         else if (j == i)
            u->m[i][j] = 1;
         else {
            u->m[i][j] = a->m[i][j] / l->m[i][i];
            for (k = 0; k < i; k++) {
               u->m[i][j] = u->m[i][j] - ((l->m[i][k] * u->m[k][j]) / l->m[i][i]);
            }
         }
      }
   }
 		return;
}




//Strassen's matrix multiplication returns null in case matrix is not 2*2
matrix * strassen_mul(const matrix *m1,const matrix* m2)
{
	if(m1->r != 2 || m1->c != 2 || m2->r != 2 || m2->c != 2)
	{
		return NULL;
	}


  matrix * result = (matrix *)malloc(sizeof(matrix));
  result->r = m1->r;
  result->c = m2->c;

  result->m = (double **)malloc(sizeof(double *)*m1->r);
    
    for(int i=0; i<m1->r; i++)
    {
    	result->m[i] = (double *)malloc(sizeof(double *)*(m2->c));
    }


  double a= (m1->m[0][0] + m1->m[1][1]) * (m2->m[0][0] + m2->m[1][1]);

  double b= (m1->m[1][0] + m1->m[1][1]) * m2->m[0][0];

  double c= m1->m[0][0] * (m2->m[0][1] - m2->m[1][1]);

  double d= m1->m[1][1] * (m2->m[1][0] - m2->m[0][0]);

  double e= (m1->m[0][0] + m1->m[0][1]) * m2->m[1][1];

  double f= (m1->m[1][0] - m1->m[0][0]) * (m2->m[0][0]+m2->m[0][1]);

  double g= (m1->m[0][1] - m1->m[1][1]) * (m2->m[1][0]+m2->m[1][1]);

 

  result->m[0][0] = a + d- e + g;

  result->m[0][1] = c + e;

  result->m[1][0] = b + d;

  result->m[1][1] = a - b + c + f;

  return result;
}


