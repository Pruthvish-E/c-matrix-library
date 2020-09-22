#include <stdlib.h>
#include "mymatrix.h"
#include <stdio.h>

int main(int argc, char const *argv[])
{
	//Matrix 1
	int r,c,n;
	printf("Enter the number of rows and columns of the matrix1 :-\n");
	scanf("%d%d",&r,&c);
	printf("Enter the elements in row major fashion :- \n");
	double * elements = (double *)malloc(sizeof(double)*(r*c));
	for(int i = 0;i<r*c;i++)
	{
		scanf("%lf",&elements[i]);
	}
	matrix * m1  = create_matrix(c,r,elements);
	displayMat(m1);

	//Matrix 2
	int r1,c1,n1;
	printf("Enter the number of rows and columns of the matrix2 :-\n");
	scanf("%d%d",&r1,&c1);
	printf("Enter the elements in row major fashion :- \n");
	double * elements1 = (double *)malloc(sizeof(double)*(r*c));
	for(int i = 0;i<r1*c1;i++)
	{
		scanf("%lf",&elements1[i]);
	}
	matrix * m2  = create_matrix(c1,r1,elements1);
	displayMat(m2);

	//Result matrix
	// printf("\nResult after multiplication :-\n");
	// matrix * result = multiply_mat(m1,m2);

	// if(result == NULL)
	// {
	// 	printf("\nIncompatible matrices\n");
	// }
	// else
	// {

	// 	displayMat(result);
	// }
	// //destroy_matrix(&result);
	// result = destroy_matrix(result);
	

	// printf("\nAfter additon :- \n");
	// result = add(m1,m2);
	// if(result)
	// {
	// 	displayMat(result);
	// }
	// else
	// {
	// 	printf("\nIncompatible matrices\n");
	// }
	// //destroy_matrix(&result);
	// result = destroy_matrix(result);
		
	// printf("\nAfter subtraction :- \n");
	// result = subtract(m1,m2);
	
	// if(result)
	// {
	// 	displayMat(result);
	// }
	// else
	// {
	// 	printf("\nIncompatible matrices\n");
	// }



	// result = destroy_matrix(result);
	// printf("\nidentity Matrix of order 4:-\n");
	// result = identity(4);
	// displayMat(result);


	// printf("\nzeroMatrix 3*2:-\n");
	// result = destroy_matrix(result);
	// result = zeroMatrix(3,2);
	// displayMat(result);

	// printf("\nafter transpose m1 :-\n");
	// result = transpose(m1);
	// displayMat(result);

	// printf("\nAfter row swap 1 and 2 :-\n");
	// row_swap(m1,1,2);
	// displayMat(m1);

	// printf("\nAfter col swap 1 and 2 :-\n");
	// col_swap(m2,1,2);
	// displayMat(m2);

	int orderDim[] = {1, 2, 3, 4, 3}; 
	matrix * l = identity(3);
	matrix * u = clone(m1);
	

	LUfactorize(m1,l,u);
	printf("\nsolutions :- L :- \n");
	displayMat(l);

	printf("\nsolutions :- U :- \n");
	displayMat(u);
	
	// for (int i = 0; i < m1->r; ++i)
	// {
	// 	printf("%lf\n",ans[i]);
	// }
	// printf("\n");
	//result3 = destroy_matrix(result3);
	//m1 = destroy_matrix(m1);
	//m2 = destroy_matrix(m2);
	return 0;
}



