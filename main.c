#include <stdlib.h>
#include "mymatrix.h"
#include <stdio.h>
#include <limits.h>
#include <math.h>

int main(int argc, char const *argv[])
{

    int ch;
    int r,c,n;
    double f;
    int a,b;
    matrix * m3;
    double * elements;
    matrix *result;
    matrix * U;
    matrix * L;
    int loop = 1;
    char proceed;
    while(loop){
        printf(" 1. create a matrix1\n 2. create a matrix2\n 3. display a matrix\n 4. multiply two matrices\n 5. add two matrices\n 6. subract two matrices\n 7. create an identity matrix\n 8. create a zero matrix\n 9. swap rows a and b of a matrix\n 10. swap columns a and b of a matrix\n 11. multiply a matrix with a scalar\n 12.Reduce row b by a factor*a times\n 13. find determinant of a matrix2\n 14.Reduce column b by a factor*a\n");
        printf(" 15. compare two matrices for equality\n 16. clone matrix1\n 17. transpose a matrix1 \n 18. generate a random matrix \n 19. augment matrix1 with matrix2 \n 20. find cofactor of a matrix1 \n 21. find adjoint of the matrix1 \n 22. find determinant of matrix1 \n 23. find inverse of matrix1 \n 24. find minimum matrix multiplications to multiply series of compatible matrices \n");
        printf(" 25. Find solutions to linear equations using guassian elimination\n 26. find eigen values of a matrix \n 27. find eigen vector of minimum eigen value of a matrix using power method \n 28. LU factorize matrix1 \n 29. apply strassens matrix multiplication \n 30. destroy matrices \n 31. abort\n");
        printf("\n\nenter choice: ");
        scanf("%d",&ch);
        switch(ch)
        {
        case 1:
            printf("enter the number of rows and columns of the matrix1: \n");
            scanf("%d%d",&r,&c);
            printf("enter the elements in row major fashion: \n");
            elements = (double *)malloc(sizeof(double)*(r*c));
            for(int i = 0;i<r*c;i++)
            {
                scanf("%lf",&elements[i]);
            }
            matrix* m1 = create_matrix(c,r,elements);
            printf("Enter p to proceed\n");
            while((proceed = getchar()) != 'p')
            {
            	continue;
            }
            break;

        case 2:
            printf("enter the number of rows and columns of the matrix2: \n");
            scanf("%d%d",&r,&c);
            printf("enter the elements in row major fashion: \n");
            elements = (double *)malloc(sizeof(double)*(r*c));
            for(int i = 0;i<r*c;i++)
            {
                scanf("%lf",&elements[i]);
            }
            matrix* m2 = create_matrix(c,r,elements);
			printf("Enter p to proceed\n");
            while((proceed = getchar()) != 'p')
            {
            	continue;
            }

            break;

        case 3: printf("1. display matrix1\n 2. display matrix2\n");
            int ch_display;
            printf("enter choice :");
            scanf("%d",&ch_display);
            switch(ch_display){
                case 1: displayMat(m1); break;
                case 2: displayMat(m2); break;
                default : printf("invalid choice\n");
            }
            printf("Enter p to proceed\n");
            while((proceed = getchar()) != 'p')
            {
            	continue;
            }
            break;

        case 4:
        		if(multiply_compat(m1,m2)){
                result = multiply_mat(m1,m2);
                
                    displayMat(result);
                }
                else{
                    printf("incompatible matrices\n");
                }
                destroy_matrix(result);
                printf("Enter p to proceed\n");
            while((proceed = getchar()) != 'p')
            {
            	continue;
            }
                break;

        case 5:
                result = add(m1,m2);
                if(result){
                    displayMat(result);
                }
                else{
                    printf("incompatible matrices\n");
                }
                destroy_matrix(result);
                printf("Enter p to proceed\n");
            while((proceed = getchar()) != 'p')
            {
            	continue;
            }
                break;

        case 6:
                result = subtract(m1,m2);
                if(result){
                    displayMat(result);
                }
                else{
                    printf("incompatible matrices\n");
                }
                destroy_matrix(result);
                printf("Enter p to proceed\n");
            while((proceed = getchar()) != 'p')
            {
            	continue;
            }
                break;

        case 7: printf("enter the order of the matrix: ");
                scanf("%d",&n);
                result = identity(n);
                displayMat(result);
                destroy_matrix(result);
               printf("Enter p to proceed\n");
            while((proceed = getchar()) != 'p')
            {
            	continue;
            }
                break;

        case 8: printf("enter the number of rows and columns of the zero matrix: ");
                scanf("%d %d",&r,&c);
                result = zeroMatrix(r,c);
                displayMat(result);
                destroy_matrix(result);
               printf("Enter p to proceed\n");
            while((proceed = getchar()) != 'p')
            {
            	continue;
            }
                break;

        case 9: printf("Enter rows a and b to swap :-");
        		
        		scanf("%d%d",&a,&b);
        		if(row_swap(m1,a,b))
        		{
        			printf("after swap :-\n");
        			displayMat(m1);
        		}
        		else printf("Invalid inputs\n");
        		printf("Enter p to proceed\n");
            while((proceed = getchar()) != 'p')
            {
            	continue;
            }
        		break;

        case 10: printf("Enter columns a and b to swap :-");
        	
        		scanf("%d%d",&a,&b);
        		if(col_swap(m1,a,b))
        		{
        			printf("after swap :-\n");
        			displayMat(m1);
        		}
        		else printf("Invalid inputs\n");
        	printf("Enter p to proceed\n");
            while((proceed = getchar()) != 'p')
            {
            	continue;
            }
        		break;


        case 11: printf("enter a scalar :");
                 scanf("%lf",&f);
                 int res = scalar_multiply(m1,f);
                 if(res){
                    displayMat(m1);
                 }
                 else{
                    printf("matrix does not exist\n");
                 }
             printf("Enter p to proceed\n");
            while((proceed = getchar()) != 'p')
            {
            	continue;
            }
                 break;

        case 12: printf("enter row a and b and factor :- ");
        		 int a,b;
        		 double factor;
        		 scanf("%d%d%lf",&a,&b,&factor);
        		 if(reduceRow(m1,b,a, factor))
        		 {
        		 	printf("After row reduction :- \n");
        		 	displayMat(m1);
        		 }
        		 else printf("Invalid inputs\n");
        	printf("Enter p to proceed\n");
            while((proceed = getchar()) != 'p')
            {
            	continue;
            }
        		 break;
       	
       	case 13:printf("Determinant of the matrix2 :-\n");
       			if(m2->r == m2->c)
       			printf("%lf\n",determinant(m2,m2->r));
       			else printf("not a square matrix\n");
       		printf("Enter p to proceed\n");
            while((proceed = getchar()) != 'p')
            {
            	continue;
            }
       			break; 
       	
       	case 14: 
       		printf("enter columnn a and b and factor :- ");
        		
        		 double factor1;
        		 scanf("%d%d%lf",&a,&b,&factor1);
        		 if(reduceCol(m1,b,a,factor1))
        		 {
        		 	printf("After column reduction :- \n");
        		 	displayMat(m1);
        		 }
        		 else printf("Invalid inputs\n");
        	printf("Enter p to proceed\n");
            while((proceed = getchar()) != 'p')
            {
            	continue;
            }
        		 break;


       	case 15:
       		if(equals(m1,m2))
       		{
       			printf("Equal\n");
       		}
       		else printf("Not equal\n");
       		printf("Enter p to proceed\n");
            while((proceed = getchar()) != 'p')
            {
            	continue;
            }
       		break;

       	case 16:
       		m3 = clone(m1);
       		printf("Matrix m3 after cloning:-\n");
       		displayMat(m3);
       		printf("Enter p to proceed\n");
            while((proceed = getchar()) != 'p')
            {
            	continue;
            }
       		break;

       	case 17:
       		printf("transpose of m1 :- \n");
       		matrix *m4 = transpose(m1);
       		displayMat(m4);
       		printf("Enter p to proceed\n");
            while((proceed = getchar()) != 'p')
            {
            	continue;
            }
       		break;

       	case 18:
       		printf("Enter number of rows and columns and max limit:-");
       		int r1,c1;
       		double limit;
       		scanf("%d%d%lf",&r1,&c1,&limit);
       		matrix * m5 = rand_matrix(r1,c1,limit);
       		printf("Generated random matrix:- \n");
       		displayMat(m5);
       		printf("Enter p to proceed\n");
            while((proceed = getchar()) != 'p')
            {
            	continue;
            }
       		break;

       	case 19:
       		m3 = aug_matrix(m1,m2);
       		if(m3)
       		{
       			printf("After augmenting :-\n");
       			displayMat(m3);
       		}
       		else
       		{
       			printf("Not compatible :- \n");
       		}
       		printf("Enter p to proceed\n");
            while((proceed = getchar()) != 'p')
            {
            	continue;
            }
       		break;

       	case 20:
       		if(m1->r != m1->c) 
       		{
       			printf("Not a square matrix\n");
       			break;
       		}
       		matrix * m7 = cofactor(m1);
       		printf("cofactor of matrix m1:- \n");
       		displayMat(m7);
       		printf("Enter p to proceed\n");
            while((proceed = getchar()) != 'p')
            {
            	continue;
            }
       		break;

       case 21:
       		if(m1->r != m1->c) 
       		{
       			printf("Not a square matrix\n");
       			break;
       		}
       		matrix * m8 = adjoint(m1);
       		printf("adjoint of matrix m1:- \n");
       		displayMat(m8);
       		printf("Enter p to proceed\n");
            while((proceed = getchar()) != 'p')
            {
            	continue;
            }
       		break;
      
      case 22:
      	printf("Determinant of the matrix1 :-\n");
       			if(m1->r == m1->c)
       			printf("%lf\n",determinant(m1,m1->r));
       			else printf("not a square matrix\n");
       			printf("Enter p to proceed\n");
            while((proceed = getchar()) != 'p')
            {
            	continue;
            }
       			break;

      case 23:
      		if(m1->r != m1->c) 
       		{
       			printf("Not a square matrix\n");
       			break;
       		}
       		matrix * m9 = inverse(m1);
       		printf("inverse of matrix1 :- \n");
       		displayMat(m9);
       		printf("Enter p to proceed\n");
            while((proceed = getchar()) != 'p')
            {
            	continue;
            }
       		break;

      case 24:
      		printf("Enter number of matrices to multiply :- \n");
      		int no;
      		scanf("%d",&no);
      		int * orderDim = (int *)malloc(sizeof(int)*(no+1));
      		printf("Enter order of matrices in sequential order :-\n");
      		for (int i = 0; i <=no; ++i)
      		{
      			scanf("%d",&orderDim[i]);
      		}
      		int ans1 = minmatrixOrder(orderDim,no);
      		printf("minimum matrix multiplications :- %d\n",ans1);
      		printf("Enter p to proceed\n");
            while((proceed = getchar()) != 'p')
            {
            	continue;
            }
      		break;

     case 25:
     		printf("Solutions to the system of linear equation :- \n");
     		double * sols = GaussianELimationSolutions(m1);
     		printf("Solutions are :- \n");
     		for (int i = 0; i < m1->c; ++i)
     		{
     			printf("%lf\n",sols[i]);
     		}
     		printf("Enter p to proceed\n");
            while((proceed = getchar()) != 'p')
            {
            	continue;
            }
     		break;

     case 26:

     		printf("Eigen values :- \n");
     		double * sols2 = eigenvalues(m1);
     		for (int i = 0; i < m1->c; ++i)
     		{
     			printf("%lf\n",sols2[i]);
     		}
     		printf("Enter p to proceed\n");
            while((proceed = getchar()) != 'p')
            {
            	continue;
            }
     		break;

     case 27:

     		printf("minimum eigen vector :-\n");
     		matrix * m100 = eigenVectorPower(m1);

     		displayMat(m100);
     		printf("Enter p to proceed\n");
            while((proceed = getchar()) != 'p')
            {
            	continue;
            }
     		break;

     case 28:

     	L = identity(m1->r);
     	U = identity(m1->r);

     	LUfactorize(m1,L,U);
     	printf("Matrix L :-\n");
     	displayMat(L);
     	printf("Matrix U :- \n");
     	displayMat(U);
     	printf("Enter p to proceed\n");
            while((proceed = getchar()) != 'p')
            {
            	continue;
            }
     	break;

     case 29:
     	printf("result of strassens matrix multiplication :- \n");
     	matrix *m200 = strassen_mul(m1,m2);
     	displayMat(m200);
     	printf("Enter p to proceed\n");
            while((proceed = getchar()) != 'p')
            {
            	continue;
            }
     	break;

     case 30:
     	destroy_matrix(m1);
     	destroy_matrix(m2);
     	printf("Done destroying matrices\n");
     	printf("Enter p to proceed\n");
            while((proceed = getchar()) != 'p')
            {
            	continue;
            }
    	break;

     case 31:
     	loop = 0;
     	break;

     default:
     	loop = 0;
     	break;

        }
    }

    return 0;

   }

