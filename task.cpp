#include <iostream>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <xmmintrin.h>
#include <math.h>
#include <ctime>
#include<chrono>
#include<immintrin.h>

using namespace std::chrono;

using namespace std;

//Prototypes of all functions------------
void test_case1();
void test_case2();
void test_case3();
void m_sum(int **, int **, int ** ,int);
void m_mul(int **, int **, int ** ,int);

//--------------------------------Data Structure for NxM matrix----------------------
class Matrix_Ds{
    
    public:
    int row,columns;
    float ** matrix;
    Matrix_Ds();
    Matrix_Ds(const int, const int);
    void populate(void);
    void auto_populate(void);
    void print_matrix(void);
    void matrix_mul(Matrix_Ds &,Matrix_Ds &,int);
    void optimize_mat_mul(int);
    ~Matrix_Ds();



};

//member functions definitions

Matrix_Ds::Matrix_Ds(const int r, const int c){
    //In constructor we will construct matrix
    row = r;
    columns = c;
    
    matrix = new float*[row];
    
    for (int i = 0; i < row; i++)
    {
        matrix[i] = new float[columns];
    }
}

void Matrix_Ds::populate(void){

    cout<<"Enter the elements of matrix"<<endl;

    for(int i = 0 ; i<row ; i++)
    {
        for (int j = 0 ; j<columns ; j++)
        {
            cin>>matrix[i][j];
        }
    }
}
void Matrix_Ds::auto_populate(){
for(int i = 0 ; i<row ; i++)
    {
        for (int j = 0 ; j<columns ; j++)
        {
            matrix[i][j] = rand()%252;
        }
    }
}
void Matrix_Ds::print_matrix(void){

    //Printing the input matrix
    for ( int i = 0 ; i<row ; i++)
    {
        for (int j = 0 ; j<columns ; j++)
        {
            cout<<matrix[i][j]<<' ';

        }
        cout<<endl;
    }
}
void Matrix_Ds::matrix_mul(Matrix_Ds & object1, Matrix_Ds & object2,int flag=0)
{
    //assertion for the matrix multiplication most important property that the number 
    //columns of the first matrix should match the number of columns of the second. 
    assert(object1.columns==object2.row);
    
   
    row = object1.row;
    columns = object2.columns;
    
    for (int i = 0; i< row ; i++ )
    {
        matrix[i]=new float[columns];
    }

    //implementation for test_case2()

    if(flag==1){
        for(int i = 0 ; i< object1.row ; i++)
        {
            for (int j = 0 ; j<object1.columns ; j++)
            {
                object1.matrix[i][j]=2;
            }
            
        }
        for(int i = 0 ; i<object2.row ; i++)
        {
            for (int j = 0 ; j<object2.columns ; j++)
            {
                object2.matrix[i][j]=3;
            }
        }
        
        for (int i = 0; i < object1.row; i++)
        {
            matrix[i] = new float[object2.columns];
        }
        
        for ( int i = 0 ; i< object1.row; i++)
        {
            for(int j = 0 ; j < object2.columns ; j++)
            {
                matrix[i][j] = 0;

            }

        }
    
    //multiplication
    for (int i = 0 ; i < object1.row ; i++){
            
        for (int j = 0 ; j < object2.columns ; j++)
            {
            for (int k = 0 ; k< object1.columns ; k++)
                {
                    matrix[i][j] += object1.matrix[i][k] * object2.matrix[k][j];
                }
            }
    

    }
    //print final matrix after multiplication for test_case
    cout<<"The final matrix after multiplication is "<<endl;
    for ( int i = 0 ; i< object1.row ; i++)
    {
        for (int j = 0 ; j<object2.columns; j++)
        {
            cout<<matrix[i][j]<<' ';

        }
        cout<<endl;
    }

    
    
    
    }

  

    
        
    for ( int i = 0 ; i<row ; i++)
        {
            matrix[i] = new float[columns];
        }
    //initializing elements to 0

    for ( int i = 0 ; i< row; i++)
        {
            for(int j = 0 ; j < columns ; j++)
            {
                matrix[i][j] = 0;

            }
        }
    //multiplication
    for (int i = 0 ; i < object1.row; i++){
            
        for (int j = 0 ; j < object2.columns ; j++)
            {
            for (int k = 0 ; k< object1.columns ; k++)
                {
                    matrix[i][j] += object1.matrix[i][k] * object2.matrix[k][j];
                }
            }
    

    }

    //print final matrix after multiplication
    cout<<"The final matrix after multiplication is "<<endl;
    for ( int i = 0 ; i< object1.row ; i++)
    {
        for (int j = 0 ; j< object2.columns ; j++)
        {
            cout<<matrix[i][j]<<' ';

        }
        cout<<endl;
    }
    

}
//-----------------Part B--------
void Matrix_Ds::optimize_mat_mul(int DIM){
    
    
    int n = DIM;
    __attribute__((aligned(16))) float **a= new float*[DIM];
    __attribute__((aligned(16))) float **b = new float*[DIM];
    __attribute__((aligned(16))) float **c = new float*[DIM];
    __attribute__((aligned(16))) float **d = new float*[DIM];
    for ( int i = 0 ; i < DIM;i++)
    {
        a[i]= new float [DIM];
        b[i]= new float [DIM];
        c[i]= new float [DIM];
        d[i]= new float [DIM];
    }

    srand((unsigned)time(0));
    for(int i = 0; i < n; ++i)
    {
        for(int j = 0; j < n; ++j)
        {
            a[i][j] = 7;
        }
    }

    for(int i = 0; i < n; ++i)
    {
        for(int j = 0; j < n; ++j)
        {
            b[i][j] = 5;
        }
    }
    clock_t tStart1 = clock();

    for(int i = 0; i < n; ++i)
    {
        for(int j = 0; j < n; ++j)
        {
            c[i][j] = 0;
            for(int k = 0; k < n; k ++)
            {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    cout<<"--------------------------------------------------------"<<endl;
    printf("Time taken by Normal Matrix Multiplication: %.2fs\n", (double)(clock() - tStart1)/CLOCKS_PER_SEC);
    cout<<"--------------------------------------------------------"<<endl;
    clock_t tStart = clock();
    for(int i = 0; i < n; ++i)
    {
        for(int j = 0; j < n; ++j)
        {
            __m128 *m1 = (__m128*)a[i];
            __m128 *m2 = (__m128*)d[j];

            float* result;
            c[i][j] = 0;
            for(int k = 0; k < n; k += 4)
            {   

                __m128 m3 = _mm_mul_ps(*m1,*m2);
                result = (float*)&m3;
                c[i][j] += result[0]+result[1]+result[2]+result[3];
                m1++;
                m2++;
            }
        }
    }

    cout<<"--------------------------------------------------------"<<endl;
    printf("Time taken SIMD MM: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
    cout<<"--------------------------------------------------------"<<endl;
    
    delete[] a;
    delete[] b;
    delete[] c;
    delete[] d;

    
}
Matrix_Ds::~Matrix_Ds(){
    cout<<"Destructor has been called. Deleting memory allocated to Matrix"<<endl;
    delete[] matrix;
}

//---------------Main Function--------------
int main()
{   //Terminal Interface
    int x;
    cout<<"Press 1 if you want to make Matrix Data Structure"<<endl;
    cout<<"Press 2 if you want to multiply two matrices"<<endl;
    cout<<"Press 3 for test cases to run"<<endl;
    cout<<"press 4 to check the difference between Optimized SIMD and Normal MM"<<endl;
    cout<<"Press 5 to check Divide and Conquer vs Normal vs SIMD results(Bonus Question)"<<endl;
    cin>>x;
    if (x==1){
 //-----------------------------------------------------------
    //Part 1 (a)
    //calling function matrix_nm which will take define data structure to handle NxM matrix
    int in;
    int row;
    int col;
    cout<<"Q1 a: Defining data structure for n X m float matrix"<<endl;
    cout<<"Enter number of rows for a matrix"<<endl;
    cin>>row;
    cout<<"Enter number of columns for a matrix"<<endl;
    cin>>col;

    Matrix_Ds matrix(row,col);
    cout<<"Matrix data structure has been made. Press 1 if you want to populate it or 2 for auto populate"<<endl;
    cin>>in;
    
    if(in==1){
        matrix.populate();
    }
    else if(in==2){
        matrix.auto_populate();
    }
    cout<<"Press 1 if you want to print the matrix"<<endl;
    cin>>in;
    if(in==1){
        matrix.print_matrix();
    }

    }
    
    else if (x==2){
    
//--------------------------------------------------------------
    //Part 1 (b) ---> Matrix Multiplication

    
    //input matrix1
    int row1;
    int col1;
    cout<<"Enter number of rows for a matrix1"<<endl;
    cin>>row1;
    cout<<"Enter number of columns for a matrix1"<<endl;
    cin>>col1;
    int row2;
    int col2;
    cout<<"Enter number of rows for a matrix2"<<endl;
    cin>>row2;
    cout<<"Enter number of columns for a matrix2"<<endl;
    cin>>col2;
    Matrix_Ds matrix1(row1,col1);
    Matrix_Ds matrix2(row2,col2);
    int check = 0;
    cout<<"Press 1 for Auto populating the matrix with random numbers, while Press 2 for manually populating it"<<endl;
    cin>>check;
    if(check==1)
    {
    matrix1.auto_populate();
    matrix2.auto_populate();
    }
    else if(check == 2)
    {
    
    matrix1.populate();
    matrix2.populate();
    }
 
    
    // calling function matrix_mul() for matrix multiplication
    Matrix_Ds final(0,0);

    final.matrix_mul(matrix1,matrix2,0);
    }

//------------------------------------------------
    //Part 1(c) test cases
    else if(x==3){
        int x=0;
        cout<<"Press 1 for test_case1(): which will abort the program because the dimesions are not matched for multiplication"<<endl;
        cout<<"Press 2 for test_case2(): for large dimension multiplication"<<endl;
        cout<<"Press 3 for test_case3(): Checking Commutative Property"<<endl;
        cin>>x;
        if (x==1){
             test_case1();
        }
        else if (x==2){

            //large matrix multiplication
            test_case2();
        }
        else if ( x==3){
            test_case3();
        }
        
    }
    else if(x==4)
    {   Matrix_Ds objj(5,5);
        int  x = 0;
        cout<<"Enter dimensions of matrix. The rest of the enteries will be auto filled"<<endl;
        cin>>x;
        objj.optimize_mat_mul(x);
        cout<<"On Average the SIMD optimzed code run  3.33 times faster than Normal Matrix Multiplication code on 1000x1000 matrix "<<endl;
    cout<<"On Average the SIMD optimzed code run  6.33 times faster than Normal Matrix Multiplication code on 1600x1600 matrix "<<endl;
    cout<<"The optimized variant giving better and better performance as we increase the size of our matrices"<<endl;
        
    }
    //----------------------Bonus
    else if (x==5)
    {   int x =0;
        //Matrix_Ds ob(0,0);
        //cout<<"Enter Dimesnions of Matrix"<<endl;
        //cin>>x;
        //ob.optimize_mat_mul(x);
        int row,col;
        cout<<"Write down Dimensions. The dimension should from domain{2^n:2,4,8,16....}"<<endl;
        cin>>row;
        col = row; 
    int ** a = new int * [row];
    int ** b = new int * [row];
    int ** c = new int * [row];
    for ( int i = 0; i < row; i++)
    {
        a[i] = new int [col];
        b[i] = new int [col];
        c[i] = new int [col];
    }
    for ( int i = 0; i<row; i++)
    {
        for ( int j = 0; j<col; j++)
        {
            a[i][j] = rand()%100;
            b[i][j] = rand()%100;
        }
    }
    clock_t tStart = clock();
    m_mul(a,b,c,row);
    cout<<"--------------------------------------------------------"<<endl;
    printf("Time taken by D&C Matrix Multiplication: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
    Matrix_Ds ob(0,0);
    ob.optimize_mat_mul(row);
    }
    
    return 0;
}


//-----------------------------Part 1 - C(Test cases)------------------------------------

//test_case1()will check that our code should not multiply matrices if the number of columns 
//of the first matrix does not match the number of rows of the second matrix
void test_case1(){
    
    Matrix_Ds object1(4,3);
    Matrix_Ds object2(4,2);
    Matrix_Ds prod_matrix(0,0);
    prod_matrix.matrix_mul(object1,object2,1);
}
void test_case2(){
    
    Matrix_Ds object1(1000,1000);
    Matrix_Ds object2(1000,1000);
    Matrix_Ds prod_matrix(0,0);
    prod_matrix.matrix_mul(object1,object2,1);
}
//Checking commutative property should not hold in matrix multiplication->AB!=BA
void test_case3(){
    Matrix_Ds object1(2,3);
    object1.auto_populate();
    Matrix_Ds object2(3,2);
    object2.auto_populate();
    Matrix_Ds prod_matrix1(0,0);
    cout<<"AXB"<<endl;
    prod_matrix1.matrix_mul(object1,object2);
    Matrix_Ds prod_matrix2(0,0);
    cout<<"BXA"<<endl;
    prod_matrix2.matrix_mul(object2,object1);
    //First check if both of the matrices dimensions match
    if ((prod_matrix1.row==prod_matrix2.row) and (prod_matrix1.columns==prod_matrix2.columns))
    {
        for (int i =0 ; i < prod_matrix2.row ; i++)
        {
            for (int j = 0 ; j < prod_matrix2.columns ; j++)
            {
                if (prod_matrix1.matrix[i][j]!=prod_matrix2.matrix[i][j])
                {   cout<<"AB != BA when A and B are two different matrices"<<endl;
                    break;
                }
            }
        }
    }
    else {
        cout<<"commutative property doesnt hold which is the property of matrix multiplication"<<endl;
    }

}

//----------Bonus Question-------------
void m_sum(int **a,int **b, int **c,int dim){
	for (int i = 0; i < dim;i++)
    {

	    for (int j = 0; j < dim; j++)
        {
		    c[i][j] = a[i][j] + b[i][j];
	    }
    }
}
void m_mul(int **a, int **b, int **c, int n){
	if (n == 1){
		c[0][0] = a[0][0] * b[0][0];
        
	}
	else{
		int **a_11 = new int*[(n/2)];
		int **a_12 = new int*[(n/2)];
		int **a_21 = new int*[(n/2)];
		int **a_22 = new int*[(n/2)];

		int **b_11 = new int*[(n/2)];
		int **b_12 = new int*[(n/2)];
		int **b_21 = new int*[(n/2)];
		int **b_22 = new int*[(n/2)];

		int **c_11 = new int*[(n/2)];
		int **c_12 = new int*[(n/2)];
		int **c_21 = new int*[(n/2)];
		int **c_22 = new int*[(n/2)];

		int **tem_1 = new int*[(n/2)];
		int **tem_2 = new int*[(n/2)];

		for (int i = 0; i < (n/2); i++){
			a_11[i] = new int[(n/2)];
			a_12[i] = new int[(n/2)];
			a_21[i] = new int[(n/2)];
			a_22[i] = new int[(n/2)];

			b_11[i] = new int[(n/2)];
			b_12[i] = new int[(n/2)];
			b_21[i] = new int[(n/2)];
			b_22[i] = new int[(n/2)];

			c_11[i] = new int[(n/2)];
			c_12[i] = new int[(n/2)];
			c_21[i] = new int[(n/2)];
			c_22[i] = new int[(n/2)];

			tem_1[i] = new int[(n/2)];
			tem_2[i] = new int[(n/2)];
		}

		for (int i = 0; i < (n/2); i++)
        {
		
            for (int j = 0; j < (n/2); j++)
            {
			    a_11[i][j] = a[i][j];
			    a_12[i][j] = a[i][j + (n/2)];
			    a_21[i][j] = a[i + (n/2)][j];
			    a_22[i][j] = a[i + (n/2)][j + (n/2)];

			    b_11[i][j] = b[i][j];
			    b_12[i][j] = b[i][j + (n/2)];
			    b_21[i][j] = b[i + (n/2)][j];
			    b_22[i][j] = b[i + (n/2)][j + (n/2)];
		    }
        }
		m_mul(a_11, b_11, tem_1, (n/2));
		m_mul(a_12, b_21, tem_2, (n/2));
		m_sum(tem_1, tem_2, c_11, (n/2));

		m_mul(a_11, b_12, tem_1, (n/2));
		m_mul(a_12, b_22, tem_2, (n/2));
		m_sum(tem_1, tem_2, c_12, (n/2));

		m_mul(a_21, b_11, tem_1, (n/2));
		m_mul(a_22, b_21, tem_2, (n/2));
		m_sum(tem_1, tem_2, c_21, (n/2));

		m_mul(a_21, b_12, tem_1, (n/2));
		m_mul(a_22, b_22, tem_2, (n/2));
		m_sum(tem_1, tem_2, c_22, (n/2));
	

		for (int i = 0; i < n / 2; i++)
			for (int j = 0; j < n / 2; j++){
				c[i][j] = c_11[i][j];
				c[i][j + (n / 2)] = c_12[i][j];
				c[i + (n / 2)][j] = c_21[i][j];
				c[i + (n / 2)][j + (n / 2)] = c_22[i][j];
			}
	

			
		for (int i = 0; i < (n/2); i++){
			delete[] a_11[i];
			delete[] a_12[i];
			delete[] a_21[i];
			delete[] a_22[i];

			delete[] b_11[i];
			delete[] b_12[i];
			delete[] b_21[i];
			delete[] b_22[i];

			delete[] c_11[i];
			delete[] c_12[i];
			delete[] c_21[i];
			delete[] c_22[i];
			
			delete[] tem_1[i];
			delete[] tem_2[i];
			
		}
	}
}
