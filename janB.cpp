#include<iostream>
#include <vector>
using namespace std;
//double laplacian



/*//function to print matrices in a nice way
//nice for small N or big monitors*/
void printMatrix(vector<double> M, int size)
{
	for(int i=0; i<size;i++)
	{
		for(int j=0; j<size;j++)
		{
			cout	<<M[i+size*j];
			if(j==size-1)
				cout<<"\n";
			else
				cout<<" ";
		}
	}
    cout<<"\n";
	
}

/*
double initialiseGrid(double* M, int size)
{
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            M[i][j]=1.0;
        }
    }
    
  
}
*/

double serial()
{
    int N = 10;
    


    vector<double> U(10,0);
    printMatrix(U,N);
    //initialiseGrid(U,N);
    U[5+N*10]=51;
    U[1+N*5]=15;
    
    printMatrix(U,N);
    
    
    





    
    return 0.0;
}


int main()
{
    serial();
    
}
