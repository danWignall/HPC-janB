#include<iostream>

//double laplacian



/*//function to print matrices in a nice way
//nice for small N or big monitors*/
void printMatrix(double **M, int size)
{
	for(int i=0; i<size;i++)
	{
		for(int j=0; j<size;j++)
		{
			std::cout	<<M[i][j];
			if(j==size-1)
				std::cout<<"\n";
			else
				std::cout<<" ";
		}
	}
    std::cout<<"\n";
	
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
    


    double U[N][N]={};
    //initialiseGrid(U,N);
    
    
    printMatrix(U[][],N);
    
    
    





    
    return 0.0;
}


int main()
{
    serial();
    
}
