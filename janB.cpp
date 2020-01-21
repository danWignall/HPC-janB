#include<iostream>
#include<cmath>
#include <vector>

using namespace std;


//prelim value for comparing doubles
double epsilon = 1e-8;
int dEqual(double a, double b)
{
    return fabs(a-b)<epsilon;
}

//checks if two matrices are equal
int matrixEqual(vector<vector<double>> A,vector<vector<double>> B)
{
    int equal=1;
    for(int i=0; i<A.size();i++)
    {
        for(int j=0; j<A.size();j++)
        {
            equal*=dEqual(A[i][j],B[i][j]);
        }
    }
    return equal;
}

//function to print matrices in a nice way
//nice for small N or big monitors
void printMatrix(vector<vector<double>> M)
{
	for(int i=0; i<M.size();i++)
	{
		for(int j=0; j<M.size();j++)
		{
			cout	<<M[i][j];
			if(j==M.size()-1)
				cout<<"\n";
			else
				cout<<" ";
		}
	}
    cout<<"\n";
}




double serial(int N,double w,double a,double b,double x1,double x2,double y1,double y2,double q1,double q2)
{
    
    vector<vector<double>> phi(N, vector<double> (N,0)); //electrostatic potential
    vector<vector<double>> rho(N, vector<double> (N,0)); //charge density
    vector<vector<double>> U(N, vector<double> (N,0)); //matrix used to iterate phi
    
    //intialize charge density
    //CHECK THIS
    rho[x1/(a*N)][y1/(b*N)]=q1/(N*N);
    rho[x2/(a*N)][y2/(a*N)]=q2/(N*N);
    
    //iterate pot until phi_new=phi_old
    while(true)
    {
        for(int i=0;i<N;i++)
        {
            for(int j=0;j<N;j++)
            {
                //calculate U which is used to update phi
                U=0.25*(phi[i+1][j]+phi[i-1][j]+phi[i][j+1]+phi[i][j-1]+N*N*rho[i][j]);
            }
        }
        

        
        //if U=phi_old, then w(U-phi_old)=0 and convergence found
        if(matrixEqual(phi,U))
        {
            break 
        } 
        else
        {
            //update phi
            for(int i=0;i<N;i++)
            {
                for(int j=0;j<N;j++)
                {
                    phi[i][j]=phi[i][j]+w*(U[i][j]-phi[i][j]);
                }
            }
            
            //boundries are earthed
            for(int i=0;i<N;i++)
            {
                phi[i][0]=0;
                phi[i][N-1]=0;
                phi[0][i]=0;
                phi[N-1][i]=0;
            }
            
        }
    }
    
    
    
    
    
    return 0.0;
}


int main()
{
    serial();
    
}
