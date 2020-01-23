#include<iostream>
#include<cmath>
#include<vector>
#include<fstream>
#include<string>

using namespace std;


//prelim value for comparing doubles
double epsilon = 1e-6;
int dEqual(double a, double b)
{
    //return fabs(a-b)<epsilon;
    if(a==b)
    {
        return 1;
    }
    return (fabs(a - b) <= ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon ))or fabs(a-b)<epsilon;
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

//outputs a file in the format for draw.py 
//and a file in the format for draw1000.py
//assumes N is a power of 10
void fileOut(vector<vector<double>> A,string name)
{
    int N=A.size();
    
    string fileName = name+".dat";
    ofstream file;
    file.open(fileName);
    
    //for N larger than 100, draw.py is very slow, so output file with N=100 approximation
    for(int i=0;i<N;i+=N/100)
    {
        for(int j=0;j<N;j+=N/100)
        {
            file<<i/(N/100)<<" "<<j/(N/100)<<" "<<A[i][j]<<"\n";
        }
    }
    
    file.close();
    
    //also output a file with all the data
    string bigFileName = ".big"+name+".dat";
    ofstream fileBig;
    fileBig.open(bigFileName);
    
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            fileBig<<i<<" "<<j<<" "<<A[i][j]<<"\n";
        }
    }
    
    fileBig.close();
}


int serial(int N,double w,double a,double b,double x1,double y1,double q1,double x2,double y2,double q2)
{
    
    vector<vector<double>> phi(N, vector<double> (N,0)); //electrostatic potential
    vector<vector<double>> rho(N, vector<double> (N,0)); //charge density
    vector<vector<double>> U(N, vector<double> (N,0)); //matrix used to iterate phi
    
    //intialize charge density
    //nearest grid point to x1 is round(x1*(N-1)/a) (similar for x2, y1, y2)
    //charge density at nearest grid point to a point charge is q/N*N
    //but rho only appears multiplied by N^2, so use q instead, and drop the N^2 later
    rho[round(x1*(N-1)/a)][round(y1*(N-1)/b)]=q1;
    rho[round(x2*(N-1)/a)][round(y2*(N-1)/b)]=q2;
    
    
    int iter =0;
    //iterate pot until phi_new=phi_old
    while(true)
    {
        //calculate U which is used to update phi
        //exclude end points as edges grounded so phi=0
        for(int i=1;i<N-1;i++)
        {
            for(int j=1;j<N-1;j++)
            {
                
                U[i][j]=0.25*(phi[i+1][j]+phi[i-1][j]+phi[i][j+1]+phi[i][j-1]+rho[i][j]);
            }
        }
        

        /*
        //if U=phi_old, then w(U-phi_old)=0 and convergence found
        //only check every 100 loops, to save time
        if(not(iter%100))
        {
            if(matrixEqual(phi,U))
            {
                cout<<"equilib reached, after "<<iter<<" Iterations.\n";
                fileOut(phi,"phi");
                break;
            } 
        }
        */
        
        int equal=1;//this will stay 1 if phi doesnt change with iteration
        //otherwise, update phi        
        for(int i=0;i<N;i++)
        {
            for(int j=0;j<N;j++)
            {
                if(dEqual(U[i][j],phi[i][j]))
                {
                    continue;
                }
                phi[i][j]=phi[i][j]+w*(U[i][j]-phi[i][j]);
                equal=0;
            }
        }
        
        if(equal)
        {
            cout<<"equilib reached, after "<<iter<<" Iterations.\n";
            fileOut(phi,"phi");
            break;
        }
        
        
        //if there is a NaN in phi, it will spread until it reaches phi[N/2][N/2]
        //if this happens, there has been a problem, stop the code
        if(isnan(phi[N/2][N/2]))
        {
            cout<<"NaN after "<<iter<<" Iterations :(\n";
            break;
        }
        
        //output some files along the way to see how phi evolves
        if(not(iter%5000))
        {
            fileOut(phi,".phi"+to_string(iter));
            cout<<iter<<"\n";
        }
        
        //iterate iter
        iter++;
    }
    return iter;
}







int main(int argc, char *argv[])
{
    int N=stoi(argv[1]);
    double w=1;
    int iter;
    
    //N, w, a, b, x1,y1,q1, x2,y2,q2
    iter=serial(N,w,1,1,0.55,0.5,1,0.45,0.5,-1);
    
    return iter;
}
