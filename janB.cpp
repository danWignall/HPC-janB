#include<iostream>
#include<cmath>
#include<vector>
#include<fstream>
#include<string>
#include<omp.h>
#include<mpi.h>

using namespace std;


//value for comparing doubles
double epsilon = 1e-8;
int dEqual(double a, double b)
{
    if(a==b)
    {
        return 1; //catch the trivial case where the doubles are literally equal
    }
    //if not, see if they are equal within tolerance
    return (fabs(a - b) <= ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon ));//or fabs(a-b)<epsilon;
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
        
        
        //
        //CALC U
        //
        
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
        
        //
        //UPDATE PHI
        //
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
        
        
        //
        //EXIT STUFF
        //
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

int omp(int N,double w,double a,double b,double x1,double y1,double q1,double x2,double y2,double q2)
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
        
        
        //
        //CALC U
        //
        
        //calculate U which is used to update phi
        //exclude end points as edges grounded so phi=0
        #pragma omp parallel for default(none) private() shared(equal) reduction()
        for(int i=1;i<N-1;i++)
        {
            for(int j=1;j<N-1;j++)
            {
                
                U[i][j]=0.25*(phi[i+1][j]+phi[i-1][j]+phi[i][j+1]+phi[i][j-1]+rho[i][j]);
            }
        }

        //
        //UPDATE PHI
        //
        int equal=1;//this will stay 1 if phi doesnt change with iteration
        //otherwise, update phi        
        
        #pragma omp parallel for default(none) private() shared(equal) reduction()
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
        
        
        //
        //EXIT STUFF
        //
        if(equal)
        {
            cout<<"equilib reached, after "<<iter<<" Iterations.\n";
            fileOut(phi,"phi");
            break;
        }


int openMP(int N,double w,double a,double b,double x1,double y1,double q1,double x2,double y2,double q2)
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
    
    //this will stay 1 iff convergence found
    int equal=1;
    
    //switch the order of for and while loops so that only incur cost of parallelisation once, this means that will need syncronisation at the end of while loop so that processors stay in same iteration
    //this introduces an issue: cant update U[i][j] unless all surrounding phi[i][j] are also updated, otherwise could be adding together phi components from different iterations
    //to get around this, loop through the grid edges first, with barriers.
    #pragma omp parallel for default(none) private() shared(equal) reduction()
    for(int i=0;i<N;i++)
    {
        
        for(int j=0;j<N;j++)
        {
            while(true)
            {
                U[i][j]=0.25*(phi[i+1][j]+phi[i-1][j]+phi[i][j+1]+phi[i][j-1]+rho[i][j]);
                
                phi[i][j]=phi[i][j]+w*(U[i][j]-phi[i][j]);
                
                #pragma omp barrier
                
                //exit stuff
            }
        }
    }
    return 0;
}

int MPI(int N,double w,double a,double b,double x1,double y1,double q1,double x2,double y2,double q2)
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
    
    int ierr;
    MPI_Init (&argc, &argv); //parallel section
    
    
    int myRank; //rank of processor
    int size; //number of processors
    MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    
    //index that each thread starts/ends on (using 1D decomposition)
    int istart = myRank * N/size;
    int iend = istart+N/size;
    
    //the edges are grounded (phi=0) so dont change them 
    if(istart==0) istart=1; 
    if(iend==N) iend=N-1;
    
    
    //loop until convergence (happens when U=phi)
    while(true)
    {
        ////////////
        //update U//
        ////////////
        
        for(int i=istart;i<iend;i++)
        {
            for(int j=1;j<N-1;j++)//the edges (j=0,N) are grounded (phi=0) so dont change them 
            {
                
                U[i][j]=0.25*(phi[i+1][j]+phi[i-1][j]+phi[i][j+1]+phi[i][j-1]+rho[i][j]);
            }
        }
        
        int equal=1;//this will stay 1 if phi doesnt change with iteration
        
        //check if U[i][j] == phi[i][j], update phi[i][j] when it doesnt
        for(int i=istart;i<iend;i++)
        {
            for(int j=1;j<N-1;j++)
            {
                if(dEqual(U[i][j],phi[i][j]))
                {
                    continue;
                }
                phi[i][j]=phi[i][j]+w*(U[i][j]-phi[i][j]);
                equal=0;
            }
        }
        
        ierr=MPI_barrier(MPI_COMM_WORLD);
            
        
        ///////////////////////////////////////////////////////////////////////////////////////////
        //sync equal. iff all equals are 1, then convergence found, output phi and exit while loop//
        ///////////////////////////////////////////////////////////////////////////////////////////

        //equal = 1 iff all processor's equal = 1. otherwise, equal = 0.
        int globalEqual;
        ierr=MPI_Allreduce(&equal,&globalEqual,1,MPI_INT,MPI_PROD,MPI_COMM_WORLD);
        equal=globalEqual;  

        if(equal)
        {
            //each processor has a different phi, with phi[i][j] = 0 except inbetween istart-1 and iend-1 (inclusive). 
            //so could use global reduce with addition, but istart-1 from one processor overlaps with iend-1 from the previous one
            
            //so set phi[istart-1][j] = 0 to fix this
            if(myRank!=0) //but dont do this for first processor
            {                
                for(int j=0;j<N;j++)
                {
                    phi[istart-1][j]=0.0;
                }
            }
            ierr=MPI_barrier(MPI_COMM_WORLD);
            
            //now no overlap between different processors phi, so reduce.
            vector<vector<double>> globalPhi(N, vector<double> (N,0));
            ierr=MPI_ALLreduce(&phi[0][0],&globalPhi[0][0],N*N,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
            for(int i=0; i<N;i++);
            {
                for(int j=0;j<N;j++)
                {
                    phi[i][j]=globalPhi[i][j];
                }
            }
            
            break;     
        }
        
        
        ///////////////////////////////////////////////////////////////////////
        //if convergence not found, sync phi (only need to sync phi at edges)//
        ///////////////////////////////////////////////////////////////////////
        
        //even thread send end (high i) of thread's phi grid, odd thread recieve 
        if(not(myRank%2))//if myRank even
        {
            if(myRank!=size-1) //last thread shouldnt send anything
            {
                //send end of my grid to thread myRank + 1
                ierr=MPI_Send(&phi[iend-1][0],N,MPI_DOUBLE,myRank+1,myRank,MPI_COMM_WORLD);
            }
        }
        else //if myRank odd
        {
            if(myRank!=0) //first thread shouldnt recieve anything
            {
                //recieve from myRank -1
                ierr=MPI_Recv(&phi[istart-1][0],N,MPI_DOUBLE,myRank-1,myRank-1,MPI_COMM_WORLD);
            }
        }
        ierr=MPI_barrier(MPI_COMM_WORLD);
        
        //odd thread send end (high i) of phi grid, even thread recieve
        if(myRank%2) // if myRank odd
        {
            if(myRank!=size-1) //last thread shouldnt send anything
            {
                //send end of my grid to thread myRank + 1
                ierr=MPI_Send(&phi[iend-1][0],N,MPI_DOUBLE,myRank+1,myRank,MPI_COMM_WORLD);
            }
        }
        else //if myRank even
        {
            
            if(myRank!=0) //first thread shouldnt recieve anything
            {
                //recieve from myRank -1
                ierr=MPI_Recv(&phi[istart-1][0],N,MPI_DOUBLE,myRank-1,myRank-1,MPI_COMM_WORLD);
            }
        }
        ierr=MPI_barrier(MPI_COMM_WORLD);
        
        
        //
        //do the same but for the start (low i) of the phi grid
        //
        
        
        //even thread send start (low i) of thread's phi grid, odd thread recieve 
        if(not(myRank%2))//if myRank even
        {
            if(myRank!=0) //first thread shouldnt send anything
            {
                //send start of my grid to thread myRank - 1
                ierr=MPI_Send(&phi[istart][0],N,MPI_DOUBLE,myRank-1,myRank,MPI_COMM_WORLD);
            }
        }
        else //if myRank odd
        {
            if(myRank!=size-1) //last thread shouldnt recieve anything
            {
                //recieve from myRank + 1
                ierr=MPI_Recv(&phi[iend][0],N,MPI_DOUBLE,myRank+1,myRank+1,MPI_COMM_WORLD);
            }
        }
        ierr=MPI_barrier(MPI_COMM_WORLD);
        
        //odd thread send start (low i) of thread's phi grid, even thread recieve 
        if(myRank%2)//if myRank odd
        {
            if(myRank!=0) //first thread shouldnt send anything
            {
                //send start of my grid to thread myRank - 1
                ierr=MPI_Send(&phi[istart][0],N,MPI_DOUBLE,myRank-1,myRank,MPI_COMM_WORLD);
            }
        }
        else //if myRank even
        {
            if(myRank!=size-1) //last thread shouldnt recieve anything
            {
                //recieve from myRank + 1
                ierr=MPI_Recv(&phi[iend][0],N,MPI_DOUBLE,myRank+1,myRank+1,MPI_COMM_WORLD);
            }
        }        
        ierr=MPI_barrier(MPI_COMM_WORLD);
        
    }
        
        
    MPI_Finalise();
    
    cout<<"equib reached.\n";
    fileOut(phi,"phiMPI");
}




int main(int argc, char *argv[])
{
    int N=100;
    double w=1;
    int iter;
    
    //N, w, a, b, x1,y1,q1, x2,y2,q2
    serial(N,w,1,1,0.5,0.5,1,0.6,0.4,2);

}
