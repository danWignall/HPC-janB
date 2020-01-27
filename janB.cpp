#include<iostream>
#include<cmath>
#include<vector>
#include<fstream>
#include<string>
#include<omp.h>
#include<mpi.h>
#include<chrono>


double epsilon = 1e-8;//value for comparing doubles
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
int matrixEqual(std::vector<std::vector<double>> A,std::vector<std::vector<double>> B)
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
void printMatrix(std::vector<std::vector<double>> M)
{
	for(int i=0; i<M.size();i++)
	{
		for(int j=0; j<M.size();j++)
		{
			std::cout	<<M[i][j];
			if(j==M.size()-1)
				std::cout<<"\n";
			else
				std::cout<<" ";
		}
	}
    std::cout<<"\n";
}
//outputs a file in the format for draw.py (assumes N is a multiple of 100) 
//and a file with the full phi data
void fileOut(std::vector<std::vector<double>> phi,std::string name)
{
    int N=phi.size();
    
    std::string fileName = name+".dat";
    std::ofstream file;
    file.open(fileName);
    
    //for N large, draw.py is very slow, so output file with N=100 approximation
    for(int i=0;i<N;i+=N/100)
    {
        for(int j=0;j<N;j+=N/100)
        {
            file<<i/(N/100)<<" "<<j/(N/100)<<" "<<phi[i][j]<<"\n";
        }
    }
    
    file.close();
    
    //also output a file with all the data
    std::string bigFileName = name+"Big"+".dat";
    std::ofstream fileBig;
    fileBig.open(bigFileName);
    
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            fileBig<<i<<" "<<j<<" "<<phi[i][j]<<"\n";
        }
    }
    
    fileBig.close();
}


double serial(int N,double w,double a,double b,double x1,double y1,double q1,double x2,double y2,double q2)
{
    //timings
    long start = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now().time_since_epoch()).count(); //wall time in ms
    long end=0;
    double tTaken; //will be used to convert end-start into a double
    double tooLong = 5; //if function takes longer than tooLong minutes, exit function and output phi.
    
    std::vector<std::vector<double>> phi(N, std::vector<double> (N,0)); //electrostatic potential
    std::vector<std::vector<double>> rho(N, std::vector<double> (N,0)); //charge density
    std::vector<std::vector<double>> U(N, std::vector<double> (N,0)); //matrix used to iterate phi
    
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
        
        //////////
        //CALC U//
        //////////
        
        //calculate U which is used to update phi
        //exclude end points as edges grounded so phi=0
        for(int i=1;i<N-1;i++)
        {
            for(int j=1;j<N-1;j++)
            {
                
                U[i][j]=0.25*(phi[i+1][j]+phi[i-1][j]+phi[i][j+1]+phi[i][j-1]+rho[i][j]);
            }
        }
        
        //////////////
        //UPDATE PHI//
        //////////////
        
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
        
        
        //////////////
        //EXIT STUFF//
        //////////////
        
        if(equal)
        {
            std::cout<<"equilibrium reached, after "<<iter<<" Iterations.\n";
            fileOut(phi,"phi");
            break;
        }
        
        /*
        //output some files along the way to see how phi evolves
        if(not(iter%5000))
        {
            fileOut(phi,".phi"+std::to_string(iter));
            std::cout<<iter<<"\n";
        }
        */
        
        //if program is taking too long, (more than tooLong minutes) output phi (check every 1000 iterations)
        if(not(iter%1000))
        {
            end=std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now().time_since_epoch()).count();
            tTaken=(end-start)/1000.0;
            if((tTaken/60)>tooLong)
            {
                std::cout<<"Serial took too long (after "<<iter<<" Iterations), Exiting.\n";
                fileOut(phi,"phi");
                break;
            }
        }

        //if there is a NaN in phi, it will spread until it reaches phi[N/2][N/2]
        //if this happens, there has been a problem, stop the program
        if(std::isnan(phi[N/2][N/2]))
        {
            std::cout<<"NaN after "<<iter<<" Iterations :(\n";
            break;
        }
        
        
        //iterate iter
        iter++;
        
    }
    
    //return time taken in seconds
    end=std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now().time_since_epoch()).count();
    tTaken=(end-start)/1000.0;
    return tTaken;
}

double parallelOMP(int N,double w,double a,double b,double x1,double y1,double q1,double x2,double y2,double q2)
{
    std::cout<<"Number of threads: "<<omp_get_num_threads()<<"\n";
    
    //timings
    double start = omp_get_wtime();
    double end=0;
    double tooLong = 5; //if function takes longer than tooLong minutes, exit function and output phi.
    
    std::vector<std::vector<double>> phi(N, std::vector<double> (N,0)); //electrostatic potential
    std::vector<std::vector<double>> rho(N, std::vector<double> (N,0)); //charge density
    std::vector<std::vector<double>> U(N, std::vector<double> (N,0)); //matrix used to iterate phi
    
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
        #pragma omp parallel for default(none) shared(phi,U,N,rho)
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
        
        #pragma omp parallel for default(none) shared(equal,phi,U,w,N)
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
            std::cout<<"equilibrium reached, after "<<iter<<" Iterations.\n";
            fileOut(phi,"phiOMP");
            break;
        }
        
        //if program is taking too long, (more than tooLong minutes) output phi (check every 1000 iterations)
        if(not(iter%1000))
        {
            end=omp_get_wtime();
            if(((end-start)/60)>tooLong)
            {
                std::cout<<"OMP took too long (after "<<iter<<" Iterations), Exiting.\n";
                fileOut(phi,"phiOMP");
                break;
            }
        }
        
        
        //if there is a NaN in phi, it will spread until it reaches phi[N/2][N/2]
        //if this happens, there has been a problem, stop the code
        if(std::isnan(phi[N/2][N/2]))
        {
            std::cout<<"NaN after "<<iter<<" Iterations :(\n";
            break;
        }
        
        iter++;
    }
    
    end=omp_get_wtime();
    return end-start;
}

double parallelMPI(int N,double w,double a,double b,double x1,double y1,double q1,double x2,double y2,double q2,int argc, char **argv)
{
    //timings
    double start = MPI_Wtime();
    double end =0;
    int tooLong = 5; //if program takes longer than tooLong minutes, exit and output phi
    int timeOutStat = 0; //set to one iff program takes longer than tooLong minutes, used to prevent program going on too long
    
    std::vector<std::vector<double>> phi(N, std::vector<double> (N,0)); //electrostatic potential
    std::vector<std::vector<double>> rho(N, std::vector<double> (N,0)); //charge density
    std::vector<std::vector<double>> U(N, std::vector<double> (N,0)); //matrix used to iterate phi
    
    
    //intialize charge density
    //nearest grid point to x1 is round(x1*(N-1)/a) (similar for x2, y1, y2)
    //charge density at nearest grid point to a point charge is q/N*N
    //but rho only appears multiplied by N^2, so use q instead, and drop the N^2 later
    rho[round(x1*(N-1)/a)][round(y1*(N-1)/b)]=q1;
    rho[round(x2*(N-1)/a)][round(y2*(N-1)/b)]=q2;
    
    int ierr;
    MPI_Init (&argc, &argv); //parallel section
    MPI_Status status;
    
    int myRank; //rank of processor
    int size; //number of processors
    MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    std::cout<<"Number of threads: "<<size<<"\n";
    
    
    //index that each thread starts/ends on (using 1D decomposition)
    int istart = myRank * N/size;
    int iend = istart+N/size;
    
    //the edges are grounded (phi=0) so dont change them 
    if(istart==0) istart=1; 
    if(iend==N) iend=N-1;
    
    
    //loop until convergence (happens when U=phi)
    int iter=0;
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
        
        ierr=MPI_Barrier(MPI_COMM_WORLD);
            
        
        ///////////////////////////////////////////////////////////////////////////////////////////
        //sync equal. iff all equals are 1, then convergence found, output phi and exit while loop//
        ///////////////////////////////////////////////////////////////////////////////////////////

        //equal = 1 iff all processor's equal = 1. otherwise, equal = 0.
        int globalEqual;
        ierr=MPI_Allreduce(&equal,&globalEqual,1,MPI_INT,MPI_PROD,MPI_COMM_WORLD);
        equal=globalEqual;  
        
        //if program is taking too long, (more than tooLong minutes) output phi (check every 1000 iterations)
        if(not(iter%1000))
        {
            ierr=MPI_Barrier(MPI_COMM_WORLD);
            end=MPI_Wtime();
            if(((end-start)/60)>tooLong)
            {
                equal=1;
                timeOutStat=1;
            }
        }
        
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
            ierr=MPI_Barrier(MPI_COMM_WORLD);
            
            //now no overlap between different processors phi, so reduce.
            std::vector<std::vector<double>> globalPhi(N, std::vector<double> (N,0));
            ierr=MPI_Allreduce(&phi[0][0],&globalPhi[0][0],N*N,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
            for(int i=0; i<N;i++)
            {
                for(int j=0;j<N;j++)
                {
                    phi[i][j]=globalPhi[i][j];
                }
            }
            
            if(timeOutStat)
            {
                std::cout<<"MPI took too long (after "<<iter<<" Iterations), Exiting\n";
            }
            else
            {
                std::cout<<"equilibrium reached, after "<<iter<<" Iterations.\n";
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
                ierr=MPI_Recv(&phi[istart-1][0],N,MPI_DOUBLE,myRank-1,myRank-1,MPI_COMM_WORLD,&status);
            }
        }
        ierr=MPI_Barrier(MPI_COMM_WORLD);
        
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
                ierr=MPI_Recv(&phi[istart-1][0],N,MPI_DOUBLE,myRank-1,myRank-1,MPI_COMM_WORLD,&status);
            }
        }
        ierr=MPI_Barrier(MPI_COMM_WORLD);
        
        
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
                ierr=MPI_Recv(&phi[iend][0],N,MPI_DOUBLE,myRank+1,myRank+1,MPI_COMM_WORLD,&status);
            }
        }
        ierr=MPI_Barrier(MPI_COMM_WORLD);
        
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
                ierr=MPI_Recv(&phi[iend][0],N,MPI_DOUBLE,myRank+1,myRank+1,MPI_COMM_WORLD,&status);
            }
        }        
        ierr=MPI_Barrier(MPI_COMM_WORLD);
    
        iter++;
    }
        
    
    MPI_Finalize();
    fileOut(phi,"phiMPI");
    
    end=MPI_Wtime();
    return end-start;
    
}




int main(int argc, char *argv[])
{
    //parameters
    int N=100;
    double w=1;
    
    double a=1;
    double b=1;
    
    double x1=0.5;
    double y1=0.5;
    double q1=1;
    
    double x2=0.6;
    double y2=0.6;
    double q2=-2;
    
    
    //timings
    double serialTime=0;
    double ompTime=0;
    double mpiTime=0;
    
    
    //serial
    std::cout<<"\nSerial\n";
    serialTime=serial(N,w,a,b,x1,y1,q1,x2,y2,q2);
    
    //omp
    std::cout<<"\nOMP\n";
    ompTime=parallelOMP(N,w,a,b,x1,y1,q1,x2,y2,q2);
    
    //mpi
    std::cout<<"\nMPI\n";    
    mpiTime=parallelMPI(N,w,a,b,x1,y1,q1,x2,y2,q2,argc,argv);
    
    std::cout   <<"\n\n"
                <<"---------------------------\n"
                <<"Serial Took: "<<serialTime<<"s\n"
                <<"OMP Took:    "<<ompTime<<"s\n"
                <<"MPI Took:    "<<mpiTime<<"s\n"
                <<"---------------------------\n";
}
