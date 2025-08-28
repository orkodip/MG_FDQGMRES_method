//2D STEADY HEAT DIFFUSION EQUATION. CELL CENTERED FINITE DIFFERENCE DISCRETIZATION USING GHOST CELLS.
//TATBE (1995) PROBLEM 2
//MULTIGRID PRECONDITIONED FDQGMRES SOLVER
#include<iostream>
#include<fstream>
#include<cmath>
#include<time.h>
#include<cstring>
#define TOL 1e-6	//convergence criteria
#define SMALL 1e-4
using namespace std;
const double WIDTH=1.0,HEIGHT=1.0;	//domain size
const int I=1024,J=1024;	//finest grid size
#include "GMG2.cpp"
#include "MG_FDQGMRES.cpp"
class HEAT:public MG_FDQGMRES
{
	double **k0,**k1,**k2,**k3,**k4;	//thermal conductivity for different grid levels
	double **TH;	//temperature field
	void gen_coeff();	//create the coefficients at the finest grid
	void res_coeff(const int I,const int J,double **k0,double **k1);	//restrict the coeffecient matrix
	void gen_mat(const int I,const int J,int *C,int *R,double *A,double **k,
		double *b=NULL);	//create the sparse coefficient matrix and if required RHS vector
	public:
			HEAT(); ~HEAT();	//memory management
			void MG_solve();	//multigrid cycle
			void write();	//file output for Tecplot360
};
HEAT::HEAT():MG_FDQGMRES(8,I*J,1000)	//MG_FDQGMRES parameterized constructor
{
	int IMAX,JMAX;
	k0=new double*[J+2];
	TH=new double*[J+1];
	for(int i=0;i<J+2;i++)
	{
		k0[i]=new double[I+2];
		if(i<J+1) TH[i]=new double[I+1];
	}
	IMAX=I/2; JMAX=J/2;
	k1=new double*[JMAX+2];
	for(int i=0;i<JMAX+2;i++)
		k1[i]=new double[IMAX+2];
	IMAX=I/4; JMAX=J/4;
	k2=new double*[JMAX+2];
	for(int i=0;i<JMAX+2;i++)
		k2[i]=new double[IMAX+2];
	IMAX=I/8; JMAX=J/8;
	k3=new double*[JMAX+2];
	for(int i=0;i<JMAX+2;i++)
		k3[i]=new double[IMAX+2];
	IMAX=I/16; JMAX=J/16;
	k4=new double*[JMAX+2];
	for(int i=0;i<JMAX+2;i++)
		k4[i]=new double[IMAX+2];
	cout<<"HEAT: MEMORY ALLOCATED"<<endl;
}
HEAT::~HEAT()
{
	int IMAX,JMAX;
	for(int i=0;i<J+2;i++)
	{
		delete[] k0[i];
		if(i<J+1) delete[] TH[i];
	}
	delete[] TH;
	delete[] k0;
	IMAX=I/2; JMAX=J/2;
	for(int i=0;i<JMAX+2;i++) delete[] k1[i];
	delete[] k1;
	IMAX=I/4; JMAX=J/4;
	for(int i=0;i<JMAX+2;i++) delete[] k2[i];
	delete[] k2;
	IMAX=I/8; JMAX=J/8;
	for(int i=0;i<JMAX+2;i++) delete[] k3[i];
	delete[] k3;
	IMAX=I/16; JMAX=J/16;
	for(int i=0;i<JMAX+2;i++) delete[] k4[i];
	delete[] k4;
	cout<<"HEAT: MEMORY RELASED"<<endl;
}
void HEAT::gen_coeff()
{
	int a=80;	//'T' shape dimension
	for(int j=0;j<=J+1;j++)	//initialization of the whole domain (including ghost nodes)
		for(int i=0;i<=I+1;i++)
			k0[j][i]=1.0;
	for(int j=a;j<=J-2*a;j++)	//bottom part of 'T'
		for(int i=5*a;i<=I-5*a;i++)
			k0[j][i]=100.0;
	for(int j=J-a;j>=J-3*a;j--)	//top part of 'T'
		for(int i=a;i<=I-a;i++)
			k0[j][i]=100.0;
}
void HEAT::res_coeff(const int I,const int J,double **k0,double **k1)
{
	for(int j=1;j<=J;j++)
		for(int i=1;i<=I;i++)
			k1[j][i]=0.25*(k0[2*j][2*i]+k0[2*j][2*i-1]+k0[2*j-1][2*i]+k0[2*j-1][2*i-1]);
}
void HEAT::gen_mat(const int I,const int J,int *C,int *R,double *A,double **k,double *b)
{
	double dx=WIDTH/I,dy=HEIGHT/J;	//grid size
	double bt=pow((dx/dy),2.0);
	double A_e,A_w,A_n,A_s;	//coefficients
	int cnt=0; R[cnt]=0;
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
		{
			if(b!=NULL)	//create the RHS vector
			{
				if(j<=J/2)
				{
					if(i<=I/2) b[(j-1)*I+i-1]=80.0*dx*dx;
					else b[(j-1)*I+i-1]=-80.0*dx*dx;
				}
				else
				{
					if(i<=I/2) b[(j-1)*I+i-1]=-80.0*dx*dx;
					else b[(j-1)*I+i-1]=80.0*dx*dx;
				}
			}
			A_w=2.0*(k[j][i]*k[j][i-1])/(k[j][i]+k[j][i-1]);	//harmonic mean
			A_e=2.0*(k[j][i]*k[j][i+1])/(k[j][i]+k[j][i+1]);
			A_s=bt*2.0*(k[j][i]*k[j-1][i])/(k[j][i]+k[j-1][i]);
			A_n=bt*2.0*(k[j][i]*k[j+1][i])/(k[j][i]+k[j+1][i]);
			if((i==1)&&(j==1))	//bottom left corner cell
			{
				C[cnt]=(j-1)*I+i-1;	//1,1 term
				A[cnt]=-(2.0*A_w+2.0*A_s+A_n+A_e);
				cnt++;
				C[cnt]=(j-1)*I+i;	//2,1 term
				A[cnt]=A_e;
				cnt++;
				C[cnt]=j*I+i-1;	//1,2 term
				A[cnt]=A_n;
				cnt++;
			}
			else if((i==1)&&(j>1)&&(j<J))	//left boundary
			{
				C[cnt]=(j-2)*I+i-1;	//1,j-1 term
				A[cnt]=A_s;
				cnt++;
				C[cnt]=(j-1)*I+i-1;	//1,j term
				A[cnt]=-(2.0*A_w+A_e+A_n+A_s);
				cnt++;
				C[cnt]=(j-1)*I+i;	//2,j term
				A[cnt]=A_e;
				cnt++;
				C[cnt]=j*I+i-1;	//1,j+1 term
				A[cnt]=A_n;
				cnt++;
			}
			else if((i==1)&&(j==J))	//top left corner cell
			{
				C[cnt]=(j-2)*I+i-1;	//1,J-1 term
				A[cnt]=A_s;
				cnt++;
				C[cnt]=(j-1)*I+i-1;	//1,J term
				A[cnt]=-(2.0*A_w+2.0*A_n+A_e+A_s);
				cnt++;
				C[cnt]=(j-1)*I+i;	//2,J term
				A[cnt]=A_e;
				cnt++;
			}
			else if((j==J)&&(i>1)&&(i<I))	//top boundary
			{
				C[cnt]=(j-2)*I+i-1;	//i,J-1 term
				A[cnt]=A_s;
				cnt++;
				C[cnt]=(j-1)*I+i-2;	//i-1,J term
				A[cnt]=A_w;
				cnt++;
				C[cnt]=(j-1)*I+i-1;	//i,J term
				A[cnt]=-(A_e+A_w+A_s+2.0*A_n);
				cnt++;
				C[cnt]=(j-1)*I+i;	//i+1,J term
				A[cnt]=A_e;
				cnt++;
			}
			else if((j==J)&&(i==I))	//top right corner cell
			{
				C[cnt]=(j-2)*I+i-1;	//I,J-1 term
				A[cnt]=A_s;
				cnt++;
				C[cnt]=(j-1)*I+i-2;	//I-1,J term
				A[cnt]=A_w;
				cnt++;
				C[cnt]=(j-1)*I+i-1;	//I,J term
				A[cnt]=-(A_s+A_w+2.0*A_n+2.0*A_e);
				cnt++;
			}
			else if((i==I)&&(j>1)&&(j<J))	//right boundary
			{
				C[cnt]=(j-2)*I+i-1;	//I,j-1 term
				A[cnt]=A_s;
				cnt++;
				C[cnt]=(j-1)*I+i-2;	//I-1,j term
				A[cnt]=A_w;
				cnt++;
				C[cnt]=(j-1)*I+i-1;	//I,j term
				A[cnt]=-(A_s+A_w+A_n+2.0*A_e);
				cnt++;
				C[cnt]=j*I+i-1;	//I,j+1 term
				A[cnt]=A_n;
				cnt++;
			}
			else if((i==I)&&(j==1))	//bottom right corner cell
			{
				C[cnt]=(j-1)*I+i-2;	//I-1,1 term
				A[cnt]=A_w;
				cnt++;
				C[cnt]=(j-1)*I+i-1;	//I,1 term
				A[cnt]=-(A_w+A_n+2.0*A_e+2.0*A_s);
				cnt++;
				C[cnt]=j*I+i-1;	//I,2 term
				A[cnt]=A_n;
				cnt++;
			}
			else if((j==1)&&(i>1)&&(i<I))	//bottom boundary
			{
				C[cnt]=(j-1)*I+i-2;	//i-1,1 term
				A[cnt]=A_w;
				cnt++;
				C[cnt]=(j-1)*I+i-1;	//i,1 term
				A[cnt]=-(A_w+A_e+A_n+2.0*A_s);
				cnt++;
				C[cnt]=(j-1)*I+i;	//i+1,1 term
				A[cnt]=A_e;
				cnt++;
				C[cnt]=j*I+i-1;	//i,2 term
				A[cnt]=A_n;
				cnt++;
			}
			else	//inner domain
			{
				C[cnt]=(j-2)*I+i-1;	//i,j-1 term
				A[cnt]=A_s;
				cnt++;
				C[cnt]=(j-1)*I+i-2;	//i-1,j term
				A[cnt]=A_w;
				cnt++;
				C[cnt]=(j-1)*I+i-1;	//i,j term
				A[cnt]=-(A_e+A_w+A_n+A_s);
				cnt++;
				C[cnt]=(j-1)*I+i;	//i+1,j term
				A[cnt]=A_e;
				cnt++;
				C[cnt]=j*I+i-1;	//i,j+1 term
				A[cnt]=A_n;
				cnt++;
			}
			R[(j-1)*I+i]=cnt;
		}
	}
}
void HEAT::MG_solve()
{
	clock_t start = clock();
	gen_coeff();
	gen_mat(I,J,C0,R0,A0,k0,b0);
	GMG2::setup(I/2,J/2,C0,R0,A0,C1,R1,A1);	//create the coarse grid coefficient matrices
	GMG2::setup(I/4,J/4,C1,R1,A1,C2,R2,A2);
	GMG2::setup(I/8,J/8,C2,R2,A2,C3,R3,A3);
	GMG2::setup(I/16,J/16,C3,R3,A3,C4,R4,A4);
	clock_t end = clock();
	cout<<"Setup time = "<<(double)(end-start)/CLOCKS_PER_SEC<<" seconds"<<endl;
	MG_FDQGMRES::solve(C0,R0,A0,X0,b0);
	for(int j=1;j<=J;j++)	//update the solution matrix
		for(int i=1;i<=I;i++)
			TH[j][i]=X0[IND(I,i,j)];
}
void HEAT::write()
{
	double dx=WIDTH/I,dy=HEIGHT/J;	//grid size
	double Xm[I+1],Ym[J+1];
	Xm[0]=Ym[0]=0.0;
	for(int i=1;i<=I;i++)
		Xm[i]=Xm[i-1]+dx;
	for(int j=1;j<=J;j++)
		Ym[j]=Ym[j-1]+dy;
	ofstream p_out("TH.dat");
	p_out<<"TITLE = \"Temperature field\""<<endl;
	p_out<<"VARIABLES = \"X\",\"Y\",\"TH\""<<endl;
	p_out<<"ZONE I="<<I+1<<", J="<<J+1<<", DATAPACKING=BLOCK, VARLOCATION=([3]=CELLCENTERED)"<<endl;
	for(int j=0;j<=J;j++)	//print X co-ordinates of mesh
	{
		for(int i=0;i<=I;i++)
			p_out<<" "<<Xm[i];
		p_out<<endl;
	}
	p_out<<endl<<endl;
	for(int j=0;j<=J;j++)	//print Y co-ordinates of mesh
	{
		for(int i=0;i<=I;i++)
			p_out<<" "<<Ym[j];
		p_out<<endl;
	}
	p_out<<endl<<endl;
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
			p_out<<" "<<TH[j][i];
		p_out<<endl;
	}
	p_out.close();
	cout<<"HEAT: WRITE SUCCESSFULL"<<endl;
}
int main()
{
	HEAT ms;
	clock_t start = clock();
	ms.MG_solve();
	clock_t end = clock();
	cout<<"Ellapsed time = "<<(double)(end-start)/CLOCKS_PER_SEC<<" seconds"<<endl;
	//ms.write();
	return 0;
}
