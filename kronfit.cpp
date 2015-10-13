#include "stdafx.h"
#include "kronecker.h"
#include "math.h"
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <stdio.h>
#include <sstream>


using namespace std;

int main(int argc, char* argv[]) 
{
	Env = TEnv(argc, argv, TNotify::StdNotify);
	Env.PrepArgs(TStr::Fmt("Kronecker graphs. build: %s, %s. Time: %s", __TIME__, __DATE__, TExeTm::GetCurTm()));
	TExeTm ExeTm;
	Try
		int nofile=0, nopara=0, fileC=0, paraC=0, fileC1=0, paraC1=0;
		ofstream outfile; outfile.open("NofNParameterised.txt", ios::app);
		ofstream outfile1; outfile1.open("NofNReal.txt", ios::app);
        ofstream outfile2; outfile2.open("ParameterVals.txt", ios::app);
		ifstream myfile("InputGraphs.txt");
		ifstream mypara("Parameters.txt");
		string line, paraline;
		if (myfile.is_open() && mypara.is_open())
		{   
			while (getline(myfile,line) )
			{   
                nofile++;
			}
			myfile.clear();
			myfile.seekg (0, ios::beg);
			outfile<<"# No. of networks considered \t"<<nofile<< endl ; 
			outfile1<<"# No. of networks considered \t"<<nofile<< endl ;
			
			while(getline(mypara,paraline))
			{   nopara++;
				
			}
			mypara.clear();
			mypara.seekg (0, ios::beg);
			fileC1=nopara;
			outfile<<"# No. of parameters considered \t"<<nopara<< endl ;
			outfile1<<"# No. of parameters considered \t"<<nopara<< endl ;
            //outfile2<<std::fixed;
            outfile2<<"FileName \t Estimated initiator Matrix(or fitted matrix) \t No. of nodes \t No. of edges \t Scaled Nodes \t Log likelihood \t Average log likelihood \t Scaled Average log likelihood"<<endl;
			while (getline(myfile,line) )
			{	fileC++;
				fileC1++;
				while(getline(mypara,paraline))
				{   paraC++;
					paraC1++;
                    std::istringstream iss(line);
                    
                    int TempNumOne=line.size();
					char Filename[100];
                    iss>> Filename;
					//for (int a1=0;a1<=TempNumOne;a1++)
					//{
					//	Filename[a1]=line[a1];
					//}

					char s1[100], i, j, s3[100], k=0;
					//string pa= argv[1]; //Use with cluster (Linux)
					string pa="InputGraphs";
                    
					
					int TempNum=pa.size();
					
					for(int a=0; a< TempNum; a++) 
					s1[a]=pa[a];
					
                    s1[TempNum]='\0';
                    printf("File Name is %s ..", Filename);
                   /* for (int a3=0; !(Filename[a3]==' ' && Filename[a3]==' ');a3++)
					{
						Filename[a3]=line[a3];
					}
                    */
				/*for(i=0; s1[i]!='\0'; ++i);    // i contains length of string s1. 
					for(j=0; Filename[j]!='"'; ++j, ++i)
					{
						s1[i]=Filename[j];
						
					}
					s1[i]='\0';
					++j;
                    printf("S1 is %s ..", s1);
			
					*/
					std::istringstream Parameters(paraline);

					int TempNumTwo=paraline.size();
					char para[100];
					//Parameters>> para;
					printf("Size %d",TempNumTwo );
					for (int a2=0;a2<TempNumTwo;a2++)
					{
						para[a2]=paraline[a2];
					}

					printf("Para is %s ..", para);
			
				
					Env = TEnv(argc, argv, TNotify::StdNotify);
					const TStr InFNm = Env.GetIfArgPrefixStr("-i:", Filename, "Input graph file (single directed edge per line)"); ///
					const TStr PFitTo = Env.GetIfArgPrefixStr("-i:", "falseLL", "Name of graph to whom whose parameters we're fitting"); ///
					TStr OutFNm = Env.GetIfArgPrefixStr("-o:", "", "Output file prefix");
					const TInt NZero = Env.GetIfArgPrefixInt("-n0:", 2, "Initiator matrix size");
					const TStr InitMtx = Env.GetIfArgPrefixStr("-m:", "0.9, 0.5; 0.5, 0.1" , "Init Gradient Descent Matrix (R=random)").GetLc(); ///
					const TStr Perm = Env.GetIfArgPrefixStr("-p:", "d", "Initial node permutation: d:Degree, r:Random, o:Order").GetLc();
					const TInt GradIter = Env.GetIfArgPrefixInt("-gi:", 3, "Gradient descent iterations");  ///
					const TFlt LrnRate = Env.GetIfArgPrefixFlt("-l:", 1e-5, "Learning rate");
					const TFlt MnStep = Env.GetIfArgPrefixFlt("-mns:", 0.005," Minimum gradient step");
					const TFlt MxStep = Env.GetIfArgPrefixFlt("-mxs:", 0.05, "Maximum gradient step");
					const TInt WarmUp =  Env.GetIfArgPrefixInt("-w:", 10000, "Samples to warm up");
					const TInt NSamples = Env.GetIfArgPrefixInt("-s:", 100000, "Samples per gradient estimation");
					//const TInt GradType = Env.GetIfArgPrefixInt(-gt:, 1, 1:Grad1, 2:Grad2);
					const bool ScaleInitMtx = Env.GetIfArgPrefixBool("-sim:", false, "Scale the initiator to match the number of edges"); ///
					const TFlt PermSwapNodeProb = Env.GetIfArgPrefixFlt("-nsp:", 1.0, "Probability of using NodeSwap (vs. EdgeSwap) MCMC proposal distribution");
					if (OutFNm.Empty()) { OutFNm = TStr::Fmt("%s-fit%d", InFNm.GetFMid().CStr(), NZero()); }
					// load graph
					PNGraph G;
					if (InFNm.GetFExt().GetLc()==".ungraph") {
						TFIn FIn(InFNm);  G=TSnap::ConvertGraph<PNGraph>(TUNGraph::Load(FIn), true); 
					}
					else if (InFNm.GetFExt().GetLc()==".ngraph") 
					{
						TFIn FIn(InFNm);  G=TNGraph::Load(FIn); 
					}
					else 
					{
						G = TSnap::LoadEdgeList<PNGraph>(InFNm, 0, 1);
					}
					// fit
					TKronMtx InitKronMtx = InitMtx=="r" ? TKronMtx::GetRndMtx(NZero, 0.1) : TKronMtx::GetMtx(InitMtx);
					InitKronMtx.Dump("INIT PARAM", true);
					TKroneckerLL KronLL(G, InitKronMtx, PermSwapNodeProb);
					if (ScaleInitMtx) 
					{
						InitKronMtx.SetForEdges(G->GetNodes(), G->GetEdges()); 
					}
					KronLL.InitLL(G, InitKronMtx);
					InitKronMtx.Dump("SCALED PARAM", true);
					KronLL.SetPerm(Perm.GetCh(0));
					double LogLike = 0;
					int noedges=G->GetEdges();
                    double nonodes=G->GetNodes();
                    int ScaledNodes=pow(2,ceil(log(nonodes)/log(2.0)));
                    double fract= nonodes/ScaledNodes;

					FILE *F = fopen(TStr::Fmt("KronFit-%s.tab", InFNm.GetFMid().CStr()).CStr(), "at"); //Individual_Results/
					fprintf(F, "Input\t%s\n", InFNm.CStr());
					TStrV ParamV; Env.GetCmLn().SplitOnAllCh(' ', ParamV);
					fprintf(F, "Command line options\n");

					FILE *FIter = fopen("KronFitIteration.txt", "at"); //Results For each Iteration In One File//
					fprintf(FIter, "Input Network \t %s \n", InFNm.CStr());
					fprintf(FIter, "Calculation Mode \t %s \n", "New");  //BaseLine for Lescovec's Method and New for Our Method
					fprintf(FIter, "Initial Initiator Matrix \t [ %s ] \n", InitMtx.CStr());
					fprintf(FIter, "Original Number of Nodes \t [ %f ] \n", nonodes);
					fprintf(FIter, "Scaled Number of Nodes \t [ %d ] \n", ScaledNodes);
					fclose(FIter);
					//if (GradType == 1) {
					LogLike = KronLL.GradDescent(GradIter, LrnRate, MnStep, MxStep, WarmUp, NSamples);
					//} else if (GradType == 2) {
					//LogLike = KronLL.GradDescent2(GradIter, LrnRate, MnStep, MxStep, WarmUp, NSamples); }
					//else{ Fail; }

					const TKronMtx& FitMtx = KronLL.GetProbMtx();
                    //mkdir("Individual_Results", 0777);
					for (int i = 0; i < ParamV.Len(); i++) 
					{
						fprintf(F, "\t%s\n", ParamV[i].CStr()+(ParamV[i][0]=='-'?1:0)); 
					}
					fprintf(F," Initial Initiator Matrix \t [ %s ] \n",InitMtx.CStr() );
					fprintf(F, "Parameters of graph are being fitted is: %s \n", PFitTo.CStr() );
					fprintf(F," Loglikelihood\t%10.2f\n", LogLike);
					fprintf(F, "Absolute error (based on expected number of edges)\t%f\n", KronLL.GetAbsErr());
					fprintf(F, "RunTime\t%g\n", ExeTm.GetSecs());
					fprintf(F, "Estimated initiator\t%s\n", FitMtx.GetMtxStr().CStr());
					fprintf(F, "Current Time \t%s\n", TExeTm::GetCurTm());
					fprintf(F, "Execution Time \t%s\n", ExeTm.GetTmStr());
					fclose(F);
					printf("\nrun time: %s (%s)\n", ExeTm.GetTmStr(), TSecTm::GetCurTm().GetTmStr().CStr());
					outfile<<std::fixed;
					outfile1<<std::fixed;
                    outfile2<<std::fixed;
                    outfile2<<InFNm.CStr()<<"\t"<<FitMtx.GetMtxStr().CStr()<<"\t"<< nonodes<<"\t"<< noedges<<"\t"<<ScaledNodes<<"\t"<<LogLike<<"\t"<<LogLike/noedges<<"\t"<<(fract*LogLike)/noedges<<endl;
                    outfile<<"P"<<paraC<<"\tN"<<fileC<<"\t"<<LogLike<<"\t"<< LogLike/noedges<<endl;
					outfile1<<paraC1<<"\t"<<fileC1<<"\t"<<LogLike<<"\t"<<LogLike/noedges<<endl;
				}
				paraC=0;
				paraC1=0;
				mypara.clear();
				mypara.seekg (0, ios::beg);

				cout<<"\\\\\\\\\\\\\\\\\\\\\ para /////////////////////";
			}
		cout<<"\\\\\\\\\\\\\\\\\\\\\ LINE /////////////////////";
		myfile.close();
		mypara.close();
		}
	else cout << "Unable to open file";
	Catch
	getchar();
	return 0;
	
}
