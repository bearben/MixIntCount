#include <fstream>
#include <iostream>
#include <chrono>
#include "unistd.h"
#include <boost/math/distributions/normal.hpp>
#include "polytope.h"

using namespace std;

typedef std::chrono::high_resolution_clock Clock;

int rseed;

void WilsonCI(double p, double z, int n, double &up, double &lo){

	double wmid = (p + z * z / (2 * n)) / (1 + z * z / n);
	double woffset = z * sqrt( (p * (1 - p) + z * z / (4 * n)) / n) / (1 + z * z / n);
	up = wmid + woffset;
	if (up > 1) up = 1;
	lo = wmid - woffset;
	if (lo < 0) lo = 0;
	
}

polyvest::polytope* ReadFile(string fileName)
{

	polyvest::polytope *p;
	
	string line;

	ifstream fileReader(fileName);
	if (fileReader.is_open())
	{
		vector<double> arr;
	
		while (getline(fileReader, line))
		{
			if (line.find('#') != string::npos)
				line.erase(line.find('#'));
				
			istringstream iss(line);
			
			double num;
			while ((iss >> num))
				arr.push_back(num);
		}
		fileReader.close();
		
		if (arr.size() < 2)
		{
			cout << "Error: Lack file content! Please Check Input File." << endl;
			exit(-1);
		}
			
		unsigned m = arr[0];
		unsigned n = arr[1] - 2;
		
		if (arr.size() != m * (n + 2) + n + 2)
		{
			cout << "Error: Incorrect number of elements! Please Check Input File" << endl;
			exit(-1);
		}
			
		p = new polyvest::polytope(m, n, rseed);
		
		unsigned index = 2;
		
		for (unsigned i = 0; i < n; i++)
		{
			p->var_flag[i] = arr[index];
			index++;
		}
		
		for (unsigned i = 0; i < m; i++) 
		{
			//P->vecOP[i] = (int)arr[index];	// 1: <=; 0: =;
			index++;		
		
			for (unsigned j = 0; j < n; j++)
			{
				p->matA(arr[index], i, j);
				index++;
			}
			
			p->vecb(arr[index], i);
			index++;
		}
		
		return p;
		
	}
	else
	{
		cout << "Error: " << fileName << " cannot open!" << endl;
		exit(-1);
	}
	
}

void print_help()
{
	cout << "Usage: ./mixIntCount [options ..] <infile>\n";
	cout << "Options:\n";
	cout << "  -h or -help\t\tPrint help.\n";
	cout << "  -a or -approx\t\tEmploy approximate method (by default).\n";
	cout << "    -seed=..\t\t  Set random seed.\n";
	cout << "    -epsilon=..\t\t  Set accuracy parameter (by default epsilon = 0.1).\n";
	cout << "    -delta=..\t\t  Set probability parameter (by default delta = 0.05).\n";
	cout << "    -maxs=..\t\t  Set the max number of samples (by default maxs = 100000).\n";
	cout << "    -vc or -vinci\t  Employ volume computation tool Vinci.\n";
	cout << "    -pv or -polyvest\t  Employ volume estimation tool PolyVest.\n";
	cout << "  -e or -exact\t\tEmploy exact method.\n";
}


int main(int argc, char **argv) 
{	
	rseed = (int)time(NULL);
	srand(rseed);

	bool	approx_method = true;
	bool	exact_method = false;
	bool	call_vinci = false;
	bool	call_polyvest = false;
	double 	epsilon = 0.1;
	double 	delta = 0.05;
	boost::math::normal dist(0.0, 1.0);
	double	z = boost::math::quantile(dist, 1.0 - delta / 2);
	int 	maxs = 100000;

	if (argc < 2)
	{
		print_help();
		exit(-1);
	}

	string fileName = argv[argc - 1];

	for (int i = 1; i < argc; i++)
	{
		string arg_str = argv[i];
		int len = arg_str.length();

		if (arg_str[0] != '-' || len < 2) continue;
		
		if (arg_str == "-h" || arg_str == "-help")
		{
			print_help();
			exit(1);
		}
		
		if (arg_str == "-a" || arg_str == "-approx")
		{
			cout << "Enable approximate method and disable exact method.\n";
			approx_method = true;
			exact_method = false;
		} else if (arg_str.find("-seed=") != string::npos && len > 6)
		{
			string val_str = arg_str.substr(6, len - 6);
			rseed = stod(val_str);
			cout << "Set random seed = " << rseed << endl;
		} else if (arg_str.find("-epsilon=") != string::npos && len > 9)
		{
			string val_str = arg_str.substr(9, len - 9);
			epsilon = stod(val_str);
			cout << "Set parameter epsilon = " << epsilon << endl;
		} else if (arg_str.find("-delta=") != string::npos && len > 7)
		{
			string val_str = arg_str.substr(7, len - 7);
			delta = stod(val_str);
			z = boost::math::quantile(dist, 1.0 - delta / 2);
			cout << "Set parameter delta = " << delta << ", z = " << z << endl;
		} else if (arg_str.find("-maxs=") != string::npos && len > 6)
		{
			string val_str = arg_str.substr(6, len - 6);
			maxs = stod(val_str);
			cout << "Set parameter maxs = " << maxs << endl;
		} else if (arg_str == "-vc" || arg_str == "-vinci")
		{
			cout << "Employ volume computation tool Vinci.\n";
			call_vinci = true;
		} else if (arg_str == "-pv" || arg_str == "-polyvest")
		{
			cout << "Employ volume estimation tool PolyVest.\n";
			call_polyvest = true;
		} else if (arg_str == "-e" || arg_str == "-exact")
		{
			cout << "Enable exact method and disable approximate method.\n";
			exact_method = true;
			approx_method = false;
		} else
		{
			cout << "Unknown option \"" << arg_str << "\".\n";
			exit(-1);
		}
	}
	
	ofstream finalOut("results.txt", ios::app);
	finalOut << fileName << ", ";
	finalOut.close();
	
	auto t1 = Clock::now();
	
	polyvest::polytope *P;
	
	P = ReadFile(fileName);
	P->Print();
	
	int m = P->m;
	int n = P->n;

	cout << "Read file \'" << fileName << "\' finished.\n";
	
	if (approx_method)
	{
	
		// copy from P then enlarge it
		polyvest::polytope *bigP = P->Clone(rseed);
		bigP->Enlarge();
	
		// copy from bigP then apply ellipsoid
		polyvest::polytope *bigQ = bigP->Clone(rseed);
		bigQ->AffineTrans();

		arma::vec point(n);
		int SampleCount = 0;
		int counterP = 0;
		int counterCSP = 0;
		int counterCSBP = 0;
		double r1, wur1, wlr1, r2, wur2, wlr2;
		r1 = wur1 = wlr1 = 1000;
		r2 = wur2 = wlr2 = 0;

		bigQ->PrepForWalk();
		for (int i = 0; i < n; i++) 
			bigQ->Walk();
		
		for (int i = 0; i < maxs; i++)
		{

			for (int j = 0; j < n; j++) 
				bigQ->Walk();
			point = bigQ->GetInvPoint(bigQ->x);
		
			if (P->isInside(point)) 
			{
				counterP++;
				counterCSP++;
				if (bigP->CUonBoundary(point)) 
					counterCSBP++;
			} else 
				if (bigP->CUonBoundary(point))
				{
					counterCSP++;
					counterCSBP++;
				}
			
			r1 = (double)counterP / (double)counterCSP;			
			r2 = (double)(counterCSP - counterCSBP) / (double)counterCSP;
			
			WilsonCI(r1, z, counterCSP, wur1, wlr1);
			WilsonCI(r2, z, counterCSP, wur2, wlr2);
		
			SampleCount = i + 1;
			
			if (wur2 / wlr1 - wlr2 / wur1 <= epsilon && wur1 / wlr1 <= 1 + epsilon)
			{
				cout << wur1 << ' ' << wlr1 << ' ' << wur2 << ' ' << wlr2 << endl;
				cout << "Stop criterion satisfied at No. " << i + 1 << " sample." << endl;
				break;
			}
		
		}
		
		double VolP = 0;
		if (call_vinci)
		{
			VolP = P->ExactVol();
			cout << "v Exa V(P):\t" << VolP << endl;
		} else if (call_polyvest)
		{
			polyvest::polytope *tP = P->Clone(rseed);
			if (tP->AffineTrans()){
				cout << "Preprocess finished.\n";
				VolP = tP->EstimateVol(0.2, 0.05, 1.0);
			}
			delete tP;
			cout << "v Est V(P):\t" << VolP << endl;
		} else
			cout << "For Vol(P) please use option \'-vinci\' and \'-polyvest\'.\n";
		
		finalOut.open("results.txt", ios::app);
		finalOut << VolP << ", ";
	
		cout << "c P:       \t" << counterP << endl;
		cout << "c CS(P):   \t" << counterCSP << endl;
		cout << "c CS(B(P)):\t" << counterCSBP << endl;
		if (counterP >= 0)
		{
			finalOut << 1 / r1 << ", " << 1 / wur1 << ", " << 1 / wlr1 << ", "
					 << r2 / r1 << ", " << wlr2 / wur1 << ", " << wur2 / wlr1 << ", ";
			cout << "b UPPER:   \t" << 1 / r1 << " * Vol(P), [" << 1 / wur1 << ", " << 1 / wlr1 << "]\n";
			cout << "b LOWER:   \t" << r2 / r1 << " * Vol(P), [" << wlr2 / wur1 << ", " << wur2 / wlr1 << "]\n";
		} else
		{
			finalOut << "inf, " << 1 / wur1 << ", " << 1 / wlr1 << ", "
					 << "nan, " << wlr2 / wur1 << ", " << wur2 / wlr1 << ", ";
			cout << "b UPPER:   \tinf" << ", [" << 1 / wur1 << ", " << 1 / wlr1 << "]\n";
			cout << "b LOWER:   \tnan" << ", [" << wlr2 / wur1 << ", " << wur2 / wlr1 << "]\n";
		}
		
		finalOut << SampleCount << ", " << counterCSP << ",";
		finalOut.close();		
		
		delete bigP;
		delete bigQ;
		
	}
	
	if (exact_method)
	{
	
		double MIP = 0;
		long vcounter = 0;
		int intcount = 0;
		int *intidx = new int[n] ();
		for (int i = 0; i < n; i++)
			if (P->var_flag[i]) 
			{
				intidx[intcount] = i;
				intcount++;
			}
		
		if (intcount == 0)
			MIP = P->ExactVol();
		else if (intcount == n)
		{
			cout << "Pure integer problem, please call counting algorithm for ILP instead.\n";
			exit(-1);
		} else
		{
		
			int assignc = 1;
			int *assign = new int[n] ();
			int *vup = new int[n] ();
			
			// init from first integer variable
			int idx0 = intidx[0];
			P->GetBounds(idx0);
			
			assign[idx0] = ceil(P->xmin[idx0]);
			vup[idx0] = floor(P->xmax[idx0]);
			
			while (assign[idx0] <= vup[idx0])
			{
				// create a sub-problem tP with the assignment
				arma::vec a(n);
				a.zeros();
				for (int i = 0; i < assignc; i++)
					a[intidx[i]] = assign[intidx[i]];
				
				polyvest::polytope *tP = new polyvest::polytope(m, n - assignc, rseed);
				
				tP->b = P->b - P->A * a;

				int counter = 0;
				int tidx = 0;
				for (int j = 0; j < n; j++)
				{
					if (j == intidx[tidx])
					{
						if (tidx < assignc - 1) tidx++;
						continue;
					}
					for (int i = 0; i < m; i++)
						tP->A(i, counter) = P->A(i, j);
					tP->var_flag[counter] = P->var_flag[j];
					counter++;
				}
				
				if (assignc == intcount)
				{
					
					//for (int i = 0; i < assignc; i++)
						//cout << assign[intidx[i]] << '\t';
					//cout << endl;
				
					if (intcount == n) 
						MIP++;
					else
					{
						MIP += tP->ExactVol();
						vcounter++;
					}
					
					while (assign[intidx[assignc - 1]] >= vup[intidx[assignc - 1]]) 
					{
						if (assignc <= 1) break;
						assignc--;
					}
					
					assign[intidx[assignc - 1]]++;
					
				} else
				{
					int tidx = intidx[assignc];
					int uidx = tidx - assignc;
					tP->GetBounds(uidx);				
					assign[tidx] = ceil(tP->xmin[uidx]);
					vup[tidx] = floor(tP->xmax[uidx]);				
					assignc++;
				}
				
				delete tP;

			}
			
		}
		
		finalOut.open("results.txt", ios::app);
		finalOut << MIP << ", " << vcounter << ", ";
		finalOut.close();
		cout << "v MI(P):   \t" << MIP << endl;
		
	}

	auto t2 = Clock::now();
	
	double avgxr = 0;
	P->GetBounds();
	for (int i = 0; i < n; i++)
		avgxr += P->xmax[i] - P->xmin[i];
	avgxr = avgxr / n;
	cout << "x Avg Range:\t" << avgxr << endl;
	
	finalOut.open("results.txt", ios::app);
	finalOut << avgxr << ", " << 
				(double)chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / 1000000000 << endl;
	finalOut.close();
	
	delete P;
	
/*
	double vol = 0;
	if (P->AffineTrans()){
		cout << "Preprocess finished.\n";
		vol = P->EstimateVol(0.5, 0.1, 1.0);
	}
*/

	//ofstream finalOut("results.txt", ios::app);
	//finalOut << fileName << ", ";
	//finalOut << (double)std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / 1000000000 << endl;
	//finalOut.close();
}





