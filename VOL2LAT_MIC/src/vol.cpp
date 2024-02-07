/*  vol.cpp
 *
 *  Copyright (C) 2016-2017 Cunjing Ge.
 *
 *  All rights reserved.
 *
 *  This file is part of VolCE.
 *  See COPYING for more information on using this software.
 */

#include "solver.h"
#include "glpk.h"
#include <limits>

bool check_all_zeros(arma::rowvec r) {
	for (arma::rowvec::iterator it = r.begin(); it != r.end(); it++) {
		if (*it != 0) return false;
	}
	return true;
}

//////////////////////////////////////////////////////////////////////
//// Initialization //////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void volce::solver::vol_init() {

	nVars = vnum_list.size();
	nFormulas = ineq_list.size();
	
	//linear constraints	
	bigA.zeros(nFormulas, nVars);
	bigb.zeros(nFormulas);
	bigop = new int[nFormulas];
	
	for (unsigned int i = 0; i < nFormulas; i++) {		
		//retrieve ineq
		volce::ineqc ie = ineq_list[i];
				
		//op: eq 0, le -10, lt -1, gt 1, ge 10
		if (ie.iseq()) bigop[i] = 0;
		else bigop[i] = -10;
		
		//terms
		for (unsigned int j = 0; j < ie.size(); j++) {
			bigA(i, ie[j].id) = ie[j].m;
		}
		bigb(i) = ie.get_const_r();
	}

}

void volce::solver::mat_init(int *bools, unsigned int nRows, std::vector<int> vars) {

	unsigned int nVars = vars.size();
	unsigned int counter = 0;

	if (wordlength > 0) {
		nRows += 2 * nVars;
		matA.set_size(nRows, nVars);
		colb.set_size(nRows);
		rowop = new int[nRows];
		
		// wordlength bounds
		for (unsigned int i = 0; i < nVars; i++) {
			matA(counter, i) = 1;
			colb(counter) = pow(2, wordlength - 1) - 1;
			rowop[counter] = -10;
			counter++;
			
			matA(counter, i) = 1;
			colb(counter) = -pow(2, wordlength - 1);
			rowop[counter] = 10;
			counter++;
		}
		
	} else {
		matA.set_size(nRows, nVars);
		colb.set_size(nRows);
		rowop = new int[nRows];		
	}
	
	// inequalities
	for(unsigned int i = 0; i < nFormulas; i++) {
		if (bools[i] < 0) continue;
			
		if ((bigop[i] == 1 && bools[i] == 1) || (bigop[i] == -10 && bools[i] == 0)) {
			// >
			rowop[counter] = 1;
		} else if ((bigop[i] == 10 && bools[i] == 1) || (bigop[i] == -1 && bools[i] == 0)) {
			// >=
			rowop[counter] = 10;			
		} else if ((bigop[i] == -1 && bools[i] == 1) || (bigop[i] == 10 && bools[i] == 0)) {
			// <
			rowop[counter] = -1;
		} else if ((bigop[i] == -10 && bools[i] == 1) || (bigop[i] == 1 && bools[i] == 0)) {
			// <=
			rowop[counter] = -10;
		} else 
			assert(bigop[i] != 0);
		
		for (unsigned int j = 0; j < nVars; j++)
			matA(counter, j) = bigA(i, vars[j]);
		colb(counter) = bigb(i);

		bool redundent = false;
		
		if (enable_ge) {
			for (unsigned int j = 0; j < counter; j++){
				if (check_all_zeros(matA.row(counter) - matA.row(j)) && colb(counter) - colb(j) == 0) {
					if (rowop[counter] == rowop[j]) {
						redundent = true;
					} else if (rowop[counter] * rowop[j] == -100) {
						// one is <=, another one is >=
						redundent = true;
						rowop[j] = 0;
					}
				}
				if (check_all_zeros(matA.row(counter) + matA.row(j)) && colb(counter) + colb(j) == 0) {
					if (rowop[counter] + rowop[j] == 0) {
						redundent = true;
					} else if (rowop[counter] * rowop[j] == 100) {
						redundent = true;
						rowop[j] = 0;
					}
				}
			}
		}
		if (!redundent) counter++;
	}
	
	if (counter < nRows) {
		matA.resize(counter, nVars);
		colb.resize(counter);
	}
	
}


//////////////////////////////////////////////////////////////////////
//// Gauss Elimination ///////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
unsigned int volce::solver::gauss_elimination() {

	if (!enable_ge) return matA.n_cols;

	unsigned int counter = 0;

	for (unsigned int eqid = 0; eqid < matA.n_rows; ) {
		if (rowop[eqid] != 0){
			eqid++;
			continue;
		}
		
		// find first non-zero element
		unsigned int nzid = 0;
		for (unsigned int i = 0; i < matA.n_cols; i++)
			if (matA(eqid, i) != 0) {
				nzid = i;
				break;
			}
			
		//std::cout << eqid << ' ' << nzid << std::endl;

		for (unsigned int i = 0; i < matA.n_rows; i++) {
			if (eqid != i) {
				double coef = matA(i, nzid) / matA(eqid, nzid);
				matA.row(i) -= coef * matA.row(eqid);
				colb(i) -= coef * colb(eqid);
			}
		}
		
		matA.shed_row(eqid);
		matA.shed_col(nzid);
		colb.shed_row(eqid);
		for (unsigned int i = eqid; i < matA.n_rows; i++){
			rowop[i] = rowop[i + 1];
		}
		
		counter++;
	}

	if (counter > 0) {
		std::cout << "Gauss Elimination: " << counter << " variables." << std::endl;
	}
	
	return matA.n_cols;

}


//////////////////////////////////////////////////////////////////////
//// Factorization ///////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
const unsigned int volce::solver::get_decided_vars(int *bools, std::vector<int> &vars){

	for (unsigned int i = 0; i < nVars; i++) {
		bool undecided = true;
		for (unsigned int j = 0; j < nFormulas; j++) {
			if (bools[j] < 0) continue;
			if (bigA(j, i) != 0) {
				undecided = false;
				break;
			}
		}
		if (!undecided) {
			vars.push_back(i);
		}
	}
	
	assert(vars.size() <= nVars);
	
	return vars.size();
}

//merge target solution into source solution
const bool volce::solver::merge_sols(int *source, int *target){

	bool need_merge = false;
	for (unsigned int i = 0; i < nFormulas; i++){
		if (source[i] < 0) continue;
		else if(source[i] == target[i]){
			need_merge = true;
			break;
		}
	}
	if (need_merge)
		for (unsigned int i = 0; i < nFormulas; i++) 
			source[i] = (target[i] == -1) ? source[i] : target[i];
	
	return need_merge;
}

const unsigned int volce::solver::factorize_bsol(int *bools, std::vector<int*> &pbools){
	std::vector<int> constraints;
	std::vector<int> vars;
	unsigned int nc, nv;
	
	for (unsigned int i = 0; i < nFormulas; i++){
		if (bools[i] < 0) continue;
		constraints.push_back(i);
	}
	
	nc = constraints.size();
	nv = get_decided_vars(bools, vars);

	for (unsigned int i = 0; i < nv; i++){
		int *pbsol = new int[nFormulas];
		for (unsigned int j = 0; j < nFormulas; j++) pbsol[j] = -1;
		for (unsigned int j = 0; j < nc; j++)
			if (bigA(constraints[j], vars[i]) != 0) 
				pbsol[constraints[j]] = bools[constraints[j]];
		pbools.push_back(pbsol);
	}

	for (unsigned int i = 0; i < pbools.size(); i++){
		unsigned int j = i + 1;
		while (j < pbools.size()){
			if (merge_sols(pbools[i], pbools[j])){
				delete[] pbools[j];
				pbools.erase(pbools.begin() + j);
				j = i + 1;
			}else
				j++;
		}
	}
	
	return pbools.size();
}


void WilsonCI(double p, double z, int n, double &up, double &lo){

	double wmid = (p + z * z / (2 * n)) / (1 + z * z / n);
	double woffset = z * sqrt( (p * (1 - p) + z * z / (4 * n)) / n) / (1 + z * z / n);
	up = wmid + woffset;
	if (up > 1) up = 1;
	lo = wmid - woffset;
	if (lo < 0) lo = 0;
	
}

//bound checking for each call
//employing linear programming
//compute the error between the volume and the number of lattices
void volce::solver::bound_computation(double &lb, double &ub) {

	double 	epsilon = 0.1;
	double 	delta = 0.05;
	boost::math::normal dist(0.0, 1.0);
	double	z = boost::math::quantile(dist, 1.0 - delta / 2);
	int 	maxs = 100000;
	
	int m = matA.n_rows;
	int n = matA.n_cols;
	
	polyvest::polytope *P = new polyvest::polytope(m, n);
	
	for (unsigned i = 0; i < n; i++)
		P->var_flag[i] = 1;	// all int
	
	for (unsigned i = 0; i < m; i++) 
	{
		//P->vecOP[i] = (int)arr[index];	// 1: <=; 0: =;	
		if (rowop[i] == -10) {
			//LE
			for (unsigned j = 0; j < n; j++)
				P->matA(matA(i, j), i, j);
			P->vecb(colb(i), i);
		} else if (rowop[i] == -1) {
			//LT
			for (unsigned j = 0; j < n; j++)
				P->matA(matA(i, j), i, j);
			P->vecb(colb(i) - 0.000001, i);
		} else if (rowop[i] == 1) {
			//GT
			for (unsigned j = 0; j < n; j++)
				P->matA(-matA(i, j), i, j);
			P->vecb(-colb(i) - 0.000001, i);			
		} else if (rowop[i] == 10) {
			//GE
			for (unsigned j = 0; j < n; j++)
				P->matA(-matA(i, j), i, j);
			P->vecb(-colb(i), i);
		} else {
			// EQ
			for (unsigned j = 0; j < n; j++)
				P->matA(matA(i, j), i, j);
			P->vecb(colb(i), i);
		}	
	}
	
	// copy from P then enlarge it
	polyvest::polytope *bigP = P->Clone();
	bigP->Enlarge();
	
	// copy from bigP then apply ellipsoid
	polyvest::polytope *bigQ = bigP->Clone();
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
		
		SampleCount++;
			
		if (wur2 / wlr1 - wlr2 / wur1 <= epsilon && wur1 / wlr1 <= 1 + epsilon)	{
			break;
		}
		
	}

	lb = wlr2 / wur1;
	ub = 1 / wlr1;
				
	delete bigP;
	delete bigQ;	
	delete P;

}

//////////////////////////////////////////////////////////////////////
//// Volume Estimation ///////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
const volce::VOL_RES_CLS volce::solver::volume_estimation_basic(int *bools, unsigned int nRows, std::vector<int> vars, 
		double epsilon, double delta, double coef) {
		
	//nVars: the number of "decided" numeric variables
	//nRows: the number of "decided" linear formulas
	//nFormulas: the number of linear formulas
	//bools[nFormulas], vars[nVars]
	mat_init(bools, nRows, vars);
	if (gauss_elimination() == 0)
		return VOL_RES_CLS(1, 1, 1);
	double lb, ub;
	lb = ub = 1;
	bound_computation(lb, ub);
	
	//bound checking
	//if (err < 0) {
		//err_unbounded_polytope();
	//}
	
	//update stats of vol calls
	stats_vol_calls++;
	stats_total_dims += matA.n_cols;
	if (matA.n_cols > stats_max_dims) stats_max_dims = matA.n_cols;

	//estimating
	polyvest::polytope p(matA.n_rows, matA.n_cols);
	
	p.msg_off = true;
	
	for(unsigned int i = 0; i < matA.n_rows; i++){

		//insert one row		
		if (rowop[i] > 0) {
			p.b(i) = -colb(i);
			for (unsigned int j = 0; j < matA.n_cols; j++)
				p.A(i, j) = -matA(i, j);
		} else if (rowop[i] < 0) {
			p.b(i) = colb(i);
			for (unsigned int j = 0; j < matA.n_cols; j++)
				p.A(i, j) = matA(i, j);
		} else assert(rowop[i] != 0);
/*	
		for (unsigned int j = 0; j < i; j++){
			if (check_all_zeros(p.A.row(i) + p.A.row(j)) &&
				p.b(i) == (-1) * p.b(j)){
				//degenerate
				//std::cout << "DEGENERATE" << std::endl;
				return VOL_RES_CLS(0, err, -err);
			}
		}
*/
	}

	if (p.AffineTrans()){
		p.EstimateVol(epsilon, delta, coef);
		return VOL_RES_CLS(p.Volume(), p.Volume() * (1 + epsilon) * ub, p.Volume() / (1 + epsilon) * lb);
	}else{
		//degenerate
		//cout << "DEGENERATE" << std::endl;
		return VOL_RES_CLS(0, ub, lb);
	}
}

const volce::VOL_RES_CLS volce::solver::volume_estimation(int *boolsol, double epsilon, double delta, double coef){
	std::vector<int> vars;
	unsigned int nRows = 0;

	// count rows	
	for (unsigned int i = 0; i < nFormulas; i++)
		if (boolsol[i] >= 0)
			nRows++;

	if (nRows == 0)
		return VOL_RES_CLS(0, 0, 0);
	
	unsigned int nVars_decided_total = get_decided_vars(boolsol, vars);

	// no factorization
	if (!enable_fact){
		unsigned int nVars_undecided = nVars - nVars_decided_total;
		
		//volume of cube consisted of undecided variables
		double cube_vol = pow(pow(2, wordlength), nVars_undecided);
		if (wordlength == 0 && nVars_undecided > 0) {
			//unbounded
			err_unbounded_polytope();
		}
		
		if (nVars_decided_total == 0)
			return VOL_RES_CLS(cube_vol, cube_vol, cube_vol);
		else if (nVars_decided_total == 1)
			return volume_computation_light(boolsol, vars.back()) * cube_vol;
		else
			return volume_estimation_basic(boolsol, nRows, vars, epsilon, delta, coef) * cube_vol;
	}
	
	volce::VOL_RES_CLS vol = VOL_RES_CLS(1, 1, 1);
	
	//factorization
	std::vector<int*> pbools;
	unsigned int nVars_tmp = 0;
	unsigned int npbools = factorize_bsol(boolsol, pbools);

	//count the bunch factorized
	if (npbools > 1) stats_fact_bunches++;
	
	//compute each subproblem
	for (unsigned int i = 0; i < npbools; i++){
		vars.clear();
		nRows = 0;
		for (unsigned int j = 0; j < nFormulas; j++){
			if (pbools[i][j] < 0) continue;
			nRows++;
		}
	
		unsigned int nVars_decided = get_decided_vars(pbools[i], vars);
		nVars_tmp += nVars_decided;
		
		if (nVars_decided == 1)
			vol = vol * volume_computation_light(pbools[i], vars.back());
		else {
			// increase coef while partitions into some pieces
			vol = vol * volume_estimation_basic(pbools[i], nRows, vars, epsilon, delta, coef * npbools);
		}
	}
	
	assert(nVars_decided_total <= nVars);
	assert(nVars_decided_total == nVars_tmp);
	
	if (wordlength == 0 && nVars - nVars_decided_total > 0) {
		//unbounded
		err_unbounded_polytope();
	} 
	
	if (nVars_decided_total == 0) {
		double v = pow(pow(2, wordlength), nVars);
		return VOL_RES_CLS(v, v, v);
	} else
		return vol * pow(pow(2, wordlength), nVars - nVars_decided_total);
}


//////////////////////////////////////////////////////////////////////
//// Volume Computation //////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
const volce::VOL_RES_CLS volce::solver::volume_computation_light(int *bools, int var){

	//update stats of vol calls
	stats_vol_calls++;
	stats_total_dims++;

	double max, min;

	if (wordlength == 0){
		max = std::numeric_limits<double>::max();
		min = -std::numeric_limits<double>::max();
	}else{
		max = pow(2, wordlength - 1) - 1;
		min = -pow(2, wordlength - 1);
	}
	for (unsigned int i = 0; i < nFormulas; i++){
		if(bools[i] < 0) continue;
		
		double v = bigb[i] / bigA(i, var);		
		int cmp = bigop[i];
		
		if (bigA(i, var) < 0) cmp = -cmp;

		if ((cmp > 0 && bools[i] == 1) || (cmp < 0 && bools[i] == 0)){
			if (v > min) min = v;
		}else{
			if (v < max) max = v;
		}
	}
	
	if (max == std::numeric_limits<double>::max() || 
		min == -std::numeric_limits<double>::max()) {
		err_unbounded_polytope();
	} 
	
	double lat = lattice_counting_light(bools, var);
	if ((max - min) < 0) return VOL_RES_CLS(lat, lat, lat);
	else {
		//std::cout << "light:" << max - min << ' ' << lat << std::endl;
		return VOL_RES_CLS(lat, lat, lat);
	}
}

const volce::VOL_RES_CLS volce::solver::volume_computation_basic(int *bools, unsigned int nRows, std::vector<int> vars){
	//nVars: the number of "decided" numeric variables
	//nRows: the number of "decided" linear formulas
	//nFormulas: the number of linear formulas
	//bools[nFormulas], vars[nVars]
	mat_init(bools, nRows, vars);
	if (gauss_elimination() == 0)
		return VOL_RES_CLS(1, 1, 1);
/*	for (unsigned int i = 0; i < matA.n_rows; i++) {
		for (unsigned int j = 0; j < matA.n_cols; j++)
			std::cout << matA(i, j) << '\t';
		std::cout << rowop[i] << '\t' << colb(i) << std::endl;
	}*/
	double lb, ub;
	lb = ub = 1;
	bound_computation(lb, ub);
	
	//bound checking
	//if (err < 0) {
		//err_unbounded_polytope();
	//}
	
	//search previous computation result for reusing
	std::vector<int> bools_vec;
	if (enable_fact) {
		for (unsigned int i = 0; i < nFormulas; i++) 
			bools_vec.push_back(bools[i]);
		std::map<std::vector<int>, double>::iterator vol_map_iter = vol_map.find(bools_vec);
		if (vol_map_iter != vol_map.end()) {
			//result exist
			stats_vol_reuses++;
			return VOL_RES_CLS(vol_map_iter->second, vol_map_iter->second * ub, vol_map_iter->second * lb);
		}
	}
	
	//update stats of vol calls
	stats_vol_calls++;
	stats_total_dims += matA.n_cols;
	if (matA.n_cols > stats_max_dims) stats_max_dims = matA.n_cols;
	
	// compute
	std::string filename = tooldir + "/vinci_input_tmp"; //.ine

	std::ofstream ofile;
	
	ofile.open(filename + ".ine");
	if (!ofile.is_open()) {
		err_open_file(filename);
	}

	//set printf format
	ofile.setf(std::ios::fixed, std::ios::floatfield);

	ofile << "H-representation" << std::endl;
	ofile << "begin" << std::endl;
	ofile << matA.n_rows << " " << matA.n_cols + 1 << " real" << std::endl;
	
	for(unsigned int i = 0; i < matA.n_rows; i++) {

		if (rowop[i] > 0) {
			ofile << -colb[i] << " ";
			for (unsigned int j = 0; j < matA.n_cols; j++)
				ofile << matA(i, j) << " ";
			ofile << std::endl;
		} else if (rowop[i] < 0) {
			ofile << colb[i] << " ";
			for (unsigned int j = 0; j < matA.n_cols; j++)
				ofile << -matA(i, j) << " ";
			ofile << std::endl;
		} else {
			assert(rowop[i] != 0);
/*
			ofile << -colb[i] << " ";
			for (unsigned int j = 0; j < matA.n_cols; j++)
				ofile << matA(i, j) << " ";
			ofile << std::endl;
			ofile << colb[i] << " ";
			for (unsigned int j = 0; j < matA.n_cols; j++)
				ofile << -matA(i, j) << " ";
			ofile << std::endl;
*/
		}
	}
	
	ofile << "end" << std::endl;
	
	ofile.close();
	
	//execute vinci
	std::string cmd = tooldir + "/vinci " + filename + " >/dev/null";
	int proc = system(cmd.c_str());
	
	//read result
	std::ifstream ifile;
	double vol = 0;
	filename = resultdir + "/vinci.result";

	ifile.open(filename);
	if (!ifile.is_open()) {
		err_open_file(filename);
	}
	
	ifile >> vol;
	
	ifile.close();	
	
	if (enable_fact) {
		//new entry
		vol_map.insert(std::pair<std::vector<int>, double>(bools_vec, vol));
	}
	return VOL_RES_CLS(vol, vol * ub, vol * lb);
	
}

const volce::VOL_RES_CLS volce::solver::volume_computation(int *boolsol){

	std::vector<int> vars;
	unsigned int nRows = 0;
	volce::VOL_RES_CLS vol = VOL_RES_CLS(1, 1, 1);
	
	// count rows
	for (unsigned int i = 0; i < nFormulas; i++)
		if (boolsol[i] >= 0) 
			nRows++;

	if (nRows == 0)
		return VOL_RES_CLS(0, 0, 0);
	
	unsigned int nVars_decided_total = get_decided_vars(boolsol, vars);
	
	// no factorization
	if (!enable_fact){
		unsigned int nVars_undecided = nVars - nVars_decided_total;
		
		//volume of cube consisted of undecided variables
		double cube_vol = pow(pow(2, wordlength), nVars_undecided);
		if (wordlength == 0 && nVars_undecided > 0) {
			//unbounded
			err_unbounded_polytope();
		}
		
		if (nVars_decided_total == 0)
			return VOL_RES_CLS(cube_vol, cube_vol, cube_vol);
		else if (nVars_decided_total == 1)
			return volume_computation_light(boolsol, vars.back()) * cube_vol;
		else
			return volume_computation_basic(boolsol, nRows, vars) * cube_vol;
	}
	
	//factorization
	std::vector<int*> pbools;
	unsigned int nVars_tmp = 0;
	unsigned int npbools = factorize_bsol(boolsol, pbools);
	
	//count the bunch factorized
	if (npbools > 1) stats_fact_bunches++;
	
	//compute each subproblem
	for (unsigned int i = 0; i < npbools; i++){
		vars.clear();
		nRows = 0;
		for (unsigned int j = 0; j < nFormulas; j++){
			if (pbools[i][j] < 0) continue;
			nRows++;
		}
	
		unsigned int nVars_decided = get_decided_vars(pbools[i], vars);
		nVars_tmp += nVars_decided;
		
		if (nVars_decided == 1)
			vol = vol * volume_computation_light(pbools[i], vars.back());
		else
			vol = vol * volume_computation_basic(pbools[i], nRows, vars);
		//std::cout << "A:" << vol.value << ' ' << vol.upper << ' ' << vol.lower << std::endl;
	}
	
	assert(nVars_decided_total <= nVars);
	assert(nVars_decided_total == nVars_tmp);
	
	if (wordlength == 0 && nVars - nVars_decided_total > 0) {
		//unbounded
		err_unbounded_polytope();
	}
	
	if (nVars_decided_total == 0) {
		double v = pow(pow(2, wordlength), nVars);
		return VOL_RES_CLS(v, v, v);
	} else
		return vol * pow(pow(2, wordlength), nVars - nVars_decided_total);
}

//////////////////////////////////////////////////////////////////////
//// Lattice Counting ////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
const double volce::solver::lattice_counting_light(int *bools, int var){

	//update stats of vol calls
	stats_vol_calls++;
	stats_total_dims++;

	double max, min;

	if (wordlength == 0){
		max = std::numeric_limits<double>::max();
		min = -std::numeric_limits<double>::max();
	}else{
		max = pow(2, wordlength - 1) - 1;
		min = -pow(2, wordlength - 1);
	}
	for (unsigned int i = 0; i < nFormulas; i++){
		if(bools[i] < 0) continue;
		
		double v = bigb[i] / bigA(i, var);		
		int cmp = bigop[i];
		
		if (bigA(i, var) < 0) cmp = -cmp;
		
		if ((cmp == 1 && bools[i] == 1) || (cmp == -10 && bools[i] == 0)){
			//GT
			if (v >= min){
				if (v == ceil(v))
					min = v + 1;
				else
					min = ceil(v);
			}
		}else if ((cmp == 1 && bools[i] == 0) || (cmp == -10 && bools[i] == 1)){
			//LE
			if (v <= max) max = floor(v);
		}else if ((cmp == 10 && bools[i] == 1) || (cmp == -1 && bools[i] == 0)){
			//GE
			if (v >= min) min = ceil(v);
		}else if ((cmp == 10 && bools[i] == 0) || (cmp == -1 && bools[i] == 1)){
			//LT
			if (v <= max){
				if (v == floor(v))
					max = v - 1;
				else
					max = floor(v);
			}
		}else if (cmp == 0){
			//EQ
			if (v >= min) min = ceil(v);
			if (v <= max) max = floor(v);
		}
	}
	
	if (max == std::numeric_limits<double>::max() ||
		min == -std::numeric_limits<double>::max()) {
		err_unbounded_polytope();
	}
	
	if ((max - min) < 0) return 0;
	else return (max - min + 1);
}

const double volce::solver::lattice_counting_basic(int *bools, unsigned int nRows, std::vector<int> vars){
	//nVars: the number of "decided" numeric variables
	//nRows: the number of "decided" linear formulas
	//nFormulas: the number of linear formulas
	//bools[nFormulas], vars[nVars]
	mat_init(bools, nRows, vars);
	if (gauss_elimination() == 0)
		return 1;	// only one solution
	double lb, ub;
	lb = ub = 1;
	bound_computation(lb, ub);
	
	//bound checking
	//if (err < 0) {
		//err_unbounded_polytope();
	//}

	//search previous counting result for reusing
	std::vector<int> bools_vec;
	if (enable_fact) {
		for (unsigned int i = 0; i < nFormulas; i++) 
			bools_vec.push_back(bools[i]);
		std::map<std::vector<int>, double>::iterator vol_map_iter = vol_map.find(bools_vec);
		if (vol_map_iter != vol_map.end()) {
			//result exist
			stats_vol_reuses++;
			return vol_map_iter->second;
		}
	}

	//update stats of vol calls
	stats_vol_calls++;
	stats_total_dims += matA.n_cols;
	if (matA.n_cols > stats_max_dims) stats_max_dims = matA.n_cols;
	
	//counting
	std::string filename = tooldir + "/latte_input_tmp"; //.ine
	
	std::ofstream ofile;
	
	ofile.open(filename);
	if (!ofile.is_open()){
		err_open_file(filename);
	}

	//set printf format
	ofile.setf(std::ios::fixed, std::ios::floatfield);
	ofile.precision(0);

	ofile << matA.n_rows << " " << matA.n_cols + 1 << std::endl;
	
	for(unsigned int i = 0; i < matA.n_rows; i++) {

		if (rowop[i] == 1){
			//GT
			ofile << (-1) * (colb[i] * 1000000 + 1)  << " ";
			for (unsigned int j = 0; j < matA.n_cols; j++)
				ofile << (-1) * matA(i, j) * 1000000 << " ";
			ofile << std::endl;
		}else if (rowop[i] == -10){
			//LE
			ofile << colb[i] * 1000000 << " ";
			for (unsigned int j = 0; j < matA.n_cols; j++)
				ofile << matA(i, j) * 1000000 << " ";
			ofile << std::endl;
		}else if (rowop[i] == 10){
			//GE
			ofile << (-1) * colb[i] * 1000000 << " ";
			for (unsigned int j = 0; j < matA.n_cols; j++)
				ofile << (-1) * matA(i, j) * 1000000 << " ";
			ofile << std::endl;
		}else if (rowop[i] == -1){
			//LT
			ofile << colb[i] * 1000000 - 1 << " ";
			for (unsigned int j = 0; j < matA.n_cols; j++)
				ofile << matA(i, j) * 1000000 << " ";
			ofile << std::endl;
		} else {
			//EQ = LE + GE
			assert(rowop[i] != 0);
/*			
			ofile << colb[i] * 1000000 << " ";
			for (unsigned int j = 0; j < matA.n_cols; j++)
				ofile << matA(i, j) * 1000000 << " ";
			ofile << std::endl;
			ofile << (-1) * colb[i] * 1000000 << " ";
			for (unsigned int j = 0; j < matA.n_cols; j++)
				ofile << (-1) * matA(i, j) * 1000000 << " ";
			ofile << std::endl;
*/
		}
	}
	
	ofile.close();

	//execute latte
	std::string cmd = tooldir + "/count " + filename + " >/dev/null 2>/dev/null";
	int proc = system(cmd.c_str());
	
	//read result
	std::ifstream ifile;
	double count = 0;
	filename = resultdir + "/numOfLatticePoints";

	ifile.open(filename);
	if (!ifile.is_open()) {
		err_open_file(filename);
	}
	
	ifile >> count;
	
	ifile.close();

	if (enable_fact) {
		//new entry
		vol_map.insert(std::pair<std::vector<int>, double>(bools_vec, count));
	}
	//std::cout << count << std::endl;
	return count;
}

const double volce::solver::lattice_counting(int *boolsol){
	std::vector<int> vars;
	unsigned int nRows = 0;
	double count = 1;
	
	// count rows
	for (unsigned int i = 0; i < nFormulas; i++)
		if (boolsol[i] >= 0)
			nRows++;

	if (nRows == 0) 
		return 0;

	unsigned int nVars_decided_total = get_decided_vars(boolsol, vars);

	// no factorization
	if (!enable_fact){
		unsigned int nVars_undecided = nVars - nVars_decided_total;
		//cout << "The number of decided variables: " << vars.size() << std::endl;
		
		//volume of cube consisted of undecided variables
		double cube_count = pow(pow(2, wordlength), nVars_undecided);
		if (wordlength == 0 && nVars_undecided > 0) {
			//unbounded
			err_unbounded_polytope();
		}
		
		if (nVars_decided_total == 0)
			return cube_count;
		else if (nVars_decided_total == 1)
			return lattice_counting_light(boolsol, vars.back()) * cube_count;
		else
			return lattice_counting_basic(boolsol, nRows, vars) * cube_count;
	}
	
	//factorization
	std::vector<int*> pbools;
	unsigned int nVars_tmp = 0;
	unsigned int npbools = factorize_bsol(boolsol, pbools);
	
	//count the bunch factorized
	if (npbools > 1) stats_fact_bunches++;
	
	//compute each subproblem
	for (unsigned int i = 0; i < npbools; i++){
		vars.clear();
		nRows = 0;
		for (unsigned int j = 0; j < nFormulas; j++){
			if (pbools[i][j] < 0) continue;
			nRows++;
		}
	
		unsigned int nVars_decided = get_decided_vars(pbools[i], vars);
		nVars_tmp += nVars_decided;
		
		if (nVars_decided == 1)
			count *= lattice_counting_light(pbools[i], vars.back());
		else
			count *= lattice_counting_basic(pbools[i], nRows, vars);
		
	}
	
	assert(nVars_decided_total <= nVars);
	assert(nVars_decided_total == nVars_tmp);
	
	if (wordlength == 0 && nVars - nVars_decided_total > 0) {
		//unbounded
		err_unbounded_polytope();
	}
	
	if (nVars_decided_total == 0)
		return pow(pow(2, wordlength), nVars);
	else
		return count * pow(pow(2, wordlength), nVars - nVars_decided_total);
}




