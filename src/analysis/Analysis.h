/// \file 		Analysis.h
/// \brief 		Calculates residual, compares analytical and numerical solutions, saves variables
/// \date 		July 11, 2016
/// \author 	Severt
/// \copyright 	<2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_ANALYSIS_ANALYSIS_H_
#define ARTSS_ANALYSIS_ANALYSIS_H_

#include "../interfaces/ISolver.h"
#include "../utility/GlobalMacrosTypes.h"
#include "Solution.h"

class Analysis {
public:
	explicit Analysis(Solution *solution);

	void Analyse(ISolver* solver, real t);
	//real* CalcL2NormMidPoint(real t, real* sum, read_ptr num_u, read_ptr num_p, read_ptr num_T);
	void CalcL2NormMidPoint(ISolver* solver, real t, real* sum);
	void CalcRMSError(real sumu, real sump, real sumT);
	bool CheckTimeStepVN(Field* u, real dt);
	bool CheckTimeStepCFL(Field* u, Field* v, Field* w, real dt);
	real SetDTwithCFL(Field* u, Field* v, Field* w);
	void SaveVariablesInFile(ISolver* solv);

private:
	real m_tol = 1e-7;

	bool CompareSolutions(read_ptr num, read_ptr ana, FieldType type, real t);
	real CalcAbsoluteSpatialError(read_ptr num, read_ptr ana);
	real CalcRelativeSpatialError(read_ptr num, read_ptr ana);
	void writeFile(const real *field, std::string filename, size_t *innerList, size_t size_innerList, size_t *boundaryList, size_t size_boundaryList, size_t *obstacleList,
                   size_t size_obstacleList);

    bool hasAnalyticSolution = false;
    Solution *m_solution;
};

#endif /* ARTSS_ANALYSIS_ANALYSIS_H_ */
