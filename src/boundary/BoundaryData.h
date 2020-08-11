/// \file 		BoundaryData.h
/// \brief 		Data class for boundary data
/// \date 		Oct 08, 2020
/// \author 	My Linh Würzburger
/// \copyright 	<2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_BOUNDARY_BOUNDARYDATA_H_
#define ARTSS_BOUNDARY_BOUNDARYDATA_H_

#include <string>
#include <vector>
#include "../utility/Utility.h"
#include "../utility/tinyxml2.h"
#include "../Field.h"

const size_t numberOfPatches = 6;
enum Patch : int {
    UNKOWN_PATCH = -1,
    FRONT = 0,
    BACK = 1,
    BOTTOM = 2,
    TOP = 3,
    LEFT = 4,
    RIGHT = 5
};

const size_t numberOfBoundaryConditions = 3;
enum BoundaryCondition : int {
    UNKOWN_CONDITION = -1,
    NEUMANN = 0,
    DIRICHLET = 1,
    PERIODIC = 2
};

class BoundaryData {
public:
    BoundaryData();

    ~BoundaryData();
    void print();


    static std::string getFieldTypeName(FieldType f);
    static std::string getBoundaryConditionName(BoundaryCondition bc);
    static std::string getPatchName(Patch p);

    static FieldType matchField(const std::string& s);
    static Patch matchPatch(const std::string& s);
    static BoundaryCondition matchBoundaryCondition(const std::string& s);

    void addBoundaryCondition(const std::vector<Patch>& patches, real value, BoundaryCondition boundaryCondition);
    BoundaryCondition getBoundaryCondition(Patch p){ return m_boundaryConditions[p];};
    real getValue(Patch p){ return m_values[p];};

    bool isEmpty(){ return !m_hasValues; };
private:

#ifndef BENCHMARKING
    std::shared_ptr<spdlog::logger> m_logger;
#endif
    BoundaryCondition* m_boundaryConditions;
    real* m_values;

    bool m_hasValues = false;
};
#endif /* ARTSS_BOUNDARY_BOUNDARYDATA_H_ */
