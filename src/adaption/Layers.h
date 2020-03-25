/// \file 		Layers.h
/// \brief 		Adaption class for initial condition with layers (layersT)
/// \date 		Dec 04, 2018
/// \author 	My Linh Würzburger
/// \copyright 	<2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_ADAPTION_LAYERS_H_
#define ARTSS_ADAPTION_LAYERS_H_


#include "../utility/GlobalMacrosTypes.h"
#include "../Field.h"
#include "../interfaces/AdaptionFunctionI.h"

class Layers:public AdaptionFunctionI {
public:
    Layers(Adaption* pAdaption, Field** fields);
    ~Layers();
    virtual bool update() ;
    virtual void applyChanges();
    virtual bool hasReduction(){return false;};
private:
    void adaptXDirection(real temperature, size_t noBufferCell);
    void adaptXDirection_serial(real temperature, size_t noBufferCell);
    void setXValues(bool start);
    size_t getExpansionSize();

    size_t m_minimal;
    size_t m_noBufferCells, m_timestep, m_timecounter, m_expansion_size;

    real m_checkValue;
    Field *m_T, *m_Ta, *m_Nu, *m_kappa, *m_gamma;
    real m_x1, m_x2, m_y1, m_y2, m_z1, m_z2;
    size_t m_nx, m_ny, m_nz;
    Adaption *m_pAdaption;
    Field *m_u, *m_v, *m_w, *m_P;
};


#endif /* ARTSS_ADAPTION_LAYERS_H_ */
