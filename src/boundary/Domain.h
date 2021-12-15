/// \file       Domain.h
/// \brief      stores index lists of domain
/// \date       Oct 01, 2019
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_BOUNDARY_DOMAIN_H_
#define ARTSS_BOUNDARY_DOMAIN_H_


#include "Obstacle.h"
#include "../utility/Utility.h"
#include "PatchObject.h"

class Domain {
 public:
    Domain(size_t *obstacle_list, size_t size_obstacle_list,
           size_t **surface_list, PatchObject &size_surface_list,
           size_t multigrid_level);

    ~Domain() = default;

    const size_t * get_domain_list() const { return m_domain_list.data(); }
    size_t get_size_domain_list() const { return m_domain_list.size(); }

    const size_t * get_inner_list() const { return m_inner_list.data(); }
    size_t get_size_inner_list() const { return m_inner_list.size(); }

    const std::vector<std::vector<size_t>> get_boundary_list() const { return m_boundary_patch_divided; }
    PatchObject * get_size_boundary_list() { return &m_size_boundary; }

    void update_lists(size_t *obstacle_list, size_t size_obstacle_list, size_t **surface_list, PatchObject &size_surface_list);
    void control(size_t size_obstacle_list, PatchObject &size_surface_list);

 private:
#ifndef BENCHMARKING
    std::shared_ptr<spdlog::logger> m_logger;
#endif
    size_t m_multigrid_level;

    std::vector<size_t> m_domain_list;

    std::vector<std::vector<size_t>> m_boundary_patch_divided;  // boundary cells, patch divided
    PatchObject m_size_boundary;

    std::vector<size_t> m_boundary_list;  // boundary cells altogether

    std::vector<size_t> m_inner_list;

    void init(size_t size_obstacle_list, PatchObject &size_surface_list);
    void inner_cells(const size_t *obstacle_list, size_t size_obstacle_list);
    void boundary_cells(size_t **surface_list, PatchObject &size_surface_list);
    void print(size_t size_obstacle_list, PatchObject &size_surface_list);
    void clear_lists();

    void joined_list();
};

#endif /* ARTSS_BOUNDARY_DOMAIN_H_ */
