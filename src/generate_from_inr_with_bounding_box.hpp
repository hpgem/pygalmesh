#ifndef GENERATE_FROM_INR_WITH_BB_HPP
#define GENERATE_FROM_INR_WITH_BB_HPP

#include <string>
#include <vector>

namespace pygalmesh {

void generate_from_inr_with_bounding_box(
    const std::string & inr_filename,
    const std::string & outfile,
    const bool lloyd = false,
    const bool odt = false,
    const bool perturb = true,
    const bool exude = true,
    const double max_edge_size_at_feature_edges = 0.0,  // std::numeric_limits<double>::max(),
    const double min_facet_angle = 0.0,
    const double max_radius_surface_delaunay_ball = 0.0,
    const double max_facet_distance = 0.0,
    const double max_circumradius_edge_ratio = 0.0,
    const double max_cell_circumradius = 0.0,
    const bool verbose = true,
    const int seed = 0
    );

void
generate_from_inr_with_subdomain_sizing_with_bounding_box(
    const std::string & inr_filename,
    const std::string & outfile,
    const double default_max_cell_circumradius,
    const std::vector<double> & max_cell_circumradiuss,
    const std::vector<int> & cell_labels,
    const bool lloyd = false,
    const bool odt = false,
    const bool perturb  = true,
    const bool exude = true,
    const double max_edge_size_at_feature_edges = 0.0,
    const double min_facet_angle = 0.0,
    const double max_radius_surface_delaunay_ball = 0.0,
    const double max_facet_distance = 0.0,
    const double max_circumradius_edge_ratio = 0.0,
    const bool verbose  = true,
    const int seed = 0
    );

} // namespace pygalmesh

#endif // GENERATE_FROM_INR_WITH_BB_HPP
