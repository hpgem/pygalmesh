#ifndef GENERATE_FROM_INR_WITH_BB_HPP
#define GENERATE_FROM_INR_WITH_BB_HPP

#include <string>
#include <vector>

#include <CGAL/Image_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_domain_with_polyline_features_3.h>
#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/Implicit_mesh_domain_3.h>




namespace pygalmesh {

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Labeled_mesh_domain_3<K> Image_domain;
typedef CGAL::Mesh_domain_with_polyline_features_3<Image_domain> Mesh_domain;

bool add_bounding_box(const CGAL::Image_3& image,
                      Mesh_domain& domain);

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

} // namespace pygalmesh

#endif // GENERATE_FROM_INR_WITH_BB_HPP
