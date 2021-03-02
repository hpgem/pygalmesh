#define CGAL_MESH_3_VERBOSE 1

#include "generate_periodic_multiple_domains.hpp"

#include <CGAL/Periodic_3_mesh_3/config.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/make_periodic_3_mesh_3.h>
#include <CGAL/optimize_periodic_3_mesh_3.h>
#include <CGAL/Periodic_3_mesh_3/IO/File_medit.h>
#include <CGAL/Periodic_3_mesh_triangulation_3.h>
#include <CGAL/Periodic_3_function_wrapper.h>
#include <CGAL/Implicit_to_labeling_function_wrapper.h>
#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/number_type_config.h> // CGAL_PI
#include <cmath>
#include <iostream>
#include <fstream>


namespace pygalmesh {

// typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
// typedef CGAL::Labeled_mesh_domain_3<K> Periodic_mesh_domain;

// // Triangulation
// typedef CGAL::Periodic_3_mesh_triangulation_3<Periodic_mesh_domain>::type Tr;
// typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;



// Kernel
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::FT                                               FT;
typedef K::Point_3                                          Point;
typedef K::Iso_cuboid_3                                     Iso_cuboid;

// Domain
typedef FT (*Function)(const Point&);

// This wrapper is needed to make 'sphere_function' periodic.
typedef CGAL::Periodic_3_function_wrapper<std::function<double(K::Point_3)>, K>      Periodic_function;
typedef CGAL::Implicit_multi_domain_to_labeling_function_wrapper<Periodic_function> Multi_domain_wrapper;

typedef CGAL::Labeled_mesh_domain_3<K>                      Periodic_mesh_domain;

// Triangulation
typedef CGAL::Periodic_3_mesh_triangulation_3<Periodic_mesh_domain>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr>                       C3t3;

// Criteria
typedef CGAL::Mesh_criteria_3<Tr>                           Periodic_mesh_criteria;

// domains
// typedef std::shared_ptr<pygalmesh::DomainBase> domain;

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;

void
generate_periodic_mesh_multiple_domains(
    const std::vector<std::shared_ptr<pygalmesh::DomainBase>> domains,
    const std::vector<std::string> vps,
    const std::string & outfile,
    const std::array<double, 6> bounding_cuboid,
    const bool lloyd,
    const bool odt,
    const bool perturb,
    const bool exude,
    const double max_edge_size_at_feature_edges,
    const double min_facet_angle,
    const double max_radius_surface_delaunay_ball,
    const double max_facet_distance,
    const double max_circumradius_edge_ratio,
    const double max_cell_circumradius,
    const int number_of_copies_in_output,
    const bool verbose,
    const int seed
    )
{
  CGAL::get_default_random() = CGAL::Random(seed);
  const int ndomains = domains.size();

  K::Iso_cuboid_3 cuboid(
      bounding_cuboid[0],
      bounding_cuboid[1],
      bounding_cuboid[2],
      bounding_cuboid[3],
      bounding_cuboid[4],
      bounding_cuboid[5]
      );


//    std::vector<std::function<double(K::Point_3)>> ds[domains.size()]; 
//    std::vector<Periodic_function> funcs;
//    for (int id=0; id<domains.size(); id++)
//    {
//        const double ds[id] = [&](K::Point_3 p){
//            return domains[id]->eval({p.x(), p.y(), p.z()});
//        }
//        funcs.push_back(Periodic_function(ds[id], cuboid));
//    }  

  std::vector<Periodic_function> funcs;
  // wrap the domains
  const auto d0 = [&](K::Point_3 p) {
    return domains[0]->eval({p.x(), p.y(), p.z()});
  };
  funcs.push_back(Periodic_function(d0, cuboid));

  const auto d1 = [&](K::Point_3 p) {
    return domains[1]->eval({p.x(), p.y(), p.z()});
  };
  funcs.push_back(Periodic_function(d1, cuboid));

  if(ndomains>2){
    const auto d2 = [&](K::Point_3 p) {
        return domains[2]->eval({p.x(), p.y(), p.z()});
    };
    funcs.push_back(Periodic_function(d2, cuboid));      
  }

  else if(ndomains>3){
    const auto d3 = [&](K::Point_3 p) {
        return domains[3]->eval({p.x(), p.y(), p.z()});
    };
    funcs.push_back(Periodic_function(d3, cuboid));      
  }
  
  else if(ndomains>4){
    const auto d4 = [&](K::Point_3 p) {
        return domains[4]->eval({p.x(), p.y(), p.z()});
    };
    funcs.push_back(Periodic_function(d4, cuboid));      
  }

  else if(ndomains>5){
    const auto d5 = [&](K::Point_3 p) {
        return domains[5]->eval({p.x(), p.y(), p.z()});
    };
    funcs.push_back(Periodic_function(d5, cuboid));      
  }

  else if(ndomains>6){
    const auto d6 = [&](K::Point_3 p) {
        return domains[6]->eval({p.x(), p.y(), p.z()});
    };
    funcs.push_back(Periodic_function(d6, cuboid));      
  }

  else if(ndomains>7){
    const auto d7 = [&](K::Point_3 p) {
        return domains[7]->eval({p.x(), p.y(), p.z()});
    };
    funcs.push_back(Periodic_function(d7, cuboid));      
  }


  Multi_domain_wrapper multi_domain_function(funcs, vps);
  Periodic_mesh_domain domain(multi_domain_function, cuboid);

  Periodic_mesh_criteria criteria(
      CGAL::parameters::edge_size=max_edge_size_at_feature_edges,
      CGAL::parameters::facet_angle=min_facet_angle,
      CGAL::parameters::facet_size=max_radius_surface_delaunay_ball,
      CGAL::parameters::facet_distance=max_facet_distance,
      CGAL::parameters::cell_radius_edge_ratio=max_circumradius_edge_ratio,
      CGAL::parameters::cell_size=max_cell_circumradius
      );

  // Mesh generation
  if (!verbose) {
    // suppress output
    std::cerr.setstate(std::ios_base::failbit);
  }
  C3t3 c3t3 = CGAL::make_periodic_3_mesh_3<C3t3>(
      domain,
      criteria,
      lloyd ? CGAL::parameters::lloyd() : CGAL::parameters::no_lloyd(),
      odt ? CGAL::parameters::odt() : CGAL::parameters::no_odt(),
      perturb ? CGAL::parameters::perturb() : CGAL::parameters::no_perturb(),
      exude ? CGAL::parameters::exude() : CGAL::parameters::no_exude()
      );
  if (!verbose) {
    std::cerr.clear();
  }

  // Output
  std::ofstream medit_file(outfile);
  CGAL::output_periodic_mesh_to_medit(medit_file, c3t3, number_of_copies_in_output);
  medit_file.close();

  return;
}

} // namespace pygalmesh
