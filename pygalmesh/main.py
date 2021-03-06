import math
import os
import tempfile

import meshio
import numpy
from _pygalmesh import (
    SizingFieldBase,
    _generate_2d,
    _generate_from_inr,
    _generate_from_inr_with_subdomain_sizing,
    _generate_from_inr_with_features,
    _generate_from_off,
    _generate_mesh,
    _generate_periodic_mesh,
    _generate_periodic_mesh_multiple_domains,
    _generate_surface_mesh,
    _remesh_surface,
)


class Wrapper(SizingFieldBase):
    def __init__(self, f):
        self.f = f
        super().__init__()

    def eval(self, x):
        return self.f(x)


def generate_mesh(
    domain,
    extra_feature_edges=None,
    bounding_sphere_radius=0.0,
    lloyd=False,
    odt=False,
    perturb=True,
    exude=True,
    max_edge_size_at_feature_edges=0.0,
    min_facet_angle=0.0,
    max_radius_surface_delaunay_ball=0.0,
    max_facet_distance=0.0,
    max_circumradius_edge_ratio=0.0,
    max_cell_circumradius=0.0,
    verbose=True,
    seed=0,
):
    """
    From <https://doc.cgal.org/latest/Mesh_3/classCGAL_1_1Mesh__criteria__3.html>:

    max_edge_size_at_feature_edges:
        a scalar field (resp. a constant) providing a space varying (resp. a uniform)
        upper bound for the lengths of curve edges. This parameter has to be set to a
        positive value when 1-dimensional features protection is used.
    min_facet_angle:
        a lower bound for the angles (in degrees) of the surface mesh facets.
    max_radius_surface_delaunay_ball:
        a scalar field (resp. a constant) describing a space varying (resp. a uniform)
        upper-bound or for the radii of the surface Delaunay balls.
    max_facet_distance:
        a scalar field (resp. a constant) describing a space varying (resp. a uniform)
        upper bound for the distance between the facet circumcenter and the center of
        its surface Delaunay ball.
    max_circumradius_edge_ratio:
        an upper bound for the radius-edge ratio of the mesh tetrahedra.
    max_cell_circumradius:
        a scalar field (resp. a constant) describing a space varying (resp. a uniform)
        upper-bound for the circumradii of the mesh tetrahedra.
    """
    extra_feature_edges = [] if extra_feature_edges is None else extra_feature_edges

    fh, outfile = tempfile.mkstemp(suffix=".mesh")
    os.close(fh)

    def _select(obj):
        if isinstance(obj, float):
            return obj, None
        assert callable(obj)
        return -1.0, Wrapper(obj)

    (
        max_edge_size_at_feature_edges_value,
        max_edge_size_at_feature_edges_field,
    ) = _select(max_edge_size_at_feature_edges)
    max_cell_circumradius_value, max_cell_circumradius_field = _select(
        max_cell_circumradius
    )
    (
        max_radius_surface_delaunay_ball_value,
        max_radius_surface_delaunay_ball_field,
    ) = _select(max_radius_surface_delaunay_ball)
    max_facet_distance_value, max_facet_distance_field = _select(max_facet_distance)

    # if feature_edges:
    #     if max_edge_size_at_feature_edges == 0.0:
    #         raise ValueError(
    #             "Need a positive max_edge_size_at_feature_edges bound if feature_edges are present."
    #         )
    # elif max_edge_size_at_feature_edges != 0.0:
    #     warnings.warn(
    #         "No feature edges. The max_edge_size_at_feature_edges argument has no effect."
    #     )

    _generate_mesh(
        domain,
        outfile,
        extra_feature_edges=extra_feature_edges,
        bounding_sphere_radius=bounding_sphere_radius,
        lloyd=lloyd,
        odt=odt,
        perturb=perturb,
        exude=exude,
        max_edge_size_at_feature_edges_value=max_edge_size_at_feature_edges_value,
        max_edge_size_at_feature_edges_field=max_edge_size_at_feature_edges_field,
        min_facet_angle=min_facet_angle,
        max_radius_surface_delaunay_ball_value=max_radius_surface_delaunay_ball_value,
        max_radius_surface_delaunay_ball_field=max_radius_surface_delaunay_ball_field,
        max_facet_distance_value=max_facet_distance_value,
        max_facet_distance_field=max_facet_distance_field,
        max_circumradius_edge_ratio=max_circumradius_edge_ratio,
        max_cell_circumradius_value=max_cell_circumradius_value,
        max_cell_circumradius_field=max_cell_circumradius_field,
        verbose=verbose,
        seed=seed,
    )

    mesh = meshio.read(outfile)
    os.remove(outfile)
    return mesh


def generate_2d(
    points,
    constraints,
    B=math.sqrt(2),
    max_edge_size=0.0,
    num_lloyd_steps=0,
):
    # some sanity checks
    points = numpy.asarray(points)
    constraints = numpy.asarray(constraints)
    assert numpy.all(constraints >= 0)
    assert numpy.all(constraints < len(points))
    # make sure there are no edges of 0 length
    edges = points[constraints[:, 0]] - points[constraints[:, 1]]
    length2 = numpy.einsum("ij,ij->i", edges, edges)
    if numpy.any(length2 < 1.0e-15):
        raise RuntimeError("Constraint of (near)-zero length.")

    points, cells = _generate_2d(
        points,
        constraints,
        B,
        max_edge_size,
        num_lloyd_steps,
    )
    return meshio.Mesh(numpy.array(points), {"triangle": numpy.array(cells)})


def generate_periodic_mesh(
    domain,
    bounding_cuboid,
    extra_feature_edges=None,
    lloyd=False,
    odt=False,
    perturb=True,
    exude=True,
    max_edge_size_at_feature_edges=0.0,
    min_facet_angle=0.0,
    max_radius_surface_delaunay_ball=0.0,
    max_facet_distance=0.0,
    max_circumradius_edge_ratio=0.0,
    max_cell_circumradius=0.0,
    number_of_copies_in_output=1,
    verbose=True,
    make_periodic=False,
    seed=0,
):

    fh, outfile = tempfile.mkstemp(suffix=".mesh")
    os.close(fh)

    assert number_of_copies_in_output in [1, 2, 4, 8]

    extra_feature_edges = [] if extra_feature_edges is None else extra_feature_edges

    _generate_periodic_mesh(
        domain,
        outfile,
        bounding_cuboid,
        extra_feature_edges=extra_feature_edges,
        lloyd=lloyd,
        odt=odt,
        perturb=perturb,
        exude=exude,
        max_edge_size_at_feature_edges=max_edge_size_at_feature_edges,
        min_facet_angle=min_facet_angle,
        max_radius_surface_delaunay_ball=max_radius_surface_delaunay_ball,
        max_facet_distance=max_facet_distance,
        max_circumradius_edge_ratio=max_circumradius_edge_ratio,
        max_cell_circumradius=max_cell_circumradius,
        number_of_copies_in_output=number_of_copies_in_output,
        verbose=verbose,
        make_periodic=make_periodic,
        seed=seed,
    )

    mesh = meshio.read(outfile)
    os.remove(outfile)
    
    return mesh


def generate_periodic_mesh_multiple_domains(
    domains,
    signs,
    bounding_cuboid,
    lloyd=False,
    odt=False,
    perturb=True,
    exude=True,
    max_edge_size_at_feature_edges=0.0,
    min_facet_angle=0.0,
    max_radius_surface_delaunay_ball=0.0,
    max_facet_distance=0.0,
    max_circumradius_edge_ratio=0.0,
    max_cell_circumradius=0.0,
    number_of_copies_in_output=1,
    verbose=True,
    outfile=None,
    seed=0,
):
    # deal with the outputfile generated by cgal
    remove_file = False
    if outfile is None:
        fh, outfile = tempfile.mkstemp(suffix=".mesh")
        os.close(fh)
        remove_file = True

    # we need more than 1 domain
    assert len(domains) > 1
    # Cgal only support up to 8 domains
    assert len(domains) < 9
    # cgal only support 1,2,4,8 copies
    assert number_of_copies_in_output in [1, 2, 4, 8]

    _generate_periodic_mesh_multiple_domains(
        domains,
        signs,
        outfile,
        bounding_cuboid,
        lloyd=lloyd,
        odt=odt,
        perturb=perturb,
        exude=exude,
        max_edge_size_at_feature_edges=max_edge_size_at_feature_edges,
        min_facet_angle=min_facet_angle,
        max_radius_surface_delaunay_ball=max_radius_surface_delaunay_ball,
        max_facet_distance=max_facet_distance,
        max_circumradius_edge_ratio=max_circumradius_edge_ratio,
        max_cell_circumradius=max_cell_circumradius,
        number_of_copies_in_output=number_of_copies_in_output,
        verbose=verbose,
        seed=seed,
    )
    
    mesh = meshio.read(outfile)

    if remove_file:
        os.remove(outfile)
    
    return mesh


def generate_surface_mesh(
    domain,
    bounding_sphere_radius=0.0,
    min_facet_angle=0.0,
    max_radius_surface_delaunay_ball=0.0,
    max_facet_distance=0.0,
    verbose=True,
    seed=0,
):
    fh, outfile = tempfile.mkstemp(suffix=".off")
    os.close(fh)

    _generate_surface_mesh(
        domain,
        outfile,
        bounding_sphere_radius=bounding_sphere_radius,
        min_facet_angle=min_facet_angle,
        max_radius_surface_delaunay_ball=max_radius_surface_delaunay_ball,
        max_facet_distance=max_facet_distance,
        verbose=verbose,
        seed=seed,
    )

    mesh = meshio.read(outfile)
    os.remove(outfile)
    return mesh


def generate_volume_mesh_from_surface_mesh(
    filename,
    lloyd=False,
    odt=False,
    perturb=True,
    exude=True,
    max_edge_size_at_feature_edges=0.0,
    min_facet_angle=0.0,
    max_radius_surface_delaunay_ball=0.0,
    max_facet_distance=0.0,
    max_circumradius_edge_ratio=0.0,
    max_cell_circumradius=0.0,
    verbose=True,
    reorient=False,
    seed=0,
):
    mesh = meshio.read(filename)

    fh, off_file = tempfile.mkstemp(suffix=".off")
    os.close(fh)
    meshio.write(off_file, mesh)

    fh, outfile = tempfile.mkstemp(suffix=".mesh")
    os.close(fh)

    _generate_from_off(
        off_file,
        outfile,
        lloyd=lloyd,
        odt=odt,
        perturb=perturb,
        exude=exude,
        max_edge_size_at_feature_edges=max_edge_size_at_feature_edges,
        min_facet_angle=min_facet_angle,
        max_radius_surface_delaunay_ball=max_radius_surface_delaunay_ball,
        max_facet_distance=max_facet_distance,
        max_circumradius_edge_ratio=max_circumradius_edge_ratio,
        max_cell_circumradius=max_cell_circumradius,
        verbose=verbose,
        reorient=reorient,
        seed=seed,
    )

    mesh = meshio.read(outfile)
    os.remove(off_file)
    os.remove(outfile)
    return mesh

def generate_from_inr_with_features(
    inr_filename,
    extra_feature_edges = None,
    bbox_feature = False,
    lloyd=False,
    odt=False,
    perturb=True,
    exude=True,
    max_edge_size_at_feature_edges=0.0,
    min_facet_angle=0.0,
    max_radius_surface_delaunay_ball=0.0,
    max_facet_distance=0.0,
    max_circumradius_edge_ratio=0.0,
    max_cell_circumradius=0.0,
    verbose=True,
    seed=0,
):

    extra_feature_edges = [] if extra_feature_edges is None else extra_feature_edges

    fh, outfile = tempfile.mkstemp(suffix=".mesh")
    os.close(fh)

    assert isinstance(max_cell_circumradius, float)

    _generate_from_inr_with_features(
        inr_filename,
        outfile,
        extra_feature_edges=extra_feature_edges,
        bbox_feature=bbox_feature,
        lloyd=lloyd,
        odt=odt,
        perturb=perturb,
        exude=exude,
        max_edge_size_at_feature_edges=max_edge_size_at_feature_edges,
        min_facet_angle=min_facet_angle,
        max_radius_surface_delaunay_ball=max_radius_surface_delaunay_ball,
        max_facet_distance=max_facet_distance,
        max_circumradius_edge_ratio=max_circumradius_edge_ratio,
        max_cell_circumradius=max_cell_circumradius,
        verbose=verbose,
        seed=seed,
    )

    mesh = meshio.read(outfile)
    os.remove(outfile)
    return mesh


def generate_from_inr(
    inr_filename,
    lloyd=False,
    odt=False,
    perturb=True,
    exude=True,
    max_edge_size_at_feature_edges=0.0,
    min_facet_angle=0.0,
    max_radius_surface_delaunay_ball=0.0,
    max_facet_distance=0.0,
    max_circumradius_edge_ratio=0.0,
    max_cell_circumradius=0.0,
    verbose=True,
    seed=0,
):
    fh, outfile = tempfile.mkstemp(suffix=".mesh")
    os.close(fh)

    if isinstance(max_cell_circumradius, float):
        
        _generate_from_inr(
            inr_filename,
            outfile,
            lloyd=lloyd,
            odt=odt,
            perturb=perturb,
            exude=exude,
            max_edge_size_at_feature_edges=max_edge_size_at_feature_edges,
            min_facet_angle=min_facet_angle,
            max_radius_surface_delaunay_ball=max_radius_surface_delaunay_ball,
            max_facet_distance=max_facet_distance,
            max_circumradius_edge_ratio=max_circumradius_edge_ratio,
            max_cell_circumradius=max_cell_circumradius,
            verbose=verbose,
            seed=seed,
        )
    else:
        assert isinstance(max_cell_circumradius, dict)
        if "default" in max_cell_circumradius.keys():
            default_max_cell_circumradius = max_cell_circumradius.pop("default")
        else:
            default_max_cell_circumradius = 0.0

        max_cell_circumradiuss = list(max_cell_circumradius.values())
        subdomain_labels = list(max_cell_circumradius.keys())


        _generate_from_inr_with_subdomain_sizing(
            inr_filename,
            outfile,
            default_max_cell_circumradius,
            max_cell_circumradiuss,
            subdomain_labels,
            lloyd=lloyd,
            odt=odt,
            perturb=perturb,
            exude=exude,
            max_edge_size_at_feature_edges=max_edge_size_at_feature_edges,
            min_facet_angle=min_facet_angle,
            max_radius_surface_delaunay_ball=max_radius_surface_delaunay_ball,
            max_facet_distance=max_facet_distance,
            max_circumradius_edge_ratio=max_circumradius_edge_ratio,
            verbose=verbose,
            seed=seed,
        )

    mesh = meshio.read(outfile)
    os.remove(outfile)
    return mesh


def remesh_surface(
    filename,
    max_edge_size_at_feature_edges=0.0,
    min_facet_angle=0.0,
    max_radius_surface_delaunay_ball=0.0,
    max_facet_distance=0.0,
    verbose=True,
    seed=0,
):
    mesh = meshio.read(filename)

    fh, off_file = tempfile.mkstemp(suffix=".off")
    os.close(fh)
    meshio.write(off_file, mesh)

    fh, outfile = tempfile.mkstemp(suffix=".off")
    os.close(fh)

    _remesh_surface(
        off_file,
        outfile,
        max_edge_size_at_feature_edges=max_edge_size_at_feature_edges,
        min_facet_angle=min_facet_angle,
        max_radius_surface_delaunay_ball=max_radius_surface_delaunay_ball,
        max_facet_distance=max_facet_distance,
        verbose=verbose,
        seed=seed,
    )

    mesh = meshio.read(outfile)
    os.remove(off_file)
    os.remove(outfile)
    return mesh


def save_inr(vol, h, fname):
    """
    Save a volume (described as a numpy array) to INR format.
    Code inspired by iso2mesh (http://iso2mesh.sf.net) by Q. Fang
    INPUTS:
    - vol: volume as numpy array
    - h: voxel sizes as list or numpy array
    - fname: filename for saving the inr file
    """
    fid = open(fname, "wb")

    btype, bitlen = {
        "uint8": ("unsigned fixed", 8),
        "uint16": ("unsigned fixed", 16),
        "float32": ("float", 32),
        "float64": ("float", 64),
    }[vol.dtype.name]

    header = (
        "#INRIMAGE-4#{8:s}\nXDIM={0:d}\nYDIM={1:d}\nZDIM={2:d}\nVDIM=1\nTYPE={3:s}\n"
        + "PIXSIZE={4:d} bits\nCPU=decm\nVX={5:f}\nVY={6:f}\nVZ={7:f}\n"
    ).format(*vol.shape, btype, bitlen, h[0], h[1], h[2], "{")

    header = header + "\n" * (256 - 4 - len(header)) + "##}\n"

    fid.write(header.encode("ascii"))
    fid.write(vol.tobytes(order="F"))


def generate_from_array(
    vol,
    h,
    extra_feature_edges = None,
    bbox_feature = False,
    lloyd=False,
    odt=False,
    perturb=True,
    exude=True,
    max_edge_size_at_feature_edges=0.0,
    min_facet_angle=0.0,
    max_radius_surface_delaunay_ball=0.0,
    max_cell_circumradius=0.0,
    max_facet_distance=0.0,
    max_circumradius_edge_ratio=0.0,
    verbose=True,
    seed=0,
):
    assert vol.dtype in ["uint8", "uint16"]
    fh, inr_filename = tempfile.mkstemp(suffix=".inr")
    os.close(fh)
    save_inr(vol, h, inr_filename)

    if (bbox_feature is False) & (extra_feature_edges is None):
        mesh = generate_from_inr(
            inr_filename,
            lloyd,
            odt,
            perturb,
            exude,
            max_edge_size_at_feature_edges,
            min_facet_angle,
            max_radius_surface_delaunay_ball,
            max_facet_distance,
            max_circumradius_edge_ratio,
            max_cell_circumradius,
            verbose,
            seed,
        )
    else :
        mesh = generate_from_inr_with_features(
            inr_filename,
            extra_feature_edges,
            bbox_feature,
            lloyd,
            odt,
            perturb,
            exude,
            max_edge_size_at_feature_edges,
            min_facet_angle,
            max_radius_surface_delaunay_ball,
            max_facet_distance,
            max_circumradius_edge_ratio,
            max_cell_circumradius,
            verbose,
            seed,
        )

    os.remove(inr_filename)
    return mesh


