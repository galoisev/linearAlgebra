#ifndef __ENTITIES_H
#define __ENTITIES_H

#include <iostream>
#include <istream>
#include <fstream>
#include <string>
#include <vector>
#include <cassert>



namespace GMSH
{


	class Entities
	{
		/*
	private:
		int _nbrePoints, _nbreLines, _nbreSurfaces, _nbreVolumes;
	
		int _ptsTag;
		double _x, _y, _z;
		int _physicalTag;*/

	public:

		int numberOfpoints, numberOfLine, numberOfSurface, numberOfVolume;

		int points_tag_ep;
		double x_ep, y_ep, z_ep;
		int physical_tag_ep;

		int line_tag_el;
		double xmin_el, ymin_el, zmin_el;
		double xmax_el, ymax_el, zmax_el;
		int physical_tag_el;//physical tag entité linéaire
		int num_boundary_nodes_el;
		int boundary_node_1_el, boundary_node_2_el, boundary_node_3_el, boundary_node_4_el;
		int* boundary_node_x_el;

		int surface_tag_es;
		double xmin_es, ymin_es, zmin_es;
		double xmax_es, ymax_es, zmax_es;
		int physical_tag_es;//physical tag entité linéaire
		int num_boundary_curves_es;
		int boundary_node_1_es, boundary_node_2_es, boundary_node_3_es, boundary_node_4_es;
		int* boundary_node_x_es;


		Entities() {};
		Entities(int points_tag_ep, double x_ep, double y_ep, double z_ep, int physical_tag_ep, int line_tag_el, double xmin_el, double ymin_el, double zmin_el, double xmax_el, double ymax_el, double zmax_el, int physical_tag_el, int num_boundary_nodes_el, int boundary_node_1_el, int boundary_node_2_el, int surface_tag_es, double xmin_es, double ymin_es, double zmin_es, double xmax_es, double ymax_es, double zmax_es, int physical_tag_es, int num_boundary_curves_es) : 
			points_tag_ep(points_tag_ep), x_ep(x_ep), y_ep(y_ep), z_ep(z_ep), physical_tag_ep(physical_tag_ep), line_tag_el(line_tag_el),
			xmin_el(xmin_el), ymin_el(ymin_el), zmin_el(zmin_el), xmax_el(xmax_el), ymax_el(ymax_el), zmax_el(zmax_el), physical_tag_el(physical_tag_el), num_boundary_nodes_el(num_boundary_nodes_el), boundary_node_1_el(boundary_node_1_el), boundary_node_2_el(boundary_node_2_el),
			surface_tag_es(surface_tag_es), xmin_es(xmin_es), ymin_es{ ymin_es }, zmin_es(zmin_es), xmax_es(xmax_es), ymax_es{ ymax_es }, zmax_es(zmax_es),
			physical_tag_es(physical_tag_es), num_boundary_curves_es( num_boundary_curves_es){};
		Entities(const Entities& entities) :
			points_tag_ep(entities.points_tag_ep), x_ep(entities.x_ep), y_ep(entities.y_ep), z_ep(entities.z_ep), physical_tag_ep(entities.physical_tag_ep), line_tag_el(entities.line_tag_el),
			xmin_el(entities.xmin_el), ymin_el(entities.ymin_el), zmin_el(entities.zmin_el), xmax_el(entities.xmax_el), ymax_el(entities.ymax_el), zmax_el(entities.zmax_el), physical_tag_el(entities.physical_tag_el), num_boundary_nodes_el(entities.num_boundary_nodes_el), boundary_node_1_el(entities.boundary_node_1_el), boundary_node_2_el(entities.boundary_node_2_el),
			surface_tag_es(entities.surface_tag_es), xmin_es(entities.xmin_es), ymin_es(entities.ymin_es), zmin_es(entities.zmin_es), xmax_es(entities.xmax_es), ymax_es(entities.ymax_es), zmax_es(entities.zmax_es),
			physical_tag_es(entities.physical_tag_es), num_boundary_curves_es(entities.num_boundary_curves_es){};
		Entities& operator=(const Entities& enti) = default;


	
		friend std::ostream& operator<<(std::ostream& f, const Entities& entities);
		friend std::istream& operator>>(std::istream& f, Entities& entities);


	};



	std::ostream& operator<<(std::ostream& f, const Entities& entities)
	{
		f << entities.numberOfpoints << " " << entities.numberOfLine << " " << entities.numberOfSurface << " " << entities.numberOfVolume << "\n";
		return f;
	}

	std::istream& operator>>(std::istream& f, Entities& entities)
	{
		f >> entities.numberOfpoints >> entities.numberOfLine >> entities.numberOfSurface >> entities.numberOfVolume;
		return f;

	}


}//end namespace GMSH

#endif // !__ENTITIES_H
