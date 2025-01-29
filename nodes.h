#ifndef __NODES_H
#define __NODES_H


#include <iostream>
#include <istream>
#include <fstream>
#include <string>
#include <vector>
#include <cassert>


namespace GMSH
{

	class Nodes
	{
	public:

		int num_entity_blocks{}, num_total_nodes{}, min_node_tag{}, max_node_tag{};
		int entity_dim{}, entity_tag{}, parametric{}, num_nodes_in_block{};
		int node_tag{};
		double x{}, y{}, z{};

		Nodes() {};
		Nodes(int nnode_tag, double xx, double yy, double zz):node_tag(nnode_tag),x(xx),y(yy),z(zz){};
		Nodes(int num_entity_blocks, int num_total_nodes, int min_node_tag, int max_node_tag, int entity_dim, int entity_tag, int parametric, int num_nodes_in_block, int node_tag, double x, double y, double z) :
			num_entity_blocks(num_entity_blocks), num_total_nodes(num_total_nodes), min_node_tag(min_node_tag), max_node_tag(max_node_tag), entity_dim(entity_dim), entity_tag(entity_tag), parametric(parametric), 
			num_nodes_in_block(num_nodes_in_block), node_tag(node_tag), x(x), y(y), z(z) {};
		Nodes(Nodes& nodes) : num_entity_blocks(nodes.num_entity_blocks), num_total_nodes(nodes.num_total_nodes), min_node_tag(nodes.min_node_tag), max_node_tag(nodes.max_node_tag), entity_dim(nodes.entity_dim), 
			entity_tag(nodes.entity_tag), parametric(nodes.parametric),	num_nodes_in_block(nodes.num_nodes_in_block), node_tag(nodes.node_tag), x(nodes.x), y(nodes.y), z(nodes.z) {};
		Nodes(Nodes a, Nodes b):x(b.x - a.x), y(b.y - a.y), z(b.z - a.z) {}
		Nodes& operator=(const Nodes& nodes) = default;

		auto  operator<=>(const Nodes& nodes)const = default;

		Nodes operator+(Nodes P) { return Nodes(num_entity_blocks, num_total_nodes, min_node_tag, max_node_tag, entity_dim, entity_tag, parametric, num_nodes_in_block, node_tag, x + P.x, y + P.y, z + P.z); }
		Nodes& operator+=(Nodes P)
		{
			x += P.x;
			y += P.y;
			z += P.z;
			return (*this);
		}
		Nodes operator-(Nodes P) { return Nodes(num_entity_blocks, num_total_nodes, min_node_tag, max_node_tag, entity_dim, entity_tag, parametric, num_nodes_in_block, node_tag, x - P.x, y - P.y, z - P.z); }
		Nodes& operator-=(Nodes P)
		{
			x -= P.x;
			y -= P.y;
			z -= P.z;
			return (*this);
		}
		double  operator,(Nodes P) { return x * P.x + y * P.y + z * P.z; }
		Nodes operator^(Nodes P) { return Nodes(num_entity_blocks, num_total_nodes, min_node_tag, max_node_tag, entity_dim, entity_tag, parametric, num_nodes_in_block, node_tag, y * P.z - z * P.y, z * P.x - x * P.z, x * P.y - y * P.x); }
		Nodes operator*(double c)
		{
			return Nodes(num_entity_blocks, num_total_nodes, min_node_tag, max_node_tag, entity_dim, entity_tag, parametric, num_nodes_in_block, node_tag, c * x, c * y, c * z);
		}
		Nodes operator/(double c)
		{
			assert(c != 0.0);
			return Nodes(num_entity_blocks, num_total_nodes, min_node_tag, max_node_tag, entity_dim, entity_tag, parametric, num_nodes_in_block, node_tag, x / c, y / c, z / c);
		}

		friend Nodes centroid(Nodes A, Nodes B, Nodes C)
		{
			Nodes G;
			G.x = (A.x + B.x + C.x) / 3;
			G.y = (A.y + B.y + C.y) / 3;
			G.z = (A.z + B.z + C.z) / 3;
			return G;
		}

		friend Nodes normal(Nodes A, Nodes B, Nodes C)
		{
			Nodes N;//normalize
			Nodes AB(A, B), AC(A, C);
			Nodes n = AB.operator^(AC);
			//assert((n.x != 0) && (n.y != 0) && (n.z != 0));
			N = n.operator*(1.0 / n.normeL2());
			return N;
		}

		double normeL2()
		{
			double s{ 0 };
			s = sqrt(x * x + y * y + z * z);
			return s;
		}

		
		friend std::ostream& operator<<(std::ostream& f, const Nodes& nodes);
		friend std::istream& operator>>(std::istream& f, Nodes& nodes);		

	};


	
	std::ostream& operator<<(std::ostream& f, const Nodes& nodes)
	{
		f << "\n";
		//f << "\tNumber of nodes: " << nodes.num_total_nodes << "\n";
		f << "\tnode number: " << nodes.node_tag << "\t\tx: " << nodes.x << "\ty: " << nodes.y << "\tz: " << nodes.z << "\n";
		f << "\n";
		return f;
	}


	std::istream& operator>>(std::istream& f, Nodes& nodes)
	{
		f >> nodes.node_tag >> nodes.x >> nodes.y >> nodes.z;
		//f >> nodes.node_tag;
		//f >> nodes.x >> nodes.y >> nodes.z;
		return f;
	}
	





}//end namespace GMSH

#endif // !__NODES_H
