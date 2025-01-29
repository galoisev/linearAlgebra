#ifndef _ELEMENTS_H
#define _ELEMENTS_H

#include "meshFormat.h"
#include "entities.h"
#include "nodes.h"



namespace GMSH
{

	typedef double dble;

	class Elements
	{
	public:
		int num_entity_blocks{}, num_total_elements{}, min_element_tag{}, max_element_tag{};
		int entity_dim{}, entity_tag{}, element_type{}, num_elements_in_block{};
		int element_tag{}, node_tag_1{}, node_tag_2{}, node_tag_3{};
		double area{ 0 }, perimeter = { 0 }, step_k = { 0 };
		Nodes CG, NG;
		Nodes* nodes[3]{};

		std::vector<Nodes*>list_nodes{};


		


		Elements() {};
		Elements(int nnum_entity_blocks, int nnum_total_elements, int mmin_element_tag, int mmax_element_tag, int eentity_dim, int eentity_tag, int eelement_type, int nnum_elements_in_block, int eelement_tag, int nnode_tag_1, int nnode_tag_2, int nnode_tag_3):
			num_entity_blocks(nnum_entity_blocks), num_total_elements(nnum_total_elements), min_element_tag(mmin_element_tag), max_element_tag(mmax_element_tag), entity_dim(eentity_dim), entity_tag(eentity_tag), element_type(eelement_type), num_elements_in_block(nnum_elements_in_block), element_tag(eelement_tag), node_tag_1(nnode_tag_1), node_tag_2(nnode_tag_2), node_tag_3(nnode_tag_3){};
		Elements(const Elements& e) :num_entity_blocks(e.num_entity_blocks), num_total_elements(e.num_total_elements), min_element_tag(e.min_element_tag), max_element_tag(e.max_element_tag), entity_dim(e.entity_dim), entity_tag(e.entity_tag), element_type(e.element_type), num_elements_in_block(e.num_elements_in_block), element_tag(e.element_tag), node_tag_1(e.node_tag_1), node_tag_2(e.node_tag_2), node_tag_3(e.node_tag_3) {};
		Elements(int element_tag, Nodes CG, double area, double perimeter) :element_tag(element_tag), CG(CG), area(area), perimeter(perimeter) {}
		Elements& operator=(const Elements& e) = default;
		auto  operator<=>(const Elements& tab)const = default;

		//~Elements();
		/*
		void allocArrays()
		{
			nodes = new Nodes * [3];
			for (int i = 0; i < num_total_elements; i++)
			{
				nodes[i] = new Nodes[num_total_elements];
			}
		}*/


		void set(Nodes* v0, int nnode_tag_1, int nnode_tag_2, int nnode_tag_3, int ir)
		{
			this->node_tag_1 = nnode_tag_1;
			this->node_tag_2 = nnode_tag_2;
			this->node_tag_3 = nnode_tag_3;

			nodes[0] = v0 + this->node_tag_1;
			//std::cout << "------->\tnodes[0] = " << nodes[0]->node_tag << ","  << nodes[0]->x << "," << nodes[0]->y << "," << nodes[0]->z << "\n";
			nodes[1] = v0 + this->node_tag_2;
			//std::cout << "------->\tnodes[1] = " << nodes[1]->node_tag << ","  << nodes[1]->x << "," << nodes[1]->y << "," << nodes[1]->z << "\n";
			nodes[2] = v0 + this->node_tag_3;
			//std::cout << "------->\tnodes[2] = " << nodes[2]->node_tag << ","  << nodes[2]->x << "," << nodes[2]->y << "," << nodes[2]->z << "\n";


			Nodes A(*nodes[0]), B(*nodes[1]), C(*nodes[2]);

			list_nodes.push_back(&A);
			list_nodes.push_back(&B);
			list_nodes.push_back(&C);

			/*
			Nodes AB(A, B);
			Nodes AC(A, C);
			Nodes BC(B, C);
			Nodes u = AB.operator^(AC);		
			Nodes prod_vect(u);
			double l_AB = AB.normeL2();
			double l_AC = AC.normeL2();
			double l_BC = BC.normeL2();
			double l_t = (l_AB + l_AC + l_BC);//périmètre du triangle t
			this->perimeter = l_t;*/

			dble perimeter_element = length(A, B, C);
			this->perimeter = perimeter_element;

			/*
			double area_calc = 0.5*u.normeL2();
			this->area = area_calc;*/

			dble area_element = surface(A, B, C);
			this->area = area_element;

			Nodes G = centroid(A, B, C);
			this->CG = G;

			Nodes NCG = normal(A, B, C);
			this->NG = NCG;

			// cas où les triangles sont pratiquement tous équilatéruax !
			double rayon_k = circumcircle(A, B, C);
			this->step_k = 2*rayon_k;

			//informations
			/*
			std::cout << this->element_tag << "\n";
			std::cout << "area= " << this->area << "," << " perimeter= " << this->perimeter << "\n";
			std::cout << "CG=" << this->element_tag << ")= " << CG.x << "," << CG.y << "," << CG.z << "\n\n";
			*/

			delete v0;

		}


		dble length(Nodes A, Nodes B, Nodes C)
		{
			Nodes AB(A, B);
			Nodes AC(A, C);
			Nodes BC(B, C);
			Nodes u = AB.operator^(AC);
			Nodes prod_vect(u);
			dble l_AB = AB.normeL2();
			dble l_AC = AC.normeL2();
			dble l_BC = BC.normeL2();
			return (l_AB + l_AC + l_BC);//périmètre du triangle tl_t
		}



		dble surface(Nodes A, Nodes B, Nodes C)
		{
			Nodes AB(A, B);
			Nodes AC(A, C);
			Nodes BC(B, C);
			Nodes u = AB.operator^(AC);
			Nodes prod_vect(u);
			return 0.5 * u.normeL2();
		}


		Nodes centroid(Nodes A, Nodes B, Nodes C)
		{
			Nodes G;
			G.x = (A.x + B.x + C.x) / 3;
			G.y = (A.y + B.y + C.y) / 3;
			G.z = (A.z + B.z + C.z) / 3;
			return G;
		}

		Nodes normal(Nodes A, Nodes B, Nodes C)
		{
			Nodes N;//normalize
			Nodes AB(A, B), AC(A, C);
			Nodes n = AB.operator^(AC);
			if ((n.x == 0) && (n.y == 0) && (n.z == 0))
			{
				std::cerr << "\tZERO DIVISION ERROR! Normal vector is zero." << "\n";
				std::system("PAUSE");
				std::exit(1);
			}				
			N = n.operator*(1.0 / n.normeL2());
			return N;
		}


		double h(Nodes& A, Nodes& B, Nodes& C);
		double circumcircle(Nodes& A, Nodes& B, Nodes& C);



		friend std::ostream& operator<<(std::ostream& f, const Elements& elements);
		friend std::istream& operator>>(std::istream& f, Elements& elements);

	};

	/*
	Elements::~Elements()
	{
		for (int i = 0; i < 3; i++)
		{
			delete[]nodes[i];
		}
		delete[]nodes;
	}
	*/



	double Elements::circumcircle(Nodes& A, Nodes& B, Nodes& C)
	{
		Nodes AB(A, B), AC(A, C), BC(B, C);
		double c = AB.normeL2(); //c: longueur du côté [A,B) du triangle (A,B,C)
		double b = AC.normeL2(); //b: longueur du côté [A,C) du triangle (A,B,C)
		double a = BC.normeL2(); //a: longueur du côté [B,C) du triangle (A,B,C)
		Nodes u = AB.operator^(AC);
		Nodes prod_vect(u);
		double area_calc = 0.5 * u.normeL2();//aire du triangle(A,B,C)
		return ((a * b * c) / (4 * area_calc));
	}

	double Elements::h(Nodes& A, Nodes& B, Nodes& C)
	{
		double step_K = { 0.0 };//pas de l'èléments K
		Nodes G = centroid(A, B, C);
		Nodes GA(G, A);
		step_K = 2.0 * GA.normeL2();//diametre
		return step_K;
	}





	std::ostream& operator<<(std::ostream& f, const Elements& elements)
	{
		//f << "\num_entity_blocks:" << elements.num_entity_blocks << "\tnum_total_elements:" << elements.num_total_elements << "\tmin_element_tag:" << elements.min_element_tag << "\tmax_element_tag:" << elements.max_element_tag << "\n";
		f << elements.element_tag << "," << elements.area << "," << elements.CG << "," << elements.NG << "," << elements.perimeter << "\n";
		return f;
	}

	std::istream& operator>>(std::istream& f, Elements elements)
	{	
		f >> elements.num_entity_blocks >> elements.num_total_elements >> elements.min_element_tag >> elements.max_element_tag ;
		/*
		elements.i++;
		f >> elements.i;
		*/

		return f;
	}

}//end namespace GMSH

#endif // ! END _ELEMENTS_H
