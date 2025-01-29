#ifndef __MESH_H
#define __MESH_H


#include "elements.h"
#include <vector>
#include <algorithm>
#include <elements.h>
#include <chrono>
#include <thread>




namespace GMSH
{

	void displayProgressBar(int progress, int total) {
		int barWidth = 50;  // Largeur de la barre de progression
		float ratio = static_cast<float>(progress) / total;
		int pos = static_cast<int>(barWidth * ratio);

		std::cout << "\t" << "[";
		for (int i = 0; i < barWidth; ++i) {
			if (i < pos) std::cout << "=";
			else if (i == pos) std::cout << ">";
			else std::cout << " ";
		}
		std::cout << "] " << int(ratio * 100.0) << " % of " << total << ".\r";
		std::cout.flush();
	}


	void show_progress_bar(int time, const std::string& message, char symbol)
	{
		std::string progress_bar;
		const double progress_level = 1.42;

		std::cout << message << "\n\n";

		for (double percentage = 0; percentage <= 100; percentage += progress_level)
		{
			progress_bar.insert(0, 1, symbol);
			std::cout << "\r [" << std::ceil(percentage) << '%' << "] " << progress_bar;
			std::this_thread::sleep_for(std::chrono::milliseconds(time));
		}
		std::cout << "\n\n";
	}


	template<class T>
	void read(std::string message, std::ifstream& file, T& classe)
	{
		std::string line{};
		while (std::getline(file, line))
		{
			if (line == message)
			{
				// Lire les trois valeurs après "$MeshFormat"
				file >> classe;
				break;  // Quitter la boucle après avoir trouvé et lu les valeurs
			}
		}
	}


	template<class T>
	void printVector(std::vector<T>& v)
	{
		for (auto&& x : v)
		{
			std::cout << x << "\n";
		}
	}


	// Function to delete pointers in the vector
	void clearVector(std::vector<Nodes*>& vptr_nodes) {
		for (auto node_ptr : vptr_nodes) {
			delete node_ptr;  // Deallocate the memory for each Nodes object
		}
		vptr_nodes.clear();  // Clear the vector
	}






	class Mesh
	{
	public:
		std::string filename{};
		int num_total_elements{};
		double hx{}, hy{};
		double area{}, step{};
		Nodes* ptr_nodes{};
		Elements* ptr_elements{};


		std::vector<Nodes*>list_ptr_nodes_cg{}, ptr_nodes_n{}, list_ptr_nodes{};
		std::vector<Nodes>vnodes{};

		std::vector<double>list_step{};
		std::vector<Nodes>list_cg_k{};

		std::vector<Elements*>list_elements{};






		Mesh() {};
		Mesh(std::string& ffilename) :filename(ffilename)
		{
			std::ifstream file(filename);
			if (!file)
			{
				std::cerr << "\n\n\n\t****** FILE NOT FOUND ! ******\n";
				std::system("PAUSE");
				exit(1);
			}
			else
			{
				std::system("cls");
				/*
				std::cout << "\tDEBUT lecture FICHIER ! " << "\n";
				std::cout << "\n\n\n";
				*/

				std::string line;
					






				//std::system("PAUSE");
				//std::cout << "\ttraitement MESH FORMAT " << "\n";
				//$MeshFormat : spécifie la version du format, ici 4.1.
				// Lire les lignes jusqu'à la ligne qui contient "$MeshFormat"
				GMSH::MeshFormat mshf;
				while (std::getline(file, line))
				{
					if (line == "$MeshFormat")
					{
						file >> mshf;
						break;  // Quitter la boucle après avoir trouvé et lu les valeurs
					}
				}
				//std::cout << mshf << "\n";


				//std::system("PAUSE");
				//std::cout << "\ttraitement ENTITIES " << "\n";
				GMSH::Entities entities{};
				GMSH::Entities* ptr_entities{};
				while (std::getline(file, line))
				{
					if (line == "$Entities")
					{
						file >> entities;
						//std::cout << entities << "\n";

						/*
						 *  nombre total d'entités de chaque type
						 * ---------------------------------------
						 *  4 points géométriques(entités de dimension 0).
						 *  4 lignes(entités de dimension 1).
						 * 	1 surface(entité de dimension 2).
						 * 	0 volumes(entités de dimension 3).
						*/


						ptr_entities = new GMSH::Entities[entities.numberOfpoints];

						for (int i = 0; i < entities.numberOfpoints; i++)
						{
							file >> ptr_entities[i].points_tag_ep >> ptr_entities[i].x_ep >> ptr_entities[i].y_ep >> ptr_entities[i].z_ep >> ptr_entities[i].physical_tag_ep;
							//std::cout << ptr_entities[i].points_tag_ep << "," << ptr_entities[i].x_ep << "," << ptr_entities[i].y_ep << "," << ptr_entities[i].z_ep << "," << ptr_entities[i].physical_tag_ep << "\n";
							//file >> entities.points_tag_ep >> entities.x_ep >> entities.y_ep >> entities.z_ep >> entities.physical_tag_ep;
							//std::cout << entities.points_tag_ep << "," << entities.x_ep << "," << entities.y_ep << "," << entities.z_ep << "," << entities.physical_tag_ep << "\n";
						}

						//int* ptr_curve_tag;
						for (int i = 0; i < entities.numberOfLine; i++)
						{
							file >> ptr_entities[i].line_tag_el >> ptr_entities[i].xmin_el >> ptr_entities[i].ymin_el >> ptr_entities[i].zmin_el >> ptr_entities[i].xmax_el >> ptr_entities[i].ymax_el >> ptr_entities[i].zmax_el >> ptr_entities[i].physical_tag_el >> ptr_entities[i].num_boundary_nodes_el >> ptr_entities[i].boundary_node_1_el >> ptr_entities[i].boundary_node_2_el; // >> ptr_entities[i].boundary_node_3_el;
							//std::cout << ptr_entities[i].line_tag_el << "," << ptr_entities[i].xmin_el << "," << ptr_entities[i].ymin_el << "," << ptr_entities[i].zmin_el << "," << ptr_entities[i].xmax_el << "," << ptr_entities[i].ymax_el << "," << ptr_entities[i].zmax_el << "," << ptr_entities[i].physical_tag_el << "," << ptr_entities[i].num_boundary_nodes_el << "," << ptr_entities[i].boundary_node_1_el << "," << ptr_entities[i].boundary_node_2_el << "\n"; // "," << ptr_entities[i].boundary_node_3_el << "\n"; // << "," << ptr_entities[i].boundary_node_3_el << "\n";
							/* TODO ....
							ptr_entities = new Entities[ptr_entities[i].num_boundary_nodes_el];
							for (int m = 0; m < ptr_entities[i].num_boundary_nodes_el; m++)
							{
								file >> ptr_entities[m].boundary_node_x_el;
							}
							*/


							//file >> entities.line_tag_el >> entities.xmin_el >> entities.ymin_el >> entities.zmin_el >> entities.xmax_el >> entities.ymax_el >> entities.zmax_el >> entities.physical_tag_el >> entities.num_boundary_nodes_el >> entities.boundary_node_1_el >> entities.boundary_node_2_el;
							//std::cout << entities.line_tag_el << "," << entities.xmin_el << "," << entities.ymin_el << "," << entities.zmin_el << "," << entities.xmax_el << "," << entities.ymax_el << "," << entities.zmax_el << "," << entities.physical_tag_el << "," << entities.num_boundary_nodes_el << "," << entities.boundary_node_1_el << "," << entities.boundary_node_2_el << "\n";
						}

						for (int i = 0; i < entities.numberOfSurface; i++)
						{
							//ptr_curve_tag = new int[entities.numberOfLine];
							if (entities.numberOfLine == 4)
							{
								file >> ptr_entities[i].surface_tag_es >> ptr_entities[i].xmin_es >> ptr_entities[i].ymin_es >> ptr_entities[i].zmin_es >> ptr_entities[i].xmax_es >> ptr_entities[i].ymax_es >> ptr_entities[i].zmax_es >> ptr_entities[i].physical_tag_es >> ptr_entities[i].num_boundary_curves_es >> ptr_entities[i].boundary_node_1_es >> ptr_entities[i].boundary_node_2_es >> ptr_entities[i].boundary_node_3_es >> ptr_entities[i].boundary_node_4_es;;
								//std::cout << ptr_entities[i].surface_tag_es << "," << ptr_entities[i].xmin_es << "," << ptr_entities[i].ymin_es << "," << ptr_entities[i].zmin_es << "," << ptr_entities[i].xmax_es << "," << ptr_entities[i].ymax_es << "," << ptr_entities[i].zmax_es << "," << ptr_entities[i].physical_tag_es << "," << ptr_entities[i].num_boundary_curves_es << "," << ptr_entities[i].boundary_node_1_es << "," << ptr_entities[i].boundary_node_2_es << "," << ptr_entities[i].boundary_node_3_es << "," << ptr_entities[i].boundary_node_4_es << "\n";

								//file >> entities.surface_tag_es >> entities.xmin_es >> entities.ymin_es >> entities.zmin_es >> entities.xmax_es >> entities.ymax_es >> entities.zmax_es >> entities.physical_tag_es >> entities.num_boundary_curves_es >> ptr_curve_tag[0] >> ptr_curve_tag[1] >> ptr_curve_tag[2] >> ptr_curve_tag[3];
								//std::cout << entities.surface_tag_es << "," << entities.xmin_es << "," << entities.ymin_es << "," << entities.zmin_es << "," << entities.xmax_es << "," << entities.ymax_es << "," << entities.zmax_es << "," << entities.physical_tag_es << "," << entities.num_boundary_curves_es << "," << ptr_curve_tag[0] << "," << ptr_curve_tag[1] << "," << ptr_curve_tag[2] << "," << ptr_curve_tag[3] << "\n";
							}

						}

						for (int i = 0; i < entities.numberOfVolume; i++)
						{
							// A FINIR...
						}

						break;  // Quitter la boucle après avoir trouvé et lu les valeurs

					}
				}

				//std::cout << "\n";





				//std::system("PAUSE");
				//std::cout << "\ttraitement NODES " << "\n";
				//$Nodes: liste des nœuds avec leurs coordonnées.


				GMSH::Nodes nodes{};
				//GMSH::Nodes* ptr_nodes{};
				std::vector<GMSH::Nodes*>v_nodes{};
				std::vector<double>vx{};


				while (std::getline(file, line))
				{
					if (line == "$Nodes")
					{
						// Lire les trois valeurs après "$MeshFormat"
						file >> nodes.num_entity_blocks >> nodes.num_total_nodes >> nodes.min_node_tag >> nodes.max_node_tag;
						///////////////std::cout << nodes.num_entity_blocks << "," << nodes.num_total_nodes << "," << nodes.min_node_tag << "," << nodes.max_node_tag << "\n";

						/*
						num_entity_blocks : Le nombre de blocs d'entités. Chaque bloc correspond à un groupe de nœuds appartenant à une entité géométrique spécifique (par exemple, une ligne, une surface).
						num_total_nodes : Le nombre total de nœuds dans le fichier.
						min_node_tag : Le plus petit identifiant de nœud.
						max_node_tag : Le plus grand identifiant de nœud.
						*/
						ptr_nodes = new Nodes[nodes.num_total_nodes];
						Nodes* P;

						int l = 1;
						int cpt = { 0 };
						for (int k = 0; k < nodes.num_total_nodes; k++)
						{



							/*
							double ratio = static_cast<double>(k) / nodes.num_total_nodes;
							double percentage = 100 * ratio;
							std::cout << k << "," << nodes.num_total_nodes << "," << percentage << "," << std::ceil(percentage) << "\n";

							//std::cout << std::ceil(percentage) << "\n";
							std::this_thread::sleep_for(std::chrono::milliseconds(50));
							std::system("cls");
							*/

							

							file >> ptr_nodes[k].entity_dim >> ptr_nodes[k].entity_tag >> ptr_nodes[k].parametric >> ptr_nodes[k].num_nodes_in_block;
							////////////std::cout << ptr_nodes[k].entity_dim << "," << ptr_nodes[k].entity_tag << "," << ptr_nodes[k].parametric << "," << ptr_nodes[k].num_nodes_in_block << "\n";
							//file >> nodes.entity_dim >> nodes.entity_tag >> nodes.parametric >> nodes.num_nodes_in_block;
							//std::cout << nodes.entity_dim << "," << nodes.entity_tag << "," << nodes.parametric << "," << nodes.num_nodes_in_block << "\n";

							/*
							entity_dim : La dimension de l'entité (0 pour un point, 1 pour une ligne, 2 pour une surface, 3 pour un volume).
							entity_tag : L'identifiant de l'entité géométrique à laquelle appartiennent les nœuds dans ce bloc.
							parametric : Indique si les coordonnées des nœuds sont paramétriques (1 si elles sont paramétriques, 0 sinon).
							num_nodes_in_block : Le nombre de nœuds dans ce bloc.
							*/


							for (int i = 0; i < ptr_nodes[k].num_entity_blocks; i++)
							{
								//file >> ptr_nodes[k].node_tag;
								//std::cout << ptr_nodes[k].node_tag << "\n";
								for (int j = 0; j < ptr_nodes[k].num_nodes_in_block; j++)
								{

								}
								cpt++;
							}


							if (ptr_nodes[k].entity_dim == 0)//pour un point
							{

								file >> ptr_nodes[k].node_tag;
								////////////std::cout << ptr_nodes[k].node_tag << "\n";

								for (int i = 0; i < ptr_nodes[k].num_nodes_in_block; i++)
								{
									file >> ptr_nodes[k].x >> ptr_nodes[k].y >> ptr_nodes[k].z;
									////////////std::cout << ptr_nodes[k].x << "," << ptr_nodes[k].y << "," << ptr_nodes[k].z << "\n";
								}

								P = new Nodes(l, ptr_nodes[k].x, ptr_nodes[k].y, ptr_nodes[k].z);
								list_ptr_nodes.push_back(P);


								l++;
							}
							else if (ptr_nodes[k].entity_dim == 1)//pour une ligne
							{
								/*
								file >> ptr_nodes[k].node_tag;
								std::cout << ptr_nodes[k].node_tag << "\n";
								*/
								for (int i = 0; i < ptr_nodes[k].num_nodes_in_block; i++)
								{
									file >> ptr_nodes[k].node_tag;
									///////////std::cout << ptr_nodes[k].node_tag << "\n";
								}
								/*
								P = new Nodes(l, ptr_nodes[k].x, ptr_nodes[k].y, ptr_nodes[k].z);
								list_ptr_nodes.push_back(P);
								//list_ptr_nodes.push_back(&ptr_nodes[k]);
								l++;*/

								for (int i = 0; i < ptr_nodes[k].num_nodes_in_block; i++)
								{
									file >> ptr_nodes[k].x >> ptr_nodes[k].y >> ptr_nodes[k].z;
									//////////////std::cout << ptr_nodes[k].x << "," << ptr_nodes[k].y << "," << ptr_nodes[k].z << "\n";
									P = new Nodes(l, ptr_nodes[k].x, ptr_nodes[k].y, ptr_nodes[k].z);
									list_ptr_nodes.push_back(P);


									l++;
								}
							}
							else if (ptr_nodes[k].entity_dim == 2)//pour une surface
							{
								for (int i = 0; i < ptr_nodes[k].num_nodes_in_block; i++)
								{
									file >> ptr_nodes[k].node_tag;
									//////////////std::cout << ptr_nodes[k].node_tag << "\n";
								}

								for (int i = 0; i < ptr_nodes[k].num_nodes_in_block; i++)
								{
									ptr_nodes[k].node_tag = l;
									file >> ptr_nodes[k].x >> ptr_nodes[k].y >> ptr_nodes[k].z;
									
									//file >> ptr_nodes[k].node_tag >> ptr_nodes[k].x >> ptr_nodes[k].y >> ptr_nodes[k].z;
									

									///////////////std::cout << ptr_nodes[k].node_tag <<  ","  << ptr_nodes[k].x << "," << ptr_nodes[k].y << "," << ptr_nodes[k].z << "\n";
									P = new Nodes(l, ptr_nodes[k].x, ptr_nodes[k].y, ptr_nodes[k].z);
									list_ptr_nodes.push_back(P);


									//progressBar
									/*
									displayProgressBar(i, ptr_nodes[k].num_nodes_in_block - 1);
									std::this_thread::sleep_for(std::chrono::milliseconds(5));  // Pause pour la simulation
								    */
									l++;
								}


								/*
								Nodes* P = new Nodes(l, ptr_nodes[k].x, ptr_nodes[k].y, ptr_nodes[k].z);
								v_nodes.push_back(P);
								list_ptr_nodes.push_back(&ptr_nodes[k]);
								l++;*/
							}
							else if (ptr_nodes[k].entity_dim == 3)//pour un volume
							{
								std::cout << "\tERROR! not yet implemented." << "\n";
							}

							//this->ptr_nodes = ptr_nodes;

							if (cpt < nodes.num_entity_blocks - 1)
							{
								//std::cout << "\t\t\t\t\tcpt= " << cpt << "\n";
								cpt++;
							}
							else
							{
								break;
							}//end if



						}//end for k



						break;  // Quitter la boucle après avoir trouvé et lu les valeurs



					}//$EndNodes

				}//end while



				//clearVector(v_nodes);

				/*
				std::cout << "\n";
				std::system("PAUSE");
				*/

				/*

				std::cout << "\t=======================\n";
				// Iterate through the array of Nodes* using a for loop
				for (int i = 0; i < nodes.num_total_nodes; ++i) {
					ptr_nodes[i].node_tag = i;
					std::cout << "Node ID: " << ptr_nodes[i].node_tag
						<< ", Coordinates: (" << ptr_nodes[i].x
						<< ", " << ptr_nodes[i].y
						<< ", " << ptr_nodes[i].z << ")\n";
					Nodes* P = new Nodes(ptr_nodes[i].node_tag, ptr_nodes[i].x, ptr_nodes[i].y, ptr_nodes[i].z);
					list_ptr_nodes.push_back(P);
				}
				std::cout << "\t=======================\n";


				*/

				//std::cout << "\n";
				//std::system("PAUSE");






				this->list_ptr_nodes = list_ptr_nodes;
				/*
				for (auto p : list_ptr_nodes)
				{
					std::cout << *p << "\n";
				}
				*/



				/*
				std::system("PAUSE");
				std::cout << "\ttraitement ELEMENTS " << "\n";
				*/
				GMSH::Elements elements{};
				//GMSH::Elements* ptr_elements{};
				//GMSH::Nodes* vertices{};

				//Nodes* current_node;


				while (std::getline(file, line))
				{
					if (line == "$Elements")
					{
						//num_entity_blocks num_total_elements min_element_tag max_element_tag
						file >> elements.num_entity_blocks >> elements.num_total_elements >> elements.min_element_tag >> elements.max_element_tag;
						////////////std::cout << elements.num_entity_blocks << "," << elements.num_total_elements << "," << elements.min_element_tag << "," << elements.max_element_tag << "\n";

						this->num_total_elements = elements.num_total_elements;
						ptr_elements = new GMSH::Elements[elements.num_total_elements];
						//std::system("PAUSE");


						/*
						* num_entity_blocks : Nombre de blocs d'entités.
						* num_total_elements : Nombre total d'éléments dans le fichier.
						* min_element_tag : Le plus petit identifiant d'élément.
						* max_element_tag : Le plus grand identifiant d'élément.
						*
						* Types d'éléments

						Le champ element_type correspond à un type d'élément géométrique, comme indiqué dans la documentation de GMSH. Quelques exemples courants :

							0 : Point (1 neud)
							1 : Ligne (2 nœuds).
							2 : Triangle (3 nœuds).
							3 : Quadrilatère (4 nœuds).
							4 : Tétraèdre (4 nœuds).
							5 : Hexaèdre (8 nœuds).
						*/
						/////////////////////////////////////////////////////////////////////////////Nodes* ptr_nodes;
						Nodes* Gk;
						double area_tot = { 0.0 };
						Nodes CG;

						

						int q = { 0 }, cpt = { 0 };
						for (int l = 0; l < elements.num_total_elements; l++)
						{
							file >> ptr_elements[l].entity_dim >> ptr_elements[l].entity_tag >> ptr_elements[l].element_type >> ptr_elements[l].num_elements_in_block;

							/*
							displayProgressBar(l, elements.num_total_elements);
							std::this_thread::sleep_for(std::chrono::milliseconds(50));  // Pause pour la simulation
							std::cout << std::endl;
							*/
							/////////////std::cout << ptr_elements[l].entity_dim << ", " << ptr_elements[l].entity_tag << "," << ptr_elements[l].element_type << "," << ptr_elements[l].num_elements_in_block << "\n";
							//file >> elements.entity_dim >> elements.entity_tag >> elements.element_type >> elements.num_elements_in_block;
							//std::cout << elements.entity_dim << ", " << elements.entity_tag << "," << elements.element_type << "," << elements.num_elements_in_block << "\n";



							for (int i = 0; i < ptr_elements[l].num_elements_in_block; i++)
							{

								//std::cout << "\tnbre elements (points, segments ou triangles): " << ptr_elements[l].num_elements_in_block << "\n";
								if (ptr_elements[l].entity_dim == 0) // 0 pour un pts
								{
									//std::system("PAUSE");
									//std::cout << "\ttraitement POINTS " << "\n";
									//file >> elements.element_tag >> elements.node_tag_1;
									//std::cout << elements.element_tag << "," << elements.node_tag_1 << "\n";
									file >> ptr_elements[l].element_tag >> ptr_elements[l].node_tag_1;
									/////////std::cout << "\t-----" << ptr_elements[l].element_tag << "," << ptr_elements[l].node_tag_1 << "\n";
								}//end if
								else if (ptr_elements[l].entity_dim == 1) // 1 pour une kigne
								{
									//std::system("PAUSE");
									//std::cout << "\ttraitement SEGMENTS " << "\n";
									//file >> elements.element_tag >> elements.node_tag_1 >> elements.node_tag_2;
									//std::cout << elements.element_tag << "," << elements.node_tag_1 << "," << elements.node_tag_2 << "\n";
									file >> ptr_elements[l].element_tag >> ptr_elements[l].node_tag_1 >> ptr_elements[l].node_tag_2;
									///////////std::cout << "\t-----" << ptr_elements[l].element_tag << "," << ptr_elements[l].node_tag_1 << "," << ptr_elements[l].node_tag_2 << "\n";
								}//end if								
								else if (ptr_elements[l].entity_dim == 2)// pour une surface
								{
									// TRIANGLES
									//std::system("PAUSE");	
									//std::cout << "\ttraitement TRIANGLES " << "\n";
									file >> ptr_elements[l].element_tag >> ptr_elements[l].node_tag_1 >> ptr_elements[l].node_tag_2 >> ptr_elements[l].node_tag_3;
									////////////std::cout << "\t-----" << ptr_elements[l].element_tag << "," << ptr_elements[l].node_tag_1 << "," << ptr_elements[l].node_tag_2 << "," << ptr_elements[l].node_tag_3 << "\n";

									
									ptr_nodes = new Nodes[list_ptr_nodes.size()];
									for (int q = 0; q < list_ptr_nodes.size(); q++)
									{
										ptr_nodes->node_tag = list_ptr_nodes[q]->node_tag;
										ptr_nodes[q].x = list_ptr_nodes[q]->x;
										ptr_nodes[q].y = list_ptr_nodes[q]->y;
										ptr_nodes[q].z = list_ptr_nodes[q]->z;
									}
									
									
									ptr_elements[l].set(ptr_nodes, ptr_elements[l].node_tag_1 - 1, ptr_elements[l].node_tag_2 - 1, ptr_elements[l].node_tag_3 - 1, 0);
									
									
									////list_area.push_back( ptr_elements[l].area );
									//std::cout << "\tsurface (" << ptr_elements[l].element_tag << ") = " << ptr_elements[l].area  << "\n";
									area_tot += ptr_elements[l].area;
									this->area = area_tot;
									

								    list_step.push_back( ptr_elements[l].step_k );


									///list_diam_k.push_back(ptr_elements[l].diam_k);

									Gk = new Nodes(ptr_elements[l].CG);
									Gk->node_tag = ptr_elements[l].element_tag;									
									list_ptr_nodes_cg.push_back( Gk );

									
									Elements elem(ptr_elements[l].element_tag, ptr_elements[l].CG, ptr_elements[l].area, ptr_elements[l].perimeter);
									list_elements.push_back(&elem);

									


								}
								else // 3 pour un volume
								{
									std::cout << "\tERROR! not yet implemented." << "\n";
								}

							}//end for


							//std::system("PAUSE");
							if (cpt < elements.num_entity_blocks - 1)
							{
								//std::cout << "\t\t\t\t\tcpt= " << cpt << "\n";
								cpt++;
								//std::system("PAUSE");
							}
							/*
							else
							{
								std::cout << "\n\t--------------------------------------------------\n";
								std::cout << "\tsurface couvertes par les triangles:" << area_tot << "\n";
								//std::cout << "\tsurface couvertes par les triangles:" << PI << "\n";
								std::cout << "\t--------------------------------------------------\n";
								std::system("PAUSE");
								break;
							}*/

						}//end for

						//delete[]Gk;
						//delete  ptr_nodes; // Free dynamically allocated array
						break;

					}//end $Elements



				}//end while
				std::cout << "\n";







			}//end else

			/*
			std::system("PAUSE");
			std::cout << "\n\n\n";
			std::cout << "\tFIN lecture FICHIER ! " << "\n";
			*/
			file.close();


		}//Mesh



		Mesh& operator=(const Mesh& msh) = default;
		auto  operator<=>(const Mesh& msh)const = default;

		//~Mesh();





	};

	/*
	Mesh::~Mesh()
	{
		delete[]ptr_elements;
	}*/


}//end namespace GMSH






#endif // !__MESH_H
