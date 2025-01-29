#ifndef __CALCULSFACTEURSDEVUE__H
#define __CALCULSFACTEURSDEVUE__H


#include <complex>
#include <vector>
#include <functional>
#include <map>
#include <ctime>
#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include <windows.h>
#include <algorithm> 
#include <thread>
#include <Windows.h>
#include <assert.h>




#include "../ToolLib/nodes.h"
#include "../ToolLib/elements.h"
#include "../ToolLib/mesh.h"
#include "../ToolLib/geo.h"
#include "../ToolLib/geo_old.h"
#include "../ToolLib/viewfactor.h"
#include "../ToolLib/array.h"
#include "../ToolLib/entities.h"
#include "../ToolLib/MatrixRotation.h"
#include "../ToolLib/meshFormat.h"
#include "../ToolLib/quantityOfHeatByRadiation.h"

#include "../ToolLib/array_class.h"
#include "array.h"





void showProgressBar_vf(int current, int total, const std::string& message) {
	int progress = (100 * current) / total;
	int barWidth = 50;
	int pos = (barWidth * current) / total;

	std::cout << "\r" << message << " [";
	for (int i = 0; i < barWidth; ++i) {
		if (i < pos) std::cout << "=";
		else if (i == pos) std::cout << ">";
		else std::cout << " ";
	}
	std::cout << "] " << progress << "% " << std::flush;
}

/*
void showProgressBar_vf(int progress, int total, std::string message1) {
	int barWidth = 50;
	float percentage = static_cast<float>(progress) / total;
	int pos = (int)(barWidth * percentage);
	std::cout << message1 << "\n\n";
	std::cout << "\t\t\t[";
	for (int i = 0; i < barWidth; ++i) {
		if (i < pos)
			std::cout << "=";
		else if (i == pos)
			std::cout << ">";
		else
			std::cout << " ";
	}
	std::cout << "] " << int(percentage * 100) << " %\r";
	std::cout.flush();
}*/


namespace GMSH
{



	template<typename T>
	class CalculsFacteursDeVue
	{
	private:
		std::vector<std::string> listeFichiersMaillage;
	public:
		CalculsFacteursDeVue() = default;//default constructor
		CalculsFacteursDeVue(const std::vector<std::string>& liste) :listeFichiersMaillage(liste)
		{
			solve();
		}
		CalculsFacteursDeVue(const CalculsFacteursDeVue& cvf) = default;
		CalculsFacteursDeVue& operator=(const CalculsFacteursDeVue& cvf) = default;//copy by assignment


		auto operator<=>(const CalculsFacteursDeVue& cfv)const = default;

    	std::vector<std::string> getListeFichierMaillage();




		/*
		linearAlgebra::Array<double> getArrayOfVF()
		{
			return this->MAT_VF;
		}


		void setArrayOfVF(linearAlgebra::Array<double>& mat)
		{
			this->MAT_VF = mat;
		}
		*/





		void solve()
		{
			if (listeFichiersMaillage.empty())
			{
				std::cerr << "Erreur : listeFichiersMaillage est vide !" << std::endl;
				return;
			}

			int numberOfFiles = listeFichiersMaillage.size() * listeFichiersMaillage.size() - listeFichiersMaillage.size();
			int nbre_ligne = listeFichiersMaillage.size();

			// Initialisation correcte de la matrice
			linearAlgebra::Array<T> MAT_VF(nbre_ligne, nbre_ligne);

			

			
			GMSH::View_Factor f12{};
			GMSH::Geo geo_i{}, geo_j{};
			double val_f12 = {};
			int i, j, cpt_mesh_number = { 0 };

			std::cout << "\n\tMATRIX of VIEW FACTORS in PROGRESS..." << "\n\n\n";
			
			


			int totalIterations = nbre_ligne * (nbre_ligne - 1);
			int currentIteration = 0;
			

			i = { 0 };
			for (auto m1 : listeFichiersMaillage)
			{
				std::string m1_msh = m1 + ".msh";
				std::string m1_geo = m1 + ".geo";

				geo_i.parsefile_geo(m1_geo);
				geo_i.surfaceNumber = i;

				GMSH::Mesh mesh1(m1_msh);
				f12.area1 = mesh1.area;
				j = { 0 };
				for (auto m2 : listeFichiersMaillage)
				{
					std::string m2_msh = m2 + ".msh";
					std::string m2_geo = m2 + ".geo";

					geo_j.parsefile_geo(m2_geo);
					geo_j.surfaceNumber = j;

					GMSH::Mesh mesh2(m2_msh);
					f12.area2 = mesh2.area;
					if (j != i)
					{
						val_f12 = f12.calculs2(mesh1.list_ptr_nodes_cg, mesh2.list_ptr_nodes_cg);
						f12.value = val_f12 * (mesh1.area / mesh1.list_elements.size()) * (mesh2.area / mesh2.list_elements.size());
						f12.emitterNumber = { i }; f12.receiverNumber = { j };
						f12.numberOfElements1 = mesh1.num_total_elements;
						f12.numberOfElements2 = mesh2.num_total_elements;
						f12.area1 = mesh1.area;
						f12.area2 = mesh2.area;

						if (i >= nbre_ligne || j >= nbre_ligne)
						{
							throw std::out_of_range("Indices i ou j hors limite !");
						}
						MAT_VF.setValue(i, j, f12.value);



						// Mise à jour de la barre de progression ici
						currentIteration++;
						showProgressBar_vf(currentIteration, totalIterations, "READING MESH FILES in PROGRESS...");
						std::this_thread::sleep_for(std::chrono::milliseconds(15));  // Pause pour la simulation
						/*
						showProgressBar_vf(cpt_mesh_number, numberOfFiles, "READING MESH FILES in PROGRESS...");
						std::this_thread::sleep_for(std::chrono::milliseconds(100));  // Pause pour la simulation
						*/		


						cpt_mesh_number++;
					}//end if j

					// Vérifications des bornes avec des assertions et exceptions
					if (j < 0 || j >= MAT_VF.columnSize())
					{
						throw std::out_of_range("Indices hors limites dans Array::operator()");
					}
					j++;

				}//end for m2
				
				 // Vérifications des bornes avec des assertions et exceptions
				if (i < 0 || i >= MAT_VF.rowSize())
				{
					throw std::out_of_range("Indices hors limites dans Array::operator()");
				}
				i++;


				std::cout << "\nTraitement termine !\n";


			}//end for m1 */


			
			for (int i = 0; i < MAT_VF.rowSize(); i++)
			{
				double sum_j = { 0.0 };
				for (int j = 0; j < MAT_VF.columnSize(); j++)
				{
					if (i != j)
					{
						sum_j += MAT_VF.getValue(i,j);
					}
				}
				double val = std::abs(1.0 - sum_j);
				MAT_VF.setValue(i, i, val); 
			}
			


			/* RENORMALISATION des FACTEURS de VUE (pour garantir la conservation de l'énergie - CF CGHATGPT) */
			/*
			for (int i = 0; i < MAT_VF.rowSize(); i++)
			{
				double value = { 0.0 };
				for (int j = 0; j < MAT_VF.columnSize(); j++)
				{
					value = MAT_VF.getValue(i, j) / MAT_VF.getSumElementsOfRow(i);
					MAT_VF.setValue(i, j, value);
				}
			}*/



			std::cout << MAT_VF << "\n";



		}//end solve

	};//end classe CalculsFacteursDeVue


}//end namespace GMSH


#endif // !__CALCULSFACTEURSDEVUE__H
