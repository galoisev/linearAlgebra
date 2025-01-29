#ifndef __QUANTITYOFHEATBYRADIATION__H
#define __QUANTITYOFHEATBYRADIATION__H



#include <iostream>
#include <assert.h>
#include <string>
#include <math.h>
#include <map>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <assert.h>
#include <stdio.h>
#include <format>



typedef double dble;




namespace GMSH
{



	class QuantityOfHeatByRadiation
	{
	public:
		dble Fij{}, Qij{};
		dble emissivity_i{}, T_i{}, Area_i{};
		dble emissivity_j{}, T_j{}, Area_j{};
		dble delta_time{ 3600.0 };

		QuantityOfHeatByRadiation(/**/) = default;//default constructor
		QuantityOfHeatByRadiation(dble& ei, dble& Ti, dble& areai, dble& Tj, dble& dt);//constructor
		QuantityOfHeatByRadiation(const QuantityOfHeatByRadiation& Q) = default;//constructor by copy
		QuantityOfHeatByRadiation& operator=(const QuantityOfHeatByRadiation& Q) = default;//constructor by assigment


		friend std::ostream& operator<<(std::ostream& f, const QuantityOfHeatByRadiation& Q);
		friend std::istream& operator>>(std::istream& f, QuantityOfHeatByRadiation& Q);

		dble radiation(dble& Fij);
		void printFile(const std::string& filename);

		void set_area_i(dble& area);//emitter
		dble get_area_i();//emitter

		void set_area_j(dble& area);//receiver
		dble get_area_j();//receiver

	};

	void QuantityOfHeatByRadiation::set_area_i(dble& area)
	{
		this->Area_i = area;
	}
	dble QuantityOfHeatByRadiation::get_area_i()
	{
		return this->Area_i;
	}
	void QuantityOfHeatByRadiation::set_area_j(dble& area)
	{
		this->Area_j = area;
	}
	dble QuantityOfHeatByRadiation::get_area_j()
	{
		return this->Area_j;
	}


	void QuantityOfHeatByRadiation::printFile(const std::string& filename)
	{
		//Sauvegarde du fichier
		std::ofstream fic(filename);
		fic << "\n";
		fic << "\t\n";
		fic << "\t\t\t*** SUMMARY *** \n";
		fic << "\t\n";
		fic << "\tArea_i[M2] = " << this->Area_i << ", emissivity_i = " << this->emissivity_i << ", T_i[K] = " << this->T_i << "\n";
		fic << "\tArea_j[M2] = " << this->Area_j << ", emissivity_j = " << this->emissivity_j << ", T_j[K] = " << this->T_j << "\n";
		fic << "\n\t\t\tRADIATION HEAT QUANTITY [J]: " << this->Qij << "\n";
		fic << "\n\n\n";
		fic.close();
	}


	std::ostream& operator<<(std::ostream& f, const QuantityOfHeatByRadiation& Q)
	{
		f << "\n";
		f << "\t\n";
		f << "\t\t\t*** SUMMARY *** \n";
		f << "\t\n";
		f << "\tArea_i[M2] = " << Q.Area_i <<  ", emissivity_i = " << Q.emissivity_i << ", T_i[K] = " << Q.T_i << "\n";
		f << "\tArea_j[M2] = " << Q.Area_j <<  ", emissivity_j = " << Q.emissivity_j << ", T_j[K] = " << Q.T_j << "\n";
		f << "\n\t\t\tRADIATION HEAT QUANTITY [J]: " << Q.Qij << "\n";
		f << "\n\n\n";
		return f;
	}

	std::istream& operator>>(std::istream& f, QuantityOfHeatByRadiation& Q)
	{
		f >> Q.Qij >> Q.Area_i >> Q.emissivity_i >> Q.T_i >> Q.Area_j >> Q.emissivity_j >> Q.T_j;
		return f;
	}


	dble QuantityOfHeatByRadiation::radiation(dble& Fij)
	{
		dble sigma = 5.67e-8; // constante de Stefan-Boltzmann en W/(m**2 * K**4)
		dble P_rad = (sigma * this->Area_i * (pow(this->T_i, 4.0) - pow(this->T_j, 4.0)) * Fij) / (1.0 / this->emissivity_i + 1.0 / this->emissivity_j - 1.0); //puissance radiative en Watts
		dble Q_ij = P_rad * this->delta_time; //quantite de chaleur echangee par radiation entre les surfaces i et j (en Joules)
		this->Qij = Q_ij;
		return this->Qij;
	}


	QuantityOfHeatByRadiation::QuantityOfHeatByRadiation(dble& ei, dble& Ti, dble& areai, dble& Tj, dble& dt) :emissivity_i(ei), T_i(Ti), Area_i(areai), T_j(Tj), delta_time(dt)
	{
		dble F_ij = this->Fij;
		dble Q_ij = radiation(F_ij);
	}


}//end GMSH namespace




#endif // ! __QUANTITYOFHEATBYRADIATION__H