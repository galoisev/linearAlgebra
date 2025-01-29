#ifndef __GEO_H
#define __GEO_H


#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <regex>
#include <map>
#include <vector>
#include <cassert>




namespace GMSH
{

	class GeoConstants
	{
	public:
		double value;         // La valeur principale
		std::string metadata; // Le "Name"

		GeoConstants() {};
		GeoConstants(double dble_val, std::string str_metadata) :value(dble_val), metadata(str_metadata) {};
		GeoConstants(const GeoConstants& g) :value(g.value), metadata(g.metadata) {};
		GeoConstants& operator=(const GeoConstants& g) = default;

		friend std::ostream& operator<<(std::ostream& f, const GeoConstants& g);
		friend std::istream& operator>>(std::istream& f, GeoConstants& g);


		void affichage();
		bool parsefile(const std::string& filename);


	};



	std::ostream& operator<<(std::ostream& f, const GeoConstants& g)
	{
		f << g.value << g.metadata << "\n";
		return f;
	}

	std::istream& operator>>(std::istream& f, GeoConstants& g)
	{
		f >> g.value >> g.metadata;
		return f;
	}



	bool GMSH::GeoConstants::parsefile(const std::string& filename)
	{
		std::ifstream file(filename);
		if (!file.is_open())
		{
			std::cerr << "\tErreur: Impossible d'ouvrir le fichier !" << "\n";
			return 1;
		}
		bool isDefineConstantSection = false;
		std::map<std::string, GeoConstants> constants;
		std::string line;
		while (std::getline(file, line))
		{
			if (line.find("DefineConstant[") != std::string::npos)
			{
				isDefineConstantSection = true;
				continue;
			}
			if (isDefineConstantSection && line.find("]") != std::string::npos)
			{
				isDefineConstantSection = false;
				break;
			}
			//parsage des lignes

			if (isDefineConstantSection)
			{
				std::string name, metadata;
				dble value{};
				//trouve la position de la valeur principale et du Name
				auto eqPos = line.find('=');
				auto namePos = line.find("Name");
				if (eqPos != std::string::npos && namePos != std::string::npos)
				{
					/*
					//Nom de la constante
					name = line.substr(0, eqPos);
					{
						(std::remove_if(name.begin(), name.end(), isspace), name.end());//Supprime les espaces
						//Valeur de la constante
						auto startValue = line.find('{', eqPos) + 1;
						auto endValue = line.find(',', startValue);
						value = std::stod(line.substr(startValue, endValue - startValue));

						//Metadonnees (Name)
						auto startName = line.find('"', namePos) + 1;
						auto endName = line.find('"', startName);
						metadata = line.substr(startName, endName - startName);

						//stocke les données dans la map
						constants[name] = { value,metadata };
					}
					*/
				}
			}
		}
	}

	

	inline void GMSH::GeoConstants::affichage()
	{
		
	}

}//end namespace GMSH

#endif // !__GEO_H


