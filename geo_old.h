#ifndef __GEO__OLD_H
#define __GEO_OLD_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <regex>
#include <map>
#include <vector>
#include <stdexcept>
#include <filesystem>

namespace GMSH
{
    typedef double dble;

	class Geo
	{
	public:
        int surfaceNumber{};
		int N;
        dble e{}, T{};

		Geo():N(-1) {};
        Geo(int surfNumb, int NN, dble ee, dble TT) :surfaceNumber{ surfNumb }, N(NN), e{ ee }, T{ TT } {};
		Geo(const Geo& geo) :surfaceNumber( geo.surfaceNumber ), N( geo.N ),e( geo.e ),T( geo.T ){};
		Geo& operator=(const Geo& geo) = default;

		friend std::ostream& operator<<(std::ostream& f, const Geo& geo);
		friend std::istream& operator>>(std::istream& f, Geo& geo);

        bool parseFile(const std::string& filename);
        void parsefile_all(std::string input);
        void parsefile_geo(const std::string& filename);


        void displayConstants() const;    
        int get_N() const { return N; }
        dble get_e() const { return e; }
        dble get_T() const { return T; }
	};


    void Geo::parsefile_geo(const std::string& filename)
    {
        std::ifstream file(filename);
        if (!file.is_open())
        {
            std::cerr << "\tFile not found !" << "\n";
            std::system("pause");
            std::exit(1);
        }
        // Charger tout le contenu dans une chaîne
        std::string content((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
        Geo::parsefile_all(content);

        file.close();
    }

    
	std::ostream& operator<<(std::ostream& f, const Geo& geo)
	{
        f << "\n\n\tsurface number: " << geo.surfaceNumber << "\n";
		f << "\tN: " << geo.get_N() << ", emissivity: " << geo.get_e() << ", Temperature (K): " << geo.get_T() << "\n\n";
		return f;
	}


	std::istream& operator>>(std::istream& f, Geo& geo)
	{
		f >> geo.N >> geo.e >> geo.T;
		return f;
	}


    void Geo::displayConstants() const {//display N
        if (N != -1) {
            std::cout << "\t\t\tN = " << N << std::endl; //<< ", e = " << e << ", T = " << T 
        }
        else {
            std::cout << "N not found in file." << std::endl;
        }
    }

    bool Geo::parseFile(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Error opening file: " << filename << std::endl;
            return false;
        }
        std::string line;
        std::smatch match;
        std::regex pattern(R"(N\s*=\s*\{(\d+))");  // Regex to match "N = {number"
        //std::regex pattern(R"(e\s*=\s*\{([\d.]+))"); 
        while (std::getline(file, line)) {
            if (std::regex_search(line, match, pattern)) {
                //e = std::stod(match[1].str());  // Convert matched value to integer
                N = std::stoi(match[1].str());
                break;
            }
        }
        file.close();
        return N != -1;  // return true if N was found, false otherwise
    }




    void Geo::parsefile_all(std::string input)
    {
        std::regex pattern(R"(([A-Za-z]+)\s*=\s*\{([\d.]+))");
        std::smatch match;

        // Itere pour trouver toutes les correspondances
        std::string::const_iterator searchStart(input.cbegin());
        while (std::regex_search(searchStart, input.cend(), match, pattern))
        {   
            if (!match[1].str().compare("N")) this->N = std::stoi(match[2].str());
            if (!match[1].str().compare("e")) this->e = std::stod(match[2].str());
            if (!match[1].str().compare("T")) this->T = std::stod(match[2].str());  
            //std::cout << match[1] << " = " << match[2] << std::endl;
            
            searchStart = match.suffix().first; //passe à la prochain correspondance
        }
    }










}//end namespace GMSH

#endif // !__GEO_OLD_H
