#ifndef __MESHFORMAT_H
#define __MESHFORMAT_H

#include <iostream>
#include <istream>
#include <fstream>
#include <string>
#include <vector>
#include <cassert>


namespace GMSH
{



	class MeshFormat
	{
	public:
		double version_number;//est un reel valant actuellement 2.1 ou 2.2 (dernieres versions du format msh).
		int file_type;//vaut 0 pour un chier texte.
		int data_size;//precise le nombre de decimales significatives d'un flottant double precision(= sizeof(double)).

		MeshFormat() {};
		MeshFormat(double vnumber, int ftype, int dsize) :version_number(vnumber), file_type(ftype), data_size(dsize) {};
		MeshFormat(const MeshFormat& mf) :version_number(mf.version_number), file_type(mf.file_type), data_size(mf.data_size) {};
		MeshFormat& operator=(const MeshFormat& mf) = default;

		friend std::ostream& operator<<(std::ostream& f, const MeshFormat& meshformat);
		friend std::istream& operator>>(std::istream& f, MeshFormat& meshformat);

	};



	std::ostream& operator<<(std::ostream& f, const MeshFormat& meshformat)
	{
		f << meshformat.version_number << " " << meshformat.file_type << " " << meshformat.data_size << "\n";
		return f;
	}

	std::istream& operator>>(std::istream& f, MeshFormat& meshformat)
	{
		f >> meshformat.version_number >> meshformat.file_type >> meshformat.data_size;
		return f;
	}

}//end namespace GMSH

#endif // !__MESHFORMAT_H
