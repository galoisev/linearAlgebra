#include "mesh.h"
#ifndef __MESH_H
#define __MESH_H



#endif // !__MESH_H

Mesh::Mesh(std::string& fileName)
{
	std::ifstream f(filename);
	if (!f)
	{
		std::cerr << "\n\tFile: " << filename << " not found !" << "\n";
		exit(1);
	}

	std::cout << "\n\tFile Reading: " << filename << " in progress. " << "\n";
	f >> nod;

}
