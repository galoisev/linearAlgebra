#ifndef _MATRICE_ROTATION_H
#define _MATRICE_ROTATION_H



//#include "../../tool/ToolLib/matrix.h"
#include <string>
#include <assert.h>
#include "../array.h"




#pragma once



using namespace std;





/* A rotation matrix is a transformation matrix that is used to perform a rotation in Euclidean space. */
template<class T>
linearAlgebra::Array<T> matrixRotation(char axis, double angle_RAD)
{
    linearAlgebra::Array<T> P(3, 3); 
        
    switch (axis)
    {
    case 'x':
        // The axis of rotation is “x”.
        P.setValue(0, 0, 1); P.setValue(0, 1, 0); P.setValue(0, 2, 0);
        P.setValue(1, 0, 0); P.setValue(1, 1, cos(angle_RAD)); P.setValue(1, 2, -sin(angle_RAD));
        P.setValue(2, 0, 0); P.setValue(2, 1, sin(angle_RAD)); P.setValue(2, 2,  cos(angle_RAD));
        break;
    case 'X':
        // The axis of rotation is “X”.
        P.setValue(0, 0, 1); P.setValue(0, 1, 0); P.setValue(0, 2, 0);
        P.setValue(1, 0, 0); P.setValue(1, 1, cos(angle_RAD)); P.setValue(1, 2, -sin(angle_RAD));
        P.setValue(2, 0, 0); P.setValue(2, 1, sin(angle_RAD)); P.setValue(2, 2, cos(angle_RAD));
        break;
    case 'y':
        // The axis of rotation is “y”.            
        P.setValue(0, 0, cos(angle_RAD)); P.setValue(0, 1, 0); P.setValue(0, 2, -sin(angle_RAD));
        P.setValue(1, 0, 0); P.setValue(1, 1, 1); P.setValue(1, 2, 0);
        P.setValue(2, 0, sin(angle_RAD)); P.setValue(2, 1, 0); P.setValue(2, 2, cos(angle_RAD));
        break;
    case 'Y':
        // The axis of rotation is “Y”.            
        P.setValue(0, 0, cos(angle_RAD)); P.setValue(0, 1, 0); P.setValue(0, 2, -sin(angle_RAD));
        P.setValue(1, 0, 0); P.setValue(1, 1, 1); P.setValue(1, 2, 0);
        P.setValue(2, 0, sin(angle_RAD)); P.setValue(2, 1, 0); P.setValue(2, 2, cos(angle_RAD));
        break;
    case 'z':
        // The axis of rotation is “z”.
        P.setValue(0, 0, cos(angle_RAD)); P.setValue(0, 1, -sin(angle_RAD)); P.setValue(0, 2, 0);
        P.setValue(1, 0, sin(angle_RAD)); P.setValue(1, 1, cos(angle_RAD)); P.setValue(1, 2, 0);
        P.setValue(2, 0, 0); P.setValue(2, 1, 0); P.setValue(2, 2, 1);
        break;
    case 'Z':
        // The axis of rotation is “Z”.
        P.setValue(0, 0, cos(angle_RAD)); P.setValue(0, 1, -sin(angle_RAD)); P.setValue(0, 2, 0);
        P.setValue(1, 0, sin(angle_RAD)); P.setValue(1, 1, cos(angle_RAD)); P.setValue(1, 2, 0);
        P.setValue(2, 0, 0); P.setValue(2, 1, 0); P.setValue(2, 2, 1);
        break;
    default:
        std::cout << "\tThe axis are either 'x', 'y' or 'z'. Type 'x' or 'X', ..." << std::endl;
        break;
    }

    return P;

}
















#endif /* _MATRICE_ROTATION_H */



