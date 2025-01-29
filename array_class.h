#ifndef _MATRIX_H
#define _MATRIX_H

#include <vector>
#include <iostream>
#include <stdexcept>
#include <sstream>
#include <fstream>
#include <limits>
#include <typeinfo>
#include <cmath> // Pour std::abs
#include <complex>
#include <format>
#include <assert.h>
#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */



typedef double dble;



template<typename T>
inline void swap(T& a, T& b)
{
    T tmp(a);
    a = b;
    b = tmp;
}

template<typename T>
inline T MAX(T& a, T& b)
{
    if (a > b)
    {
        return a;
    }
    else
    {
        return b;
    }
}

template<typename T>
inline T MIN(T& a, T& b)
{
    if (a < b)
    {
        return a;
    }
    else
    {
        return b;
    }
}

template<typename T>
inline T ABS(T& a)
{
    if (a < 0)
    {
        return -a;
    }
    else
    {
        return a;
    }
}



template<typename T>
class Matrix
{
private:
    size_t _rows, _cols;
    std::vector<T> x;

public: 
    
    Matrix(/**/) :Matrix<T>(3, 3, 0) {}// Matrix default constructor
    Matrix(size_t rows, size_t cols, T initValue = T())//Matrix constructor
        : _rows(rows), _cols(cols), x(_rows* _cols, initValue) {}

    //Matrix(const Matrix<T>& mat) :_rows(mat._rows), _cols(mat._cols), x(_rows* _cols, 0) {}; //Matrix copy constructor

    /*
    Matrix<T>& operator=(const Matrix<T>& mat)// Matrix copy assignment
    {
        if (this != &mat) {  // Vérification de l'auto-affectation
            this->_rows = mat._rows;
            this->_cols = mat._cols;
            this->x = mat.x; // Copier les données (selon votre structure)
        }
        return *this; // Retourner l'objet courant
    }*/

    Matrix<T>& operator=(const Matrix<T>& mat) = default;//copy by assignment


    bool operator==(const Matrix<T>& mat);// compare les matrices élément par élément
    auto operator<=>(const Matrix<T>& tab)const = default;

        
    Matrix<T>& operator+=(const Matrix<T>& mat);
    Matrix<T> operator+(const Matrix<T>& mat);
    Matrix<T> operator+(Matrix<T>& mat);    
    Matrix<T>& operator-=(const Matrix<T>& mat);
    Matrix<T> operator-(const Matrix<T>& mat);    
    Matrix<T> operator-(Matrix<T>& mat);    
    Matrix<T> operator*=(const T& c);
    Matrix<T> operator*(Matrix<T>& mat);
    Matrix<T> operator*(T& c);
    Matrix<T> operator/=(const T& c);
    Matrix<T> operator/(T& c);
    Matrix<T> operator/(const T& c)const;
    
    
    //Matrix<T> operator/(Matrix<T>& vec);

    T& operator()(size_t i, size_t j);


    Matrix<T> transpose();
    T trace();

    Matrix<T> eye();
    Matrix<T> ones();
    Matrix<T> zeros();

    friend Matrix<T> operator-(Matrix<T> mat);
    friend Matrix<T> operator*(const Matrix<T>& mat1, const Matrix<T>& mat2);
    friend Matrix<T> operator*(long c, const Matrix<T>& mat2);
    friend Matrix<T> operator*(const Matrix<T>& mat2, long c);
    //friend Matrix<T> operator*(Matrix<T>& mat2, double c);


    void setValue(size_t i, size_t j, T val);
    void setValue2(size_t k, T val);
    T getValue(size_t i, size_t j);


    bool isSymetric();
    bool isTridiagonal();
    bool isDiagonal();
    bool isLowerTriangular();//todo
    bool isUpperTriangular();//todo
    bool isHessenberg();//todo
    bool isToeplitz();//todo

    bool isDiagonalDominant(size_t& i);

    Matrix<T> renormalization();

    Matrix<T> transformation(char& axis, dble& angle_RAD);
    Matrix<T> translation(char& axis, dble& angle_RAD, dble& tx, dble& ty, dble& tz);

    void print(std::string message);

    // Obtenir les dimensions
    size_t getRows() const { return this->_rows; }
    size_t getCols() const { return this->_cols; }

    void setSizeRows(size_t m) { _rows = m; }
    void setSizeCols(size_t n) { _cols = n; }

  
    friend std::ostream& operator<<(std::ostream& f, Matrix<T>& mat)
    {
        f << "\n";
        for (size_t i = 0; i < mat.getRows(); i++)
        {
            f << "\t\t";
            for (size_t j = 0; j < mat.getCols(); j++)
            {
                f << mat(i,j) << " ";
            }
            f << "\n";
        }
        f << "\n";
        return f;
    }


    //linear algebra
    Matrix<T> rowDeletion(Matrix<T>& mat, int rowNumber, int columnNumber);
    T expo(int n);
    Matrix<T> coMatrix();    
    T det(Matrix<T>& mat1);
    T det();
    Matrix<T>inv();
    

    Matrix<T> reshape(int _newRows, int _newCols);
    Matrix<T> random(); 

    T norm_1();//standard norm 1 of the table
    T norm_2();//standard norm 2 (Frobenius) of the table
    T sumOfLineElement(int& lineNumber);
    T sumOfColumnElements(int& columnNumber);
    T getSumElementsOfColumn(int& col);
    T getSumElementsOfRow(int& row);

    T getElementMaxOf(char choice, int number);//find most important values ​​from rows.


    //operator special si Matrix<T> est une matrice de dim(m,1) ou dim(1,n)
    T normL2();
    T& operator()(size_t i);
    T operator,(Matrix<T>& vec);//dot product 
    Matrix<T> operator^(Matrix<T>& vec);//cross product (dim vector vec is equal to 3.)



};





/* MATRIX<T> */

/*
template<typename T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T>& mat)
{
    return (*this);
}*/


template<typename T>
bool Matrix<T>::isDiagonalDominant(size_t& i)
{
    bool test = { true };//on suppose que la matrice A est a diagonal dominante !
    T sum_j{ 0.0 };
    for (size_t j = 0; j < this->getCols(); j++)
    {
        if (i != j)
        {
            sum_j += std::abs(this->getValue(i, j));
        }        
    }
    if (std::abs(this->getValue(i, i)) <= sum_j)
        test = false;
    return test;
}


template<typename T>
Matrix<T> Matrix<T>:: operator*(T& c)
{
    Matrix<T> temp(*this);
    for (size_t i = 0; i < this->_rows; i++)
    {
        for (size_t j = 0; j < this->_cols; j++)
        {
            dble val = this->getValue(i, j);
            temp.setValue(i, j, val * c);
        }
    }
    return temp;
}

/*
template<typename T>
Matrix<T> Matrix<T>::operator/(Matrix<T>& vec)
{
    assert(this->getRows() == vec.getRows());
    assert(this->getCols() == vec.getCols());
    Matrix<T>temp(vec.getRows(), vec.getCols());
    for (size_t i = 0; i < this->getRows(); i++)
    {
        assert(vec.getValue(i, 0) != 0);
        T val = this->getValue(i, 0) / vec.getValue(i, 0);
        temp.setValue(i, 0, val);
    }
    return temp;
}*/

template<typename T>
T Matrix<T>::getSumElementsOfColumn(int& col)
{
    T sum = { 0.0 };
    for (int i = 0; i < this->getRows(); i++)
    {
        sum += this->getValue(i, col);
    }
    return sum;
}

template<typename T>
T Matrix<T>::getSumElementsOfRow(int& row)
{
    T sum = T(0);
    for (int j = 0; j < this->getCols(); j++)
    {
        sum += this->getValue(row, j);
    }
    return sum;
}


template<typename T>
T Matrix<T>::getElementMaxOf(char choice, int number)//find most important values ​​from rows('r') or columns('c').
{
    // Initialisation avec la plus petite valeur possible pour T
    T val = T(0); // std::numeric_limits<T>::lowest();
    switch (choice) {
    case 'r':// Rechercher dans une ligne
        // code block
        if (number < 0 || number >= this->getCols()) {
            std::cerr << "L'indice de colonne est invalide.\n";
            return val;
        }
        for (int j = 0; j < this->getCols(); j++) {
            T currentValue = this->getValue(number, j);
            if (std::abs(currentValue) > std::abs(val)) {
                val = currentValue;
            }
        }
        break;
    case 'c':// Rechercher dans une colonne
        // code block
        if (number < 0 || number >= this->getRows()) {
            std::cerr << "L'indice de ligne est invalide.\n";
            return val;
        }
        for (int i = 0; i < this->getRows(); i++) {
            T currentValue = this->getValue(i, number);
            if (std::abs(currentValue) > std::abs(val)) {
                val = currentValue;
            }
        }
        break;
    default:
        // code block
        std::cerr << "Le premier parametre doit etre 'r' (ligne) ou 'c' (colonne).\n";
    } 
    return val;    
}



template<typename T>
T Matrix<T>::normL2()
{
    assert(this->getCols() == 1);
    T sum{ 0 };
    for (size_t i = 0; i < this->getRows(); i++)
    {
        sum += this->getValue(i, 0) * this->getValue(i, 0);
    }
    return sqrt(sum);
}


template<typename T>
T Matrix<T>::norm_2()
{
    T norm_l2{ 0 };
    T sum_i{ 0 };
    for (int i = 0; i < this->getRows(); i++)
    {
        T sum_j{ 0 };
        for (int j = 0; j < this->getCols(); j++)
        {
            T squarematcoefij = this->getValue(i, j) * this->getValue(i, j);
            sum_j += squarematcoefij;
        }
        sum_i += sum_j;
    }
    norm_l2 = sum_i;
    return sqrt(norm_l2);
}


template<typename T>
T Matrix<T>::norm_1()
{
    std::vector<T> vmatcoef;
    T norm_l1{ 0 };
    for (int j = 0; j < this->getCols(); j++)
    {
        T sum{ 0 };
        for (int i = 0; i < this->getRows(); i++)
        {
            T matcoefij = this->getValue(i, j);
            sum += ope::ABS(matcoefij);
        }
        vmatcoef.push_back(sum);
    }
    norm_l1 = ope::maxVector(vmatcoef);
    return norm_l1;
}



template<typename T>
Matrix<T> Matrix<T>::renormalization()
{
    Matrix<T>tmp(this->getRows(), this->getCols(), T(0));
    for (int i = 0; i < this->getRows(); i++)
    {
        dble value = { 0.0 }; 
        for (int j = 0; j < this->getCols(); j++)
        {
            value = dble( (this->getValue(i, j)) / (this->getElementMaxOf('r',i)));
            tmp.setValue(i, j, value);
        }
    }
    return tmp;
}


template<typename T>
void Matrix<T>::setValue2(size_t k, T val)
{
    if (k >= this->_rows*this->_cols) {
        throw std::out_of_range("Index hors limite !");
    }
    x[k] = val;
}


template<typename T>
Matrix<T> Matrix<T>::reshape(int _newRows, int _newCols)
{
    if (_newRows*_newCols != _rows * _cols) {
        throw std::invalid_argument("Nombre d'éléments incompatible avec les dimensions.");
    }
    Matrix<T>tmp(_newRows, _newCols, 0);
    for (size_t k = 0; k < _rows * _cols; k++)
    {
        tmp.setValue2(k, this->x[k]);
    }
    return tmp;
}



template<typename T>
Matrix<T> Matrix<T>::random()
{
    Matrix<T> R(_rows, _cols, 0);
    for (size_t i = 0; i < _rows; i++)
    {
        for (size_t j = 0; j < _cols; j++)
        {
            int val  = rand() % 100;
            R.setValue(i, j, val);
        }
    }
    return R;
}


template<typename T>
bool Matrix<T>::isDiagonal()
{
    assert(_rows == _cols);
    bool test = false;
    size_t cpt = { 0 };
    for (size_t i = 0; i < _rows; i++)
    {
        for (size_t j = 0; j < _cols; j++)
        {
            if ((this->getValue(i, i) != 0) && (this->getValue(i, j) == 0))
                cpt++;
        }
    }
    if (cpt == _rows)
        test = true;
    return test;
}



template<typename T>
Matrix<T> Matrix<T>::translation(char& axis, dble& angle_RAD, dble& tx, dble& ty, dble& tz)
{
    assert((_rows == 4) && (_cols == 4));
    Matrix<T> P(4, 4);

    switch (axis)
    {
    case 'x':
        // The axis of rotation is “x”.
        P.setValue(0, 0, 1); P.setValue(0, 1, 0); P.setValue(0, 2, 0); P.setValue(0, 3, tx);
        P.setValue(1, 0, 0); P.setValue(1, 1, cos(angle_RAD)); P.setValue(1, 2, -sin(angle_RAD)); P.setValue(1, 3, ty);
        P.setValue(2, 0, 0); P.setValue(2, 1, sin(angle_RAD)); P.setValue(2, 2, cos(angle_RAD)); P.setValue(2, 3, tz);
        P.setValue(3, 0, 0); P.setValue(3, 1, 0); P.setValue(3, 2, 0); P.setValue(3, 3, 1);
        break;
    case 'X':
        // The axis of rotation is “X”.
        P.setValue(0, 0, 1); P.setValue(0, 1, 0); P.setValue(0, 2, 0); P.setValue(0, 3, tx);
        P.setValue(1, 0, 0); P.setValue(1, 1, cos(angle_RAD)); P.setValue(1, 2, -sin(angle_RAD)); P.setValue(1, 3, ty);
        P.setValue(2, 0, 0); P.setValue(2, 1, sin(angle_RAD)); P.setValue(2, 2, cos(angle_RAD)); P.setValue(2, 3, tz);
        P.setValue(3, 0, 0); P.setValue(3, 1, 0); P.setValue(3, 2, 0); P.setValue(3, 3, 1);
        break;
    case 'y':
        // The axis of rotation is “y”.            
        P.setValue(0, 0, cos(angle_RAD)); P.setValue(0, 1, 0); P.setValue(0, 2, -sin(angle_RAD)); P.setValue(0, 3, tx);
        P.setValue(1, 0, 0); P.setValue(1, 1, 1); P.setValue(1, 2, 0); P.setValue(1, 3, ty); P.setValue(1, 3, ty);
        P.setValue(2, 0, sin(angle_RAD)); P.setValue(2, 1, 0); P.setValue(2, 2, cos(angle_RAD)); P.setValue(2, 3, tz);
        P.setValue(3, 0, 0); P.setValue(3, 1, 0); P.setValue(3, 2, 0); P.setValue(3, 3, 1);
        break;
    case 'Y':
        // The axis of rotation is “Y”.            
        P.setValue(0, 0, cos(angle_RAD)); P.setValue(0, 1, 0); P.setValue(0, 2, -sin(angle_RAD)); P.setValue(0, 3, tx);
        P.setValue(1, 0, 0); P.setValue(1, 1, 1); P.setValue(1, 2, 0); P.setValue(1, 3, ty);
        P.setValue(2, 0, sin(angle_RAD)); P.setValue(2, 1, 0); P.setValue(2, 2, cos(angle_RAD)); P.setValue(2, 3, tz);
        P.setValue(3, 0, 0); P.setValue(3, 1, 0); P.setValue(3, 2, 0); P.setValue(3, 3, 1);
        break;
    case 'z':
        // The axis of rotation is “z”.
        P.setValue(0, 0, cos(angle_RAD)); P.setValue(0, 1, -sin(angle_RAD)); P.setValue(0, 2, 0); P.setValue(0, 3, tx);
        P.setValue(1, 0, sin(angle_RAD)); P.setValue(1, 1, cos(angle_RAD)); P.setValue(1, 2, 0); P.setValue(1, 3, ty);
        P.setValue(2, 0, 0); P.setValue(2, 1, 0); P.setValue(2, 2, 1); P.setValue(2, 3, tz);
        P.setValue(3, 0, 0); P.setValue(3, 1, 0); P.setValue(3, 2, 0); P.setValue(3, 3, 1);
        break;
    case 'Z':
        // The axis of rotation is “Z”.
        P.setValue(0, 0, cos(angle_RAD)); P.setValue(0, 1, -sin(angle_RAD)); P.setValue(0, 2, 0); P.setValue(0, 3, tx);
        P.setValue(1, 0, sin(angle_RAD)); P.setValue(1, 1, cos(angle_RAD)); P.setValue(1, 2, 0); P.setValue(1, 3, ty);
        P.setValue(2, 0, 0); P.setValue(2, 1, 0); P.setValue(2, 2, 1); P.setValue(2, 3, tz);
        P.setValue(3, 0, 0); P.setValue(3, 1, 0); P.setValue(3, 2, 0); P.setValue(3, 3, 1);
        break;
    case 'xyz':
        angle_RAD = 0.0;
        P.setValue(0, 0, cos(angle_RAD)); P.setValue(0, 1, -sin(angle_RAD)); P.setValue(0, 2, 0); P.setValue(0, 3, tx);
        P.setValue(1, 0, sin(angle_RAD)); P.setValue(1, 1, cos(angle_RAD)); P.setValue(1, 2, 0);
        P.setValue(2, 0, 0); P.setValue(2, 1, 0); P.setValue(2, 2, 1); P.setValue(2, 3, tz);
        P.setValue(3, 0, 0); P.setValue(3, 1, 0); P.setValue(3, 2, 0); P.setValue(3, 3, 1);
        break;
    case 'XYZ':
        angle_RAD = 0.0;
        P.setValue(0, 0, cos(angle_RAD)); P.setValue(0, 1, -sin(angle_RAD)); P.setValue(0, 2, 0); P.setValue(0, 3, tx);
        P.setValue(1, 0, sin(angle_RAD)); P.setValue(1, 1, cos(angle_RAD)); P.setValue(1, 2, 0);
        P.setValue(2, 0, 0); P.setValue(2, 1, 0); P.setValue(2, 2, 1); P.setValue(2, 3, tz);
        P.setValue(3, 0, 0); P.setValue(3, 1, 0); P.setValue(3, 2, 0); P.setValue(3, 3, 1);
        break;

    default:
        std::cout << "\tThe axis are either 'x', 'y' or 'z'. Type 'x' or 'X', ..." << std::endl;
        break;
    }

    return P;
}


template<typename T>
Matrix<T> Matrix<T>::transformation(char& axis, dble& angle_RAD)
{
    assert((_rows == 3) && (_cols == 3));
    Matrix<T> P(3, 3);

    switch (axis)
    {
    case 'x':
        // The axis of rotation is “x”.
        P.setValue(0, 0, 1); P.setValue(0, 1, 0); P.setValue(0, 2, 0);
        P.setValue(1, 0, 0); P.setValue(1, 1, cos(angle_RAD)); P.setValue(1, 2, -sin(angle_RAD));
        P.setValue(2, 0, 0); P.setValue(2, 1, sin(angle_RAD)); P.setValue(2, 2, cos(angle_RAD));
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
    case 'xyz':
        angle_RAD = 0.0;
        P.setValue(0, 0, cos(angle_RAD)); P.setValue(0, 1, -sin(angle_RAD)); P.setValue(0, 2, 0);
        P.setValue(1, 0, sin(angle_RAD)); P.setValue(1, 1, cos(angle_RAD)); P.setValue(1, 2, 0);
        P.setValue(2, 0, 0); P.setValue(2, 1, 0); P.setValue(2, 2, 1);
        break;
    case 'XYZ':
        angle_RAD = 0.0;
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







template<typename T>
Matrix<T> Matrix<T>::rowDeletion(Matrix<T>& mat, int rowNumber, int columnNumber)
{
    assert(mat.getCols() == mat.getRows());
    Matrix<T> dest(mat.getRows() - 1, mat.getCols() - 1, 0);
    size_t l = { 0 }, c;
    for (size_t ii = 0; ii < mat.getRows(); ii++)
    {
        if (ii != rowNumber)
        {
            c = { 0 };
            for (size_t jj = 0; jj < mat.getCols(); jj++)
            {
                if (jj != columnNumber)
                {
                    T val = mat.getValue(ii, jj);
                    dest.setValue(l, c, val);
                    c++;
                }
            }
            l++;
        }
    }
    return dest;
}

template<typename T>
T Matrix<T>::expo(int n)
{
    if ((n % 2) != 0) //n est pair
        return 1;
    else
        return -1;
}

template<typename T>
Matrix<T> Matrix<T>::coMatrix()
{
    Matrix<T> mat2(_rows, _cols, 0);
    Matrix<T> tmp(_rows, _cols, 0);
    assert(_rows == _cols);
    if (_rows == 1)
    {
        T val = 1.0;
        tmp.setValue(0, 0, val);
    }
    else
    {
        for (size_t i = 0; i < _rows; i++)
        {
            for (size_t j = 0; j < _cols; j++)
            {
                mat2 = this->rowDeletion(*this, i, j);
                T val = this->expo(i + j) * mat2.det();
                tmp.setValue(i, j, val);
            }
        }
    }
    return tmp;
}




template<typename T>
Matrix<T> Matrix<T>::inv()
{
    Matrix<T> tmp(_rows, _cols, 0);
    assert(this->det() != 0);
    tmp = this->coMatrix().transpose() / this->det();
    return tmp;
}

template<typename T>
T Matrix<T>::det(Matrix<T>& mat1)
{
    //recursivité
    assert(mat1.getCols() == mat1.getRows());
    Matrix<T> mat2(mat1.getRows(), mat1.getCols(), 0);
    T x = { 0 };
    if (mat1.getRows() == 1)
        return mat1.getValue(0, 0);
    else
    {
        for (int i = 0; i < mat1.getRows(); i++)
        {
            mat2 = mat2.rowDeletion(mat1, i, 0);//suivant les lignes
            //x += std::pow(-1.0, i) * mat1.getValue(i, 0) * det(mat2);
            x += mat1.expo(i) * mat1.getValue(i, 0) * det(mat2);
        }
        return x;
    }
}

template<typename T>
T Matrix<T>::det()
{
    return det(*this);
}













template<typename T>
bool Matrix<T>::isSymetric()
{
    bool test = false;
    Matrix<T> AT = this->transpose();
    if (this->operator==(AT))
        test = true;    
    return test;
}


template<typename T>
bool Matrix<T>::operator==(const Matrix<T>& mat)// compare les matrices élément par élément
{
    bool isequivalent = true;
    for (size_t i = 0; i < _rows; i++)
    {
        for (size_t j = 0; j < _cols; j++)
        {
            size_t k = i * _cols + j;
            if (this->getValue(i, j) != mat.x[k])
            {
                isequivalent = false;
                break;  // On peut arrêter ici si les matrices diffèrent
            }
        }
        if (!isequivalent) break;
    }
    return isequivalent;
}


template<typename T>
Matrix<T>& Matrix<T>::operator+=(const Matrix<T>& mat)
{
    for (size_t i = 0; i < mat.getRows(); i++)
    {
        for (size_t j = 0; j < mat.getCols(); j++)
        {
            if (i >= this->_rows || j >= this->_cols)
            {
                throw std::out_of_range("Index hors limite !");
            }
            size_t k = i * this->_cols + j;
            x[k] += mat.x[k];
        }
    }
    return (*this);
}


template<typename T>
Matrix<T>& Matrix<T>::operator-=(const Matrix<T>& mat)
{
    for (size_t i = 0; i < mat.getRows(); i++)
    {
        for (size_t j = 0; j < mat.getCols(); j++)
        {
            if (i >= _rows || j >= this->_cols)
            {
                throw std::out_of_range("Index hors limite !");
            }
            size_t k = i * this->_cols + j;
            x[k] -= mat.x[k];
        }
    }
    return (*this);
}


template<typename T>
T& Matrix<T>::operator()(size_t i, size_t j)
{
    if (i >= this->_rows || j >= this->_cols) {
        throw std::out_of_range("Index hors limite !");
    }
    size_t k = i * this->_cols + j;
    return x[k];
}



template<typename T>
Matrix<T> Matrix<T>::operator+(const Matrix<T>& mat)
{
    assert(mat.getCols() == this->getCols());
    assert(mat.getRows() == this->getRows());
    Matrix<T> temp(this->getRows(), this->getCols(), 0);
    for (size_t i = 0; i < mat.getRows(); i++)
    {
        for (size_t j = 0; j < mat.getCols(); j++)
        {
            size_t k = i * mat.getCols() + j;
            temp.x[k] = x[k] + mat.x[k];
        }
    }
    return temp;
}

template<typename T>
Matrix<T> Matrix<T>::operator+(Matrix<T>& mat)
{
    assert(mat.getRows() == this->getRows());
    assert(mat.getCols() == this->getCols());
    Matrix<T> temp(this->getRows(), this->getCols(), 0);
    for (size_t i = 0; i < mat.getRows(); i++)
    {
        for (size_t j = 0; j < mat.getCols(); j++)
        {
            dble val = this->getValue(i, j) + mat.getValue(i, j);
            temp.setValue(i, j, val);
        }
    }
    return temp;
}





template<typename T>
Matrix<T> Matrix<T>:: operator-(const Matrix<T>& mat)
{
    assert(mat.getRows() == this->getRows());
    assert(mat.getCols() == this->getCols());
    Matrix<T> temp(this->getRows(), this->getCols(), 0);
    for (size_t i = 0; i < mat.getRows(); i++)
    {
        for (size_t j = 0; j < mat.getCols(); j++)
        {
            size_t k = i * this->_cols + j;
            temp.x[k] = this->x[k] - mat.x[k];
        }
    }
    return temp;
}


template<typename T>
Matrix<T> Matrix<T>:: operator-(Matrix<T>& mat)
{
    assert(mat.getRows() == this->getRows());
    assert(mat.getCols() == this->getCols());
    Matrix<T> temp(this->getRows(), this->getCols(), 0);
    for (size_t i = 0; i < mat.getRows(); i++)
    {
        for (size_t j = 0; j < mat.getCols(); j++)
        {
            dble val = this->getValue(i, j) - mat.getValue(i, j);
            temp.setValue(i, j, val);
        }
    }
    return temp;
}



template<typename T>
std::ostream& operator<<(std::ostream& f, Matrix<T>& mat)
{
    f << "\n";
    for (size_t i = 0; i < mat.getRows(); i++)
    {
        f << "\t\t";
        for (size_t j = 0; j < mat.getCols(); j++)
        {
            f << mat(i, j) << " ";
        }
        f << "\n";
    }
    f << "\n";
    return f;
}



template<typename T>
Matrix<T> Matrix<T>:: operator*(Matrix<T>& mat)
{
    assert(mat.getRows() == this->getCols());
    Matrix<T> temp(this->getRows(), mat.getCols(), 0);
    for (size_t i = 0; i < this->getRows(); i++)
    {
        for (size_t j = 0; j < mat.getCols(); j++)
        {
            double val{ 0 };
            for (size_t k = 0; k < this->getCols(); k++)
            {
                val += this->getValue(i, k) * mat.getValue(k, j);
            }
            temp.setValue(i, j, val);
        }
    }
    return temp;
}


template<typename T>
Matrix<T> Matrix<T>:: operator/(const T& c)const
{
    assert(c != 0);
    Matrix<T> temp(*this);
    for (size_t i = 0; i < this->_rows; i++)
    {
        for (size_t j = 0; j < this->_cols; j++)
        {
            size_t k = i * this->_cols + j;
            //dble val = this->getValue(i, j);
            //val = val / c;
            //temp.setValue(i, j, val);
            temp.x[k] = x[k] / c;
        }
    }
    return temp;
}


template<typename T>
Matrix<T> Matrix<T>:: operator/(T& c)
{
    assert(c != 0);
    Matrix<T> temp(*this);
    for (size_t i = 0; i < this->_rows; i++)
    {
        for (size_t j = 0; j < this->_cols; j++)
        {
            dble val = this->getValue(i, j);
            //val = val / c;
            temp.setValue(i, j, val/c);
        }
    }
    return temp;
}


template<typename T>
Matrix<T> Matrix<T>::operator/=(const T& c)
{
    assert(c != 0);
    for (size_t i = 0; i < _rows; i++)
    {
        for (size_t j = 0; j < _cols; j++)
        {
            size_t k = i * this->_cols + j;
            this->x[k] /= c;
        }
    }
    return (*this);
}

template<typename T>
Matrix<T> Matrix<T>::operator*=(const T& c)
{
    assert(c != 0);
    for (size_t i = 0; i < _rows; i++)
    {
        for (size_t j = 0; j < _cols; j++)
        {
            size_t k = i * this->_cols + j;
            this->x[k] *= c;
        }
    }
    return (*this);
}




template<typename T>
Matrix<T> Matrix<T>:: eye()
{
    Matrix<T> eye(this->_rows, this->_cols, 0);
    for (size_t i = 0; i < this->_rows; i++)
    {
        for (size_t j = 0; j < this->_cols; j++)
        {
            if (i >= this->_rows || j >= this->_cols)
            {
                throw std::out_of_range("Index hors limite !");
            }
            if (i == j)
            {
                eye(i, j) = 1.0;
            }
            else
            {
                eye(i, j) = 0.0;
            }
        }
    }
    return eye;
}

template<typename T>
Matrix<T> Matrix<T>:: ones()
{
    Matrix<T> ones(this->_rows, this->_cols, 1);
    return ones;
}

template<typename T>
Matrix<T> Matrix<T>:: zeros()
{
    Matrix<T> temp(this->_rows, this->_cols, 0);
    return temp;
}







template<typename T>
T Matrix<T>::getValue(size_t i, size_t j)
{
    if (i >= this->_rows)
    {
        throw std::out_of_range("Index hors limite !");
    }
    size_t k = i * this->_cols + j;
    return x[k];
}



template <typename T>
void Matrix<T>::setValue(size_t i, size_t j, T val)
{
    assert((i >= 0) && (i < _rows));
    assert((j >= 0) && (j < _cols));
    size_t k = i * _cols + j;
    x[k] = val; // Accès sécurisé
}


template<typename T>
T Matrix<T>::trace()
{
    T s{ 0 };
    for (size_t i = 0; i < this->_rows; i++)
    {
        if (i >= this->_rows)
        {
            throw std::out_of_range("Index hors limite !");
        }
        s += this->getValue(i, i);
    }
    return s;
}

template<typename T>
Matrix<T> Matrix<T>:: transpose()
{
    Matrix<T> temp(this->_cols, this->_rows, 0);
    for (size_t i = 0; i < this->_rows; i++)
    {
        for (size_t j = 0; j < this->_cols; j++)
        {
            temp(j, i) = (*this)(i, j);
        }
    }
    return temp;
}


template<typename T>
void Matrix<T>::print(std::string message)
{
    std::cout << "\n";
    std::cout << "\n\n\t" << message << "\n";
    for (size_t i = 0; i < this->_rows; i++)
    {
        for (size_t j = 0; j < this->_cols; j++)
        {
            if (this->getValue(i, j) < 0.0)
            {
                std::cout << "\t" << this->getValue(i, j) << " ";
            }
            else
            {
                std::cout << "\t " << this->getValue(i, j) << " ";
            }

        }
        std::cout << "\n";
    }
    std::cout << "\n";
}













template<typename T>
Matrix<T> operator*(const Matrix<T>& mat1, const Matrix<T>& mat2)
{
    Matrix< T> temp(mat1._rows, mat2._rows, 0);
    for (size_t i = 0; i < mat1.getRows(); i++)
    {
        for (size_t j = 0; j < mat2.getCols(); j++)
        {
            T sum = { 0 };
            for (size_t k = 0; k < mat1.getCols(); k++)
            {
                sum += mat1(i,k) * mat2(k,j);
            }
            temp(i,j) = sum;
        }
    }
    return temp;
}

template<typename T>
Matrix<T> operator*(T c, Matrix<T>& mat2)
{
    Matrix<T> temp(mat2.getRows(), mat2.getCols(), 0);
    for (size_t i = 0; i < mat2.getRows(); i++)
    {
        for (size_t j = 0; j < mat2.getCols(); j++)
        {
            temp(i, j) = c * mat2(i, j);
        }
    }
    return temp;
}

template<typename T>
Matrix<T> operator*(Matrix<T>& mat2, T c)
{
    Matrix<T> temp(mat2.getRows(), mat2.getCols(), 0);
    for (size_t i = 0; i < mat2.getRows(); i++)
    {
        for (size_t j = 0; j < mat2.getCols(); j++)
        {
            temp(i, j) = mat2(i, j) * c;
        }
    }
    return temp;
}


/*

template<typename T>
Matrix<T> operator*(Matrix<T>& mat2, double c)
{
    Matrix<T> temp(mat2.getRows(), mat2.getCols(), 0);
    for (size_t i = 0; i < mat2.getRows(); i++)
    {
        for (size_t j = 0; j < mat2.getCols(); j++)
        {
            double val = this->getValue(i, j) * c;
            temp.setValue(i, j, val);
        }
    }
    return temp;
}

*/








/*
* ici: PARTIE SPECIALE
* //operator special si Matrix<T> est une matrice de dim(m,1) ou dim(1,n)
*/

template<typename T>
T& Matrix<T>::operator()(size_t i)
{
    assert(this->getCols() == 1);
    assert(i < this->getRows());
    return this->getValue(i, 0);
}

template<typename T>
T Matrix<T>::operator,(Matrix<T>& vec)//dot product 
{
    assert(this->getRows() == vec.getRows());
    assert(this->getCols() == vec.getCols() == 1);
    T val = { 0 };
    for (size_t i = 0; i < vec.getRows(); i++)
    {
        val += this->getValue(i, 0) * vec.getValue(i, 0);
    }
    return val;
}

template<typename T>
Matrix<T> Matrix<T>::operator^(Matrix<T>& vec)//cross product (dim vector vec is equal to 3.)
{
    assert(this->getRows() == vec.getRows() == 3);
    assert(this->getCols() == vec.getCols() == 1);
    Matrix<T>tmp(3, 1, 0);
    tmp.x[0] = this->x[1] * vec.x[2] - this->x[2] * vec.x[1];
    tmp.x[1] = this->x[2] * vec.x[0] - this->x[0] * vec.x[2];
    tmp.x[2] = this->x[0] * vec.x[1] - this->x[1] * vec.x[0];
    return tmp;
}
















#endif /* _MATRIX_H */


