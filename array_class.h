#ifndef _MATRIX_H
#define _MATRIX_H

#include <vector>
#include <iostream>
#include <stdexcept>
#include <sstream>
#include <fstream>
#include <assert.h>







template<typename T>
class Matrix
{
private:
    size_t _rows, _cols;
    std::vector<T> x;

public: 
    // Constructeur par défaut
    Matrix(size_t rows, size_t cols, T initValue = T())
        : _rows(rows), _cols(cols), x(_rows* _cols, initValue) {}


    Matrix<T>& operator=(const Matrix<T>& mat)
    {
        if (this != &mat) {  // Vérification de l'auto-affectation
            this->_rows = mat._rows;
            this->_cols = mat._cols;
            this->x = mat.x; // Copier les données (selon votre structure)
        }
        return *this; // Retourner l'objet courant
    }

    T& operator()(size_t i, size_t j);
    bool operator==(const Matrix<T>& mat);// compare les matrices élément par élément

    Matrix<T> operator+(const Matrix<T>& mat);
    Matrix<T>& operator+=(const Matrix<T>& mat);
    Matrix<T> operator-(const Matrix<T>& mat);
    Matrix<T>& operator-=(const Matrix<T>& mat);
    
    Matrix<T> operator*(Matrix<T>& mat);
    Matrix<T> operator/(const T& c); 

    Matrix<T> transpose();
    T trace();

    Matrix<T> eye();
    Matrix<T> ones();
    Matrix<T> zeros();

    friend Matrix<T> operator*(const Matrix<T>& mat1, const Matrix<T>& mat2);
    friend Matrix<T> operator*(long c, const Matrix<T>& mat2);
    friend Matrix<T> operator*(const Matrix<T>& mat2, long c);


    void setValue(size_t i, size_t j, T val);
    T getValue(size_t i, size_t j);


    bool isSymetric();
    bool isTridiagonal();

    void print(std::string message);

    // Obtenir les dimensions
    size_t getRows() const { return this->_rows; }
    size_t getCols() const { return this->_cols; }


  
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



};



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
    Matrix<T> temp(this->getRows(), this->getCols(), 0);
    for (size_t i = 0; i < mat.getRows(); i++)
    {
        for (size_t j = 0; j < mat.getCols(); j++)
        {
            if (i >= this->_rows || j >= this->_cols)
            {
                throw std::out_of_range("Index hors limite !");
            }
            size_t k = i * _rows + j;
            temp.x[k] = x[k] + mat.x[k];
        }
    }
    return temp;
}



template<typename T>
Matrix<T> Matrix<T>:: operator-(const Matrix<T>& mat)
{
    Matrix<T> temp(this->getRows(), this->getCols(), 0);
    for (size_t i = 0; i < mat.getRows(); i++)
    {
        for (size_t j = 0; j < mat.getCols(); j++)
        {
            if (i >= this->_rows || j >= this->_cols)
            {
                throw std::out_of_range("Index hors limite !");
            }
            size_t k = i * this->_cols + j;
            temp.x[k] = x[k] - mat.x[k];
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
    Matrix<T> temp(this->getRows(), this->getCols(), 0);
    for (size_t i = 0; i < this->getRows(); i++)
    {
        for (size_t j = 0; j < mat.getCols(); j++)
        {
            temp(i, j) = 0.0;
            for (size_t k = 0; k < this->getCols(); k++)
            {
                temp(i, j) += this->getValue(i, k) * mat(k, j);
            }
        }
    }
    return temp;
}


template<typename T>
Matrix<T> Matrix<T>:: operator/(const T& c)
{
    if (c == 0)
    {
        std::cerr << "\tZeroDivisionError ! " << "\n";
        exit(-1);
    }
    Matrix<T> temp(*this);
    for (size_t i = 0; i < this->_rows; i++)
    {
        for (size_t j = 0; j < this->_cols; j++)
        {
            if (i >= _rows || j >= this->_cols)
            {
                throw std::out_of_range("Index hors limite !");
            }
            size_t k = i * this->_cols + j;
            temp.x[k] = x[k] / c;
        }
    }
    return temp;
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
    if (i >= this->_rows || j >= this->_cols) {
        throw std::out_of_range("Index hors limite !");
    }
    size_t k = i * _cols + j;
    return x[k];
}



template<typename T>
void Matrix<T>::setValue(size_t i, size_t j, T val)
{
    if (i >= this->_rows || j >= this->_cols) {
        throw std::out_of_range("Index hors limite !");
    }
    size_t k = i * this->_cols + j;
    x[k] = val;
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










#endif /* _MATRIX_H */


