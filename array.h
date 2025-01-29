#ifndef __ARRAY_H__
#define __ARRAY_H__



#pragma once



#include <iostream>
#include <sstream>
#include <fstream>
#include <math.h>
#include <vector>
#include <complex>
#include <format>
#include <stdexcept>
#include <cassert>






namespace ope
{
	template<class T>
	inline void swap(T& a, T& b)
	{
		T tmp(a);
		a = b;
		b = tmp;
	}

	template<class T>
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

	template<class T>
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

	template<class T>
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

	/*
	template<class T>
	inline T MIN(T& a, T& b) { a > b ? b : a; }

	template<class T>
	inline T ABS(T& a) { a > 0 ? a : -a; }
	*/

	template<class T>
	inline void print(T& a, const char* msg)
	{
		std::cout << "\n";
		std::cout << "\t" << msg << a << std::endl;
		std::cout << "\n";
	}

	template<class T>
	inline T maxVector(std::vector<T>& v)
	{
		T maxvalue = { -999999 };
		for (auto& x : v)
		{
			if (MAX(x, maxvalue) == x)
			{
				maxvalue = x;
			}
		}
		return maxvalue;
	}

}//end ope

namespace linearAlgebra
{

	using namespace ope;

	template<class T>
	class Array
	{
	private:
		void allocArrays();
		void freeMemoryArrays();
	public:
		int nx{}, ny{};
		T** x{};

		Array() = default; //default constructor
		Array(int nnx, int nny);//constructor
		~Array();//destructor
		Array(const Array<T>& tab);//copy constructor
		Array& operator=(const Array<T>& tab);//copy by assignment


		void setValue(int i, int j, T val);
		void setValue2(int i, int j, T& value);
		T getValue(int i, int j);
		T getValue2(int& i, int& j);

		Array operator  +(const Array<T>& tab);
		Array& operator+=(const Array<T>& tab);
		Array operator  -(const Array<T>& tab);
		Array& operator-=(const Array<T>& tab);
		Array operator  *(const Array<T>& tab);
		Array& operator*=(const Array<T>& tab);
		Array operator  /(const T& c);
		Array& operator/=(const T& c);
		Array operator  *(const T& c);
		auto  operator<=>(const Array<T>& tab)const = default;

		T& operator()(int i, int j);

		T& operator* ()const { return *x; }
		T* operator->()const { return  x; }

		T operator,(const linearAlgebra::Array<T>& tab);


		template<class U>
		friend std::ostream& operator<<(std::ostream& f, const linearAlgebra::Array<U>& tab);
		template<class U>
		friend std::istream& operator>>(std::istream& f, linearAlgebra::Array<U>& tab);




		friend linearAlgebra::Array<T> operator*(const linearAlgebra::Array<T>& tab1, const linearAlgebra::Array<T>& tab2);
		friend linearAlgebra::Array<T> operator*(T c, const linearAlgebra::Array<T>& tab);
		friend linearAlgebra::Array<T> operator*(const linearAlgebra::Array<T>& tab, T c);

		



		void printing(const char* message);

		T getSumElementsOfColumn(int& col);
		T getSumElementsOfRow(int& row);

		void setRowNumber(int& nrow) { this->nx = nrow; }
		void setColNumber(int& ncol) { this->ny = ncol; }

		int rowSize() { return nx; }
		int columnSize() { return ny; }
		int getSize() { return nx * ny; }

		//T sum(const int& i , const Array<T>& tab); A FAIRE...
		bool all(const linearAlgebra::Array<T>& tab);

		linearAlgebra::Array<T> ones();
		linearAlgebra::Array<T> zeros();
		linearAlgebra::Array<T> eye();

		T trace();

		linearAlgebra::Array<T> transformation(char& axis, double& angle_RAD);
		linearAlgebra::Array<T> translation(char& axis, double& angle_RAD, double& tx, double& ty, double& tz);


		linearAlgebra::Array<T> transpose();


		void printFile(const std::string& message, const std::string& filename);

		bool isLowerTriangular();
		bool isUpperTriangular();
		bool isTridiagonal();
		bool isDiagonal();
		bool isSymetric();
		bool isHessenberg();//todo
		bool isToeplitz();//todo
		bool isDiagonalDominant(size_t& i);//todo





		T norm_1();//standard norm 1 of the table
		T norm_2();//standard norm 2 (Frobenius) of the table

		T sumOfLineElement(int& lineNumber);
		T sumOfColumnElements(int& columnNumber);

		Array<T> renormalization();

		Array<T> reshape(int _newRows, int _newCols);//TODO/brief


		//inverse et determinant d'une matrice
		Array<T> rowDeletion(Array<T>& mat, int rowNumber, int columnNumber);
		T expo(int n);
		Array<T> coMatrix();
		T det(Array<T>& mat1);
		T det();
		Array<T>inv();


	};





	template<typename T>
	Array<T> Array<T>::rowDeletion(Array<T>& mat, int rowNumber, int columnNumber)
	{
		assert(mat.rowSize() == mat.columnSize());
		Array<T> dest(mat.rowSize() - 1, mat.columnSize()- 1, 0);
		size_t l = { 0 }, c;
		for (size_t ii = 0; ii < mat.rowSize(); ii++)
		{
			if (ii != rowNumber)
			{
				c = { 0 };
				for (size_t jj = 0; jj < mat.columnSize(); jj++)
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
	T Array<T>::expo(int n)
	{
		if ((n % 2) != 0) //n est pair
			return 1;
		else
			return -1;
	}

	template<typename T>
	Array<T> Array<T>::coMatrix()
	{
		Array<T> mat2(this->rowSize(), this->columnSize(), 0);
		Array<T> tmp(this->rowSize(), this->columnSize(), 0);
		assert(this->rowSize() == this->columnSize());
		if (this->rowSize() == 1)
		{
			T val = 1.0;
			tmp.setValue(0, 0, val);
		}
		else
		{
			for (size_t i = 0; i < this->rowSize(); i++)
			{
				for (size_t j = 0; j < this->columnSize(); j++)
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
	Array<T> Array<T>::inv()
	{
		Array<T> tmp(this->rowSize(), this->columnSize(), 0);
		assert(this->det() != 0);
		tmp = this->coMatrix().transpose() / this->det();
		return tmp;
	}

	template<typename T>
	T Array<T>::det(Array<T>& mat1)
	{
		//recursivité
		assert(mat1.columnSize() == mat1.rowSize());
		Array<T> mat2(mat1.rowSize(), mat1.columnSize(), 0);
		T x = { 0 };
		if (mat1.rowSize() == 1)
			return mat1.getValue(0, 0);
		else
		{
			for (int i = 0; i < mat1.rowSize(); i++)
			{
				mat2 = mat2.rowDeletion(mat1, i, 0);//suivant les lignes
				//x += std::pow(-1.0, i) * mat1.getValue(i, 0) * det(mat2);
				x += mat1.expo(i) * mat1.getValue(i, 0) * det(mat2);
			}
			return x;
		}
	}

	template<typename T>
	T Array<T>::det()
	{
		return det(*this);
	}


















	////////



	template<typename T>
	void Array<T>::freeMemoryArrays()
	{
		if (x)
		{
			for (int i = 0; i < nx; i++)
			{
				delete[] x[i];
			}
			delete[] x;
		}
		x = nullptr;
	}


	template<typename T>
	Array<T> Array<T>::reshape(int _newRows, int _newCols)
	{
		if (_newRows * _newCols != nx * ny) {
			throw std::invalid_argument("Nombre d'éléments incompatible avec les dimensions.");
		}
		Array<T>tmp(_newRows, _newCols, 0);
		for (size_t k = 0; k < nx *ny; k++)
		{
			tmp.setValue2(k, this->x[k]);
		}
		return tmp;
	}


	template<class T>
	Array<T> Array<T>::renormalization()
	{
		for (int i = 0; i < this->nx; i++)
		{
			double value = { 0.0 };
			for (int j = 0; j < this->ny; j++)
			{
				value = this->getValue(i, j) / this->getSumElementsOfRow(i);
				this->setValue(i, j, value);
			}
		}
	}





	template<class T>
	T linearAlgebra::Array<T>::sumOfLineElement(int& lineNumber)
	{
		double s = { 0 };
		for (int j = 0; j < ny; j++)
		{
			s += this->getValue(lineNumber, j);
		}
		return s;
	}


	template<class T>
	T linearAlgebra::Array<T>::sumOfColumnElements(int& columnNumber)
	{
		double s = { 0 };
		for (int i = 0; i < nx; i++)
		{
			s += this->getValue(i, columnNumber);
		}
		return s;
	}





	template<class T>
	Array<T> linearAlgebra::Array<T>::operator*(const T& c)
	{
		Array<T>tmp;
		for (int i = 0; i < this->nx; i++)
		{
			for (int j = 0; j < this->ny; j++)
			{
				tmp.x[i][j] = c * this->x[i][j];
			}
		}
		return tmp;
	}

	// à tester !!!
	template<class T>
	bool linearAlgebra::Array<T>::isDiagonal()
	{
		bool isdiagonal = { true };
		for (int i = 0; i < this->nx; i++)
		{
			int cpt_j = { 0 };
			for (int j = 0; j < this->ny; j++)
			{
				if ((i != j) && (this->getValue(i, j) != 0))
				{
					isdiagonal = { false };
					break;
				}
				else
				{
					continue;
				}
			}
		}
		return isdiagonal;
	}


	template<typename T>
	T linearAlgebra::Array<T>::norm_2()
	{
		T norm_l2{ 0 };
		T sum_i{ 0 };
		for (int i = 0; i < this->rowSize(); i++)
		{
			T sum_j{ 0 };
			for (int j = 0; j < this->columnSize(); j++)
			{
				T squarematcoefij = this->getValue2(i, j) * this->getValue2(i, j);
				sum_j += squarematcoefij;
			}
			sum_i += sum_j;
		}
		norm_l2 = sum_i;
		return sqrt(norm_l2);
	}


	template<typename T>
	T linearAlgebra::Array<T>::norm_1()
	{
		std::vector<T> vmatcoef;
		T norm_l1{ 0 };
		for (int j = 0; j < this->columnSize(); j++)
		{
			T sum{ 0 };
			for (int i = 0; i < this->rowSize(); i++)
			{
				T matcoefij = this->getValue2(i, j);
				sum += ope::ABS(matcoefij);
			}
			vmatcoef.push_back(sum);
		}		
		norm_l1 = ope::maxVector(vmatcoef);
		return norm_l1;		
	}

	template<typename T>
	bool linearAlgebra::Array<T>::isTridiagonal()
	{
		assert(this->nx() == this->ny());
		int NUMBER_OF_ZERO_ELEMENTS = this->rowSize()* this->rowSize() - (this->rowSize() + 2 * (this->rowSize() - 1));
		bool istridiagonal = { false };
		int cpt{ 0 };
		for (int i = 0; i < this->rowSize(); i++)
		{
			for (int j = 0; j < this->columnSize(); j++)
			{
				if ((i >= j + 2) || (j >= i + 2))
				{
					if (this->getValue2(i,j) == 0.0)
					{
						cpt += 1;
					}
				}
			}
		}
		if (cpt == NUMBER_OF_ZERO_ELEMENTS)
			istridiagonal = true;
		return istridiagonal;
	}







	template<class T>
	bool linearAlgebra::Array<T>::isUpperTriangular()
	{
		assert(this->nx() == this->ny());
		bool isUpper = { false };		
		if (this->transpose().isLowerTriangular())
			isUpper = true;
		/*
		int MAX_VALUE = (this->getSize() - this->rowSize()) / 2;
		int cpt{ 0 };
		for (int i = 0; i < this->nx; i++)
		{
			for (int j = 0; j < this->ny; j++)
			{
				if (i > j)
				{
					if (this->getValue2(i, j) == 0.0)
					{
						cpt += 1;
					}
				}
			}
		}
		if (cpt == MAX_VALUE)
			isUpper = true;
			*/
		return isUpper;
	}

	template<class T>
	bool linearAlgebra::Array<T>::isLowerTriangular()
	{
		assert(this->nx() == this->ny());
		bool isLower = { false };
		int MAX_VALUE = (this->getSize() - this->rowSize()) / 2;
		int cpt{ 0 };
		for (int i = 0; i < this->nx; i++)
		{
			for (int j = 0; j < this->ny; j++)
			{
				if (j > i)
				{
					if (this->getValue2(i, j) == 0.0)
					{
						cpt += 1;
					}
				}
			}
		}
		if (cpt == MAX_VALUE)
			isLower = true;
		return isLower;
	}



	template<class T>
	void linearAlgebra::Array<T>::printFile(const std::string& message, const std::string& filename)
	{
		std::ofstream fic(filename);
		fic << "\n";
		fic << "\t" << message << "\n";
		fic << "\n";
		for (int i = 0; i < nx; i++)
		{
			for (int j = 0; j < ny; j++)
			{
				fic << std::format( "\t{0:.5e}", this->getValue2(i, j) ) << " ";
				//fic << "\t" << this->getValue2(i,j) << " ";
			}
			fic << "\n";
		}
		fic << "\n";
		fic.close();
	}


	template<class T>
	T linearAlgebra::Array<T>::operator,(const linearAlgebra::Array<T>& tab)
	{
		assert((nx > 1) || (ny > 1));
		T value = {0.0};
		if (this->ny == 1)
		{
			for (int i = 0; i < nx; i++)
			{
				value += this->x[i][0] * tab.x[i][0];
			}
		}
		else if (this->nx == 1)
		{
			for (int j = 0; j < ny; j++)
			{
				value += this->x[0][j] * tab.x[0][j];
			}
		}
		else
		{
			std::cerr << "\tOupsss ! The dot product must have one of the two indices equal to one !";
		}

		return value;
	}


	template<class T>
	T linearAlgebra::Array<T>::getSumElementsOfColumn(int& col)
	{
		T sum = { 0.0 };
		for (int i = 0; i < nx; i++)
		{
			sum += this->getValue(i, col);
		}
		return sum;
	}

	template<class T>
	T linearAlgebra::Array<T>::getSumElementsOfRow(int& row)
	{
		T sum = { 0.0 };
		for (int j = 0; j < ny; j++)
		{
			sum += this->getValue(row, j);
		}
		return sum;
	}





	template<class T>
	void linearAlgebra::Array<T>::setValue2(int i, int j, T& value)
	{
		this->x[i][j] = value;
	}


	template<class T>
	std::ostream& operator<<(std::ostream& f, const linearAlgebra::Array<T>& tab)
	{
		f << "\n";
		for (int i = 0; i < tab.nx; i++)
		{
			f << "\t\t";
			for (int j = 0; j < tab.ny; j++)
			{
				f << std::format("\t{0:.2f}", tab.x[i][j]) << " ";//avec precision dans l'affichage des tremes de la matrice.
				//f << "\t" << tab.x[i][j] << " ";
			}
			f << "\n";
			
		}
		f << "\n";
		return f;
	}

	template<class T>
	std::istream& operator>>(std::istream& f, linearAlgebra::Array<T>& tab)
	{
		for (int i = 0; i < tab.nx; i++)
		{
			for (int j = 0; j < tab.ny; j++)
			{
				f >> tab.x[i][j] ;
			}
		}
		return f;
	}


	template<class T>
	linearAlgebra::Array<T> operator*(const linearAlgebra::Array<T>& tab1, const linearAlgebra::Array<T>tab2)
	{
		linearAlgebra::Array<T>prod(tab1.nx, tab2.ny);
		for (int i = 0; i < prod.nx; i++)
		{
			for (int j = 0; j < prod.ny; j++)
			{
				for (int k = 0; k < tab1.ny; k++)
				{
					prod.x[i][j] += tab1.x[i][k] * tab2.x[k][j];
				}				
			}
		}
		return prod;
	}

	template<class T>
	linearAlgebra::Array<T> operator*(T c, const linearAlgebra::Array<T>& tab)
	{
		linearAlgebra::Array<T>tmp(tab);
		for (int i = 0; i < tmp.nx; i++)
		{
			for (int j = 0; j < tmp.ny; j++)
			{
				tmp.x[i][j] = c * tab[i][j];
			}
		}
		return tab;
	}

	template<class T>
	linearAlgebra::Array<T> operator*(const linearAlgebra::Array<T>& tab, T c)
	{
		return c * tab;
	}








	template<class T>
	T linearAlgebra::Array<T>::trace()
	{
		T s = { 0 };
		for (int i = 0; i < nx; i++)
		{
			s += x[i][i];
		}
		return s;
	}


	template<class T>
	linearAlgebra::Array<T> linearAlgebra::Array<T>::transpose()
	{
		linearAlgebra::Array<T> tmp(*this);
		for (int i = 0; i < nx; i++)
		{
			for (int j = 0; j < ny; j++)
			{
				if (j > i)
				{
					ope::swap(tmp.x[i][j] , tmp.x[j][i]);
					//this->setValue(i, j, tmp.x[i][j]);
					//this->setValue(j, i, tmp.x[j][i]);
				}
			}
		}
		return tmp;
	}



	template<class T>
	linearAlgebra::Array<T> linearAlgebra::Array<T>::translation(char& axis, double& angle_RAD, double& tx, double& ty, double& tz)
	{
		linearAlgebra::Array<T> P(4, 4);

		switch (axis)
		{
		case 'x':
			// The axis of rotation is “x”.
			P.setValue(0, 0, 1); P.setValue(0, 1, 0); P.setValue(0, 2, 0); P.setValue(0,3,tx);
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

	template<class T>
	linearAlgebra::Array<T> linearAlgebra::Array<T>::transformation(char& axis, double& angle_RAD)
	{
		linearAlgebra::Array<T> P(3, 3);

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


	template<class T>
	linearAlgebra::Array<T> linearAlgebra::Array<T>::ones()
	{
		linearAlgebra::Array<T> tmp(nx,ny);
		for (int i = 0; i < nx; i++)
		{
			for (int j = 0; j < ny; j++)
			{
				tmp.x[i][j] = 1.0;
			}
		}
		return (tmp);
	}

	template<class T>
	linearAlgebra::Array<T> linearAlgebra::Array<T>::zeros()
	{
		linearAlgebra::Array<T>tmp(nx,ny);
		for (int i = 0; i < nx; i++)
		{
			for (int j = 0; j < ny; j++)
			{
				tmp.x[i][j] = 0.0;
			}
		}
		return (tmp);
	}

	template<class T>
	linearAlgebra::Array<T> linearAlgebra::Array<T>::eye()
	{
		Array<T>tmp(nx,ny);
		for (int i = 0; i < nx; i++)
		{
			T value{ 1.0 };
			tmp.setValue2(i, i, value);
		}
		return (tmp);
	}




	template<class T>
	linearAlgebra::Array<T>& linearAlgebra::Array<T>::operator*=(const linearAlgebra::Array<T>& tab)
	{
		assert(this->ny == tab.nx);

		for (int i = 0; i < this->nx; i++)
		{
			for (int j = 0; j < tab.ny; j++)
			{
				for (int k = 0; k < tab.nx; k++)
				{
					x[i][j] += x[i][k] * tab.x[k][j];
				}
			}
		}
		return(*this);
	}



	
	template<class T>
	linearAlgebra::Array<T> Array<T>::operator*(const Array<T>& tab)
	{
		//assert(tab.nx = ny);
		Array<T> prod(nx,tab.ny);
		for (int i = 0; i < nx; i++)
		{
			for (int j = 0; j < tab.ny; j++)
			{
				T sum(0.0);
				for (int k = 0; k < tab.nx; k++)
				{
					sum += x[i][k] * tab.x[k][j];
				}
				prod.x[i][j] = sum;
			}
		}
		return prod;
	}







	template<class T>
	linearAlgebra::Array<T> linearAlgebra::Array<T>::operator-(const linearAlgebra::Array<T>& tab)
	{
		assert(this->rowSize() == tab.nx);
		assert(this->columnSize()== tab.ny);
		linearAlgebra::Array<T>tmp(*this);
		for (int i = 0; i < tab.nx; i++)
		{
			for (int j = 0; j < tab.ny; j++)
			{
				tmp.x[i][j] = x[i][j] - tab.x[i][j];
			}
		}
		return tmp;
	}

	template<class T>
	linearAlgebra::Array<T>& linearAlgebra::Array<T>::operator-=(const linearAlgebra::Array<T>& tab)
	{
		assert(this->rowSize() == tab.nx);
		assert(this->columnSize() == tab.ny);
		for (int i = 0; i < tab.nx; i++)
		{
			for (int j = 0; j < tab.ny; j++)
			{
				x[i][j] -= tab.x[i][j];
			}
		}
		return (*this);
	}

	template<class T>
	linearAlgebra::Array<T> linearAlgebra::Array<T>::operator+(const linearAlgebra::Array<T>& tab)
	{
		assert(this->rowSize() == tab.nx);
		assert(this->columnSize() == tab.ny);
		linearAlgebra::Array<T>tmp(*this);
		for (int i = 0; i < tab.nx; i++)
		{
			for (int j = 0; j < tab.ny; j++)
			{
				tmp.x[i][j] = x[i][j] + tab.x[i][j];
			}
		}
		return tmp;
	}

	template<class T>
	linearAlgebra::Array<T>& linearAlgebra::Array<T>::operator+=(const linearAlgebra::Array<T>& tab)
	{
		assert(this->rowSize() == tab.nx);
		assert(this->columnSize() == tab.ny);
		for (int i = 0; i < tab.nx; i++)
		{
			for (int j = 0; j < tab.ny; j++)
			{
				x[i][j] += tab.x[i][j];
			}
		}
		return (*this);
	}

	template<class T>
	linearAlgebra::Array<T>& linearAlgebra::Array<T>::operator/=(const T& c)
	{
		assert(c > 0.0 || c < 0.0);

		for (int i = 0; i < this->rowSize(); i++)
		{
			for (int j = 0; j < this->columnSize(); j++)
			{
				x[i][j] /= c;
			}
		}

		return (*this);
	}

	template<class T>
	linearAlgebra::Array<T> linearAlgebra::Array<T>::operator/(const T& c)
	{
		assert((c > 0) || (c < 0));
		linearAlgebra::Array<T>tmp(*this);
		for (int i = 0; i < this->rowSize(); i++)
		{
			for (int j = 0; j < this->columnSize(); j++)
			{
				tmp.x[i][j] = x[i][j] / c;
			}
		}
		return tmp;
	}




    template<class T>
	T& linearAlgebra::Array<T>::operator()(int i, int j)
	{
		// Vérifications des bornes avec des assertions et exceptions
		if (i < 0 || i >= this->rowSize() || j < 0 || j >= this->columnSize()) {
			throw std::out_of_range("Indices hors limites dans Array::operator()");
		}
		return x[i][j];
	}






	template<class T>
	inline void linearAlgebra::Array<T>::allocArrays()
	{
		x = new T * [nx];
		for (int i = 0; i < nx; i++)
		{
			x[i] = new T[ny];
		}
	}

	template<class T>
	inline void linearAlgebra::Array<T>::setValue(int i, int j, T val)
	{
		// Vérifications des bornes avec des assertions et exceptions
		if (i < 0 || i >= this->rowSize() || j < 0 || j >= this->columnSize())
		{
			throw std::out_of_range("Indices hors limites dans Array::operator()");
		}
		x[i][j] = val;
	}

	template<class T>
	inline T linearAlgebra::Array<T>::getValue(int i, int j)
	{
		// Vérifications des bornes avec des assertions et exceptions
		if (i < 0 || i >= this->rowSize() || j < 0 || j >= this->columnSize())
		{
			throw std::out_of_range("Indices hors limites dans Array::operator()");
		}
		return x[i][j];
	}

	template<class T>
	inline T linearAlgebra::Array<T>::getValue2(int& i, int& j)
	{
		// Vérifications des bornes avec des assertions et exceptions
		if (i < 0 || i >= this->rowSize() || j < 0 || j >= this->columnSize())
		{
			throw std::out_of_range("Indices hors limites dans Array::operator()");
		}
		return x[i][j];
	}



	template<class T>
	inline void linearAlgebra::Array<T>::printing(const char* message)
	{
		std::cout << "\n";
		std::cout << "\t" << message << "\n";
		for (int i = 0; i < this->rowSize(); i++)
		{
			for (int j = 0; j < this->columnSize(); j++)
			{
				std::cout << "\t" << x[i][j] << " ";
			}
			std::cout << "\n";
		}
		std::cout << "\n";
	}

	template<class T>
	linearAlgebra::Array<T>::Array<T>(int nnx, int nny) :nx(nnx), ny(nny)
	{
		allocArrays();
		for (int i = 0; i < this->rowSize(); i++)
		{
			for (int j = 0; j < this->columnSize(); j++)
			{
				x[i][j] = 0.0;
			}
		}
	}

	template<class T>
	linearAlgebra::Array<T>::Array<T>(const linearAlgebra::Array<T>& tab) :nx(tab.nx), ny(tab.ny)
	{
		allocArrays();
		for (int i = 0; i < nx; i++)
		{
			for (int j = 0; j < ny; j++)
			{
				x[i][j] = tab.x[i][j];
			}
		}

	}


	template<class T>
	linearAlgebra::Array<T>::~Array()
	{
		for (int i = 0; i < this->rowSize(); i++)
		{
			delete [] x[i];
		}
		delete[]x;
	}


	template<class T>
	linearAlgebra::Array<T>& linearAlgebra::Array<T>::operator=(const linearAlgebra::Array<T>& tab)
	{
		if (this == &tab)
		{
			return(*this);
		}
		else
		{
			if (nx != tab.nx || ny != tab.ny)
			{
				this->~Array();
				nx = tab.nx;
				ny = tab.ny;
				allocArrays();
			}
			for (int i = 0; i < nx; i++)
			{
				for (int j = 0; j < ny; j++)
				{
					x[i][j] = tab.x[i][j];
				}
			}
			return (*this);
		}

	}











}//end namespace linearAlgebra







#endif //__ARRAY_H__


