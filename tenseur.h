#ifndef _TENSOR_H
#define _TENSOR_H



#include <iostream>
#include <vector>
#include <stdexcept>
#include <cassert>


namespace LINALG
{



	template<typename T>
	class Tensor
	{
	public:
		virtual void allocationMemory() = 0;
		virtual ~Tensor() {} //destructeur vituel nécessaire

	};//end class cirtual Tensor



	template<typename T>
	class Vector :public Tensor<T>
	{
	private:
		size_t size;
		T* data;
	public:
		Vector() :size(0), data(nullptr) {}
		Vector(const T& value, size_t n) :data(new T[n]), size(n)
		{
			for (size_t i = 0; i < size; i++)
			{
				data[i] = value;
			}
		}
		~Vector()
		{
			delete[]data;
		}
		Vector(const Vector<T>& v) :size(v.size), data(new T[v.size])
		{
			for (size_t i = 0; i < size; i++)
			{
				data[i] = v.data[i];
			}
		}
		Vector& operator=(const Vector<T>& v)
		{
			if (this != &v)
			{
				delete[]data;//liberation ancienne memoire
				size = v.size;
				data = new T[v.size];
				for (size_t i = 0; i < size; i++)
				{
					data[i] = v.data[i];
				}
			}
			return (*this);				
		}

		T& operator[](size_t i)
		{
			assert(i > 0 && i < size);
			return data[i];
		}


		Vector operator+(const Vector<T>& v);
		Vector operator-(const Vector<T>& v);
		Vector operator*(const Vector<T>& v);
		Vector operator^(const Vector<T>& v);//cross product


		Vector operator*()const { return *data; }
		Vector* operator->() const { return *data; }

		auto operator<=>(const Vector<T>& vec)const = default;
		

		inline void setValue(int i, T val)
		{
			assert(i >= 0 && i < this->getSize());
			this->data[i] = val;

		}




		virtual void allocationMemory()const
		{
			if (size > 0)
			{
				data = new T * [size];
			}			
		}

		size_t getSize()
		{
			return size;
		}


	};//end class Vector





	template<typename T>
	Vector<T> Vector<T>::operator^(const Vector<T>& v)
	{
		assert(v.size == this->getSize() == 3);
		Vector<T>temp;
		temp.data[0] = this->data[1] * v.data[2] - this->data[2] * v.data[1];
		temp.data[1] = this->data[2] * v.data[0] - this->data[0] * v.data[2];
		temp.data[0] = this->data[0] * v.data[1] - this->data[1] * v.data[0];
		return temp;
	}


	/*************** END CLASS VECTOR ******************/


	template<typename T>
	class Matrix :public Tensor<T>
	{
	private:
		int rowSize, columnSize;
		T** x;
	public:
		Matrix();
		Matrix(const Matrix<T>& mat);
		~Matrix()
		{
			for (int i = 0; i < rowSize; i++)
			{
				delete[]x[i];
			}
			delete[]x;
		}

		virtual void allocationMemory()const
		{
			x = new T * [rowSize];
			for (int i = 0; i < rowSize; i++)
			{
				x[i] = new T[columnSize];
			}
		}


		auto operator<=>(const Matrix<T>& mat)const = default;
	};


	template<typename T>
	Matrix<T>::Matrix(const Matrix<T>& mat):rowSize(mat.rowSize),columnSize(mat.columnSize)
	{
		allocationMemory();
		for (int i = 0; i < rowSize; i++)
		{
			for (int j = 0; j < columnSize; j++)
			{
				x[i][j] = mat.x[i][j];
			}
		}
	}








}//end namespace LINALG



#endif // !_TENSOR_H
