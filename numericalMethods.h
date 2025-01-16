#ifndef _NUMERICAL_METHODS_H


#include "array_class.h"
#include <typeinfo>


namespace methods
{

	namespace systLin
	{

		namespace direct
		{

			template<typename T>
			class numAnalysis
			{
			public:
				virtual void solve() = 0;
			};

			template<typename T>
			class Pivot:public numAnalysis<T>
			{
			private:
			public:
				Matrix<T> A;
				Vector<T> b, x;
				Pivot(Matrix<T>& AA, Vector<T>& bb) :A(AA), b(bb)
				{
					std::cout << "\n\tLa solution x du système A*x = b est: " << "\n";
					solve();					
				}//end pivot
				Pivot& operator=(const Pivot& pivot) = default;
				virtual void solve()
				{
					int m = A.getRows();
					int n = A.getCols();
					assert(m == b.getRows());
					Vector<T>x(n, 0);

					Matrix<T> mat_A1(m, n, 0);//init
					mat_A1 = A;
					Vector<T> vec_b1(n, 0);//init
					vec_b1 = b;

					Matrix<T> mat_A2(m, n, 0);
					Vector<T> vec_b2(m, 0);
					for (size_t k = 0; k < m; k++)
					{
						for (size_t i = k + 1; i < m; i++)
						{
							for (size_t j = k + 1; j < n; j++)
							{
								dble val_A2 = mat_A1.getValue(i, j) - (mat_A1.getValue(i, k) / mat_A1.getValue(k, k)) * mat_A1.getValue(k, j);
								mat_A2.setValue(i, j, val_A2);
							}
							mat_A2.setValue(i, k, 0);

							dble val_b2 = vec_b1.getValue(i) - (mat_A1.getValue(i, k) / mat_A1.getValue(k, k)) * vec_b1.getValue(k);
							vec_b2.setValue(i, val_b2);
						}

						//les valeurs de la matrice A(k+1) sont les memes que celles de la matrice A(k) 
						//les valeurs du vecteur    b(k+1) sont les memes que celles du vecteur    b(k)
						for (size_t i = 0; i < k + 1; i++)
						{
							for (size_t j = k; j < n; j++)
							{
								dble val_A1 = mat_A1.getValue(i, j);
								mat_A2.setValue(i, j, val_A1);
							}
							dble val_b1 = vec_b1.getValue(i);
							vec_b2.setValue(i, val_b1);
						}

						/*MàJ*/
						mat_A1 = mat_A2;
						vec_b1 = vec_b2;
					}

					//resolution du systeme en remontant les equations...
					for (int i = m - 1; i >= 0; i--)
					{
						dble sum_j = { 0 };
						for (int j = n - 1; j > i; j--)
						{
							sum_j += mat_A1.getValue(i, j) * x.getValue(j);
						}
						dble val = (vec_b1.getValue(i) - sum_j) / mat_A1.getValue(i, i);
						x.setValue(i, val);
					}
					std::cout << x << "\n";
				}
			};


			class LU
			{
			private:
			public:
			};


		}
		namespace project
		{

		}
		namespace iter
		{

		}

	}






	namespace solveEquations
	{

		
		double g(double& x)
		{
			return x * x * x - x - 2;//à modifier selon...
		}


		class numAnalysis
		{
		public:
			virtual void solve() = 0;
		};


		class Dichotomie:public numAnalysis
		{
		private:
		public:
			double a{}, b{}, tolerance{1.e-3};
			int max_iter{100};
			Dichotomie(double aa, double bb, double eps, int maxiter) :a(aa), b(bb), tolerance(eps), max_iter(maxiter)
			{
				std::cout << "\n\tLa solution approchee x de l equation g(x) = 0  est: " << "\n";
				solve();
			}//end dichotomie

			Dichotomie(const Dichotomie& dico) :a(dico.a), b(dico.b), tolerance(dico.tolerance), max_iter(dico.max_iter) {}
			Dichotomie& operator=(const Dichotomie& dico) = default;


			virtual void solve()
			{
				if (solveEquations::g(this->a) * solveEquations::g(this->b) >= 0.0)
				{
					std::cerr << "\tLa fonction doit changer de signe sur [a, b]\n";
					std::system("pause");
					std::exit(1);
				}
				else
				{
					int iter = { 0 };
					double c{};
					do 
					{
						c = (this->a + this->b) / 2.0;
						if ((std::abs(solveEquations::g(c)) < tolerance) || (std::abs(this->b - this->a) < tolerance))
						{
							std::cout << "\t" << c << "\n";
							break;
						}
						if (solveEquations::g(this->a) * solveEquations::g(c) < 0.0)
						{
							this->b = c; //la racine se trouve dans [a, c]
						}
						else
						{
							this->a = c; //la racine se trouve dans [c, b]
						}
						if (iter > this->max_iter)
						{
							std::cout << "\tLa methode ne CV pas avec ce nombre d'iteration max !\n\n";
							std::system("pause");
							std::exit(1);
						}
	
						iter++;						
					} while (iter < this->max_iter);//end while					
				}//end else
			}//end algo

		};



		class NewtonRaphson
		{
		private:
		public:
		};

	}



}




Matrix<dble> solve(Matrix<dble>& A, Matrix<dble>& b)
{
	int m = A.getRows();
	int n = A.getCols();

	Matrix<dble> mat_A1(m, n, 0);//init
	mat_A1 = A;
	Matrix<dble> vec_b1(m, 1, 0);//init
	vec_b1 = b;

	Matrix<dble>x(n, 1, 0);

	Matrix<dble> mat_A2(m, n, 0);
	Matrix<dble> vec_b2(m, 1, 0);
	for (size_t k = 0; k < m ; k++)
	{
		for (size_t i = k + 1; i < m; i++)
		{	
			for (size_t j = k + 1; j < n; j++)
			{
				dble val_A2 = mat_A1.getValue(i, j) - (mat_A1.getValue(i, k) / mat_A1.getValue(k, k)) * mat_A1.getValue(k, j);
				mat_A2.setValue(i, j, val_A2);				
			}
			mat_A2.setValue(i, k, 0);

			dble val_b2 = vec_b1.getValue(i, 0) - (mat_A1.getValue(i, k) / mat_A1.getValue(k, k)) * vec_b1.getValue(k, 0);
			vec_b2.setValue(i, 0, val_b2);
		}

		//les valeurs de la matrice A(k+1) sont les memes que celles de la matrice A(k) 
		//les valeurs du vecteur    b(k+1) sont les memes que celles du vecteur    b(k)
		for (size_t i = 0; i < k + 1; i++)
		{
			for (size_t j = k; j < n; j++)
			{
				dble val_A1 = mat_A1.getValue(i, j);
				mat_A2.setValue(i, j, val_A1);
			}
			dble val_b1= vec_b1.getValue(i, 0);
			vec_b2.setValue(i, 0, val_b1);
		}

		/*MàJ*/
		mat_A1 = mat_A2;
		vec_b1 = vec_b2;
	}

	//resolution du systeme en remontant les equations...
	for (int i = m - 1; i >= 0; i--)
	{
		dble sum_j = { 0 };
		for (int j = n - 1; j > i; j--)
		{
			sum_j += mat_A1.getValue(i, j) * x.getValue(j, 0);
		}
		dble val = (vec_b1.getValue(i, 0) - sum_j) / mat_A1.getValue(i, i);
		x.setValue(i, 0, val);
	}

	return x;
}







namespace LINALG
{

	template<typename T>
	class SystLin // class linear system
	{
	private:
		std::string method;
	public:			
		SystLin() :method(" ") {}
		SystLin(std::string& new_method) :method(new_method) {}
		virtual std::string getMethod()const { return method; }
		//virtual Matrix<T> solve(Matrix<T>& A, Matrix<T>& b)

	};






	template<typename T>
	class MethodDirect:public SystLin<T>
	{
	private:
		int _rows, _cols;
		Matrix<T>A;
	public:
		MethodDirect() {};
		Vector<T> pivotDeGauss(Matrix<T>&A, Vector<T>& b);
		Vector<T> crout();//TODO
		Vector<T> cholesky();//TODO
	};

	template<typename T>
	Vector<T> MethodDirect<T>::pivotDeGauss(Matrix<T>& A, Vector<T>& b)
	{
		int m = A.getRows();
		int n = A.getCols();
		assert(m == b.getRows());
		Vector<T>x(n, 0);


		Matrix<T> mat_A1(m, n, 0);//init
		mat_A1 = A;
		Vector<T> vec_b1(n, 0);//init
		vec_b1 = b;

		
		Matrix<T> mat_A2(m, n, 0);
		Vector<T> vec_b2(m, 0);
		for (size_t k = 0; k < m; k++)
		{
			for (size_t i = k + 1; i < m; i++)
			{
				for (size_t j = k + 1; j < n; j++)
				{
					dble val_A2 = mat_A1.getValue(i, j) - (mat_A1.getValue(i, k) / mat_A1.getValue(k, k)) * mat_A1.getValue(k, j);
					mat_A2.setValue(i, j, val_A2);
				}
				mat_A2.setValue(i, k, 0);

				dble val_b2 = vec_b1.getValue(i) - (mat_A1.getValue(i, k) / mat_A1.getValue(k, k)) * vec_b1.getValue(k);
				vec_b2.setValue(i, val_b2);
			}

			//les valeurs de la matrice A(k+1) sont les memes que celles de la matrice A(k) 
			//les valeurs du vecteur    b(k+1) sont les memes que celles du vecteur    b(k)
			for (size_t i = 0; i < k + 1; i++)
			{
				for (size_t j = k; j < n; j++)
				{
					dble val_A1 = mat_A1.getValue(i, j);
					mat_A2.setValue(i, j, val_A1);
				}
				dble val_b1 = vec_b1.getValue(i);
				vec_b2.setValue(i, val_b1);
			}

			/*MàJ*/
			mat_A1 = mat_A2;
			vec_b1 = vec_b2;
		}

		//resolution du systeme en remontant les equations...
		for (int i = m - 1; i >= 0; i--)
		{
			dble sum_j = { 0 };
			for (int j = n - 1; j > i; j--)
			{
				sum_j += mat_A1.getValue(i, j) * x.getValue(j);
			}
			dble val = (vec_b1.getValue(i) - sum_j) / mat_A1.getValue(i, i);
			x.setValue(i, val);
		}
		return x;
	}


	template<typename T>
	class MethodIterativ :public SystLin<T>
	{
	private:
	public:
		MethodIterativ() {};
		void relaxation();//TODO

	};


	template<typename T>
	class MethodProjectiv :public SystLin<T>
	{
	private:
	public:
		MethodProjectiv() {};
		void gradientConjugue();//TODO

	};





}//end namespace LINALG




#endif // !_NUMERICAL_METHODS_H
