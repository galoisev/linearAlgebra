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
			class NumAnalysis
			{
			public:
				virtual void solve() = 0;
			}; //end class NumAnalysis

			template<typename T>
			class Pivot:public NumAnalysis<T>
			{
			private:
			public:
				Matrix<T> A;
				Matrix<T> b, x;
				Pivot(Matrix<T>& AA, Matrix<T>& bb) :A(AA), b(bb)
				{
					std::cout << "\n\tLa solution x du système A*x = b est: " << "\n";
					solve();					
				}//end pivot
				Pivot& operator=(const Pivot& pivot) = default;
				virtual void solve()
				{
					int m = A.getRows();
					int n = A.getCols();
					Matrix<T>x(n, 1, 0);

					Matrix<T> mat_A1(m, n, 0);//init
					mat_A1 = A;
					Matrix<T> vec_b1(n, 1, 0);//init
					vec_b1 = b;

					Matrix<T> mat_A2(m, n, 0);
					Matrix<T> vec_b2(m, 1, 0);
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
							dble val_b1 = vec_b1.getValue(i, 0);
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
							sum_j += mat_A1.getValue(i, j) * x.getValue(j,0);
						}
						dble val = (vec_b1.getValue(i,0) - sum_j) / mat_A1.getValue(i, i);
						x.setValue(i, 0, val);
					}
					std::cout << x << "\n";
				}
			};


			class LU
			{
			private:
			public:
			};//end class LU

		}//end method directe





		namespace project
		{
			template<typename T>
			class NumAnalysis
			{
			public:
				virtual void solve() = 0;
			}; //end class NumAnalysis
			template<typename T>
			class GC :public NumAnalysis<T>
			{
			private:
			public:
				Matrix<T> A;
				Matrix<T> b, x;
				GC(Matrix<T>& AA, Matrix<T>& bb) :A(AA), b(bb)
				{
					std::cout << "\n\tLa solution x du système A*x = b est: " << "\n";
					solve();
				}
			};

		}//end methode projective




		namespace iter
		{

		}

	}//end namespace systLin






	namespace solveEquations
	{

		
		dble g(double& x)
		{
			return x * x * x - x - 2;//à modifier selon...
		}


		/* les fonctions */
		double f1(dble& x, dble& y)
		{
			return (x - y * y + x * exp(y) - 2.0);
		}
		double f2(dble& x, dble& y)
		{
			return (y * exp(y) + x * x * x - 1.0);
		}

		/* leurs derivees */
		double a22(dble& x, dble& y)//f11 = f1,x
		{
			return (1 + exp(y));
		}
		double a12(dble& x, dble& y)//f12 = f1,y
		{
			return (2 * y - x * exp(y));
		}
		double a21(dble& x, dble& y)//f21 = f2,x
		{
			return (-3 * x * x);
		}
		double a11(dble& x, dble& y)//f22 = f2,y
		{
			return ((y + 1) * exp(y));
		}



		template<typename T>
		class NumAnalysis
		{
		public:
			virtual void solve() = 0;
		};//end class NumAnalysis


		template<typename T>
		class NewtonRaphson :public NumAnalysis<T>
		{
		private:
		public:
			double tolerance{ 1.e-9 };
			int max_iter{ 100 };
			Matrix<T> A;
			Matrix<T> b, x;
			NewtonRaphson(Matrix<T>& AA, Matrix<T>& bb, Matrix<T>& xx) :A(AA), b(bb),x(xx)
			{
				std::cout << "\n\tLa solution approchee x de l equation f(x) = 0  est: " << "\n";
				solve();
			}
			NewtonRaphson(const NewtonRaphson& nr) :A(nr.A), b(nr.b),x(nr.x) {}
			virtual void solve()
			{
				/*
				* A: derivee ou matrice jacobienne
				* delta: determinant(A)
				*/				
				Matrix<T> A1(2, 2, 0), J(2, 2, 0);//matrice
				Matrix<T> X1(2, 1, 0), X2(2, 1, 0), X(2, 1, 0), B(2, 1, 0), C(2, 1, 0);//vecteurs
				A1 = A; B = b; X1 = x;

				dble xi = X1.getValue(0, 0);
				dble yi = X1.getValue(1, 0);			

				J = A1.inv() / (1. / A1.det());				
				dble DELTA = J.det();
				DELTA = DELTA * (-1);

				//init
				T ecart{1.e12};
				int iter = { 0 };	
				do
				{
					X2 = X1 - A1 * B / DELTA;

					if (iter > this->max_iter)
					{
						std::cout << "\tnombre iteration max est depasse !" << "\n";
						break;
					}

					xi = X2.getValue(0, 0);
					yi = X2.getValue(1, 0);
					B.setValue(0, 0, f1(xi, yi));
					B.setValue(1, 0, f2(xi, yi));
					A1.setValue(0, 0, a11(xi, yi)); A1.setValue(0, 1, a12(xi, yi));
					A1.setValue(1, 0, a21(xi, yi)); A1.setValue(1, 1, a22(xi, yi));
					J = A1.inv() / (1. / A1.det());
					DELTA = J.det();
					DELTA = DELTA * (-1);

					X = X2 - X1;
					ecart = X.normL2();
					
					X1 = X2;
					X2 = X2.zeros();
					iter++;

				} while ((ecart > this->tolerance) && (iter < max_iter));
				
				X1.print("X = ");
			}//end solve
			
		};//end class NewtonRaphson




		template<typename T>
		class Dichotomie:public NumAnalysis<T>
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

		};//end class Dichotomie





	}//end namespace solveEquations

}//end namespace methods






#endif // !_NUMERICAL_METHODS_H
