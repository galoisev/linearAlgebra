#ifndef _NUMERICAL_METHODS_H


#include "array_class.h"
#include <typeinfo>
#include <map>


namespace methods
{



	namespace eigenValues
	{
		template<typename T>
		class SubMethods
		{
		public:
			virtual void solve() = 0;
		};




		template<typename T>
		class PuissancesIterees:public SubMethods<T>
		{
		private:
			double tolerance = { 1e-9 };
		public:
			Matrix<T> A;
			Matrix<T> x;

			PuissancesIterees(Matrix<T>& AA, Matrix<T> xx) :A(AA),x(xx)
			{
				std::cout << "\n\tRecherche des valeurs propres de la matrice: " << "\n";
				solve();
			}	


			double getTolerance()const
			{
				return tolerance;
			}


			virtual void solve()
			{
				int m = A.getRows();
				int n = A.getCols();
				assert(m == n);//A est une matrice carrée !

				//1 
				Matrix<T>V(n, 1, 1);//initialisation du vecteur V1.
				V = x;
				std::cout << "\n\n\tVecteur arbitraire: \n";
				std::cout << V << "\n\n";
				//2
				Matrix<T>U1(m, 1, 0), W1(m, 1, 0);

				//recherche de la composante le plus grande du vecteur U1.
				double valMax{ 1e-9 };				
				double k1 = { 1.e12 }, k2 = { 0. };
				

				double delta_k = std::abs(k2 - k1);
				int iter = 0;
				while (delta_k > this->getTolerance())
				{
					U1 = A * V;
					valMax = 1e-9;
					for (int i = 0; i < m; i++)
					{
						if ( std::abs( U1.getValue(i, 0) ) > valMax )
						{
							valMax = U1.getValue(i, 0);
						}
					}
					k1 = valMax;

					//normalisation du vecteur U1.
					W1 = U1 / k1;

					//std::cout << W1 << "\n\n";
					delta_k = std::abs(k2 - k1);
					if (iter > 1000)
					{
						std::cerr << "\tPAS de CV: nombre iteration max depassee !\n";
						std::system("pause");
						std::exit(1);
					}						

					//MàJ
					iter++;
					k2 = k1;
					k1 = 1.e12;
					V = W1;

				}//end while


				/*
				template<class T, class U>
				void results2(const std::map<T, U>&mresults, std::string & filename)
				{
					std::ofstream fic(filename);
					for (const auto& [key, value] : mresults)
					{
						fic << key << ";" << value << "\n";
					}
					fic.close();
					std::cout << "\tWriting (summary of view factor calculations) to the file is COMPLETE !   " << filename << std::endl;
				}

				void print_map(std::string comment, const std::map<T, U>& m)
				{
					std::cout << comment << "\n";
					// Iterate using C++17 facilities
					for (const auto& [key, value] : m)
						std::cout << '[' << key << "] = " << value << "; ";
				}*/
				// C++11 alternative:
//  for (const auto& n : m)
//      std::cout << n.first << " = " << n.second << "; ";
//
// C++98 alternative:
//  for (std::map<std::string, int>::const_iterator it = m.begin(); it != m.end(); ++it)
//      std::cout << it->first << " = " << it->second << "; ";

				
				double eigenValue = k2;
				std::cout << "\n\n\titer : " << iter << "\n";
				std::cout << "\tvaleur propre: " << eigenValue << "\n";
				std::cout << "\tvecteur propre associe:\n";
				std::cout << W1 << "\n\n";

				
				std::map<double, Matrix<double>> m_dmat = {};
				std::map<int, std::map<double, Matrix<double>>> mm_dmat = {};
				m_dmat[eigenValue] = W1;
				mm_dmat[1] = m_dmat;

				for (auto& [num, m] : mm_dmat)
				{
					std::cout << "\t" << num << "\n";
					for (auto& [x, w] : m)
					{
						std::cout << "\tvalueur propre: " << x << ", vecteur propre: " << w << "\n";
					}
					std::cout << "\n";
					//std::cout << key << "," << val << "\n";
				}
				//std::map<int, std::map<double,Matrix<double>>> map_map_vecPropre = {};

                /*
				map_map_vecPropre[1] = W1;
				for (auto& [key, m] : mm)
				{
					for (auto& [val_prop, vec_prop] : map_val_vec)
					{
						std::cout <<  "\t" << key << " : " << ", val. propre: " << val_prop << ", vect. propre associé: " << vec_prop << "\n";
					}
					
				}
				*/
				/*
				for (const auto& [id, x] : map_vecPropre)
				{
					std::cout << "\valeurPropre n°: " << id << " valeur_propre: \n";
					std::cout << x << "\n";
				}*/
				

			}




		};

	}//end namespace eigenValues







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
					std::cout << "\n\tLa solution x du systeme A*x = b est: " << "\n";
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
					for (int k = 0; k < m; k++)
					{
						for (int i = k + 1; i < m; i++)
						{
							for (int j = k + 1; j < n; j++)
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
						for (int i = 0; i < k + 1; i++)
						{
							for (int j = k; j < n; j++)
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





			template<typename T>
			class LU
			{
			private:
			public:
				Matrix<T> A;
				Matrix<T> b, x;
				LU(Matrix<T>& AA, Matrix<T>& bb, Matrix<T>& xx) : A(AA), b(bb), x(xx)
				{
					std::cout << "\n\tLa solution x du systeme A*x = b est: " << "\n";
					solve2();
				}//end LU
				LU& operator=(const LU<T>& lu) = default;
				void solve()
				{

					int m = A.getRows();
					int n = A.getCols();
					assert(m == n);//matrice carree
					Matrix<T>x(n, 1, 0);

					Matrix<T> L(m, n, 0), U(m, n, 0);//init

					//initialisation
					L = L.eye();
					L.print("L=");
					U = A;
					U.print("U=");

					//factorisation
					for (int i = 0; i < m; i++)
					{
						for (int j = 0; j < i; j++)
						{
							double sum_k = 0;
							for (int k = 0; k < i; k++)
							{
								sum_k += L.getValue(i, k) * U.getValue(k, j);
							}
							double val_L = (A.getValue(i, j) - sum_k) / U.getValue(i, i);
							L.setValue(i, j, val_L);
						}
					}
					L.print("L=");

					//màj de U
					for (int i = 0; i < m; i++)
					{
						for (int j = 0; j > i-1; j++)
						{
							double sum_k = 0;
							for (int k = 0; k < i; k++)
							{
								sum_k += L.getValue(i, k) * U.getValue(k, j);
							}
							double val_U = A.getValue(i, j) - sum_k;
							U.setValue(i, j, val_U);
						}
					}
					U.print("U=");



					x.print("x=");

				}




				void solve2()
				{
					int m = A.getRows();
					int n = A.getCols();
					assert(m == n);//matrice carree
					Matrix<T>x(n, 1, 0);

					Matrix<T> L(m, n, 0), U(m, n, 0);//init

					for (int j = 0; j < n; j++)
					{
						double val = A.getValue(0, j);
						U.setValue(0, j, val);
					}

					for (int i = 0; i < m; i++)
					{
						double val = A.getValue(i, 0) / U.getValue(0, 0);
						L.setValue(i, 0, val);
					}

					for (int i = 0; i < m; i++)
					{
						L.setValue(i, i, 1);
					}

					for (int i = 0; i < m; i++)
					{
						for (int j = i + 1; j < n; j++)
						{
							L.setValue(i, j, 0);
						}
					}
					L.print("L=");

					for (int i = 0; i < m; i++)
					{
						for (int j = 0; j < i; j++)
						{
							U.setValue(i, j, 0);
						}
					}

					for (int i = 0; i < m; i++)
					{
						for (int j = 0; j < n; j++)
						{
							double sum = 0;
							for (int k = 0; k < i; k++)
							{
								sum += L.getValue(i, k) * U.getValue(k, j);
							}
							double val_U = (A.getValue(i, j) - sum);
							U.setValue(i, j, val_U);
						}
					}
					
					for (int i = 0; i < m; i++)
					{
						for (int j = i; j < n; j++)
						{
							double sum = 0;
							for (int k = 0; k < j; k++)
							{
								sum += L.getValue(i, k) * U.getValue(k, j);
							}
							double val_L = (A.getValue(i, j) - sum) / U.getValue(j, j);
							L.setValue(i, j, val_L);
						}
					}
					
					// on "remonte" les x 
					Matrix<T>y = b;		
					Matrix<T> z(n, 1, 0);


					double val_z0 = y.getValue(0, 0) / L.getValue(0, 0);
					z.setValue(0, 0, val_z0);
					z.setValue(0, 0, y.getValue(0, 0));

					double val_L = 0, val_U = 0, som_j;
					for (int i = 1; i < m; i++)
					{
						som_j = 0;
						for (int j = 0; j < i; j++)
						{
							som_j += L.getValue(i, j) * z.getValue(j, 0);
						}
						val_L = y.getValue(i, 0) - som_j;
						z.setValue(i, 0, val_L);
					}

					double val_xn = z.getValue(n - 1, 0) / U.getValue(n - 1, n - 1);
					x.setValue(n - 1, 0, val_xn);
					
					for (int p = n - 2; p >= 0; p--)
					{
						double sum_q = 0;
						for (int q = p + 1; q < n; q++)
						{
							sum_q += U.getValue(p, q) * x.getValue(q, 0);
						}
						double val = (z.getValue(p, 0) - sum_q) / U.getValue(p, p);
						x.setValue(p, 0, val);
					}					

					x.print("x=");
					
				}


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
			class GradientConjugue :public NumAnalysis<T>
			{
			private:
			public:
				double tolerance{ 1.e-9 };
				int max_iter{ 100 };
				Matrix<T> A;
				Matrix<T> b, x;
				GradientConjugue(Matrix<T>& AA, Matrix<T>& bb, Matrix<T> xx) :A(AA), b(bb), x(xx)
				{					
					std::cout << "\n\tLa solution x du systeme A*x = b est: " << "\n";
					solve();
				}
				virtual void solve()
				{
					Matrix<T>x0(2, 1, 0), x1(2, 1, 0), r0(2, 1, 0), r1(2, 1, 0), p0(2, 1, 0), p1(2, 1, 0), q0(2, 1, 0), q1(2, 1, 0);
					T alpha0{ 0.0 }, alpha1{ 0.0 }, beta0{ 0.0 }, beta1{ 0.0 }, gamma0{ 0.0 }, mu0{ 0.0 }, gamma1{ 0.0 }, mu1{ 0.0 };
					int m = A.getRows();
					int n = A.getCols();
					Matrix<T>x(n, 1, 0);
					if (!A.isSymetric())
					{
						std::cerr << "\tThe matrix A must be both positive definite and symmetric !" << "\n";
						std::system("pause");
						std::exit(1);
					}

					x = this->x;
					A = this->A;
					b = this->b;


					x0 = x;					
					r0 = b - A * x0;
					p0 = r0;

					T ecart{ 1.e12 };
					int iter = { 0 };
					do
					{
						gamma0 = r0.operator,(r0);
						q0 = A * p0;

						alpha1 = (r0, r0) * (1. / (p0, q0));
						x1 = x0 + p0 * alpha1;

						r1 = r0 - (A * p0) * alpha1;
						beta1 = (r1, r1) * (1. / (r0, r0));

						p1 = r1 + p0 * beta1;
						if (iter > this->max_iter)
						{
							std::cout << "\tnombre iteration max est depasse !" << "\n";
							break;
						}

						ecart = r0.normL2();

						/*MaJ*/						
						x0 = x1;
						r0 = b - A * x0;
						p0 = r0;
						iter++;

					} while ((ecart > this->tolerance) && (iter < max_iter));
					
					x = x0;
					x.print("x=");
				}
			};

		}//end methode projective




		namespace iter
		{
			template<typename T>
			class NumAnalysis
			{
			public:
				virtual void solve() = 0;
			}; //end class NumAnalysis

			template<typename T>
			class Jacobi :public NumAnalysis<T>
			{
			private:
			public:
				double tolerance{ 1.e-9 };
				int max_iter{ 100 }, iter;
				Matrix<T> A;
				Matrix<T> b, x, x1, x2, x21;
				Jacobi(Matrix<T>& AA, Matrix<T>& bb, Matrix<T> xx) :A(AA), b(bb), x(xx)
				{
					/*
					for (int i = 0; i < A.getRows(); i++)
					{
						if (!A.isDiagonalDominant(i))
						{
							std::cerr << "\tThe matrix A must be diagonal dominant !" << "\n";
							std::system("pause");
							std::exit(1);
						}
					}*/
					std::cout << "\n\tLa solution x du systeme A*x = b est: " << "\n";
					solve();
				}
				virtual void solve()
				{					
					int m = A.getRows();
					int n = A.getCols();
					Matrix<T> x1(n,1,0), x2(n,1,0), x21(n,1,0);

					x1 = this->x;
					double delta = 1e9;
					iter = 0;
					while (delta > this->tolerance)
					{

						for (int i = 0; i < m; i++)
						{
							double sum_j = { 0 };
							for (int j = 0; j < n; j++)
							{
								if (j != i)
								{
									sum_j += this->A.getValue(i, j) * x1.getValue(j, 0);
								}								
							}
							double val = (this->b.getValue(i, 0) - sum_j) / this->A.getValue(i, i);
							x2.setValue(i, 0, val);
						}

						if (iter > this->max_iter)
						{
							std::cout << "\tnombre iteration max est depasse !" << "\n";
							break;
						}
						x21 = x2 - x1;
						/*MàJ*/
						delta = x21.normL2();
						x1 = x2;
						iter++;	

					}//end while

					//std::cout << iter << "\n";
					x1.print("x=");

				}//end solve
			};//end class Jacobi

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
