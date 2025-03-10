

#if _WIN32 || _WIN64
# define WINDOWS_MODE
#endif
#ifdef WINDOWS_MODE
//#include <windows.h>
#include <GeometryDll.cpp>
#define PREFIX __cdecl
#else
#define PREFIX __attribute__((__cdecl__))
#include <dlfcn.h>
#include <GeometryDll.h>
#endif


/*
#include "../nodes.h"
#include "../elements.h"
#include "../mesh.h"
#include "../geo.h"
#include "../geo_old.h"
#include "../viewfactor.h"
#include "../array.h"
#include "../entities.h"
#include "../MatrixRotation.h"
#include "../meshFormat.h"
#include "../quantityOfHeatByRadiation.h"
*/

#include "../calculsFacteursDeVue.h"
#include "../array.h"

#include "../array_class.h"
#include "../numericalMethods.h"

#include "tenseur.h"




#include <complex>
#include <vector>
#include <functional>
#include <map>


#include <ctime>
#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include <windows.h>

#include <algorithm> // Pour std::find

#include <thread>
#include <Windows.h>





void showProgressBar(int progress, int total, std::string message1) {
	int barWidth = 50;
	float percentage = static_cast<float>(progress) / total;
	int pos = (int)(barWidth * percentage);
	std::cout << message1 << "\n\n";
	std::cout << "\t\t\t[";
	for (int i = 0; i < barWidth; ++i) {
		if (i < pos)
			std::cout << "=";
		else if (i == pos)
			std::cout << ">";
		else
			std::cout << " ";
	}
	std::cout << "] " << int(percentage * 100) << " %\r";
	std::cout.flush();
}










void PrintSurfaceValues(double* surfaces)
{
	for (int i = 0; i < 10; i++)
	{
		std::cout << to_string(surfaces[i]) << std::endl;
	}
}




template<class T>
T matProduct(std::vector<T>& vH)
{
	T prod = vH.operator[](0);
	for (int i = 1; i < vH.size(); ++i)
	{
		prod = prod * vH.operator[](i);
		//prod.operator*=( v );
	}
	return prod;
}


template<class T>
void results(std::vector<T>& vresults, std::string& filename)
{
	//Sauvegarde du fichier
	std::ofstream fic(filename);
	for (auto&& x : vresults)
	{
		fic << "\tfactor of view = " << x << "\n";
	}
	fic.close();
	std::cout << "\t\tWriting (summary of view factor calculations) to the file is COMPLETE ! \n\n\n\n" << filename << std::endl;
}

template<class T, class U>
void results2(const std::map<T, U>& mresults, std::string& filename)
{
	std::ofstream fic(filename);
	for (const auto& [key, value] : mresults)
	{
		fic << key << ";" << value << "\n";
	}
	fic.close();
	std::cout << "\tWriting (summary of view factor calculations) to the file is COMPLETE !   " << filename << std::endl;
}


template<class T>
void printFile(T& data, std::string& filename, double& z)
{
	//Sauvegarde du fichier
	std::ofstream fic(filename);
	fic << "\t\t *** Matrice des facteurs de vue pour le reservoir TANK 2. *** \n\n";
	fic << " ";
	fic << "\t\t\t Niveau haut du liquide = " << z << " m.\n";
	fic << " ";
	fic << data << "\n";
	fic.close();
	//std::cout << "\n\n\n\n\tWriting (summary of view factor calculations) to the file is COMPLETE ! \n\n\n\n" << filename << std::endl;
}


/*
template<class T>
void print_vect(std::string comment, const std::vector<T>& v)
{
	std::cout << comment << "\n";
	for (auto&& value : v)
	{
		std::cout << "\t" << value << "\n";
	}
}*/


template<class T, class U>
void print_map(std::string comment, const std::map<T, U>& m)
{
	std::cout << comment << "\n";
	// Iterate using C++17 facilities
	for (const auto& [key, value] : m)
		std::cout << '[' << key << "] = " << value << "; ";

	// C++11 alternative:
	//  for (const auto& n : m)
	//      std::cout << n.first << " = " << n.second << "; ";
	//
	// C++98 alternative:
	//  for (std::map<std::string, int>::const_iterator it = m.begin(); it != m.end(); ++it)
	//      std::cout << it->first << " = " << it->second << "; ";

	std::cout << '\n';
}


std::complex<double> g1(double& w)
{
	std::complex<double> j = { (0,1) };
	std::complex<double> z = 1.0 / (1.0 + j * w);
	return z;
	//return 293.0;
}

std::complex<double> F(double& w, double& t)
{
	std::complex<double> j = { (0,1) };
	return exp(j * w * t) * g1(w);
}

std::complex<double> rectangle_method(std::complex<double>(*func)(double& w, double& t), std::vector<double>& omega, double& t)
{
	int n = omega.size();
	double hw = omega[1] - omega[0];
	std::complex<double> som = { 0.0 };
	for (auto&& w : omega)
	{
		som += hw * func(w, t);
	}
	std::complex<double> fw_t = som / (2 * PI);
	return fw_t;
}

std::vector<std::complex<double>> inverse_fourier_transform_rectangles(std::complex<double>(*func)(double& w, double& t), std::vector<double>& omega, std::vector<double>& t_values)
{
	std::vector<std::complex<double>> ft_values;
	for (auto&& t : t_values)
	{
		std::complex<double> ft = rectangle_method(func, omega, t);
		ft_values.push_back(ft);
	}
	return ft_values;
}




/* LAPLACE TRANSFORM et LAPLACE TRANSFORM INVERSE */
double f_t(double& t)
{
	return exp(-t);
}

double integrand_Laplace(double(*func)(double& t), double& s, double& t)
{
	return func(t) * exp(-s * t);
}

double rectangle_method_Laplace(double(*func)(double& s, double& t), std::vector<double>& t_values, double& s)
{
	double total_area(0.0);
	int n = t_values.size();
	double h = static_cast<double>(t_values[1] - t_values[0]);
	h = h / (n - 1);
	for (auto&& t : t_values)
	{
		double height = func(s, t);
		total_area += h * height;
	}
	return total_area;
}

double Laplace_transform(double(*func)(double& s, double& t), std::vector<double>& t_values, double& s)
{
	double Fs = rectangle_method_Laplace(func, t_values, s);
	return Fs;
}








/*
double g(double& x)
{
	return x * x - x - 1;
}*/

/*
double dichotomy_method(double& a, double& b, double& precision, int& max_iter)
{
	if (g(a) * g(b) >= 0.)
	{
		std::cerr << "\tOupsss! Value error: f(a) et f(b) must have opposite signs." << "\n";
	}

	double sol{}, m;
	int iter = { 0 };
	while ((b - a) / 2 > precision)
	{
		if (iter > max_iter)
		{
			std::cerr << "Maximum iterations reached !" << "\n";
			break;
		}
		else
		{
			m = (a + b) / 2;
			if (g(m) == 0)
			{
				sol = m;
			}
			else if (g(a) * g(b) < 0)
			{
				b = m;
			}
			else
			{
				a = m;
			}
			iter += 1;
		}
	}
	return sol;
}*/


int readFile(std::string& filename)
{
	// Ouvrir le fichier en lecture
	std::ifstream file(filename);
	if (!file.is_open())
	{
		std::cerr << "Erreur: Impossible d'ouvrit le fichier " + filename << "\n";
		exit(2);
	}
	// Lecture des données depuis le fichier dqns l'objet classe
	std::string line;
	int n = { 0 };
	while (std::getline(file, line))
	{
		//std::cout << line << "\n";
		n++;
	}

	// Fermeture du fichier
	file.close();

	return n;
}

template<class T>
bool exist(std::vector<T>& tmp, T& tmp_target)
{
	bool test = false;
	for (auto&& t : tmp)
	{
		if (t == tmp_target)
		{
			test = true;
			break;
		}
	}
	return test;
}


template<class T>
bool contains(const std::vector<T>& vec, T element) {
	return std::find(vec.begin(), vec.end(), element) != vec.end();
}

template<class T>
void addElementWithoutDuplicates(std::vector<T>& vec, T element) {
	if (!contains(vec, element)) {
		vec.push_back(element);
	}
}

/*
dble get_surface(point<dble>& A, point<dble>& B, point<dble>& C) {
	dble area = {};
	geometry::vecteur<dble> AB(A, B), AC(A, C);
	geometry::vecteur<dble> u = AB ^ AC;
	area = u.normL2();
	return area;
}*/


std::vector<int> decomposition(const std::string& s)
{
	std::vector<int>listPlanesNumber;
	for (size_t i = 0; i < s.length(); ++i)
	{
		char ch = s[i];
		int planeNumber = ch - '0';
		listPlanesNumber.push_back(planeNumber);
	}
	return listPlanesNumber;
}

/*
void algo(const std::string& filename, composite_polygone<dble>& poly_i, composite_polygone<dble>& poly_j)
{
	dble area_i{}, area_j{};
	dble emissivity_i{}, emissivity_j{};
	dble T_i{}, T_j{};

	int size_mesh_by_default{ 50 };

	area_i = 2.0 * poly_i.get_surface();
	area_j = 2.0 * poly_j.get_surface();
	THERMAL::View_factor3<dble> vfij(poly_i, poly_j);
}
*/

template<class T>
std::vector<T> solve_2ndDegreeEquation(T& a, T& b, T& c)
{
	std::vector<T> sol{};
	T delta = b * b - 4 * a * c;
	T m1 = (-b + sqrt(delta)) / (2 * a);
	sol.push_back(m1);
	if (m1 != 0)
	{
		T m2 = c / m1;
		sol.push_back(m2);
	}
	return sol;
}



template<class T>
bool isOrthogonal(geometry::vecteur<T>& u, geometry::vecteur<T>& v)
{
	bool test = { false };
	T w = u.operator,(v);//dot product
	//std::cout << w << "\n";
	/*geometry::vecteur<T>& w = u.operator^(v);//cross product
	if ( (w.x() == 0) && (w.y() == 0) && (w.z() == 0) ) test = false;*/
	if (w == 0.0)
	{
		test = true;
	}
	return test;
}



double getMinElementList(std::vector<double>& vd)
{
	double D = { 9999999.0 };
	for (auto&& d : vd)
	{
		if (d < D)
		{
			D = d;
		}
	}
	return D;
}


double distancePointToPlan(point<double>& M, vector<point<double>>& vpoint)
{
	double D = {};
	vector<double> vddistance;
	double d_PM = { 9999999.0 };
	for (auto&& P : vpoint)
	{
		d_PM = distance(P, M);
		vddistance.push_back(d_PM);
	}
	D = getMinElementList(vddistance);
	return D;
}


template<class T, class U>
void results(const std::map<T, U>& mresults, const std::string& filename)
{
	std::ofstream fic(filename);
	for (const auto& [key, value] : mresults)
	{
		fic << key << ";" << value << "\n";
	}
	fic.close();
	std::cout << "\t\\t\tWriting (summary of view factor calculations) to the file is COMPLETE !   " << filename << std::endl;
}


template<class T>
void summary(std::string& filename)
{
	//Sauvegarde du fichier
	std::ofstream fic(filename);
	fic << "\n";
	fic << "\t\n";
	fic << "\t\t\t*** SUMMARY *** \n";
	fic << "\t\n";
	//fic << "\tAngle_1 (RAD): " << this->angle1_rad << ", emitter_plan_number: " << this->emitter_plan_number << "\tAngle_2 (RAD): " << this->angle2_rad << ", receiver_plan_number: " << this->receiver_plan_number << "\n\n";
	//fic << "\tSize mesh: " << this->size_mesh << ", Global axis ref.: " << this->axis << "\n\n";
	fic << "\tArea of emitter [M2]:  " << this->area_1 << ", Area of receiver [M2]: " << this->area_2 << ",  FACTOR of VIEW: " << this->value << "\n";
	fic << "\n\n\n";
	fic.close();
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
double a11(dble& x, dble& y)//f22 = f2,y
{
	return ((y + 1) * exp(y));
}
double a12(dble& x, dble& y)//f12 = f1,y
{
	return (2 * y - x * exp(y));
}
double a21(dble& x, dble& y)//f21 = f2,x
{
	return (-3 * x * x);
}
double a22(dble& x, dble& y)//f11 = f1,x
{
	return (1 + exp(y));
}



















int	main(void)
{


	std::system("cls");

	typedef double dble;

	typedef linearAlgebra::Array<dble> matd_viewFactor;


	int i{}, j{};
	dble vf_12 = {}, vfij_upper = {};
	dble coef_ij__MAT_VF{};
	dble coef_ij_Q{};
	GMSH::QuantityOfHeatByRadiation Q_rad{};
	std::string stri{}, strj{};
	string pathname = {};
	string dirname = {};
	string filename{}, fic{}, ext{};
	vector<point<dble>>	vpts1{}, vpts2{}, vpts3{};
	vector<vector<point<dble>>> vvpts{};
	dble area_1 = {}, area_2 = {};

	char axis{};
	int size_mesh_by_default = { 15 };//by default
	int ref_plane_number = {};

	dble H = {};
	dble emissivity_i{ 0.8 }, T_i{ 133.0 }, emissivity_j{ 0.9 }, T_j{ 165.0 }, delta_time{ 3600.0 };//by default
	map< int, vector<point<dble>> > mvpoint = {};
	vector<vector<point<dble>>> vvpoint = {};
	vector<point<dble>>vpoint = {};
	map< int, dble> marea = {};
	std::vector< geometry::composite_polygone<dble>> vpoly{};
	map<int, dble> mtemp = {}, memis = {};
	map <int, dble> area_final{};

	/* COMMENT DEFINIR (REMPLIR) LES POINTS P0, P1 et P2 ?
	* 1) Chaque plan i doit avoir pour base-locale (xi,yi)
	* 2) l'indice 1 (2) correspond au point se trouvant sur l'axe xi (yi)
	* 3) les COORDONNEES de chaque points sont données dans la base GLOBALE !!!
	*/

	std::vector< geometry::composite_polygone<dble>> v_poly;
	vector<vector<point<dble>>> vvP{};
	std::vector<point<dble>> v_pts_plan, vP{}, VT{};
	geometry::point<dble> points_plan, P{};
	std::vector<int> listePlanesNumber{};
	std::vector<dble> vd_fij{};
	dble vfij_k{};
	int nbre_ligne{};
	int nbre_plan{};
	std::map<std::string, std::string> mfile;
	std::vector<std::string>list_mesh;
	std::string information = {};
	std::string version = {};
	std::string date = {};
	int B1{}, C1{}, C2{}, C3{}, C4{}, HAUT{}, L{};
	int N1 = {};
	dble e1{}, T1{};
	int N2 = {};
	dble e2{}, T2{};
	dble hx1{}, hy1{}, hx2{}, hy2{};
	dble val_f12 = { 0.0 };
	dble z{};

	GMSH::View_Factor f12{};
	GMSH::Geo* geo_old = new GMSH::Geo[2];




	/*************************** DEBUT *******************************/

	std::system("cls");
	std::system("PAUSE");


	std::cout << "\n\n\n\ttenseur (vecteur, matrice,...) !\n";

	typedef std::vector<LINALG::Tensor<double>*> typ_vect_dble;
	typ_vect_dble* ptr_x;

	ptr_x = new typ_vect_dble[10];

	delete[] ptr_x;
	ptr_x = 0;








	LINALG::Vector<double>* ptr_Vector;
	//ptr_Vector = new LINALG::Vector<double>[10];
	//double* ptr_x = new double[10];




	std::system("PAUSE");
	std::system("cls");





	
	enum METHNUM {PIVOT=1, DICO=2, NEWTONRAPH=3, GC=4, JACOBI=5, LU=6, PUISS_ITER=100, DEFAULT};
	switch (PUISS_ITER)
	{
		case 1:
		{
			std::cout << "\n\tmethode du pivot." << "\n";
			Matrix<dble> A1(2, 2, 0), inv_A1(2, 2, 0);
			A1.setValue(0, 0, 2); A1.setValue(0, 1,-1);
			A1.setValue(1, 0,-1); A1.setValue(1, 1, 1);
			A1.print("A=");
			//inv_A1 = A1.inv();
			//std::cout << inv_A1<< "\n";
			Matrix<dble> b1(2, 1, 0);
			b1.setValue(0, 0, 1);
			b1.setValue(1, 0, 0);
			b1.print("b=");
			methods::systLin::direct::Pivot<dble> Piv(A1, b1);
			break;
		}		
		case 2:
		{
			std::cout << "\n\tmethode de la dichotomie." << "\n";
			methods::solveEquations::Dichotomie<dble> Dico(1.0, 2.0, 1e-6, 100);
			break;
		}		
		case 3:
		{		
			std::cout << "\n\tmethode de Newton-Raphson." << "\n";
			Matrix<dble> A(2, 2, 0), B(2, 1, 0), X(2, 1, 0);
			dble xi = {  2.0 };
			dble yi = { -2.0 };//on suppose uniforme
			B.setValue(0, 0, f1(xi, yi));
			B.setValue(1, 0, f2(xi, yi));
			B.print("B = ");
			A.setValue(0, 0, a11(xi, yi)); A.setValue(0, 1, a12(xi, yi));
			A.setValue(1, 0, a21(xi, yi)); A.setValue(1, 1, a22(xi, yi));
			A.print("A = ");
			X.setValue(0, 0, 2);
			X.setValue(1, 0, -2);
			X.print("X0 = ");
			std::system("pause");
			methods::solveEquations::NewtonRaphson<dble> NewtonRaph(A, B, X);
			break;
		}
		case 4:
		{
			std::cout << "\n\tmethode du gradient conjugue." << "\n";
			Matrix<dble> A(2, 2, 0), B(2, 1, 0), X(2, 1, 0);
			A.setValue(0, 0, 2); A.setValue(0, 1, -1);
			A.setValue(1, 0, -1); A.setValue(1, 1, 1);
			A.print("A=");
			B.setValue(0, 0, 1);
			B.setValue(1, 0, 0);
			B.print("B=");
			X.setValue(0, 0, 0);
			X.setValue(1, 0, 0);
			X.print("X0 = ");
			methods::systLin::project::GradientConjugue<dble> GC(A, B, X);
			//methods::solveEquations::Dichotomie<dble> Dico(1.0, 2.0, 1e-6, 100);
			break;
		}
		case 5:
		{
			std::cout << "\n\tmethode de Jacobi." << "\n";
			Matrix<dble> A(3, 3, 0), B(3, 1, 0), X(3, 1, 0);
			A.setValue(0, 0, 4); A.setValue(0, 1, 2); A.setValue(0, 2, 1);
			A.setValue(1, 0, -1); A.setValue(1, 1, 2); A.setValue(1, 2, 0);
			A.setValue(2, 0, 2); A.setValue(2, 1, 1); A.setValue(2, 2, 4);
			A.print("A=");
			B.setValue(0, 0, 4);
			B.setValue(1, 0, 2);
			B.setValue(2, 0, 9);
			B.print("B=");
			X.setValue(0, 0, 0);
			X.setValue(1, 0, 0);
			X.setValue(2, 0, 0);
			X.print("X0 = ");
			methods::systLin::iter::Jacobi<dble> jacob(A, B, X);
			//methods::solveEquations::Dichotomie<dble> Dico(1.0, 2.0, 1e-6, 100);
			break;
		}
		case 6:
		{
			std::cout << "\n\tmethode PA = LU." << "\n";
			/*
			Matrix<dble> A(2, 2, 0), B(2, 1, 0), X(2, 1, 0);
			A.setValue(0, 0, 2); A.setValue(0, 1, -1);
			A.setValue(1, 0, -1); A.setValue(1, 1, 1);
			A.print("A=");
			B.setValue(0, 0, 1);
			B.setValue(1, 0, 0);
			B.print("B=");
			X.setValue(0, 0, 0);
			X.setValue(1, 0, 0);
			X.print("X0 = ");			
			*/			
			Matrix<dble> A(3, 3, 0), B(3, 1, 0), X(3, 1, 0);
			A.setValue(0, 0, 1); A.setValue(0, 1, 1); A.setValue(0, 2, 0);
			A.setValue(1, 0, 0); A.setValue(1, 1, 2); A.setValue(1, 2, 0);
			A.setValue(2, 0, 0); A.setValue(2, 1,-1); A.setValue(2, 2, 3);
			A.print("A=");
			B.setValue(0, 0, 4);
			B.setValue(1, 0, 2);
			B.setValue(2, 0, 9);
			B.print("B=");
			X.setValue(0, 0, 0);
			X.setValue(1, 0, 0);
			X.setValue(2, 0, 0);
			X.print("X0 = ");						
			methods::systLin::direct::LU<dble> lu(A, B, X);
			break;
		}
		case 100:
		{
			std::cout << "\n\tpuissances iterees: " << "\n";
			Matrix<dble> A(3, 3, 0), X(3, 1, 0);
			A.setValue(0, 0, 1); A.setValue(0, 1, 1); A.setValue(0, 2, 0);
			A.setValue(1, 0, 0); A.setValue(1, 1, 2); A.setValue(1, 2, 0);
			A.setValue(2, 0, 0); A.setValue(2, 1, -1); A.setValue(2, 2, 3);
			A.print("A=");
			X.setValue(0, 0, 0);
			X.setValue(1, 0, 1);
			X.setValue(2, 0, 0);
			X.print("X=");
			methods::eigenValues::PuissancesIterees<dble> puis_iter(A, X);

			break;
		}
		default:
		{
			std::cout << "\taucune methode numerique connue !!!" << "\n";
			break;
		}
	}






	/* testing namespace SYSLIN */
	//methods::solveEquations::Dichotomie Dico(1.0, 2.0, 1e-6, 100);

	std::cout << "\n\n\t -----" << "\n\n";



	




	std::cout << "\tFIN   test classe du fichier array_class.h " << "\n\n\n";


	std::system("PAUSE");
	std::system("cls");
	/***************************  FIN  *******************************/







	pathname = { "C:/Users/jab/dev/cpp/projets/sloshing/data/geometries_GMSH/surfaces_reservoir_tank2/" };
	std::vector<std::string> list_maillages;
	std::string message = { "\t\t\t'view factor' computations: IN PROGRESS..." };
	std::string folder;
	int option = { 4 };


	switch (option)
	{
	case 1:
		z = { 29.125 - 8.794 };
		folder = { "cas1" };
		break;
	case 12:
		z = { 29.125 - 7.694 };
		folder = { "cas12" };
		break;
	case 2:
		z = { 29.125 - 6.596 };
		folder = { "cas2" };
		break;
	case 23:
		z = { 29.125 - 5.496 };
		folder = { "cas23" };
		break;
	case 3:
		z = { 29.125 - 4.397 };
		folder = { "cas3" };
		break;
	case 34:
		z = { 29.125 - 3.297 };
		folder = { "cas34" };
		break;
	case 4:
		z = { 29.125 - 2.198 };
		folder = { "cas4" };
		break;
	case 45:
		z = { 29.125 - 1.099 };
		folder = { "cas45" };
		break;
	case 5:
		z = { 29.125 - 0.0005 };
		folder = { "cas5" };
		break;
	case 6:
		z = { 29.125 - 0.0 };
		folder = { "cas6" };
		break;
	case 7:
		z = { 29.125 - 4.397 };
		folder = { "cas_spe" };
		break;
	default:
		std::cout << "Option invalid !" << std::endl;
		break;
	}


	std::cout << "\n\t\t*** INFORMATIONS ***" << "\n\n";
	std::cout << "\t\t - fluid level in tank  : " << z << " m." << std::endl;
	std::cout << "\t\t - tank fill percentage : " << (z / 29.125) * 100 << " %." << "\n\n\n";
	std::this_thread::sleep_for(std::chrono::milliseconds(1500));





	pathname = pathname + "/" + folder + "/";
	//list_maillages = { pathname + "surf0", pathname + "surf1" };
	list_maillages = { pathname + "surf0", pathname + "surf1", pathname + "surf2", pathname + "surf3", pathname + "surf4", pathname + "surf5" };
	size_t numberOfFiles = list_maillages.size() * list_maillages.size() - list_maillages.size();


	//std::cout << "\n\tMATRIX of VIEW FACTORS is DONE !" << "\n\n\n";
	/* AFFICHAGE de la MATRICE des FACTEURS de VUE aprés RENORMALISATION */
	/* ECRITURE de la MATRICE des FACTEURS de VUE dans un FICHIER */	



	std::cout << "\tPour informations: \n";
	std::cout << "\t-----------------  \n\n";
	std::cout << "Si vous souhaitez ameliorer la precision des valeurs sur les coefficients de la matrice des facteurs de vue, " << "\n";
	std::cout << "il vous faudra modifier la valeur du nombre de triangles dans le maillage des surfaces (Voir parametres dans les fichiers .geo pour GMSH)." << "\n\n\n";
	std::system("PAUSE");
	std::system("cls");
	
	

	
	std::cout << "\n\n";
	std::cout << "\t\t\t***********************************\n";  
	std::cout << "\t\t\t        VIEW FACTOR MATRIX (N)     \n";
	std::cout << "\t\t\t***********************************\n";
    GMSH::CalculsFacteursDeVue<double> cfv(list_maillages);



	
	fic = { "View_Factor_Matrix_" }, ext = { ".txt" };
	filename = pathname + fic + folder + ext;
	//printFile(MAT_VF2, filename, z);

	std::cout << "\n\n\t\tWRITING VIEW FACTOR MATRIX TO FILE IS DONE !\n\n";
	std::this_thread::sleep_for(std::chrono::milliseconds(2000));//pause de 2000 ms
	



	std::system("PAUSE");
	std::system("cls");





	return 0;




}

