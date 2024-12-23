

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



#include "array_class.h"



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
	int pos = barWidth * percentage;
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









double g(double& x)
{
	return x * x - x - 1;
}

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
}


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




	std::cout << "\tDEBUT test classe du fichier array_class.h " << "\n";


	typedef Matrix<dble> matd22;
	matd22 A(2,2,0), B(2,2,0), C(2,2,0), D(2,2,0), G(2,2,0), Q(2,2,23.14);
	A.setValue(0, 0, 1.0); A.setValue(0, 1, 4.0);
	A.setValue(1, 0, 3.0); A.setValue(1, 1, 2.0);
	A.print("A=");



	B = A;
	B.print("B = ");


	std:: cout << "\tB est-elle symetrique ? " << B.isSymetric() << "\n";

	C = A + B;
	C.print("C = ");





	matd22 I(2,2,0), J(2,2,0), K(2,2,0);
	I = I.eye();
	I.print("I");

	std::cout << "\tI est-elle symetrique ? " << I.isSymetric() << "\n";

	J = J.ones();
	J.print("J");
	std::cout << "\tJ est-elle symetrique ? " << J.isSymetric() << "\n";



	K = J - I/2;
	K += I;
	K.print("K");


	
	G = A.operator*(A);
	G.print("G = A*A");


	
	Q.print("Q");



	std::cout << "\ttrace(G) = " <<  G.trace() << "\n";

	/*
	typedef Matrix<dble> matd33;
	matd33 matH(3,3,0);

	std::cout << "\n\tH = " << "\n";
	//std::cout << matH << "\n";


	typedef Matrix<dble> matd32;
	matd32 AmatH(2,3,0);
	AmatH.print("A*matH = ");
	*/



	// Multiplier la matrice par un scalaire
	int scalar = 2;
	Matrix<int> MAT_A(2,2,1), MAT_B(2,2,2);
	std::cout << MAT_A << "\n";
	std::cout << MAT_B << "\n";
	
	
	

	Matrix<int> result(2, 2, 7);
	result.setValue(0, 1, -8);
	result = scalar* result;
	result.print("result = ");
	result = result * 6;
	std::cout << result << "\n";

	/*
	Matrix<int> MAT_C = MAT_A * MAT_B;
	MAT_C.print("MAT_C = ");*/

	Matrix<int> MAT_D(3,2,6); 
	MAT_D.setValue(2, 0, -5);
	MAT_D.setValue(1, 1, -11);
	MAT_D.print("MAT_D = ");

	MAT_D = MAT_D.transpose();
	MAT_D.print("tMAT_D = ");

	Matrix<int> MAT_E(2,7,2); MAT_E.print("MAT_E = ");



	
	
	/////////Matrix<int> MAT_F = MAT_D * MAT_E;
	/*
	MAT_E = MAT_D * MAT_D;
	MAT_F.print("MAT_C = ");*/
	

	std::cout << "\tFIN   test classe du fichier array_class.h " << "\n\n\n";


	std::system("PAUSE");
	std::system("cls");
	/***************************  FIN  *******************************/




	return 0;




}

