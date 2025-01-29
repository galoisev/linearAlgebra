#ifndef __VIEWFACTOR_H
#define __VIEWFACTOR_H


#include "nodes.h"
#include "mesh.h"

#include <iostream>
#include <istream>
#include <fstream>
#include <string>
#include <vector>
#include <cassert>


namespace GMSH
{



    typedef double dble;
    //#define PI 3.14159265
    //constexpr double PI = 3.14159265358979;


    class View_Factor //nouvelle version dans le cas d'un maillage fait avec GMSH
    {
    public:
        //std::vector<Point>vP1, vP2;
        dble area1{}, area2{};
        int numberOfElements1{}, numberOfElements2{};
        int emitterNumber{}, receiverNumber{};
        dble value{0.0};
        Mesh msh;
        Elements elem1, elem2;
        std::vector<Elements>list_elem1, list_elem2;
        Nodes N1, N2;


        std::vector<GMSH::Nodes>vp1, vp2;
        GMSH::Nodes* ptr_node_plan1, * ptr_node_plan2;


        View_Factor() {}
        View_Factor(const View_Factor& vf) = default;
        View_Factor& operator=(const View_Factor& vf) = default;

        View_Factor operator*(double c) {}


        dble coaxialParallelDisks(dble& L, dble& Ri, dble& Rj); //special case
        dble equalRectangularPlates(dble& Wi, dble& Wj, dble& H);
        dble rectangleToRectangle(dble& H_ij, std::vector<dble>& ksi, std::vector<dble>& eta, std::vector<dble>& x, std::vector<dble>& y);


        /*
        friend dble coef_ij(GMSH::Nodes& P1, GMSH::Nodes& P2)
        {
            GMSH::Nodes P1P2(P1, P2), P2P1(P2, P1);
            dble rij = P1P2.normeL2();

            GMSH::Nodes N1(11, 0.0, 0.0, 1.0), N2(22, 1.0, 0.0, 0.0);//TODO

            dble cij = (P1P2, N1) * (P2P1, N2) / (rij * rij * rij * rij) / PI;
            return cij;
        }*/

        void printFile(const std::string& filename);

        dble f(Nodes* p1, Nodes* N1, Nodes* p2, Nodes* N2)
        {
            GMSH::Nodes P1P2(*p1, *p2), P2P1(*p2, *p1);
            dble rij = P1P2.normeL2();
            double integrand = (P1P2, *N1) * (P2P1, *N2) * (1.0 / (rij * rij * rij * rij)) * (1.0 / PI);
            return integrand;
        }


        dble calculs(dble hx1, dble hy1, dble hx2, dble hy2, std::vector<GMSH::Nodes*> vp1, std::vector<GMSH::Nodes*> vp2)
        {
            GMSH::Nodes X(0, 1E6, 0.0, 0.0);


            GMSH::Nodes O1(int("O1"), vp1.operator[](0)->x, vp1.operator[](0)->y, vp1.operator[](0)->z);
            //std::cout << O1 << "\n";
            GMSH::Nodes P1(int("P1"), vp1.operator[](1)->x, vp1.operator[](1)->y, vp1.operator[](1)->z);
            //std::cout << P1 << "\n";
            GMSH::Nodes Q1(int("Q1"), vp1.operator[](3)->x, vp1.operator[](3)->y, vp1.operator[](3)->z);
            //std::cout << Q1 << "\n";

            //////////////////////////////////////////list_elem1.push_back(elem1);

            GMSH::Nodes O2(int("O2"), vp2.operator[](0)->x, vp2.operator[](0)->y, vp2.operator[](0)->z);
            //std::cout << O2 << "\n";
            GMSH::Nodes P2(int("P2"), vp2.operator[](1)->x, vp2.operator[](1)->y, vp2.operator[](1)->z);
            //std::cout << P2 << "\n";
            GMSH::Nodes Q2(int("Q2"), vp2.operator[](3)->x, vp2.operator[](3)->y, vp2.operator[](3)->z);
            //std::cout << Q2 << "\n";



            GMSH::Nodes u1(O1, X), v1(O1, P1), u2(O2, X), v2(O2, P2);
            GMSH::Nodes O1P1(O1, P1), O1Q1(O1, Q1), O2P2(O2, P2), O2Q2(O2, Q2);

            /*
            GMSH::Nodes N1 = O1P1.operator^(O1Q1);
            N1 = N1.operator*(N1.normeL2());//normalisation du vecteur normal N1
            */

            Nodes N1 = normal(O1, P1, Q1);
            this->N1 = N1;
            //std::cout << N1 << "\n";

            /*
            GMSH::Nodes N2 = O2P2.operator^(O2Q2);
            N2 = N2.operator*(N2.normeL2());//normalisation du vecteur normal N2
            */

            Nodes N2 = normal(O2, P2, Q2);
            this->N2 = N2;
            //std::cout << N2 << "\n";

            /*
            this->angle1_rad = acos((u1.operator,(v1)) / (u1.normeL2() * v1.normeL2()));
            this->angle2_rad = acos((u2.operator,(v2)) / (u2.normeL2() * v2.normeL2()));
            */        
            /*
            this->angle1_rad = acos((N1.operator,(O1P1)) / (N1.normeL2() * O1P1.normeL2()));
            this->angle2_rad = acos((N2.operator,(O2P2)) / (N1.normeL2() * O2P2.normeL2()));
            */

            dble sum_i = { 0.0 };
            for (auto&& p1 : vp1)//p1 correspond au point (xi,yi,zi)
            {
                dble sum_j = { 0.0 };
                for (auto&& p2 : vp2)//p2 correspond au point (xj,yj,zj)
                {
                    GMSH::Nodes P1P2(*p1, *p2), P2P1(*p2, *p1);

                    dble rij = P1P2.normeL2();                   
                    if ( (this->area1 != 0.0) && (rij != 0.0) )
                    {    
                        sum_j += abs(f(p1, &N1, p2, &N2));
                    }
                }
                sum_i += sum_j;
            }
            if (this->area1 > 0.0)
            {
                this->value = sum_i;
            }
            return this->value * (1.0 / this->area1);
        }




        dble calculs2(std::vector<GMSH::Nodes*> vp1, std::vector<GMSH::Nodes*> vp2)
        {
            GMSH::Nodes X(0, 1E6, 0.0, 0.0);

            GMSH::Nodes O1(int("O1"), vp1.operator[](0)->x, vp1.operator[](0)->y, vp1.operator[](0)->z);
            GMSH::Nodes P1(int("P1"), vp1.operator[](1)->x, vp1.operator[](1)->y, vp1.operator[](1)->z);
            GMSH::Nodes Q1(int("Q1"), vp1.operator[](3)->x, vp1.operator[](3)->y, vp1.operator[](3)->z);

            GMSH::Nodes O2(int("O2"), vp2.operator[](0)->x, vp2.operator[](0)->y, vp2.operator[](0)->z);
            GMSH::Nodes P2(int("P2"), vp2.operator[](1)->x, vp2.operator[](1)->y, vp2.operator[](1)->z);
            GMSH::Nodes Q2(int("Q2"), vp2.operator[](3)->x, vp2.operator[](3)->y, vp2.operator[](3)->z);


            GMSH::Nodes u1(O1, X), v1(O1, P1), u2(O2, X), v2(O2, P2);
            GMSH::Nodes O1P1(O1, P1), O1Q1(O1, Q1), O2P2(O2, P2), O2Q2(O2, Q2);

            Nodes N1 = normal(O1, P1, Q1);
            this->N1 = N1;

            Nodes N2 = normal(O2, P2, Q2);
            this->N2 = N2;

            /*
            this->angle1_rad = acos((u1.operator,(v1)) / (u1.normeL2() * v1.normeL2()));
            this->angle2_rad = acos((u2.operator,(v2)) / (u2.normeL2() * v2.normeL2()));
            */
            /*
            this->angle1_rad = acos((O1P1.operator,(O1Q1)) / (O1P1.normeL2() * O1Q1.normeL2()));
            this->angle2_rad = acos((O2P2.operator,(O2Q2)) / (O2P2.normeL2() * O2Q2.normeL2()));
            */

            //dble mean_area_i = { 0.022 }, mean_area_j = { 0.022 };
            dble sum_i = { 0.0 };
            for (auto&& p1 : vp1)//p1 correspond au point (xi,yi,zi)
            {
                dble sum_j = { 0.0 };
                for (auto&& p2 : vp2)//p2 correspond au point (xj,yj,zj)
                {
                    GMSH::Nodes P1P2(*p1, *p2), P2P1(*p2, *p1);

                    dble rij = P1P2.normeL2();
                    if ((this->area1 != 0.0) && (rij != 0.0))
                    {
                        sum_j += std::abs(f(p1, &N1, p2, &N2));
                    }
                }
                sum_i += sum_j;
            }
            if (this->area1 > 0.0)
            {
                this->value = sum_i* (1.0 / this->area1);
            }
            return this->value;
        }



        friend std::ostream& operator<<(std::ostream& f, const View_Factor& vf);
        friend std::istream& operator>>(std::istream& f, View_Factor& vf);


    };






    dble View_Factor::coaxialParallelDisks(dble& L, dble& Ri, dble& Rj)
    {
        // L : separation distance
        // Ri: Radius of Emitter
        // Rj: Radius of Receiver
        dble Fij = {};
        dble ri = Ri / L, rj = Rj / L;
        dble S = 1.0 + ((1.0 + pow(rj, 2.0)) / pow(ri, 2.0));
        Fij = 0.5 * (S - pow(pow(S, 2.0) - 4.0 * pow((Rj / Ri), 2.0), 0.5));
        return Fij;
    }



    dble View_Factor::equalRectangularPlates(dble& Wi, dble& Wj, dble& H)
    {
        dble fij{};
        dble x1{}, y1{}, x{}, y{};
        x = Wi / H;
        y = Wj / H;
        x1 = sqrt(1 + pow(x, 2.0));
        y1 = sqrt(1 + pow(y, 2.0));
        fij = (1. / (PI * x * y)) * (log((pow(x1 * y1, 2.0) / (pow(x1, 2.0) + pow(y1, 2.0) - 1.0))) + 2 * x * (y1 * atan(x / y1) - atan(x)) + 2 * y * (x1 * atan(y / x1) - atan(y)));
        return fij;
    }


    dble View_Factor::rectangleToRectangle(dble& H_ij, std::vector<dble>& ksi, std::vector<dble>& eta, std::vector<dble>& x, std::vector<dble>& y)
    {
        dble u, v, p, q, z{ H_ij }, B{ 0.0 };
        dble fij{};
        dble A1 = (x.operator[](1) - x.operator[](0)) * (y.operator[](1) - y.operator[](0));
        dble sum{ 0.0 };
        for (int i = 0; i < x.size(); i++)
        {
            for (int j = 0; j < y.size(); j++)
            {
                for (int k = 0; k < eta.size(); k++)
                {
                    for (int l = 0; l < ksi.size(); l++)
                    {
                        u = x.operator[](i) - ksi.operator[](l);
                        v = y.operator[](j) - eta.operator[](k);
                        p = sqrt(u * u + z * z);
                        q = sqrt(v * v + z * z);
                        B = v * p * atan(v / p) + u * q * atan(u / q) - 0.5 * z * z * log(u * u + v * v + z * z);
                        sum += pow(-1.0, (i + j + k + l)) * B;
                    }
                }
            }
        }
        fij = (1. / (2 * PI * A1)) * sum;
        return fij;
    }


    std::ostream& operator<<(std::ostream& f, const View_Factor& vf)
    {
        f << "\n";
        f << "\t\n";
        f << "\t\t\t*** SUMMARY *** \n";
        f << "\t\n";
        //fic << "\tAngle_1 (RAD): " << this->angle1_rad << ", emitter_plan_number: " << this->emitter_plan_number << "\tAngle_2 (RAD): " << this->angle2_rad << ", receiver_plan_number: " << this->receiver_plan_number << "\n\n";
        //fic << "\tSize mesh: " << this->size_mesh << ", Global axis ref.: " << this->axis << "\n\n";
        f << "\tArea of emitter  [M2] (" << vf.emitterNumber  << ") = " << vf.area1 << ", Number of elements in the emitter plane:  " << vf.numberOfElements1 << "\n";
        f << "\tArea of receiver [M2] (" << vf.receiverNumber << ") = " << vf.area2 << ", Number of elements in the emitter plane:  " << vf.numberOfElements2 << "\n";
        f << "\n\t\t\tFACTOR of VIEW:  " << vf.value << "\n";
        f << "\n\n\n";
        return f;
    }

    std::istream& operator>>(std::istream& f, View_Factor& vf)
    {
        f >> vf.value;
        return f;
    }

    void View_Factor::printFile(const std::string& filename)
    {
        //Sauvegarde du fichier
        std::ofstream fic(filename);
        fic << "\n";
        fic << "\t\n";
        fic << "\t\t\t*** SUMMARY *** \n";
        fic << "\t\n";
        //fic << "\tAngle_1 (RAD): " << this->angle1_rad << ", emitter_plan_number: " << this->emitter_plan_number << "\tAngle_2 (RAD): " << this->angle2_rad << ", receiver_plan_number: " << this->receiver_plan_number << "\n\n";
        //fic << "\tSize mesh: " << this->size_mesh << ", Global axis ref.: " << this->axis << "\n\n";
        fic << "\tArea of emitter  [M2] (" << this->emitterNumber << ") = " << this->area1 << ", Number of elements in the emitter plane:  " << this->numberOfElements1 << "\n";
        fic << "\tArea of receiver [M2] (" << this->receiverNumber << ") = " << this->area2 << ", Number of elements in the emitter plane:  " << this->numberOfElements2 << "\n";
        fic << "\n\t\t\tFACTOR of VIEW:  " << this->value << "\n";
        fic << "\n\n\n";
        fic.close();
    }







}//end namespace GMSH

#endif // __VIEWFACTOR_H
