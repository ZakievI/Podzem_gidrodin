#include <iostream>
#include "Zadanie#1.h"
#include "Zadanie#2.h"
#include <vector>
   
int main()
{
	std::string mesh = "mesh.dat";
	Zadanie_2(0, 1e+6, mesh); //MKV(R, r, p_left, p_right, N_mesh)
	//MKR(1, 0.1, 10, 100, 100);
	//Zadanie_1(1e+6,0, mesh); //MKV(R, r, p_left, p_right, N_mesh)
	
}
