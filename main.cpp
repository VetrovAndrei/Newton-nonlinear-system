#include "solve.h"

void main()
{
	   std::ifstream size ("size.txt");
       std::ifstream X("X.txt");
       solve mat(size, X);
	   try
	   {
		   mat.snu(false,true);
	   }
	   catch(int error)
	   {
		   if (error == 1)
				std::cout << "bad matrix" << std::endl;
		   if (error == 2)
			   std::cout << "bad point" << std::endl;
	   }
       system("pause");
}