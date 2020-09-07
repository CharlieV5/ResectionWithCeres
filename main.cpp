#include "Resection.h"
#include "ResectionWithDistortion.h"
#include "ResectionResidualWithDistortion.h"


int main(int argc, char** argv)
{
	Resection r1;
	r1.doResection();

	ResectionWithDistortion r2;
	r2.doResection();

	getchar();
	return 0;
}

