/*
 * main.cpp
 *
 *  Created on: May 11, 2018
 *      Author: Alex Luschan <alexander.luschan@gmail.com>
 */

#include "combustionengine.h"

int main(int argc, char **argv) {
	CombustionEngine* pCE = CombustionEngine::getInst();
	double w = 60;
	double thrPos = 100;
	while(true){
		pCE->run(w, thrPos);
	}
    return 0;
}


