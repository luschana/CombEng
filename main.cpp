/*
 * main.cpp
 *
 *  Created on: May 11, 2018
 *      Author: Alex Luschan <alexander.luschan@gmail.com>
 */

#include "combustionengine.h"

void testPt(bool *bArray, GasComponent **pgc){
	for (int i = 0; i < 3; i++) {
		if(bArray[i]){
			bArray[i] = !bArray[i];
		}else{
			bArray[i] = !bArray[i];
		}
		pgc[i]->setCombustionStarted(bArray[i]);
	}
}

int main(int argc, char **argv) {
	CombustionEngine* pCE = CombustionEngine::getInst();
	double w = 60;
	double thrPos = 0.0;
	while(true){
		pCE->run(w, thrPos);
	}
	/*bool bVals[3] = {true, false, true};
	GasComponent * pgc[3];
	for (int i = 0; i < 3; i++) {
		pgc[i] = new GasComponent();
	}
	testPt(bVals, pgc);*/
    return 0;
}


