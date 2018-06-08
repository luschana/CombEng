/* Copyright 2003-2004 The MathWorks, Inc. */

// *******************************************************************
// **** To build this mex function use: mex mexengine.cpp ****
// *******************************************************************

#define S_FUNCTION_LEVEL 2
#define S_FUNCTION_NAME mexengine

// Need to include simstruc.h for the definition of the SimStruct and
// its associated macro definitions.
#include "simstruc.h"
#include "combustionengine.h"


#define IS_PARAM_DOUBLE(pVal) (mxIsNumeric(pVal) && !mxIsLogical(pVal) &&\
!mxIsEmpty(pVal) && !mxIsSparse(pVal) && !mxIsComplex(pVal) && mxIsDouble(pVal))


// Function: mdlInitializeSizes ===============================================
// Abstract:
//    The sizes information is used by Simulink to determine the S-function
//    block's characteristics (number of inputs, outputs, states, etc.).
static void mdlInitializeSizes(SimStruct *S){
    // No expected parameters
    ssSetNumSFcnParams(S, 0);

    // Parameter mismatch will be reported by Simulink
    if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) {
        return;
    }

    // Specify I/O
    int nPortCntIn = 7;
    int nPortCntOut = 32;
    if (!ssSetNumInputPorts(S, 1)) return;
    ssSetInputPortWidth(S, 0, nPortCntIn);
    ssSetInputPortDirectFeedThrough(S, 0, 1);
    if (!ssSetNumOutputPorts(S,1)) return;
    ssSetOutputPortWidth(S, 0, nPortCntOut);

    ssSetNumSampleTimes(S, 1);

    // Reserve place for C++ object
    ssSetNumPWork(S, 1);

    ssSetSimStateCompliance(S, USE_CUSTOM_SIM_STATE);

    //ssSetOptions(S, SS_OPTION_WORKS_WITH_CODE_REUSE | SS_OPTION_EXCEPTION_FREE_CODE);
    ssSetOptions(S, SS_OPTION_EXCEPTION_FREE_CODE);

}


// Function: mdlInitializeSampleTimes =========================================
// Abstract:
//   This function is used to specify the sample time(s) for your
//   S-function. You must register the same number of sample times as
//   specified in ssSetNumSampleTimes.
static void mdlInitializeSampleTimes(SimStruct *S)
{
    ssSetSampleTime(S, 0, INHERITED_SAMPLE_TIME);
    ssSetOffsetTime(S, 0, 0.0);
    ssSetModelReferenceSampleTimeDefaultInheritance(S); 
}

// Function: mdlStart =======================================================
// Abstract:
//   This function is called once at start of model execution. If you
//   have states that should be initialized once, this is the place
//   to do it.
#define MDL_START
static void mdlStart(SimStruct *S)
{
    // Store new C++ object in the pointers vector
    CombustionEngine *engine  = CombustionEngine::getInst();
    ssGetPWork(S)[0] = engine;
}

// Function: mdlOutputs =======================================================
// Abstract:
//   In this function, you compute the outputs of your S-function
//   block.
static void mdlOutputs(SimStruct *S, int_T tid){
    // Retrieve C++ object from the pointers vector
    CombustionEngine *engine = static_cast<CombustionEngine *>(ssGetPWork(S)[0]);
    // Get data addresses of I/O
    InputRealPtrsType  u = ssGetInputPortRealSignalPtrs(S,0);
               real_T *y = ssGetOutputPortRealSignal(S, 0);

    // Call AddTo method and return peak value
	engine->setPhiSpark(*u[2]); // [degree] if(sparkAngle >= -40.0 && sparkAngle <= 20.0)
	engine->setPhiInjection(*u[3]);
	engine->setTempCoolingWater(*u[4]);
	/*
	 * [0::3]  M, dotQ, phi, freq
	 */
	y[0] = 0.0;
    y[1] = 0.0;
    y[5] = engine->getFuelConsumption();
    y[6] = 0.0;
	y[9] = 0.0;
	y[3] = engine->getPhi(); // get phi from prev. time step
    for (int i = 0; i < 10; i++) {
		engine->run(*u[0], *u[1]); // [rad/s] , [-] 0..1
		y[0] += engine->getMShaft();
        y[1] += engine->getH_Cooling(); //deltaH
        y[6] += 0.0;
        y[9] += 0.0;
	}
    y[0]/=10.0;
    y[1]/=(Ts*10.0); //dotQ
	y[2] = engine->getPhi();
	if((y[3] > y[2])&&(y[3]-y[2] > 11)){ // phi reset!
		y[3] -= 4.0*M_PI;
	}
	y[3]=(y[2]-y[3])/(2*M_PI)*10000; //frequency of rotation

    y[5] = (engine->getFuelConsumption() - y[5])/(Ts*10.0); // delta to dot

    y[4] = 0.0;
	for (int i = 0; i < 6; i++) {
		if (i < Ncyl){
			y[4] += engine->getCyl_T(i);
			y[19+i] = engine->getp_Cyl(i);
			y[26+i] = engine->getT_Cyl(i);
		}else{
			y[19+i] = 0.0;
			y[26+i] = 0.0;
		}
	}
	y[4] /= Ncyl;

    y[6] = engine->getIntake().getP();
	y[7] = engine->getIntake().getT();
	y[8] = engine->getExhaust().getP();
	y[9] = engine->getExhaust().getT();

	y[10] = engine->getCyl1().getValveIn()->getGasComponent()->getP();
	y[11] = engine->getCyl1().getValveIn()->getGasComponent()->getT();
	y[12] = engine->getCyl1().getValveIn()->getGasComponent()->getMols();

	y[13] = engine->getCyl1().getGasComponent().getP();
	y[14] = engine->getCyl1().getGasComponent().getT();
	y[15] = engine->getCyl1().getGasComponent().getMols();

	y[16] = engine->getCyl1().getValveOut()->getGasComponent()->getP();
	y[17] = engine->getCyl1().getValveOut()->getGasComponent()->getT();
	y[18] = engine->getCyl1().getValveOut()->getGasComponent()->getMols();

/*	y[9] = engine->getExhaust().getNu()[defs::O2];
	y[10] = engine->getExhaust().getNu()[defs::H2O];
	y[11] = engine->getExhaust().getNu()[defs::CO2];
	y[12] = engine->getExhaust().getNu()[defs::CO];
	y[13] = engine->getExhaust().getNu()[defs::H2];
	y[14] = engine->getExhaust().getNu()[defs::C];
	y[15] = engine->getExhaust().getNu()[defs::Fuel];
	// p1..6
	y[18] = engine->getCyl1().getValveIn()->getGasComponent()->getP();
	y[19] = engine->getCyl1().getValveIn()->getGasComponent()->getT();
	y[20] = engine->getCyl1().getValveIn()->getGasComponent()->getMols();

	y[22] = Environment::getInst()->getAmbientAir()->getMols();
	y[23] = Environment::getInst()->getExhaustGas()->getMols();
	// T1..6
	y[24] = engine->getCyl1().getGasComponent().getP();
	y[25] = engine->getCyl1().getGasComponent().getT();
	y[26] = engine->getCyl1().getGasComponent().getMols();

	y[27] = engine->getCyl1().getValveOut()->getGasComponent()->getP();
	y[28] = engine->getCyl1().getValveOut()->getGasComponent()->getT();
	y[29] = engine->getCyl1().getValveOut()->getGasComponent()->getMols();

	y[30] = engine->getH_Cooling()/Ts;
	y[31] = Environment::getInst()->getExhaustGas()->getP();
*/
}

/* Define to indicate that this S-Function has the mdlG[S]etSimState methods */
#define MDL_SIM_STATE

/* Function: mdlGetSimState =====================================================
 * Abstract:
 */
#ifdef ARTELAB
	//RT - version
	static double mdlGetSimState(SimStruct* S) {
	   return 0.0;
	}
#else
	//nonRT - version
	static mxArray* mdlGetSimState(SimStruct* S) {
	   // Retrieve C++ object from the pointers vector
	   return mxCreateDoubleScalar(0.0);
	}
#endif


/* Function: mdlGetSimState =====================================================
 * Abstract:
 *
 */
static void mdlSetSimState(SimStruct* S, const mxArray* ma)
{
    // Retrieve C++ object from the pointers vector
    CombustionEngine *engine = static_cast<CombustionEngine*>(ssGetPWork(S)[0]);    
	engine->run(mxGetPr(ma)[0], mxGetPr(ma)[0]);
}

// Function: mdlTerminate =====================================================
// Abstract:
//   In this function, you should perform any actions that are necessary
//   at the termination of a simulation.  For example, if memory was
//   allocated in mdlStart, this is the place to free it.
static void mdlTerminate(SimStruct *S)
{
    // Retrieve and destroy C++ object
    CombustionEngine *engine = static_cast<CombustionEngine*>(ssGetPWork(S)[0]);
    delete engine;
}


// Required S-function trailer
#ifdef  MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "simulink.c"      /* MEX-file interface mechanism */
#else
#include "cg_sfun.h"       /* Code generation registration function */
#endif
