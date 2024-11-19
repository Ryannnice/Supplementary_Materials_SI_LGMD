/*
 * Description: LGMD-Plus
 * Filename: LGMD_model.h
 * Author: Qinbing Fu
 * Date: 2019 Aug
*/

#ifndef __LGMD_PC_H
#define __LGMD_PC_H


#include <stdint.h>


#define Layer_Width 99
#define Layer_Height 72


typedef struct
{
  uint8_t P_hat[2][Layer_Height][Layer_Width];
	uint8_t E_hat[2][Layer_Height][Layer_Width];
	uint8_t SI[Layer_Height][Layer_Width];
	uint8_t S[Layer_Height][Layer_Width];	//Summation Layer
}LGMDstruct_Layer;

typedef struct
{
	
	uint8_t Tffi;	//FFI threshold
	uint8_t nsp;
	uint8_t nts;
	uint8_t n; //compute the number of the spike
	uint8_t flag;  // the flag of the continuous spike
	uint8_t tau_P;
	uint8_t tau_E;
	float Threshold;
	float interval;
	float beta1;
	float beta2;
	float beta3;
	float alpha1;
	float alpha2;
	float Taz;
	float Ts;
	float Kspi;
}LGMDstruct_Params;


typedef struct
{
	uint8_t SPIKE;
	float MP;	//membrane potential
	float FFI_out[2];	//FFI output
	float SMP;	//sigmoid membrane potential
	float LGMD_out[2];	//LGMDs output
}LGMDstruct_Result;

typedef struct
{
	LGMDstruct_Layer Layers;
	LGMDstruct_Params Params;
	LGMDstruct_Result Results;
	int16_t* pImg;
	int8_t* pDiffImg;
}LGMDType;


#endif

