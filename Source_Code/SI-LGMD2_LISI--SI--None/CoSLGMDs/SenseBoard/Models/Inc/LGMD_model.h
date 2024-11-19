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
#define H_Start 15
#define W_Start 15


typedef struct
{
	uint8_t ON[2][Layer_Height][Layer_Width];	  //ON Cells
	uint8_t OFF[2][Layer_Height][Layer_Width];	//OFF Cells
	//uint8_t SI_ON[Layer_Height][Layer_Width];	  //SI layer
	//uint8_t SI_OFF[Layer_Height][Layer_Width];
	uint8_t ON_hat[Layer_Height][Layer_Width];
	uint8_t OFF_hat[Layer_Height][Layer_Width];
}LGMDstruct_Layer;

typedef struct
{
	uint8_t Tpm;            // threshold in PM pathway
	uint8_t Cw;             //constant to calculate w
	uint8_t CdeTde;         //decay threshold in G layer
	uint8_t flag;
	uint8_t n;
	uint8_t a4;             //scale coefficient in spiking
	uint8_t clip;
	float Delatc;           //small real number 
	float a1;               //coefficient in half-wave rectifying Eq
	float interval;         //the time interval between the consecutive frames of digital signals
	
	//the delay of L_inhibition
	float on_center_tau;    //delay time of the center in ON channel
	float on_near_tau;      //delay time of the near in ON channel
	float on_diag_tau;      //delay time of the diag in ON channel
	float off_center_tau;   //delay time of the center in OFF channel
	float off_near_tau;     //delay time of the near in OFF channel
	float off_diag_tau;     //delay time of the diag in OFF channel
  float on_delay_center;  //delay coefficient of the center in ON channel
  float on_delay_near;    //delay coefficoent of the near in ON channel
	float on_delay_diag;    //delay coefficient of the diag in ON channel
	float off_delay_center; //delay coefficient of the center in OFF channel
	float off_delay_near;   //delay coefficient of the near in OFF channel
	float off_delay_diag;   //delay coefficient of the diag in OFF channel
	float on_tau;
	float off_tau;
	float on_delay;
	float off_delay;
	float beta1;
	float beta2;
	float beta3;
	float beta4;
	float ffi_tau;          //latency in PM pathway
	float a2;               //delay coefficient in Adaptive Inhibition Mechanism
	float Won_base;         //bias baseline in ON channel
	float Woff_base;        //bias baseline in OFF channel
	float o1;               //term coefficient in S layer
	float o2;               //term coefficient in S layer
	float o3;               //term coefficient in S layer
	float Kspi;             //coefficient in sigmoid function
	float Tspi;             //spiking threshold
	float Tsfa;             //small threshold in SFA
	float sfa_tau;          //time constant in SFA
	float a3;               //a coefficient that indicates adaption rate to visual movements
	float w1;               //time-varying local biase to control the intensity of inhibitory flows in ON channel
	float w2;               //time-varying local biase to control the intensity of inhibitory flows in OFF channel
	float r1;
	float r2;
}LGMDstruct_Params; 

typedef struct
{
	uint8_t SPIKE;
	float MP;	//membrane potential
	float FFI_out[2];	//FFI output
	float SMP[2];	//sigmoid membrane potential
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

