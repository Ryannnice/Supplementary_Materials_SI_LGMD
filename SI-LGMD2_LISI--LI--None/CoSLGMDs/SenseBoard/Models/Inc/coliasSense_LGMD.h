/*
 * Description: LGMD-Plus
 * File: coliasSense_LGMD.h
 * Author: Qinbing Fu
 * Date: 2019 Aug
*/

#ifndef __CoSLGMD_H
#define __CoSLGMD_H

#include "stm32f4xx_hal.h"
#include "coliasSense_def.h"
#include "LGMD_model.h"
#include "coliasSense_BIO.h"
#include "coliasSense_Camera.h"
#include <stdint.h>

#define Image_Width 99
#define Image_Height 72

typedef struct{
	uint32_t processCount;
	uint32_t processCountLast;
	uint32_t* hFrameCount;
	LGMDType* Model;
	uint8_t Enable;
	uint8_t status;
	uint8_t processRate;
	uint8_t currentImage;
	uint8_t currentDiffImage;
	uint8_t AGC_enable_period;
}LGMD_pControlTypedef;

//attention
/******************************************************/
extern uint16_t Image[3][Image_Height][Image_Width];
extern int8_t Diff_Image[2][Image_Height][Image_Width];
extern LGMD_pControlTypedef hLGMD;
/******************************************************/


uint8_t LGMD_Param_Init(LGMD_pControlTypedef* hLGMD);
uint8_t Calc_LGMDs_OutOfChannel(LGMD_pControlTypedef* hLGMD);
uint8_t Fast_LGMD(LGMD_pControlTypedef* hLGMD);
uint8_t Decision_making(LGMD_pControlTypedef* hLGMD,uint8_t allow_motion);
uint8_t LGMD_demo(LGMD_pControlTypedef* hLGMD);


#endif