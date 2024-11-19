#ifndef __CoSCam_H
#define __CoSCam_H

#include "stm32f4xx_hal.h"
#include "coliasSense_def.h"
#include "coliasSense_BIO.h"
#include "coliasSense_interfaces.h"

#define ID_OV7670 0x73
#define OV7670_FPS 30
#define A_OV7670W 0x42
#define A_OV7670R 0x43

typedef struct 
{
	CoS_CCSTypeDef CC;		//common controller
	uint32_t cam_DMACount;
	uint32_t cam_hFrameCountount;
}CoS_CAM_TypeDef;

typedef struct{
	CoS_CAM_TypeDef* Instance;
	CoS_sI2C_HandleTypeDef* hSSCB;
	DMA_HandleTypeDef* hdma_dcmi;
	DCMI_HandleTypeDef* hdcmi;
	void* DMA_Base;
	uint32_t DMA_Length;
}CoS_CAM_InitTypeDef;

typedef struct 
{
	CoS_CAM_TypeDef* Instance; //instance address
	CoS_sI2C_HandleTypeDef* hSSCB;
	DCMI_HandleTypeDef* hdcmi;
	CoS_LockTypeDef Lock; 	//Locking object
	__IO CoS_StatusTypeDef iState; //init state
}CoS_CAM_HandleTypeDef;




CoS_StatusTypeDef CoS_CAM_Init(CoS_CAM_HandleTypeDef *hObj,CoS_CAM_InitTypeDef* hInit);
CoS_StatusTypeDef CoS_CAM_cameraReset(CoS_CAM_HandleTypeDef *hObj);
CoS_StatusTypeDef CoS_CAM_regReset(CoS_CAM_HandleTypeDef *hObj);
uint8_t CoS_CAM_ReadID(CoS_CAM_HandleTypeDef *hObj);
CoS_StatusTypeDef CoS_CAM_Start(CoS_CAM_HandleTypeDef *hObj);
CoS_StatusTypeDef CoS_CAM_DMASet(CoS_CAM_HandleTypeDef *hObj,uint32_t dmabase,uint32_t length);
CoS_StatusTypeDef CoS_CAM_Stop(CoS_CAM_HandleTypeDef *hObj);
void CoS_camera_StateReport(CoS_CAM_HandleTypeDef *hObj,u8 trans_mode_blocking);
#endif
