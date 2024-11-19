#ifndef __CoSEXBD_H
#define __CoSEXBD_H

#include "stm32f4xx_hal.h"
#include "coliasSense_def.h"

typedef enum 
{
	Exbd_none =0x00u,
	Exbd_USB = 0x01u,
	Exbd_BT = 0x02u,
	Exbd_FPV = 0x03u,
	Exbd_SWARM = 0x04u
}CoS_EXBD_NUMTypeDef;

typedef enum 
{
	ExPD_VALUE = 0x00u,
	ExPD_ON =0x01u,
	ExPD_OFF = 0x02u,
	ExPD_TOGGLE = 0x03u
}CoS_EXBD_PDCONTROL;


typedef struct{
	CoS_CCSTypeDef CC;
	CoS_EXBD_NUMTypeDef current_extend;
}CoS_EXBD_TypeDef;


typedef struct 
{
	CoS_EXBD_TypeDef * Instance;
	CoS_LockTypeDef Lock; 
	__IO CoS_StatusTypeDef iState;
}CoS_EXBD_HandleTypeDef;

typedef struct{
	CoS_EXBD_TypeDef* Instance;

}CoS_EXBD_InitTypeDef;

CoS_StatusTypeDef CoS_EXBD_Init(CoS_EXBD_HandleTypeDef *hObj,
								CoS_EXBD_InitTypeDef* Init);
CoS_StatusTypeDef CoS_EXBD_PDControl(CoS_EXBD_HandleTypeDef *hObj,u8 pins,CoS_EXBD_PDCONTROL ctrl);
#endif