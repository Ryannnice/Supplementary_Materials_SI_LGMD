#ifndef __CoSCLED_H
#define __CoSCLED_H

#include "stm32f4xx_hal.h"
#include "coliasSense_def.h"


#define CLED_ChannelB TIM_CHANNEL_1
#define CLED_ChannelG TIM_CHANNEL_2
#define CLED_ChannelR TIM_CHANNEL_3

typedef enum
{
	channel_Red=TIM_CHANNEL_3,
	channel_Green=TIM_CHANNEL_2,
	channel_Blue=TIM_CHANNEL_1,
}CoS_CLED_ChannelTypeDef;

typedef enum
{
	CLED_OFF=0U,
	CLED_Random=1U,
	CLED_Breath=2U,
	CLED_R=3U,
	CLED_G=4U,
	CLED_B=5U,
	CLED_Freeze=6U
}CLED_ModeType;

typedef struct
{
	u16 A;
	u16 shift;
	float freq;
	float phy;
	u16 val;
}CLED_ColorStructType;

typedef struct
{
	CoS_CCSTypeDef CC;		//common controller
	CLED_ColorStructType Rs;
	CLED_ColorStructType Gs;
	CLED_ColorStructType Bs;
}CoS_CLED_TypeDef;

typedef struct 
{
	CoS_CLED_TypeDef * Instance; //CLED  instance address
	CoS_LockTypeDef Lock; 	//Locking object
	__IO CoS_StatusTypeDef iState; //init state
	TIM_HandleTypeDef* htim1;
	RNG_HandleTypeDef* hrng;
}CoS_CLED_HandleTypeDef;



typedef struct{
	CoS_CLED_TypeDef* Instance;
	RNG_HandleTypeDef* hrng;
	TIM_HandleTypeDef* htim1;
}CoS_CLED_InitTypeDef;

//CoS_CLED_TypeDef CLED1;
//CoS_CCSTypeDef CC_CLED1;
//CoS_CLED_HandleTypeDef hCLED1;





CoS_StatusTypeDef CoS_CLED_Init(CoS_CLED_HandleTypeDef *hObj,CoS_CLED_InitTypeDef* Init);
void CoS_CLED_StateReport(CoS_CLED_HandleTypeDef *hObj);
CoS_StatusTypeDef CoS_CLED_DeInit(CoS_CLED_HandleTypeDef *hObj);
void CoS_CLED_ErrorCallback(CoS_CLED_HandleTypeDef *hObj);
void CoS_LED_Update(CoS_CLED_HandleTypeDef *hObj);
CoS_StatusTypeDef CoS_CLED_Enable(CoS_CLED_HandleTypeDef* hObj);
CoS_StatusTypeDef CoS_CLED_Disable(CoS_CLED_HandleTypeDef* hObj);
void CoS_CLED_setRGB(u16 R, u16 G, u16 B);
#endif
