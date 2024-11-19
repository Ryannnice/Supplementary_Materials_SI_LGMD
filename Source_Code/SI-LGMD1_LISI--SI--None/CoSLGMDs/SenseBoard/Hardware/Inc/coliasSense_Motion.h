#ifndef __COLIASSENCE_MOTION_H
#define __COLIASSENCE_MOTION_H

#include "stm32f4xx_hal.h"
#include "coliasSense_def.h"
#include "coliasSense_BIO.h"
#include "delay.h"

#define CHECK_ID Enable // Enable or DisEnable
#define BASICCHIP_ID 0x78
#define BASIC_WHOAMI 0x34

#pragma pack(4)
typedef enum
{
	DisEnable	= 0,
	Enable		= 1
} ColiasBasicEN_Type;
//Head byte format
#define	WriteWithoutACK	 (0x78)
#define	ACKWithError	 (0x7A)
#define	ACKOK			 (0x7B)
#define	WriteWithACK	 (0x7C)
#define	ReadWithACK		 (0x7D)

//Check ID
typedef struct
{
	uint8_t whoAmI;
	uint8_t ChipID;
} ColiasBasicID_TypeDef;

#define TCRT_Left_IT			(0x80)
#define TCRT_Middle_IT			(0x40)
#define TCRT_Right_IT			(0x20)
#define LightSensor_Left_IT		(0x10)
#define LightSensor_Right_IT	(0x08)
#define Infrared_IT				(0x04)

#define QueueStep_Finish_IT		(0x80)
#define Queue_Empty_IT			(0x40)
#define MotionQueue_Overflow_IT (0x20)
#define CustomMotion_Finish_IT	(0x10)

//interrupt events 
typedef struct
{
	uint8_t SensorIT;
	uint8_t MotionIT;
} IT_Events_TypeDef;


typedef struct 
{
	ColiasBasicEN_Type IT_EN;
	ColiasBasicEN_Type ReservedFinshIT_EN;
	ColiasBasicEN_Type DataTransIT_EN;
	ColiasBasicEN_Type DataTransACK_EN;
} ColiasBasicIT_TypeDef;



//Enum for LED
typedef enum
{
	BlueLED		= 0x10,
	RedLED		= 0x11,
	WhiteLED	= 0x12
} ColiasBasicLED_Type;

//LED
typedef struct
{
	ColiasBasicLED_Type LEDx;
	ColiasBasicEN_Type EN;
	ColiasBasicEN_Type ON;
	ColiasBasicEN_Type Flash;
	ColiasBasicEN_Type State;
	uint8_t FlashPeriod;
} ColiasBasicLED_TypeDef;


//Enum for TCRT
typedef enum 
{
	TCRT_Left	= 0x20,
	TCRT_Middle = 0x24,
	TCRT_Right	= 0x28
} TCRT_Type;

//TCRT
typedef struct
{
	TCRT_Type		TCRTx;
	ColiasBasicEN_Type	EN;
	ColiasBasicEN_Type	IT_EN;
	ColiasBasicEN_Type  Indicator_LED_EN;
	ColiasBasicLED_Type	Indicator_LED;
	ColiasBasicEN_Type	Indicator_LEDStatus;
	uint16_t ThresholdForIT;
	uint16_t ADC_Value;

} ColiasBasicTCRT_TypeDef;

//Enum for Light sensor
typedef enum
{
	Light_Left	= 0x80,
	Light_Right = 0x84
} LightSensor_Type;

//Light sensor
typedef struct
{
	LightSensor_Type Light_x;
	ColiasBasicEN_Type	EN;
	ColiasBasicEN_Type	IT_EN;
	ColiasBasicEN_Type	Indicator_LED_EN;
	ColiasBasicLED_Type	Indicator_LED;
	ColiasBasicEN_Type	Indicator_LEDStatus;
	uint16_t ThresholdForIT;
	uint16_t ADC_Value;

} ColiasBasicLightSensor_TypeDef;

//Infrared
typedef struct
{
	ColiasBasicEN_Type	EN;
	ColiasBasicEN_Type	IT_EN;
	uint16_t ThresholdForIT;
} ColiasBasicInfrared_TypeDef;

//Enum for Default Motion Mode
typedef enum
{
	Stop	= 0,
	Front	= 1,
	Avoid	= 2
} DefaultMotionMode_Type;

//Enum for Sequence direction
typedef enum
{
	TakeFromBottom	= 0,
	TakeFromTop		= 1
} SequenceDirection_Type;

//Enum for Sequence Motion Type
typedef enum SenquenceMotionID
{
	ColiasStop	= 0,							//stop
	ColiasFront	= 1,							//go forward
	ColiasLeft	= 2,							//go left
	ColiasRight	= 3,							//go right
	ColiasBack	= 4,							//go back
	ColiasNew1	= 5,							//quick turn left
	ColiasNew2	= 6,							//quick turn right
	ColiasNew3	= 7,							//got forward slowly
	DefaultMotion = 8
} SenquenceMotionID_Type;

typedef enum WheelDirection
{
	RollFront	= 0,
	RollBack	= 1
} WheelDirection_Type;


//Motion
typedef struct
{
	ColiasBasicEN_Type	EN;
	ColiasBasicEN_Type	IT_EN;
	DefaultMotionMode_Type DefaultMotionMode;
	ColiasBasicEN_Type	StepIT_EN;
	ColiasBasicEN_Type QueueEmptyIT_EN;
	ColiasBasicEN_Type QueueOverflowIT_EN;
	SequenceDirection_Type SequenceDirection;
	ColiasBasicEN_Type F_motion;
	SenquenceMotionID_Type CurrentMotionID;
	uint8_t	EmptyQueueSlotsCount;
	uint8_t MotionPeriod;
	uint8_t Bias;
	ColiasBasicEN_Type BiasBypass;
	ColiasBasicEN_Type QueueAbort;
	ColiasBasicEN_Type QueuePause;

} ColiasBasicMotion_TypeDef;

//Sequence Motion Define
typedef struct
{
	SenquenceMotionID_Type	MotionID;
	uint8_t					LeftWheelPower;
	uint8_t					RightWheelPower;
	int8_t					PeriodOffset;
	WheelDirection_Type		LeftwheelDirection;
	WheelDirection_Type		RightwheelDirection;
	ColiasBasicEN_Type		BiasBypass;
	
} ColiasBasicMotionConf_TypeDef;

typedef enum 
{
	FromQ0 = 0xc0,
	FromQ1 = 0xc1,
	FromQ10 = 0xca
} QueueWritePosition_Type;

typedef struct
{
	SenquenceMotionID_Type QueueMotionIDBuffer[10];
	uint8_t QueueAmountBuffer[10];
	uint8_t AmountOfWirte;
	QueueWritePosition_Type WritePosition;


}ColiasBasicMotionQueue_TypeDef;

//Customize Motion Define
typedef struct
{
	ColiasBasicEN_Type	EN;
	ColiasBasicEN_Type	CustomMotionFinishIT_EN;
	uint8_t				LeftWheelPower;
	uint8_t				RightWheelPower;
	uint16_t			Period;
	WheelDirection_Type LeftwheelDirection;
	WheelDirection_Type RightwheelDirection;
	ColiasBasicEN_Type	BiasBypass;
} ColiasBasicCustomizeMotion_TypeDef;




typedef struct
{
	ColiasBasicID_TypeDef ID;
	ColiasBasicIT_TypeDef IT;
	ColiasBasicLED_TypeDef LED1;
	ColiasBasicLED_TypeDef LED2;
	ColiasBasicLED_TypeDef LED3;
	ColiasBasicTCRT_TypeDef TCRT_Left;
	ColiasBasicTCRT_TypeDef TCRT_Middle;
	ColiasBasicTCRT_TypeDef TCRT_Right;
	ColiasBasicLightSensor_TypeDef Light_Left;
	ColiasBasicLightSensor_TypeDef Light_Right;
	ColiasBasicInfrared_TypeDef Infrared;
	ColiasBasicMotion_TypeDef Motion;


} ColiasBasicInit_TypeDef;

typedef struct
{
	uint8_t right_it;
	uint8_t left_it;
	uint8_t middle_it;
}Colias_Bumper_Inter;

typedef struct
{
	CoS_CCSTypeDef CC;		//common controller
	CoS_UART_HandleTypeDef* hHUART;
	ColiasBasicEN_Type	ACK;
	IT_Events_TypeDef IT_Events;
	uint8_t timecount;
	uint8_t MotionFinishFlag;
	Colias_Bumper_Inter bumper_inter;
	
}CoS_Motion_TypeDef;

typedef struct
{
	CoS_Motion_TypeDef* Instance;
	CoS_UART_HandleTypeDef* hHUART;
	
	ColiasBasicEN_Type	ACK;
	ColiasBasicID_TypeDef ID;
	ColiasBasicLED_TypeDef LED1;
	ColiasBasicLED_TypeDef LED2;
	ColiasBasicLED_TypeDef LED3;
	ColiasBasicTCRT_TypeDef TCRT_left;
	ColiasBasicTCRT_TypeDef TCRT_middle;
	ColiasBasicTCRT_TypeDef TCRT_right;
	ColiasBasicLightSensor_TypeDef LightSensor_left;
	ColiasBasicLightSensor_TypeDef LightSensor_right;
	ColiasBasicInfrared_TypeDef Infrared;
	ColiasBasicMotion_TypeDef Motion;
	ColiasBasicMotionConf_TypeDef Conf_DefaultMotion;
	ColiasBasicMotionConf_TypeDef Conf_Front;
	ColiasBasicMotionConf_TypeDef Conf_Left;
	ColiasBasicMotionConf_TypeDef Conf_Right;
	ColiasBasicMotionConf_TypeDef Conf_Back;
	ColiasBasicMotionConf_TypeDef Conf_New1;
	ColiasBasicMotionConf_TypeDef Conf_New2;
	ColiasBasicMotionConf_TypeDef Conf_New3;

}CoS_Motion_InitTypeDef;


typedef struct
{
	CoS_Motion_TypeDef* Instance; //instance address
	CoS_LockTypeDef Lock; 	//Locking object
	__IO CoS_StatusTypeDef iState; //init state
	//CoS_UART_HandleTypeDef* hHUART;

}CoS_Motion_HandleTypeDef;


CoS_StatusTypeDef ColiasBasic_Init(ColiasBasicInit_TypeDef *ColiasBasicInit_struct);
CoS_StatusTypeDef ColiasBasicID_Check(CoS_Motion_TypeDef *hObj, ColiasBasicID_TypeDef *ID);
CoS_StatusTypeDef ColiasBasicIT_Enable(CoS_Motion_TypeDef *hObj);
CoS_StatusTypeDef ColiasBasicIT_DisEnable(CoS_Motion_TypeDef *hObj);
CoS_StatusTypeDef ColiasBasicITEvents_Check(CoS_Motion_TypeDef *hObj, IT_Events_TypeDef *IT_Events);
CoS_StatusTypeDef ColiasBasicITEvents_Clean(CoS_Motion_TypeDef *hObj, IT_Events_TypeDef *IT_Events);
CoS_StatusTypeDef ColiasBasicLED_Config(CoS_Motion_TypeDef *hObj, ColiasBasicLED_TypeDef *LEDConfig);
CoS_StatusTypeDef ColiasBasicTCRT_Config(CoS_Motion_TypeDef *hObj, ColiasBasicTCRT_TypeDef *TCRTConfig);
uint16_t ColiasBasicTCRTADC_Get(CoS_Motion_TypeDef *hObj, uint8_t X);
CoS_StatusTypeDef ColiasBasicLightS_Config(CoS_Motion_TypeDef *hObj, ColiasBasicLightSensor_TypeDef *LightS_Config);
CoS_StatusTypeDef ColiasBasicLightSADC_Get(CoS_Motion_TypeDef *hObj, ColiasBasicLightSensor_TypeDef *LightS_Config);
CoS_StatusTypeDef ColiasBasicInfrared_Config(CoS_Motion_TypeDef *hObj, ColiasBasicInfrared_TypeDef *Infrared_Config);
CoS_StatusTypeDef ColiasBasicMotionControler_Config(CoS_Motion_TypeDef *hObj, ColiasBasicMotion_TypeDef *Motion_Config);
CoS_StatusTypeDef ColiasBasicMotion_Config(CoS_Motion_TypeDef *hObj, ColiasBasicMotionConf_TypeDef *Sequence_Config);
CoS_StatusTypeDef ColiasBasic_WriteMotionQueue(CoS_Motion_TypeDef *hObj, ColiasBasicMotionQueue_TypeDef *WriteMotionQueue);
CoS_StatusTypeDef ColiasBasic_RunCustomizeMotion(CoS_Motion_TypeDef *hObj, ColiasBasicCustomizeMotion_TypeDef *CustomizeMotion);

///////////////////////
CoS_StatusTypeDef ColiasBasicLED_ReadConfig(CoS_Motion_TypeDef *hObj, ColiasBasicLED_TypeDef *LEDConfig);
CoS_StatusTypeDef CoS_MotionParameter_Init(CoS_Motion_TypeDef *hObj, CoS_Motion_InitTypeDef* hInit);


CoS_StatusTypeDef WriteBytesWithACKToColiasBasic(CoS_UART_HandleTypeDef *hObj, uint8_t *pData, uint8_t Size);
CoS_StatusTypeDef WriteBytesWithoutACKToColiasBasic(CoS_UART_HandleTypeDef *hObj, uint8_t *pData, uint8_t Size);
CoS_StatusTypeDef ReadBytesWithACKFromColiasBasic(CoS_UART_HandleTypeDef *hObj, uint8_t Address, uint8_t * pData, uint8_t Size);

CoS_StatusTypeDef CoS_Motion_Init(CoS_Motion_HandleTypeDef *hObj, CoS_Motion_InitTypeDef* hInit);
CoS_StatusTypeDef MotionTask(CoS_Motion_TypeDef *hObj);
void CoS_MotionEXIT_Init(void);
u8 fast_motion_control(CoS_Motion_TypeDef *hObj,u8 motionID,u8 mode);
CoS_StatusTypeDef BitSetColiasBasicReg(CoS_UART_HandleTypeDef *hObj, uint8_t Address, uint8_t mask,uint8_t mode);
//CoS_StatusTypeDef ColiasBasicTCRTADC_Get(CoS_Motion_HandleTypeDef *hObj, ColiasBasicTCRT_TypeDef *TCRTConfig);
#endif // !__CORRESPOND_H
