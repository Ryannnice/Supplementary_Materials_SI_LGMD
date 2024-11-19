


#include "coliasSense_board.h"
#include "coliasSense_Motion.h"
#include "delay.h"
#include <string.h>
#include <stdlib.h>
static CoS_Motion_HandleTypeDef* hHMotion;

static CoS_StatusTypeDef motion_reg_init(CoS_Motion_HandleTypeDef *hObj, CoS_Motion_InitTypeDef* hInit)
{
	hObj->Instance->ACK = Enable;
	hObj->Instance->hHUART = hInit->hHUART;
	return 0;
}

CoS_StatusTypeDef CoS_Motion_Init(CoS_Motion_HandleTypeDef *hObj, CoS_Motion_InitTypeDef* hInit)
{
	CoS_StatusTypeDef status;

	/* Check the handle allocation */
	if (hObj == NULL ||
		hInit->hHUART == NULL)
		return CoS_STATE_ERROR;
	__CoS_LOCK(hObj);
	if (hObj->iState == CoS_STATE_UNSET)
	{
		/* Allocate lock resource and initialize it */
		hObj->Lock = CoS_UNLOCKED;
		hObj->Instance = hInit->Instance;
	}
	/* Change peripheral state */
	hObj->iState = CoS_STATE_BUSY;
	if (hInit->hHUART->iState!=CoS_STATE_OK)//must init related UART before 
	{
		printf("need to init related UART before initing motion\r\n");
		__CoS_UNLOCK(hObj);
		return CoS_STATE_ERROR;
	}
	/* Init the low level hardware */
	motion_reg_init(hObj,hInit);
//	status = ColiasBasicID_Check(hObj->ID);
	CoS_MotionEXIT_Init();
	status = CoS_MotionParameter_Init(hObj->Instance, hInit);
	if (status==CoS_STATE_OK)
	{
		/* Initialize the Instance state */
		hObj->iState = CoS_STATE_OK;
		__CoS_UNLOCK(hObj);
		/* Return function status */
		return status;
	}
	else
	{
		hObj->iState = CoS_STATE_UNSET;
		__CoS_UNLOCK(hObj);
		/* Return function status */
		return status;
	}
}



CoS_StatusTypeDef CoS_MotionParameter_Init(CoS_Motion_TypeDef *hObj, CoS_Motion_InitTypeDef* hInit)
{

	CoS_StatusTypeDef status;
	u8 trialmax=20;
	u8 trial=0;
	//uint8_t RxBuffer[140]={0};

	

	ColiasBasicID_TypeDef ID;
	IT_Events_TypeDef IT_Events;
	ColiasBasicLED_TypeDef MotionLED;
	ColiasBasicTCRT_TypeDef TCRT_CC;
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

#ifndef CHECK_ID

	status = ColiasBasicID_Check(hObj->ID);
	if (status != CoS_STATE_OK)
		return status;


	if (ID.ChipID != BASICCHIP_ID || ID.whoAmI != BASIC_WHOAMI)
		return CoS_STATE_ERROR;

#endif 

	hObj->timecount = 0;
	hObj->MotionFinishFlag = 1;

	//HAL_Delay(2000);


	memset(&MotionLED, 0, sizeof(MotionLED));
	MotionLED.LEDx = BlueLED;
	MotionLED.EN = Enable;
	//status = ColiasBasicLED_Config(hObj, &(MotionLED1));
	do
	{
		printf("{  Motion: Configuring LED1...\r\n");
		LED1_Toggle;
		status = ColiasBasicLED_Config(hObj, &(MotionLED));
		trial++;
	}while (status != CoS_STATE_OK &&trial<trialmax);

	memset(&MotionLED, 0, sizeof(MotionLED));
	MotionLED.LEDx = RedLED;
	MotionLED.EN = Enable;
	MotionLED.ON = Enable;
	do
	{
		printf("{  Motion: Configuring LED2...\r\n");
		LED1_Toggle;
		status = ColiasBasicLED_Config(hObj, &(MotionLED));
		trial++;
	}while (status != CoS_STATE_OK &&trial<trialmax);

	memset(&MotionLED, 0, sizeof(MotionLED));
	MotionLED.LEDx = WhiteLED;
	MotionLED.EN = Enable;
	MotionLED.Flash = Enable;
	MotionLED.FlashPeriod = 5;
	do
	
	{
		printf("{  Motion: Configuring LED3...\r\n");
		LED1_Toggle;
		status = ColiasBasicLED_Config(hObj, &(MotionLED));
		trial++;
	}while (status != CoS_STATE_OK &&trial<trialmax);

	memset(&TCRT_CC, 0, sizeof(TCRT_CC));
	TCRT_CC.TCRTx = TCRT_Left;
	TCRT_CC.EN = Enable;
	TCRT_CC.IT_EN = Enable;
	TCRT_CC.Indicator_LED_EN = Enable;
	TCRT_CC.Indicator_LED = BlueLED;
	TCRT_CC.ThresholdForIT = 0x300;
	do
	{
		printf("{  Motion: Configuring TCRT_CC...\r\n");
		LED1_Toggle;
		status = ColiasBasicTCRT_Config(hObj, &(TCRT_CC));
		trial++;
	}while (status != CoS_STATE_OK &&trial<trialmax);

	memset(&TCRT_CC, 0, sizeof(TCRT_CC));
	TCRT_CC.TCRTx = TCRT_Middle;
	TCRT_CC.EN = Enable;
	TCRT_CC.IT_EN = Enable;
	TCRT_CC.Indicator_LED_EN = Enable;
	TCRT_CC.Indicator_LED = BlueLED;
	TCRT_CC.ThresholdForIT = 0x300;
	do
	{
		printf("{  Motion: Configuring TCRT_CC...\r\n");
		LED1_Toggle;
		status = ColiasBasicTCRT_Config(hObj, &(TCRT_CC));
		trial++;
	}while (status != CoS_STATE_OK &&trial<trialmax);

	memset(&TCRT_CC, 0, sizeof(TCRT_CC));
	TCRT_CC.TCRTx = TCRT_Right;
	TCRT_CC.EN = Enable;
	TCRT_CC.IT_EN = Enable;
	TCRT_CC.Indicator_LED_EN = Enable;
	TCRT_CC.Indicator_LED = BlueLED;
	TCRT_CC.ThresholdForIT = 0x300;
	do
	{
		printf("{  Motion: Configuring TCRT_CC...\r\n");
		LED1_Toggle;
		status = ColiasBasicTCRT_Config(hObj, &(TCRT_CC));
		trial++;
	}while (status != CoS_STATE_OK &&trial<trialmax);

	memset(&LightSensor_left, 0, sizeof(LightSensor_left));
	LightSensor_left.Light_x = Light_Left;
	LightSensor_left.EN = Enable;
	LightSensor_left.IT_EN = Enable;
	LightSensor_left.Indicator_LED = RedLED;
	LightSensor_left.ThresholdForIT = 0x0FFF;
	do
	{
		printf("{  Motion: Configuring LightSensor_left...\r\n");
		LED1_Toggle;
		status = ColiasBasicLightS_Config(hObj, &(LightSensor_left));
		trial++;
	}while (status != CoS_STATE_OK &&trial<trialmax);

	memset(&LightSensor_right, 0, sizeof(LightSensor_right));
	LightSensor_right.Light_x = Light_Right;
	LightSensor_right.EN = Enable;
	LightSensor_right.IT_EN = Enable;
	LightSensor_right.Indicator_LED = RedLED;
	LightSensor_right.ThresholdForIT = 0x0FFF;
	do
	{
		printf("{  Motion: Configuring LightSensor_right...\r\n");
		LED1_Toggle;
		status = ColiasBasicLightS_Config(hObj, &(LightSensor_right));
		trial++;
	}while (status != CoS_STATE_OK &&trial<trialmax);

	memset(&Infrared, 0, sizeof(Infrared));
	Infrared.EN = Enable;
	Infrared.IT_EN = Enable;
	Infrared.ThresholdForIT = 0x0700;
	do
	{
		printf("{  Motion: Configuring Infrared...\r\n");
		LED1_Toggle;
		status = ColiasBasicInfrared_Config(hObj, &(Infrared));
		trial++;
	}while (status != CoS_STATE_OK &&trial<trialmax);

	memset(&Motion, 0, sizeof(Motion));
	Motion.EN = Enable;
	/**********************************/
	Motion.DefaultMotionMode = Stop;
	/**********************************/
	Motion.IT_EN = Enable;
	Motion.StepIT_EN = Enable;
	Motion.QueueEmptyIT_EN = Enable;
	Motion.QueueOverflowIT_EN = Enable;
	Motion.SequenceDirection = TakeFromBottom;
	Motion.MotionPeriod = 0x64;// unit: 100ms
	Motion.Bias = 128;
	do
	{
		printf("{  Motion: Configuring Motion settings...\r\n");
		LED1_Toggle;
		status = ColiasBasicMotionControler_Config(hObj, &(Motion));
		trial++;
	}while (status != CoS_STATE_OK &&trial<trialmax);
	/**********************************/
	memset(&Conf_DefaultMotion, 0, sizeof(Conf_DefaultMotion));
	Conf_DefaultMotion.MotionID = DefaultMotion;
	Conf_DefaultMotion.LeftWheelPower = 77; //   30/256 battery power
	Conf_DefaultMotion.RightWheelPower = 80;
	Conf_DefaultMotion.BiasBypass = DisEnable;
	Conf_DefaultMotion.LeftwheelDirection = RollFront;
	Conf_DefaultMotion.RightwheelDirection = RollFront;
	/**********************************/
	do
	{
		printf("{  Motion: Configuring Conf_DefaultMotion...\r\n");
		LED1_Toggle;
		status = ColiasBasicMotion_Config(hObj, &(Conf_DefaultMotion));
		trial++;
	}while (status != CoS_STATE_OK &&trial<trialmax);

	memset(&Conf_Front, 0, sizeof(Conf_Front));
	Conf_Front.MotionID = ColiasFront;
	Conf_Front.LeftWheelPower = 77; //   80/256 battery power
	Conf_Front.RightWheelPower = 80;
	Conf_Front.PeriodOffset = 250;
	Conf_Front.BiasBypass = DisEnable;
	Conf_Front.LeftwheelDirection = RollFront;
	Conf_Front.RightwheelDirection = RollFront;
	do
	{
		printf("{  Motion: Configuring Conf_Front...\r\n");
		LED1_Toggle;
		status = ColiasBasicMotion_Config(hObj, &(Conf_Front));
		trial++;
	}while (status != CoS_STATE_OK &&trial<trialmax);

	memset(&Conf_Left, 0, sizeof(Conf_Left));
	Conf_Left.MotionID = ColiasLeft;
	Conf_Left.LeftWheelPower = 90; //   90/256 battery power
	Conf_Left.RightWheelPower = 60;
	Conf_Left.PeriodOffset = 250;
	Conf_Left.BiasBypass = DisEnable;
	Conf_Left.LeftwheelDirection = RollBack;
	Conf_Left.RightwheelDirection = RollFront;
	do
	{
		printf("{  Motion: Configuring Conf_Left...\r\n");
		LED1_Toggle;
		status = ColiasBasicMotion_Config(hObj, &(Conf_Left));
		trial++;
	}while (status != CoS_STATE_OK &&trial<trialmax);

	memset(&Conf_Right, 0, sizeof(Conf_Right));
	Conf_Right.MotionID = ColiasRight;
	Conf_Right.LeftWheelPower = 60; //   30/256 battery power
	Conf_Right.RightWheelPower = 90;
	Conf_Right.PeriodOffset = 250;
	Conf_Right.BiasBypass = DisEnable;
	Conf_Right.LeftwheelDirection = RollFront;
	Conf_Right.RightwheelDirection = RollBack;
	do
	{
		printf("{  Motion: Configuring Conf_Right...\r\n");
		LED1_Toggle;
		status = ColiasBasicMotion_Config(hObj, &(Conf_Right));
		trial++;
	}while (status != CoS_STATE_OK &&trial<trialmax);

	memset(&Conf_Back, 0, sizeof(Conf_Back));
	Conf_Back.MotionID = ColiasBack;
	Conf_Back.LeftWheelPower = 50; //   50/256 battery power
	Conf_Back.RightWheelPower = 50;
	Conf_Back.PeriodOffset = 100;
	Conf_Back.BiasBypass = DisEnable;
	Conf_Back.LeftwheelDirection = RollBack;
	Conf_Back.RightwheelDirection = RollBack;
	do
	{
		printf("{  Motion: Configuring Conf_Back...\r\n");
		LED1_Toggle;
		status = ColiasBasicMotion_Config(hObj, &(Conf_Back));
		trial++;
	}while (status != CoS_STATE_OK &&trial<trialmax);

	memset(&Conf_New1, 0, sizeof(Conf_New1));
	Conf_New1.MotionID = ColiasNew1;
	Conf_New1.LeftWheelPower = 90; //   150/256 battery power
	Conf_New1.RightWheelPower = 10;
	Conf_New1.PeriodOffset = 200;
	Conf_New1.BiasBypass = DisEnable;
	Conf_New1.LeftwheelDirection = RollBack;
	Conf_New1.RightwheelDirection = RollBack;
	do
	{
		printf("{  Motion: Configuring Motion Conf_New1...\r\n");
		LED1_Toggle;
		status = ColiasBasicMotion_Config(hObj, &(Conf_New1));
		trial++;
	}while (status != CoS_STATE_OK &&trial<trialmax);

	memset(&Conf_New2, 0, sizeof(Conf_New2));
	Conf_New2.MotionID = ColiasNew2;
	Conf_New2.LeftWheelPower = 60; //   30/256 battery power
	Conf_New2.RightWheelPower = 160;
	Conf_New2.PeriodOffset = 200;
	Conf_New2.BiasBypass = DisEnable;
	Conf_New2.LeftwheelDirection = RollFront;
	Conf_New2.RightwheelDirection = RollBack;
	do
	{
		printf("{  Motion: Configuring Motion Conf_New2...\r\n");
		LED1_Toggle;
		status = ColiasBasicMotion_Config(hObj, &(Conf_New2));
		trial++;
	}while (status != CoS_STATE_OK &&trial<trialmax);

	memset(&Conf_New3, 0, sizeof(Conf_New3));
	Conf_New3.MotionID = ColiasNew3;
	Conf_New3.LeftWheelPower = 27; //   20/256 battery power
	Conf_New3.RightWheelPower = 30;
	Conf_New3.PeriodOffset = 0;
	Conf_New3.BiasBypass = DisEnable;
	Conf_New3.LeftwheelDirection = RollFront;
	Conf_New3.RightwheelDirection = RollFront;
	do
	{
		printf("{  Motion: Configuring Motion Conf_New3...\r\n");
		LED1_Toggle;
		status = ColiasBasicMotion_Config(hObj, &(Conf_New3));
		trial++;
	}while (status != CoS_STATE_OK &&trial<trialmax);
	if (trial>=trialmax)
	{
		return CoS_STATE_ERROR;
	}

	//status = ReadBytesWithACKFromColiasBasic(hObj->hHUART, 0x00, RxBuffer, 50);


	//status = ReadBytesWithACKFromColiasBasic(hObj->hHUART, 0x80, (RxBuffer+50), 70);
	status = ColiasBasicIT_Enable(hObj);
	while (status != CoS_STATE_OK)
	{
		LED1_Toggle;
		status = ColiasBasicIT_Enable(hObj);
	}

	HAL_NVIC_EnableIRQ(EXTI15_10_IRQn);
	ColiasBasicITEvents_Check((hCoS->hHCos_Motion->Instance), &(hCoS->hHCos_Motion->Instance->IT_Events));
	//Bumper interruption flag init
	//hObj->bumper_inter.left_it=0;
	//hObj->bumper_inter.right_it=0;
	//hObj->bumper_inter.middle_it=0;

	return CoS_STATE_OK;


}






CoS_StatusTypeDef ColiasBasicID_Check(CoS_Motion_TypeDef *hObj, ColiasBasicID_TypeDef *ID)
{
	CoS_StatusTypeDef status;
	uint8_t RxBuffer[10] = { 0 };

	status = ReadBytesWithACKFromColiasBasic(hObj->hHUART, 0x00, RxBuffer, 2);
	if (status != CoS_STATE_OK)
		return	status;

	ID->whoAmI = RxBuffer[1];
	ID->ChipID = RxBuffer[2];

	return CoS_STATE_OK;
}


CoS_StatusTypeDef ColiasBasicIT_Enable(CoS_Motion_TypeDef *hObj)
{
	uint8_t TxBuffer[2];
	TxBuffer[0] = 0x02;
	TxBuffer[1] =0x80;
	CoS_StatusTypeDef status;

	if (hObj->ACK)
	{
		status = WriteBytesWithACKToColiasBasic(hObj->hHUART, TxBuffer, 2);
	}
	else
	{
		status = WriteBytesWithoutACKToColiasBasic(hObj->hHUART, TxBuffer, 2);
	}
	if (status != CoS_STATE_OK)
		return	status;

	return CoS_STATE_OK;
}

CoS_StatusTypeDef ColiasBasicIT_DisEnable(CoS_Motion_TypeDef *hObj)
{
	uint8_t TxBuffer[2] ;
	TxBuffer[0] = 0x02;
	TxBuffer[1] =0x00;
	CoS_StatusTypeDef status;
	if (hObj->ACK)
	{
		status = WriteBytesWithACKToColiasBasic(hObj->hHUART, TxBuffer, 2);
	}
	else
	{
		status = WriteBytesWithoutACKToColiasBasic(hObj->hHUART, TxBuffer, 2);
	}
	if (status != CoS_STATE_OK)
		return	status;

	return CoS_STATE_OK;
}



CoS_StatusTypeDef ColiasBasicITEvents_Check(CoS_Motion_TypeDef *hObj, IT_Events_TypeDef *IT_Events)
{
	uint8_t	RxBuffer[10] = { 0 };
	CoS_StatusTypeDef status;

	status = ReadBytesWithACKFromColiasBasic(hObj->hHUART, 0x03, RxBuffer, 2);
	if (status != CoS_STATE_OK)
		return status;

	IT_Events->SensorIT = RxBuffer[1];
	IT_Events->MotionIT = RxBuffer[2];

	return CoS_STATE_OK;
}

CoS_StatusTypeDef ColiasBasicITEvents_Clean(CoS_Motion_TypeDef *hObj, IT_Events_TypeDef *IT_Events)
{
	uint8_t	TxBuffer[3] = { 0 };
	CoS_StatusTypeDef status;

	TxBuffer[0] = 0x03;
	if (hObj->ACK)
	{
		status = WriteBytesWithACKToColiasBasic(hObj->hHUART, TxBuffer, 3);
	}
	else
	{
		status = WriteBytesWithoutACKToColiasBasic(hObj->hHUART, TxBuffer, 3);
	}
	if (status != CoS_STATE_OK)
		return	status;

	return CoS_STATE_OK;
}


CoS_StatusTypeDef ColiasBasicLED_Config(CoS_Motion_TypeDef *hObj, ColiasBasicLED_TypeDef *LEDConfig)
{
	uint8_t TxBuffer[2] = { 0 };
	CoS_StatusTypeDef status;

	TxBuffer[0] = LEDConfig->LEDx;

	if (LEDConfig->EN == Enable)
	{
		TxBuffer[1] |= 0x80;
	}
	else
	{
		TxBuffer[1] &= ~0x80;
	}
	if (LEDConfig->ON == Enable)
	{
		TxBuffer[1] |= 0x40;
	}
	else
	{
		TxBuffer[1] &= ~0x40;
	}

	if (LEDConfig->Flash == Enable)
	{
		TxBuffer[1] |= 0x20;
	}
	else
	{
		TxBuffer[1] &= ~0x20;
	}

	if (LEDConfig->FlashPeriod >= 0x0f)
	{
		TxBuffer[1] |= 0x0f;
	}
	else
	{
		TxBuffer[1] &= ~0x0f;
		TxBuffer[1] |= ((LEDConfig->FlashPeriod) & 0x0f);
	}

	if (hObj->ACK)
	{
		status = WriteBytesWithACKToColiasBasic(hObj->hHUART, TxBuffer, 2);
	}
	else
	{
		status = WriteBytesWithoutACKToColiasBasic(hObj->hHUART, TxBuffer, 2);
	}
	if (status != CoS_STATE_OK)
		return	status;


	return CoS_STATE_OK;
}

CoS_StatusTypeDef ColiasBasicLED_ReadConfig(CoS_Motion_TypeDef *hObj, ColiasBasicLED_TypeDef *LEDConfig)
{
	uint8_t RxBuffer[10] = { 0 };
	CoS_StatusTypeDef status;

	status = ReadBytesWithACKFromColiasBasic(hObj->hHUART, LEDConfig->LEDx, RxBuffer, 1);
	if (status != CoS_STATE_OK)
		return	status;

	if (RxBuffer[1] | 0x80)
	{
		LEDConfig->EN = Enable;
	}
	else
	{
		LEDConfig->EN = DisEnable;
	}

	if (RxBuffer[1] | 0x40)
	{
		LEDConfig->ON = Enable;
	}
	else
	{
		LEDConfig->ON = DisEnable;
	}

	if (RxBuffer[1] | 0x20)
	{
		LEDConfig->Flash = Enable;
	}
	else
	{
		LEDConfig->Flash = DisEnable;
	}

	LEDConfig->FlashPeriod = RxBuffer[1] & 0x0f;

	return CoS_STATE_OK;

}

CoS_StatusTypeDef ColiasBasicTCRT_Config(CoS_Motion_TypeDef *hObj, ColiasBasicTCRT_TypeDef *TCRTConfig)
{
	uint8_t TxBuffer[3] = { 0 };
	CoS_StatusTypeDef status;

	TxBuffer[0] = TCRTConfig->TCRTx;
	if (TCRTConfig->EN == Enable)
	{
		TxBuffer[1] |= 0x80;
	}


	if (TCRTConfig->IT_EN == Enable)
	{
		TxBuffer[1] |= 0x40;
	}
	if (TCRTConfig->Indicator_LED_EN == Enable)
	{
		TxBuffer[1] |= 0x20;
	}

	if (TCRTConfig->Indicator_LED == RedLED)
	{
		TxBuffer[1] |= 0x10;
	}

	//Limit the threshold into 0x0fff
		if (TCRTConfig->ThresholdForIT >= 0x03ff)
		{
			TxBuffer[0] |= 0x03;
			TxBuffer[1] = 0xff;
		}
		else
		{
			//TxBuffer[1] &= ~0x0f;
			TxBuffer[0] |= ((TCRTConfig->ThresholdForIT) & 0x0003);
			TxBuffer[1] = (TCRTConfig->ThresholdForIT) >> 2;
		}
	//
	if (hObj->ACK)
		return WriteBytesWithACKToColiasBasic(hObj->hHUART, TxBuffer, 3);
	else
		return WriteBytesWithoutACKToColiasBasic(hObj->hHUART, TxBuffer, 3);
}

uint16_t ColiasBasicTCRTADC_Get(CoS_Motion_TypeDef *hObj,uint8_t X)
{
	uint8_t	RxBuffer[10] = { 0 };
	uint16_t ADC_Value = 0;
	CoS_StatusTypeDef status;
	switch (X)
	{
	case 1:
		status = ReadBytesWithACKFromColiasBasic(hObj->hHUART, 0x22, RxBuffer, 2);
		if (status != CoS_STATE_OK)
			return	status;
		
		
		
		break;
	case 2:
		status = ReadBytesWithACKFromColiasBasic(hObj->hHUART, 0x26, RxBuffer, 2);
		if (status != CoS_STATE_OK)
			return	status;
		break;
	case 3:
		status = ReadBytesWithACKFromColiasBasic(hObj->hHUART, 0x2A, RxBuffer, 2);
		if (status != CoS_STATE_OK)
			return	status;
		break;
	default:
		return CoS_STATE_ERROR;
	}




	//TCRTConfig->ADC_Value = RxBuffer[1];
	ADC_Value = (uint16_t)(((uint16_t)RxBuffer[1]) << 2) | (RxBuffer[2] & 0x03);
	return ADC_Value;
	
	return CoS_STATE_OK;

}

CoS_StatusTypeDef ColiasBasicLightS_Config(CoS_Motion_TypeDef *hObj, ColiasBasicLightSensor_TypeDef *LightS_Config)
{
	uint8_t TxBuffer[3] = { 0 };
	CoS_StatusTypeDef status;

	TxBuffer[0] = LightS_Config->Light_x;

	if (LightS_Config->EN == Enable)
	{
		TxBuffer[1] |= 0x80;
	}

	if (LightS_Config->IT_EN == Enable)
	{
		TxBuffer[1] |= 0x40;
	}

	if (LightS_Config->Indicator_LED_EN == Enable)
	{
		TxBuffer[1] |= 0x20;
	}

	if (LightS_Config->Indicator_LED == RedLED)
	{
		TxBuffer[1] |= 0x10;
	}
	//Limit the threshold into 0x0fff
	if (LightS_Config->ThresholdForIT >= 0x0fff)
	{

		TxBuffer[1] |= 0x0f;

		TxBuffer[2] = 0xff;
	}
	else
	{
		TxBuffer[1] |= ((LightS_Config->ThresholdForIT) & 0x000f);

		TxBuffer[2] = (LightS_Config->ThresholdForIT) >> 4;

	}
	if (hObj->ACK)
		return WriteBytesWithACKToColiasBasic(hObj->hHUART, TxBuffer, 3);
	else
		return WriteBytesWithoutACKToColiasBasic(hObj->hHUART, TxBuffer, 3);
}


CoS_StatusTypeDef ColiasBasicLightSADC_Get(CoS_Motion_TypeDef *hObj, ColiasBasicLightSensor_TypeDef *LightS_Config)
{

	uint8_t	RxBuffer[10] = { 0 };
	CoS_StatusTypeDef status;


	switch (LightS_Config->Light_x)
	{
	case Light_Left:		status = ReadBytesWithACKFromColiasBasic(hObj->hHUART, Light_Left, RxBuffer, 2);
		if (status != CoS_STATE_OK)
			return	status;
		break;
	case Light_Right:	status = ReadBytesWithACKFromColiasBasic(hObj->hHUART, Light_Right, RxBuffer, 2);
		if (status != CoS_STATE_OK)
			return	status;
		break;
	default:
		return CoS_STATE_ERROR;
	}


	LightS_Config->ADC_Value = RxBuffer[0];
	LightS_Config->ADC_Value = ((LightS_Config->ADC_Value) << 4) | (RxBuffer[1] & 0x0f);

	return CoS_STATE_OK;

}

CoS_StatusTypeDef ColiasBasicInfrared_Config(CoS_Motion_TypeDef *hObj, ColiasBasicInfrared_TypeDef *Infrared_Config)
{
	uint8_t TxBuffer[3] = { 0 };
	CoS_StatusTypeDef status;
	TxBuffer[0] = 0x88;
	if (Infrared_Config->EN == Enable)
	{
		TxBuffer[1] |= 0x80;
	}

	if (Infrared_Config->IT_EN == Enable)
	{
		TxBuffer[1] |= 0x40;
	}

	//Limit the threshold into 0x0fff
	if (Infrared_Config->ThresholdForIT >= 0x0fff)
	{

		TxBuffer[1] |= 0x0f;

		TxBuffer[2] = 0xff;
	}
	else
	{
		TxBuffer[1] |= ((Infrared_Config->ThresholdForIT) & 0x000f);

		TxBuffer[2] = (Infrared_Config->ThresholdForIT) >> 4;

	}
	if (hObj->ACK)
		return WriteBytesWithACKToColiasBasic(hObj->hHUART, TxBuffer, 2);
	else
		return WriteBytesWithoutACKToColiasBasic(hObj->hHUART, TxBuffer, 2);


}

CoS_StatusTypeDef ColiasBasicMotionControler_Config(CoS_Motion_TypeDef *hObj, ColiasBasicMotion_TypeDef *Motion_Config)
{
	uint8_t TxBuffer[5] = { 0 };
	CoS_StatusTypeDef status;

	TxBuffer[0] = 0x90;

	if (Motion_Config->EN == Enable)
	{
		TxBuffer[1] |= 0x80;
	}

	switch (Motion_Config->DefaultMotionMode)
	{
	case Stop:
		break;
	case Front:		TxBuffer[1] |= 0x20;
		break;
	case Avoid:		TxBuffer[1] |= 0x40;
		break;
	default:
		return CoS_STATE_ERROR;
	}

	if (Motion_Config->IT_EN == Enable)
	{
		TxBuffer[1] |= 0x10;
	}

	if (Motion_Config->StepIT_EN == Enable)
	{
		TxBuffer[2] |= 0x80;
	}

	if (Motion_Config->QueueEmptyIT_EN == Enable)
	{
		TxBuffer[2] |= 0x40;
	}
	if (Motion_Config->QueueOverflowIT_EN == Enable)
	{
		TxBuffer[2] |= 0x20;
	}
	if (Motion_Config->SequenceDirection == TakeFromTop)
	{
		TxBuffer[2] |= 0x10;
	}

	if (Motion_Config->QueueAbort == Enable)
	{
		TxBuffer[2] |= 0x08;
	}
	if (Motion_Config->QueuePause == Enable)
	{
		TxBuffer[2] |= 0x04;
	}

	TxBuffer[3] = Motion_Config->MotionPeriod;
	TxBuffer[4] = Motion_Config->Bias;

	if (hObj->ACK)
		return WriteBytesWithACKToColiasBasic(hObj->hHUART, TxBuffer, 5);
	else
		return WriteBytesWithoutACKToColiasBasic(hObj->hHUART, TxBuffer, 5);
}

CoS_StatusTypeDef ColiasBasicMotion_Config(CoS_Motion_TypeDef *hObj, ColiasBasicMotionConf_TypeDef *Sequence_Config)
{
	uint8_t TxBuffer[4] = { 0 };
	CoS_StatusTypeDef status;


	switch (Sequence_Config->MotionID)
	{
	case ColiasFront:
		TxBuffer[0] = 0xa3;
		break;
	case ColiasLeft:
		TxBuffer[0] = 0xa6;
		break;
	case ColiasRight:
		TxBuffer[0] = 0xa9;
		break;
	case ColiasBack:
		TxBuffer[0] = 0xac;
		break;
	case ColiasNew1:
		TxBuffer[0] = 0xAf;
		break;
	case ColiasNew2:
		TxBuffer[0] = 0xb2;
		break;
	case ColiasNew3:
		TxBuffer[0] = 0xb5;
		break;
	case DefaultMotion:
		TxBuffer[0] = 0xa0;
		break;
	default:
		return CoS_STATE_ERROR;
	}

	TxBuffer[1] = Sequence_Config->LeftWheelPower;
	TxBuffer[2] = Sequence_Config->RightWheelPower;
	if (Sequence_Config->BiasBypass == Enable)
	{
		TxBuffer[3] |= 0x04;
	}

	if (Sequence_Config->LeftwheelDirection == RollBack)
	{
		TxBuffer[3] |= 0x02;
	}

	if (Sequence_Config->RightwheelDirection == RollBack)
	{
		TxBuffer[3] |= 0x01;
	}

	if (Sequence_Config->PeriodOffset < 0)
	{
		TxBuffer[3] |= 0x80;
	}

	if (abs(Sequence_Config->PeriodOffset) >= 16)
	{
		TxBuffer[3] |= 0x78;
	}
	else
	{
		TxBuffer[3] |= (((Sequence_Config->PeriodOffset) << 3) & 0x78);
	}


	if (hObj->ACK)
		return WriteBytesWithACKToColiasBasic(hObj->hHUART, TxBuffer, 4);
	else
		return WriteBytesWithoutACKToColiasBasic(hObj->hHUART, TxBuffer, 4);
}


CoS_StatusTypeDef ColiasBasic_WriteMotionQueue(CoS_Motion_TypeDef *hObj, ColiasBasicMotionQueue_TypeDef *WriteMotionQueue)
{
	uint8_t TxBuffer[11] = { 0 };
	CoS_StatusTypeDef status;
	if (WriteMotionQueue->WritePosition != NULL)
	{
		TxBuffer[0] = WriteMotionQueue->WritePosition;
	}
	else
	{
		TxBuffer[0] = FromQ1;
	}


	for (uint8_t i = 0; i < WriteMotionQueue->AmountOfWirte; i++)
	{
		TxBuffer[i + 1] = (WriteMotionQueue->QueueMotionIDBuffer[i]) << 5;
		
		if(WriteMotionQueue->QueueAmountBuffer[i] > 0x1f)
		{
			TxBuffer[i + 1] |= 0x1F;
		}
		else
		{
			TxBuffer[i + 1] |= (WriteMotionQueue->QueueAmountBuffer[i] & 0x1F);
		}
		//printf("%2x\r\n", TxBuffer[i+1]);
	}
	//printf("%2x %2x \r\n", TxBuffer[0], TxBuffer[1]);
	
	if (hObj->ACK)
		return WriteBytesWithACKToColiasBasic(hObj->hHUART, TxBuffer, (WriteMotionQueue->AmountOfWirte) + 1);
	else
		return WriteBytesWithoutACKToColiasBasic(hObj->hHUART, TxBuffer, (WriteMotionQueue->AmountOfWirte) + 1);
	
}


CoS_StatusTypeDef ColiasBasic_RunCustomizeMotion(CoS_Motion_TypeDef *hObj, ColiasBasicCustomizeMotion_TypeDef *CustomizeMotion)
{
	uint8_t TxBuffer[6] = { 0 };
	CoS_StatusTypeDef status;

	TxBuffer[0] = 0xD0;

	if (CustomizeMotion->EN == Enable)
	{
		TxBuffer[1] |= 0x80;
	}

	if (CustomizeMotion->CustomMotionFinishIT_EN == Enable)
	{
		TxBuffer[1] = 0x40;
	}
	if (CustomizeMotion->BiasBypass == Enable)
	{
		TxBuffer[1] |= 0x04;
	}

	if (CustomizeMotion->LeftwheelDirection == RollBack)
	{
		TxBuffer[1] |= 0x02;
	}

	if (CustomizeMotion->RightwheelDirection == RollBack)
	{
		TxBuffer[1] |= 0x01;
	}

	TxBuffer[2] = CustomizeMotion->LeftWheelPower;
	TxBuffer[3] = CustomizeMotion->RightWheelPower;
	TxBuffer[4] = (uint8_t)((CustomizeMotion->Period) >> 8);
	TxBuffer[5] = (uint8_t)CustomizeMotion->Period;

	if (hObj->ACK)
		return WriteBytesWithACKToColiasBasic(hObj->hHUART, TxBuffer, 5);
	else
		return WriteBytesWithoutACKToColiasBasic(hObj->hHUART, TxBuffer, 5);
}

//ACK have not been written
CoS_StatusTypeDef WriteBytesWithACKToColiasBasic(CoS_UART_HandleTypeDef *hObj, uint8_t *pTxData, uint8_t Size)
{

	CoS_StatusTypeDef status;
	uint8_t RxBuffer[10] = { 0 };
	uint8_t timecount = 0;
	CoS_UART_TypeDef* Instance = hObj->Instance;

	Instance->escape1 = WriteWithACK;
	status = CoS_UART_sendFrame(hObj, pTxData, Size);
	if (status != CoS_STATE_OK)
	{
		return	status;
	}

	Instance->Rx_BufSize = 10;
	Instance->Rx_Base = RxBuffer;
	status = CoS_UART_RxRST(hObj);

	UART_transmit_run(hObj);



	while (!(timecount > 20 || Instance->RxOK))
	{
		timecount++;
		delay_short(100000);
	}
	
	if (RxBuffer[0] == 0x7B)
	{
		return CoS_STATE_OK;
	}
	else if (RxBuffer[0] == 0x7A)
	{
		return CoS_STATE_ERROR;
	}
	if(timecount>=20)
		return CoS_STATE_TIMEOUT;
	else
		return CoS_STATE_ERROR;
}

//ACK have not been written
CoS_StatusTypeDef WriteBytesWithoutACKToColiasBasic(CoS_UART_HandleTypeDef *hObj, uint8_t *pTxData, uint8_t Size)
{
	hObj->Instance->escape1 = WriteWithoutACK;
	return CoS_UART_sendFrame(hObj, pTxData, Size);
}

CoS_StatusTypeDef ReadBytesWithACKFromColiasBasic(CoS_UART_HandleTypeDef *hObj, uint8_t Address, uint8_t * pRxData, uint8_t Size)
{
	CoS_StatusTypeDef status;
	uint8_t TxBuffer[2] = { 0 };
	uint8_t timecount = 0;
	CoS_UART_TypeDef* Instance = hObj->Instance;


	TxBuffer[0] = Address;
	TxBuffer[1] = Size;
	Instance->escape1 = ReadWithACK;
	status = CoS_UART_sendFrame(hObj, TxBuffer, 2);
	if (status != CoS_STATE_OK)
	{
		return status;
	}
	

	Instance->Rx_BufSize = Size + 5;
	Instance->Rx_Base = pRxData;
	status = CoS_UART_RxRST(hObj);
	
	UART_transmit_run(hObj);

	while (!(timecount > 20 || Instance->RxOK))
	{
		//LED2_Toggle;
		timecount++;
		delay_short(10000);
	}
	if (timecount>=20)
		return CoS_STATE_TIMEOUT;
	if (pRxData[0] == 0x7B)
	{
		return CoS_STATE_OK;
	}
	else if (pRxData[0] == 0x7A)
	{
		return CoS_STATE_ERROR;
	}

	return CoS_STATE_ERROR;
}
//Set Colias Basic regs bitwise
/*
Address: reg address
mask: bits to be modified
mode:
    0: reset masked bits
    1: set masked bits
    2: toggle maksed bits 
*/
CoS_StatusTypeDef BitSetColiasBasicReg(CoS_UART_HandleTypeDef *hObj, uint8_t Address, uint8_t mask,uint8_t mode)
{
 u8 data,pStr[10]={0},sta;
	sta=ReadBytesWithACKFromColiasBasic(hObj,Address,&pStr[0],1);
	if (sta==CoS_STATE_OK)
	{
		data=pStr[1];
		if (mode==0)
			pStr[1]=data&mask;
		else if (mode==1)
			pStr[1]=data|mask;
		else if (mode==2)
			pStr[1]=data^mask;
		else
			return CoS_STATE_ERROR;
		pStr[0]=Address;
		return WriteBytesWithACKToColiasBasic(hObj,&pStr[0],2);
	}
	return sta;
}


void CoS_MotionEXIT_Init(void)
{
	GPIO_InitTypeDef GPIO_InitStruct;
	GPIO_InitStruct.Pin = GP1_Pin;
	GPIO_InitStruct.Mode = GPIO_MODE_IT_FALLING;
	GPIO_InitStruct.Pull = GPIO_NOPULL;
	GPIO_InitStruct.Speed = GPIO_SPEED_FREQ_MEDIUM;
	HAL_GPIO_Init(GP1_GPIO_Port, &GPIO_InitStruct);

	/* EXTI interrupt init*/
	HAL_NVIC_SetPriority(EXTI15_10_IRQn, 7, 0);
	HAL_NVIC_EnableIRQ(EXTI15_10_IRQn);
}

CoS_StatusTypeDef MotionTask(CoS_Motion_TypeDef *hObj)
{
	CoS_StatusTypeDef status;

	ColiasBasicTCRT_TypeDef TCRT;
	ColiasBasicMotionQueue_TypeDef WriteMotionQueue;
	memset(&WriteMotionQueue, 0, sizeof(WriteMotionQueue));
	ColiasBasicCustomizeMotion_TypeDef CustomMotion;
	memset(&CustomMotion, 0, sizeof(CustomMotion));
	//uint8_t RxBuffer[0xff] = {0};

	
	status = ColiasBasicITEvents_Check(hObj, &hObj->IT_Events);
	if (status != CoS_STATE_OK)
		return status;
	
	if (hObj->IT_Events.MotionIT & Queue_Empty_IT)
	{
		hObj->MotionFinishFlag = 1;
	}

	if (hObj->IT_Events.SensorIT & TCRT_Left_IT)
	{
		if (hObj->MotionFinishFlag == 1)
		{
			WriteMotionQueue.QueueMotionIDBuffer[0] = ColiasRight;
			WriteMotionQueue.QueueAmountBuffer[0] = 1;
			WriteMotionQueue.AmountOfWirte = 1;
			WriteMotionQueue.WritePosition = FromQ1;
			status = ColiasBasic_WriteMotionQueue(hObj, &WriteMotionQueue);
			while (status != CoS_STATE_OK)
			{
				LED1_Toggle;
				status = ColiasBasic_WriteMotionQueue(hObj, &WriteMotionQueue);
			}
			hObj->MotionFinishFlag = 0;
			LED2_Toggle;
		}
		else
		{
			(hObj->timecount)++;
		}
	}


	if (hObj->IT_Events.SensorIT & TCRT_Right_IT)
	{
		if (hObj->MotionFinishFlag == 1)
		{
			WriteMotionQueue.QueueMotionIDBuffer[0] = ColiasLeft;
			WriteMotionQueue.QueueAmountBuffer[0] = 1;
			WriteMotionQueue.AmountOfWirte = 1;
			WriteMotionQueue.WritePosition = FromQ1;
			status = ColiasBasic_WriteMotionQueue(hObj, &WriteMotionQueue);

			while (status != CoS_STATE_OK)
			{
				LED1_Toggle;
				status = ColiasBasic_WriteMotionQueue(hObj, &WriteMotionQueue);
			}
			hObj->MotionFinishFlag = 0;
			LED2_Toggle;
		}
		else
		{
			(hObj->timecount)++;
		}
	}


	if (hObj->IT_Events.SensorIT & TCRT_Middle_IT)
	{
		if (hObj->MotionFinishFlag == 1)
		{
			WriteMotionQueue.QueueMotionIDBuffer[0] = ColiasRight;
			WriteMotionQueue.QueueAmountBuffer[0] = 8;
			WriteMotionQueue.AmountOfWirte = 1;
			WriteMotionQueue.WritePosition = FromQ1;
			status = ColiasBasic_WriteMotionQueue(hObj, &WriteMotionQueue);

			while (status != CoS_STATE_OK)
			{
				LED1_Toggle;
				status = ColiasBasic_WriteMotionQueue(hObj, &WriteMotionQueue);
			}
			hObj->MotionFinishFlag = 0;
			LED2_Toggle;
		}
		else
		{
			(hObj->timecount)++;
		}
	}

	

	if (hObj->timecount > 50)
	{
		hObj->MotionFinishFlag = 1;
		hObj->timecount = 0;
	}


	return CoS_STATE_OK;
}


/*
For Qinbing Attention
define your own motion quere as following pattern
u8 motionStr_NEW[4]={0*,0x23**,0x63,0***}; 
*leave first byte for writing address, as it would affect preempting feature in baisc board
** start queue from second byte. Refer to the register description file fore detailed information
***string end in zero.
*/
u8 motionStr_Front[3]={0,0x21,0};	//forward: 1
u8 motionStr_Left_turn[4]={0,0x02,0x48,0};	//left: 8 and stop: 3
u8 motionStr_Right_turn[4]={0,0x02,0x68,0};	//right: 9 and stop: 3
u8 motionStr_Left_back[4]={0,0x03,0xa9,0};
u8 motionStr_Right_back[4]={0,0x03,0xc9,0};
u8 motionStr_back[3]={0,0x81,0};	//back: 1
u8 motionStr_stop[3]={0,0x03,0};	//stop: 3
u8 motionStr_avoid_left[4]={0,0x03,0xa9,0};	//long left 3 + backward(stop)
u8 motionStr_avoid_right[4]={0,0x03,0xc9,0};	//long right 3 + backward(stop)
u8 motionStr_long_stop[3]={0,0x0f,0};				//long stop: 16
u8 motionStr_long_back[3]={0,0x83,0};				//long back: 3

u8 motionStr_new_Q[5]={0,0xa6,0xc6,0x88,0};

/*
For Qinbing Attention
chage the charecter in cases for your own motion.
*/
u8 fast_motion_control(CoS_Motion_TypeDef *hObj,u8 motionID,u8 mode)
{
	u8 addr;
	if (mode)
		addr=0xc1; //norlam queue (1)
	else
		addr=0xc0;//preempted queue (0)
	u8* pStr;
	switch (motionID)
	{
		case 'F':
		{
			pStr=motionStr_Front;
			break;
		}
		case 'L':
		{
			pStr=motionStr_Left_turn;
			//pStr=motionStr_new_Q;
			break;
		}
		case 'R':
		{
			pStr=motionStr_Right_turn;
			//pStr=motionStr_new_Q;
			break;
		}
		case 'S':
		{
			pStr=motionStr_stop;
			break;
		}
		case 'A':
		{
			pStr=motionStr_avoid_left;
			break;
		}
		case 'D':
		{
			pStr=motionStr_avoid_right;
			break;
		}
		case 'O':
		{
			//pStr=motionStr_long_stop;
			pStr=motionStr_Left_back;
			break;
		}
		case 'B':
		{
			//pStr=motionStr_long_back;
			pStr=motionStr_Right_back;
			break;
		}
	}
	*pStr=addr;
	return WriteBytesWithACKToColiasBasic(hCoS->hHUART4,pStr,strlen(pStr));
}


