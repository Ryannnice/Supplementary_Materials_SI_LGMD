#include "coliasSense_CLED.h"
#include "coliasSense_board.h"
#include "delay.h"
#include "math.h"
#include "arm_math.h"

static u8* mode_nameStr[]={"OFF","Random","Breath","R","G","B","FROZEN"};
static CoS_CLED_HandleTypeDef* hHCLED1;
//static CoS_CCSTypeDef CC_CLED1;

//static TIM_HandleTypeDef*  htim1;
//static RNG_HandleTypeDef* hrng;


//CoS_CLED_TypeDef* const pCLED1=&CLED1;
//CoS_CCSTypeDef* const pCC_CLED1=&CC_CLED1;

CLED_ColorStructType const defaultcolor=
	{.A=32760,.shift=32760,.freq=6,.phy=0};

static void CoS_CLED_GPIO_Init(void)
{
	GPIO_InitTypeDef GPIO_InitStruct;
	__HAL_RCC_GPIOE_CLK_ENABLE();
	GPIO_InitStruct.Pin = LED4_Pin|LED5_Pin|LED6_Pin;
	GPIO_InitStruct.Mode = GPIO_MODE_AF_PP;
	GPIO_InitStruct.Pull = GPIO_PULLUP;
	GPIO_InitStruct.Speed = GPIO_SPEED_FREQ_MEDIUM;
	GPIO_InitStruct.Alternate = GPIO_AF1_TIM1;
	HAL_GPIO_Init(LED4_GPIO_Port, &GPIO_InitStruct);
	
	GPIO_InitStruct.Pin = LED5_Pin;
	HAL_GPIO_Init(LED5_GPIO_Port, &GPIO_InitStruct);
	GPIO_InitStruct.Pin = LED6_Pin;
	HAL_GPIO_Init(LED6_GPIO_Port, &GPIO_InitStruct);
}


/* TIM1 init function */
static void CoS_TIM1_Init(TIM_HandleTypeDef* htim1)
{
	TIM_ClockConfigTypeDef sClockSourceConfig;
	TIM_MasterConfigTypeDef sMasterConfig;
	TIM_BreakDeadTimeConfigTypeDef sBreakDeadTimeConfig;
	TIM_OC_InitTypeDef sConfigOC;
	__TIM1_CLK_ENABLE();
	htim1->Instance = TIM1;
	htim1->Init.Prescaler = 1;
	htim1->Init.CounterMode = TIM_COUNTERMODE_UP;
	htim1->Init.Period = 65535;
	htim1->Init.ClockDivision = TIM_CLOCKDIVISION_DIV1;
	htim1->Init.RepetitionCounter = 0;
	if (HAL_TIM_Base_Init(htim1) != HAL_OK)
		Error_Handler();
	sClockSourceConfig.ClockSource = TIM_CLOCKSOURCE_INTERNAL;
	if (HAL_TIM_ConfigClockSource(htim1, &sClockSourceConfig) != HAL_OK)
		Error_Handler();
	if (HAL_TIM_OC_Init(htim1) != HAL_OK)
		Error_Handler();
	sMasterConfig.MasterOutputTrigger = TIM_TRGO_ENABLE;
	sMasterConfig.MasterSlaveMode = TIM_MASTERSLAVEMODE_DISABLE;
	if (HAL_TIMEx_MasterConfigSynchronization(htim1, &sMasterConfig) != HAL_OK)
		Error_Handler();
	sBreakDeadTimeConfig.OffStateRunMode = TIM_OSSR_DISABLE;
	sBreakDeadTimeConfig.OffStateIDLEMode = TIM_OSSI_DISABLE;
	sBreakDeadTimeConfig.LockLevel = TIM_LOCKLEVEL_OFF;
	sBreakDeadTimeConfig.DeadTime = 0;
	sBreakDeadTimeConfig.BreakState = TIM_BREAK_DISABLE;
	sBreakDeadTimeConfig.BreakPolarity = TIM_BREAKPOLARITY_HIGH;
	sBreakDeadTimeConfig.AutomaticOutput = TIM_AUTOMATICOUTPUT_DISABLE;
	if (HAL_TIMEx_ConfigBreakDeadTime(htim1, &sBreakDeadTimeConfig) != HAL_OK)
		Error_Handler();
	sConfigOC.OCMode = TIM_OCMODE_PWM1;
	sConfigOC.Pulse = 0;
	sConfigOC.OCPolarity = TIM_OCPOLARITY_LOW;
	sConfigOC.OCNPolarity = TIM_OCNPOLARITY_HIGH;
	sConfigOC.OCFastMode = TIM_OCFAST_DISABLE;
	sConfigOC.OCIdleState = TIM_OCIDLESTATE_RESET;
	sConfigOC.OCNIdleState = TIM_OCNIDLESTATE_RESET;
	if (HAL_TIM_OC_ConfigChannel(htim1, &sConfigOC, TIM_CHANNEL_1) != HAL_OK)
		Error_Handler();
	sConfigOC.Pulse = 0;
	if (HAL_TIM_OC_ConfigChannel(htim1, &sConfigOC, TIM_CHANNEL_2) != HAL_OK)
		Error_Handler();
	sConfigOC.Pulse = 0;
	if (HAL_TIM_OC_ConfigChannel(htim1, &sConfigOC, TIM_CHANNEL_3) != HAL_OK)
		Error_Handler();
	HAL_TIM_Base_MspInit(htim1);
	//HAL_TIM_MspPostInit(&htim1);
	//HAL_TIM_Base_Start(&htim1);
	//

}


static void CoS_RNG_Init(RNG_HandleTypeDef* hrng)
{
	__HAL_RCC_RNG_CLK_ENABLE();
	hrng->Instance=RNG;
	HAL_RNG_Init(hrng);
}	



	
static void CLEDCC_Init(CoS_CLED_TypeDef* hObj)
{
	
	hObj->CC.mode=CLED_Breath;
	hObj->Rs=defaultcolor;
	hObj->Gs=defaultcolor;
	hObj->Bs=defaultcolor;
	
	hObj->Gs.A=15000;
	hObj->Gs.shift=15000;
	hObj->Gs.freq=8;
	hObj->Gs.phy=1.5708;
	hObj->Bs.freq=10;
	hObj->Bs.phy=2.618;
}

CoS_StatusTypeDef CoS_CLED_Init(CoS_CLED_HandleTypeDef *hObj,
								CoS_CLED_InitTypeDef* Init)
{
	/* Check the handle allocation */
	if(hObj == NULL)
	{
		return CoS_STATE_ERROR;
	}
	if( Init->hrng == NULL ||
		Init->htim1 == NULL ||
		Init->Instance==NULL)
	{
		return CoS_STATE_ERROR;
	}
	__CoS_LOCK(hObj);
	if(hObj->iState == CoS_STATE_UNSET)
	{  
		/* Allocate lock resource and initialize it */
		hObj->Lock = CoS_UNLOCKED;
		hObj->Instance=Init->Instance;
		hObj->hrng=Init->hrng;
		hObj->htim1=Init->htim1;
		hHCLED1=hObj;
		
	}
	/* Change peripheral state */
	hObj->iState = CoS_STATE_BUSY;
	/* Init the low level hardware */
	CoS_CLED_GPIO_Init();
	CoS_RNG_Init(hObj->hrng);
	CoS_TIM1_Init(hObj->htim1);
	CLEDCC_Init(hObj->Instance);
	if (INIT_INFO_printf)
		printf("{}Multi Color LED Initialize OK.\r\n");
	/* Initialize the Instance state */
	hObj->iState = CoS_STATE_OK;
	CoS_CLED_StateReport(hObj);
	__CoS_UNLOCK(hObj);
	/* Return function status */
	return CoS_STATE_OK;
}
void CoS_CLED_StateReport(CoS_CLED_HandleTypeDef *hObj)
{
	BIO_IO_TypdeDef HPin;
	printf("{  Multi Color LED state:%d}\r\n",hObj->iState);
	if (hObj->iState!=CoS_STATE_OK)
		return;
	printf("{  mode: %s}\r\n",mode_nameStr[hObj->Instance->CC.mode]);
	printf("{  pin LED4:");
	HPin.Pin=LED4_Pin;
	HPin.Port=LED4_GPIO_Port;
	CoS_Pin_StateReport(&HPin,1);
	HPin.Pin=LED5_Pin;
	HPin.Port=LED5_GPIO_Port;
	printf("{  pin LED5:");
	CoS_Pin_StateReport(&HPin,1);
	HPin.Pin=LED6_Pin;
	HPin.Port=LED6_GPIO_Port;
	printf("{  pin LED6:");
	CoS_Pin_StateReport(&HPin,1);

}
CoS_StatusTypeDef CoS_CLED_DeInit(CoS_CLED_HandleTypeDef *hObj)
{
	/* Check the handle allocation */
	if(hObj == NULL)
	{
		return CoS_STATE_ERROR;
	}
	__CoS_LOCK(hObj);
	if(hObj->iState == CoS_STATE_UNSET)
	{  
		/* Allocate lock resource and initialize it */
		hObj->Lock = CoS_UNLOCKED;
		/* Init the low level hardware */
		//CoS_LowInit(handle);
	}
	/* Change peripheral state */
	hObj->iState = CoS_STATE_BUSY;
	/* Enable the Peripheral */
	//CoS_LED_IO_Init();
	/* Initialize the Instance state */
	hObj->iState = CoS_STATE_OK;
	__CoS_UNLOCK(hObj);
	/* Return function status */
	return CoS_STATE_OK;
}

CoS_StatusTypeDef CoS_CLED_Enable(CoS_CLED_HandleTypeDef* hObj)
{
		/* Check the handle allocation */
	if(hObj == NULL)
	{
		return CoS_STATE_ERROR;
	}
	hObj->Instance->CC.enable=1;
	HAL_TIM_Base_Start(hObj->htim1);
	HAL_TIM_PWM_Start(hObj->htim1,CLED_ChannelB);
	HAL_TIM_PWM_Start(hObj->htim1,CLED_ChannelG);
	HAL_TIM_PWM_Start(hObj->htim1,CLED_ChannelR);
	
	return CoS_STATE_OK;
}

CoS_StatusTypeDef CoS_CLED_Disable(CoS_CLED_HandleTypeDef* hObj)
{
			/* Check the handle allocation */
	if(hObj == NULL)
	{
		return CoS_STATE_ERROR;
	}
	hObj->Instance->CC.enable=0;
	return CoS_STATE_OK;
}

void CoS_CLED_setRGB(u16 R, u16 G, u16 B)
{
	__HAL_TIM_SET_COMPARE(hHCLED1->htim1,CLED_ChannelB,B);
	__HAL_TIM_SET_COMPARE(hHCLED1->htim1,CLED_ChannelG,G);
	__HAL_TIM_SET_COMPARE(hHCLED1->htim1,CLED_ChannelR,R);
}

static void CLED_RandomGo(void)
{
	u16 t[3];
	u32 rd;
	HAL_RNG_GenerateRandomNumber(hHCLED1->hrng,&rd);
	t[0]=rd>>16;
	t[1]=rd&0x7ffe;
	t[2]=(rd&0xff)<<8;
	if (t[0]<t[1])
		if (t[0]<t[2])
			t[0]=0;
		else
			t[2]=0;
	else if (t[1]>t[2])
		t[2]=0;
	else
		t[1]=0;
	//SB_CLED_setRGB(rd>>16,rd&0xfffe,((rd&0xff)<<8)+((rd&0xff0000)>>16));
	hHCLED1->Instance->Rs.val=t[0];
	hHCLED1->Instance->Gs.val=t[1];
	hHCLED1->Instance->Bs.val=t[2];
	CoS_CLED_setRGB(t[0],t[1],t[2]);
}

u16 CoS_CLED_getV(CoS_CLED_ChannelTypeDef channel)
{
	return __HAL_TIM_GET_COMPARE(hHCLED1->htim1,channel);
}

static u16 CLED_ColorCalc(CoS_CLED_ChannelTypeDef channel)
{
	CLED_ColorStructType* pcc;
	u16 val;
	float f1,f2;
	if (channel==CLED_ChannelR)
		pcc=&(hHCLED1->Instance->Rs);
	else if (channel==CLED_ChannelG)
		pcc=&(hHCLED1->Instance->Gs);
	else if (channel==CLED_ChannelB)
		pcc=&(hHCLED1->Instance->Bs);
	else
		return 0;
	//f1=arm_cos_f32((float)(TIM3TickCount)/pcc->freq+pcc->phy);
	f2=((float)(*pTIM3Tick))/pcc->freq+pcc->phy;
	f1=arm_cos_f32(f2);
	val=(u16)(f1*pcc->A+pcc->shift);
	pcc->val=val;
	return val;
	
}

void CoS_LED_Update(CoS_CLED_HandleTypeDef *hObj)
{
	u16 t[3]={0};
	CoS_CLED_TypeDef* h=hObj->Instance;
	switch (h->CC.mode)
	{
		case CLED_Random:
		{
			CLED_RandomGo();
			h->CC.mode=CLED_Freeze;
			return;
		}
		case CLED_Breath:
		{
			CLED_ColorCalc(channel_Red);
			CLED_ColorCalc(channel_Green);
			CLED_ColorCalc(channel_Blue);
			break;
		}
		case CLED_R:
		{
			h->Rs.val=65535;
			h->Gs.val=0;
			h->Bs.val=0;
			break;
		}
		case CLED_G:
		{
			h->Rs.val=0;
			h->Gs.val=65535;
			h->Bs.val=0;
			break;
		}
		case CLED_B:
		{
			h->Rs.val=0;
			h->Gs.val=0;
			h->Bs.val=65535;
			break;
		}
		case CLED_Freeze:
			return;
		default:
		{
			h->Rs.val=0;
			h->Gs.val=0;
			h->Bs.val=0;
			break;
		}
	}
	CoS_CLED_setRGB(h->Rs.val,h->Gs.val,h->Bs.val);
}

