/**
  ******************************************************************************
  * @file    coliasSense_BIO.c
  * @author  CIL, University of Lincoln
  * @version V1.5.0
  * @date    28-December-2016
  * @brief   Colias Sense board Basic IO module driver.
  *          This file provides firmware functions to manage the following 
  *          functionalities of the Basic IO functions:
  *           + Initialize functions
  *           + LED control functions
  *           + Key state read and interrupt control
  *           + GPIO IO control functions
  *               
  ==============================================================================
                        ##### How to use this driver #####
  ==============================================================================
  [..]
    The BIO driver can be used as follows:
    
    (#) Declare a BIO_IO_TypdeDef instance structure.
    (#) Declare a CoS_BIO_HandleTypeDef handle structure.
    (#) Add a field of the pointer to CoS_BIO_HandleTypeDef in the SenseBoard_HandleTypeDef structure.
    (#) Point the handle to the correct CoS_BIO_HandleTypeDef location.
    (#) Initialize the BIO resources by implementing the CoS_BIO_Init() API:
        (##) Declare a CoS_BIO_InitTypeDef structure which contains the required resources.
            (+++) The location of BIO_IO_TypdeDef instance.
            (+++) The location of TIM_HandleTypeDef htim2, which is used as tic toc timer.
            (+++) The location of TIM_HandleTypeDef htim3, which is used for general background tasks calling. 
            (+++) The location of TIM_HandleTypeDef htim5, which is used for high priority background tasks calling. 
            (+++) The IRQn number of  general background tasks calling (htim2).
        (##) Define the background task calling function CoS_TIM3_IRQ_Task()
        
     [..] 
       (@) The background tasks will be managed by  BackgroundTask_Enable() and 
           BackgroundTask_Disable().
          
     [..] 
       (@) The instance BIO is singleton. There are two global pointer variable defined in the 
           coliasSense_BIO.c file
          

  ******************************************************************************
  * @attention
  *
  * <h2><center>&copy; COPYRIGHT(c) </center></h2>
  */ 

#include "coliasSense_BIO.h"
#include "coliasSense_board.h"

/****use pointers to the controller below to increase efficiency,
but should be careful. 

The object BIO is SINGLETON

************/
#ifdef Singleton_BIO
CoS_BIO_TypeDef* hBIO1;
uint64_t* pTIM3Tick;
#endif


static void CoS_TIM2_Init(TIM_HandleTypeDef* hObj);
static void CoS_TIM3_Init(TIM_HandleTypeDef* hObj);
static void CoS_BIO_LED_Init(void);
static void CoS_BIO_Keys_Init(void);


void BackgroundTask_Enable(void) 
{
	HAL_NVIC_EnableIRQ(hBIO1->BackgroundIRQn);
}
void BackgroundTask_Disable(void) 
{
	HAL_NVIC_DisableIRQ(hBIO1->BackgroundIRQn);
}



 __weak void Error_Handler(void)
{
	while(1) 
	{
	}
}
static void CoS_SystemClock_Config(void)
{

  RCC_OscInitTypeDef RCC_OscInitStruct;
  RCC_ClkInitTypeDef RCC_ClkInitStruct;
  RCC_PeriphCLKInitTypeDef PeriphClkInitStruct;

  __HAL_RCC_PWR_CLK_ENABLE();

  __HAL_PWR_VOLTAGESCALING_CONFIG(PWR_REGULATOR_VOLTAGE_SCALE1);

  RCC_OscInitStruct.OscillatorType = RCC_OSCILLATORTYPE_HSE;
  RCC_OscInitStruct.HSEState = RCC_HSE_ON;
  RCC_OscInitStruct.PLL.PLLState = RCC_PLL_ON;
  RCC_OscInitStruct.PLL.PLLSource = RCC_PLLSOURCE_HSE;
  RCC_OscInitStruct.PLL.PLLM = 25;
  RCC_OscInitStruct.PLL.PLLN = 360;
  RCC_OscInitStruct.PLL.PLLP = RCC_PLLP_DIV2;
  RCC_OscInitStruct.PLL.PLLQ = 4;
  if (HAL_RCC_OscConfig(&RCC_OscInitStruct) != HAL_OK)
  {
    Error_Handler();
  }

  if (HAL_PWREx_EnableOverDrive() != HAL_OK)
  {
    Error_Handler();
  }

  RCC_ClkInitStruct.ClockType = RCC_CLOCKTYPE_HCLK|RCC_CLOCKTYPE_SYSCLK
                              |RCC_CLOCKTYPE_PCLK1|RCC_CLOCKTYPE_PCLK2;
  RCC_ClkInitStruct.SYSCLKSource = RCC_SYSCLKSOURCE_PLLCLK;
  RCC_ClkInitStruct.AHBCLKDivider = RCC_SYSCLK_DIV1;
  RCC_ClkInitStruct.APB1CLKDivider = RCC_HCLK_DIV4;
  RCC_ClkInitStruct.APB2CLKDivider = RCC_HCLK_DIV2;
  if (HAL_RCC_ClockConfig(&RCC_ClkInitStruct, FLASH_LATENCY_5) != HAL_OK)
  {
    Error_Handler();
  }

  PeriphClkInitStruct.PeriphClockSelection = RCC_PERIPHCLK_I2S;
  PeriphClkInitStruct.PLLI2S.PLLI2SN = 213;
  PeriphClkInitStruct.PLLI2S.PLLI2SR = 4;
  if (HAL_RCCEx_PeriphCLKConfig(&PeriphClkInitStruct) != HAL_OK)
  {
    Error_Handler();
  }

  HAL_RCC_MCOConfig(RCC_MCO1, RCC_MCO1SOURCE_HSE, RCC_MCODIV_1);

  HAL_SYSTICK_Config(HAL_RCC_GetHCLKFreq()/1000);

  HAL_SYSTICK_CLKSourceConfig(SYSTICK_CLKSOURCE_HCLK);

  /* SysTick_IRQn interrupt configuration */
  HAL_NVIC_SetPriority(SysTick_IRQn, 0, 1);
}

static void HAL_NVIC_Init(void)
{

	HAL_NVIC_SetPriorityGrouping(NVIC_PRIORITYGROUP_3);

	/* System interrupt init*/

	HAL_NVIC_SetPriority(MemoryManagement_IRQn, 0, 0);
	HAL_NVIC_SetPriority(BusFault_IRQn, 0, 0);
	HAL_NVIC_SetPriority(UsageFault_IRQn, 0, 0);
	HAL_NVIC_SetPriority(SVCall_IRQn, 0, 0);
	HAL_NVIC_SetPriority(DebugMonitor_IRQn, 0, 0);
	HAL_NVIC_SetPriority(PendSV_IRQn, 0, 0);
	HAL_NVIC_SetPriority(SysTick_IRQn, 0, 1);
}

static void CoS_GPIO_Init(void)
{
	//TODO: this is the default IO pin config
	GPIO_InitTypeDef GPIO_InitStruct;
	GPIO_InitTypeDef GPIO_InitStructIn;
	/* GPIO Ports Clock Enable */
	__HAL_RCC_GPIOA_CLK_ENABLE();
	__HAL_RCC_GPIOB_CLK_ENABLE();
	__HAL_RCC_GPIOC_CLK_ENABLE();
	__HAL_RCC_GPIOD_CLK_ENABLE();
	__HAL_RCC_GPIOE_CLK_ENABLE();
	__HAL_RCC_GPIOH_CLK_ENABLE();
	GPIO_InitStruct.Mode = GPIO_MODE_OUTPUT_PP;
	GPIO_InitStruct.Pull = GPIO_PULLUP;
	GPIO_InitStruct.Speed = GPIO_SPEED_FREQ_MEDIUM;
	GPIO_InitStructIn.Mode = GPIO_MODE_INPUT;
	GPIO_InitStructIn.Pull = GPIO_NOPULL;
	GPIO_InitStructIn.Speed = GPIO_SPEED_FREQ_MEDIUM;
	GPIO_InitStruct.Pin = LED1_Pin;//output Pullup
	HAL_GPIO_Init(LED1_GPIO_Port, &GPIO_InitStruct);
	
	GPIO_InitStruct.Pin = LED2_Pin;//output Pullup
	HAL_GPIO_Init(LED2_GPIO_Port, &GPIO_InitStruct);
	
	GPIO_InitStruct.Pin = LED3_Pin;//output Pullup
	HAL_GPIO_Init(LED3_GPIO_Port, &GPIO_InitStruct);
	
	GPIO_InitStruct.Pin = LED4_Pin;//output Pullup
	HAL_GPIO_Init(LED4_GPIO_Port, &GPIO_InitStruct);
	
	GPIO_InitStruct.Pin = LED5_Pin;//output Pullup
	HAL_GPIO_Init(LED5_GPIO_Port, &GPIO_InitStruct);
	
	GPIO_InitStructIn.Pin = S1_Pin;//input nopull
	HAL_GPIO_Init(S1_GPIO_Port, &GPIO_InitStructIn);
	
	GPIO_InitStructIn.Pin = S2_Pin;//input nopull
	HAL_GPIO_Init(S2_GPIO_Port, &GPIO_InitStructIn);
	
	GPIO_InitStructIn.Pin = GP1_Pin;//input nopull
	HAL_GPIO_Init(GP1_GPIO_Port, &GPIO_InitStructIn);
	
	GPIO_InitStructIn.Pin = GP2_Pin;//input nopull
	HAL_GPIO_Init(GP2_GPIO_Port, &GPIO_InitStructIn);
	
	GPIO_InitStructIn.Pin = FLASH_CS_Pin;//input nopull
	HAL_GPIO_Init(FLASH_CS_GPIO_Port, &GPIO_InitStructIn);
	
	GPIO_InitStructIn.Pin = FLASH_D2_Pin;//input nopull
	HAL_GPIO_Init(FLASH_D2_GPIO_Port, &GPIO_InitStructIn);
	
	GPIO_InitStructIn.Pin = FLASH_D3_Pin;//input nopull
	HAL_GPIO_Init(FLASH_D3_GPIO_Port, &GPIO_InitStructIn);
	
	GPIO_InitStructIn.Pin = eSPI_P2_Pin;//input nopull
	HAL_GPIO_Init(eSPI_P2_GPIO_Port, &GPIO_InitStructIn);
	
	GPIO_InitStructIn.Pin = eSPI_P3_Pin;//input nopull
	HAL_GPIO_Init(eSPI_P3_GPIO_Port, &GPIO_InitStructIn);
	
	GPIO_InitStructIn.Pin = eSPI_MOSI_Pin;//input nopull
	HAL_GPIO_Init(eSPI_MOSI_GPIO_Port, &GPIO_InitStructIn);
	
	GPIO_InitStructIn.Pin = sI2C_SCL_Pin;//input nopull
	HAL_GPIO_Init(sI2C_SCL_GPIO_Port, &GPIO_InitStructIn);
	
	GPIO_InitStructIn.Pin = sI2C_SDA_Pin;//input nopull
	HAL_GPIO_Init(sI2C_SDA_GPIO_Port, &GPIO_InitStructIn);
	
	GPIO_InitStructIn.Pin = ACC_INT_Pin;//input nopull
	HAL_GPIO_Init(ACC_INT_GPIO_Port, &GPIO_InitStructIn);
	
	GPIO_InitStructIn.Pin = ACC_AD0_Pin;//input nopull
	HAL_GPIO_Init(ACC_AD0_GPIO_Port, &GPIO_InitStructIn);
	
	GPIO_InitStructIn.Pin = ACC_nCS_Pin;//input nopull
	HAL_GPIO_Init(ACC_nCS_GPIO_Port, &GPIO_InitStructIn);
	
	GPIO_InitStructIn.Pin = eSPI_P1_Pin;//input nopull
	HAL_GPIO_Init(eSPI_P1_GPIO_Port, &GPIO_InitStructIn);
	
	GPIO_InitStructIn.Pin = eSPI_SCK_Pin;//input nopull
	HAL_GPIO_Init(eSPI_SCK_GPIO_Port, &GPIO_InitStructIn);
	
	GPIO_InitStructIn.Pin = eSPI_MISO_Pin;//input nopull
	HAL_GPIO_Init(eSPI_MISO_GPIO_Port, &GPIO_InitStructIn);
	
	GPIO_InitStructIn.Pin = MCO_CAM_Pin;//input nopull
	HAL_GPIO_Init(MCO_CAM_GPIO_Port, &GPIO_InitStructIn);
	
	GPIO_InitStructIn.Pin = SCCB_SDA_Pin;//input nopull
	HAL_GPIO_Init(SCCB_SDA_GPIO_Port, &GPIO_InitStructIn);
	
	GPIO_InitStructIn.Pin = SCCB_SCL_Pin;//input nopull
	HAL_GPIO_Init(SCCB_SCL_GPIO_Port, &GPIO_InitStructIn);
}

/* TIM2 init function */
//TIM2 is a global timer, used for tic and toc
//each tic is 0.1us
static void CoS_TIM2_Init(TIM_HandleTypeDef* hObj)
{
	TIM_ClockConfigTypeDef sClockSourceConfig;
	TIM_MasterConfigTypeDef sMasterConfig;
	__TIM2_CLK_ENABLE();
	hObj->Instance = TIM2;
	hObj->Init.Prescaler = 8;
	hObj->Init.CounterMode = TIM_COUNTERMODE_UP;
	hObj->Init.Period = 0xffffffff;
	hObj->Init.ClockDivision = TIM_CLOCKDIVISION_DIV1;
	if (HAL_TIM_Base_Init(hObj) != HAL_OK)
		Error_Handler();
	sClockSourceConfig.ClockSource = TIM_CLOCKSOURCE_INTERNAL;
	if (HAL_TIM_ConfigClockSource(hObj, &sClockSourceConfig) != HAL_OK)
		Error_Handler();
	sMasterConfig.MasterOutputTrigger = TIM_TRGO_RESET;
	sMasterConfig.MasterSlaveMode = TIM_MASTERSLAVEMODE_DISABLE;
	if (HAL_TIMEx_MasterConfigSynchronization(hObj, &sMasterConfig) != HAL_OK)
		Error_Handler();
	HAL_TIM_Base_Start(hObj);
}


/* TIM3 init function */
//TIM3 is a timer used for background task calls (heartbeat)
//each interrupt is one heartbeat in us
static void CoS_TIM3_Init(TIM_HandleTypeDef* hObj)
{
	TIM_ClockConfigTypeDef sClockSourceConfig;
	TIM_MasterConfigTypeDef sMasterConfig;
	__TIM3_CLK_ENABLE();
	hObj->Instance = TIM3;
	hObj->Init.Prescaler = 89;
	hObj->Init.CounterMode = TIM_COUNTERMODE_UP;
	hObj->Init.Period = Heartbeat_us;
	hObj->Init.ClockDivision = TIM_CLOCKDIVISION_DIV1;
	if (HAL_TIM_Base_Init(hObj) != HAL_OK)
		Error_Handler();
	sClockSourceConfig.ClockSource = TIM_CLOCKSOURCE_INTERNAL;
	if (HAL_TIM_ConfigClockSource(hObj, &sClockSourceConfig) != HAL_OK)
		Error_Handler();
	sMasterConfig.MasterOutputTrigger = TIM_TRGO_RESET;
	sMasterConfig.MasterSlaveMode = TIM_MASTERSLAVEMODE_DISABLE;
	if (HAL_TIMEx_MasterConfigSynchronization(hObj, &sMasterConfig) != HAL_OK)
		Error_Handler();
	HAL_TIM_Base_Start(hObj);
	HAL_NVIC_SetPriority(hBIO1->BackgroundIRQn, 7, 1);//lowest priority
	__HAL_TIM_ENABLE_IT(hObj,TIM_IT_UPDATE);
	BackgroundTask_Disable();	

}



static void CoS_BIO_LED_Init(void)
{
	GPIO_InitTypeDef GPIO_InitStruct;
	__HAL_RCC_GPIOA_CLK_ENABLE();
	__HAL_RCC_GPIOC_CLK_ENABLE();
	__HAL_RCC_GPIOD_CLK_ENABLE();
	GPIO_InitStruct.Pin = LED1_Pin;
	GPIO_InitStruct.Mode = GPIO_MODE_OUTPUT_PP;
	GPIO_InitStruct.Pull = GPIO_PULLUP;
	GPIO_InitStruct.Speed = GPIO_SPEED_FREQ_MEDIUM;
	HAL_GPIO_Init(LED1_GPIO_Port, &GPIO_InitStruct);
	
	GPIO_InitStruct.Pin = LED2_Pin;
	HAL_GPIO_Init(LED2_GPIO_Port, &GPIO_InitStruct);
	
	GPIO_InitStruct.Pin = LED3_Pin;
	HAL_GPIO_Init(LED3_GPIO_Port, &GPIO_InitStruct);
}

static void CoS_BIO_Keys_Init(void)
{
	GPIO_InitTypeDef GPIO_InitStruct;
	GPIO_InitStruct.Pin = S1_Pin;
	GPIO_InitStruct.Mode = GPIO_MODE_IT_RISING_FALLING;
	GPIO_InitStruct.Pull = GPIO_PULLUP;
	GPIO_InitStruct.Speed = GPIO_SPEED_FREQ_MEDIUM;
	HAL_GPIO_Init(S1_GPIO_Port, &GPIO_InitStruct);
	
	GPIO_InitStruct.Pin = S2_Pin;
	HAL_GPIO_Init(S2_GPIO_Port, &GPIO_InitStruct);
	
	/* EXTI interrupt init*/
	HAL_NVIC_SetPriority(EXTI0_IRQn, 7, 0);
	HAL_NVIC_SetPriority(EXTI1_IRQn, 7, 0);
	HAL_NVIC_EnableIRQ(EXTI0_IRQn);
	HAL_NVIC_EnableIRQ(EXTI1_IRQn);	
}


/**
  * @brief  Initializes the Basic IO by default setting
  * @param  hObj: pointer to a CoS_BIO_HandleTypeDef structure that contains
  *                the configuration information for BIO.
  * @retval CoS status
  */
CoS_StatusTypeDef CoS_BIO_Init(CoS_BIO_HandleTypeDef *hObj,CoS_BIO_InitTypeDef* Init)
{
	/* Check the handle allocation */
	if(hObj == NULL)
	{
		return CoS_STATE_ERROR;
	}
	if( Init->htim2 == NULL ||
		Init->htim3 == NULL ||
		Init->htim5 == NULL || 
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
		hObj->htim2=Init->htim2;
		hObj->htim3=Init->htim3;
		hObj->htim5=Init->htim5;
		hBIO1=hObj->Instance;
		pTIM3Tick=&(hObj->Instance->TIM3TickCount);
	}
	/* Change peripheral state */
	hObj->iState = CoS_STATE_BUSY;
	/* Init the low level hardware */
		CoS_SystemClock_Config();
	if (SysTick_Config(SystemCoreClock / 1000))//1ms
		while (1);	
	HAL_NVIC_Init();
	
	hObj->Instance->keys[0].Port=S1_GPIO_Port;
	hObj->Instance->keys[0].Pin=S1_Pin;
	hObj->Instance->keys[1].Port=S2_GPIO_Port;
	hObj->Instance->keys[1].Pin=S2_Pin;
	hObj->Instance->LEDs[0].Port=LED1_GPIO_Port;
	hObj->Instance->LEDs[0].Pin=LED1_Pin;
	hObj->Instance->LEDs[1].Port=LED2_GPIO_Port;
	hObj->Instance->LEDs[1].Pin=LED2_Pin;
	hObj->Instance->LEDs[2].Port=LED3_GPIO_Port;
	hObj->Instance->LEDs[2].Pin=LED3_Pin;
	hObj->Instance->BackgroundIRQn=Init->BackgroundIRQn;
	hObj->Instance->presstime_threshold=5000;
	CoS_GPIO_Init();
	CoS_TIM2_Init(hObj->htim2);
	CoS_TIM3_Init(hObj->htim3);
	CoS_BIO_LED_Init();
	CoS_BIO_Keys_Init();
	CoS_BIO_Toggle(&(hBIO1->LEDs[LED1]));
	LED1_ON;
	/* Initialize the Instance state */
	hObj->iState = CoS_STATE_OK;
	__CoS_UNLOCK(hObj);
	/* Return function status */
	return CoS_STATE_OK;
}

CoS_StatusTypeDef CoS_BIO_DeInit(CoS_BIO_HandleTypeDef *hObj)
{
	/* Check the handle allocation */
	if(hObj == NULL)
	{
		return CoS_STATE_ERROR;
	}
	hObj->iState = CoS_STATE_BUSY;

	/*
	TODO:Deinit_IO(void);
	*/
	/* Initialize the Instance state */
	hObj->iState = CoS_STATE_UNSET;
	__CoS_UNLOCK(hObj);
	/* Return function status */
	return CoS_STATE_OK;
}


void CoS_BIO_write(BIO_IO_TypdeDef* hIO,GPIO_PinState sta)
{
	HAL_GPIO_WritePin(hIO->Port, hIO->Pin, (GPIO_PinState)sta); 
}

void CoS_BIO_Toggle(BIO_IO_TypdeDef* hIO)
{
	HAL_GPIO_TogglePin(hIO->Port, hIO->Pin); 
}

GPIO_PinState CoS_BIO_read(BIO_IO_TypdeDef* hIO)
{
	return (GPIO_PinState)HAL_GPIO_ReadPin(hIO->Port, hIO->Pin); 
}

CoS_StatusTypeDef CoS_BIO_enable(CoS_BIO_HandleTypeDef* hObj)
{
		/* Check the handle allocation */
	if(hObj == NULL)
	{
		return CoS_STATE_ERROR;
	}
	hObj->Instance->CC.enable=1;
	return CoS_STATE_OK;
	
}

CoS_StatusTypeDef CoS_BIO_disable(CoS_BIO_HandleTypeDef* hObj)
{
		/* Check the handle allocation */
	if(hObj == NULL)
	{
		return CoS_STATE_ERROR;
	}
	hObj->Instance->CC.enable=0;
	return CoS_STATE_OK;
	
}


CoS_StatusTypeDef CoS_BIO_control(CoS_BIO_HandleTypeDef* hObj)
{
		/* Check the handle allocation */
	if(hObj == NULL)
	{
		return CoS_STATE_ERROR;
	}
	hObj->iState = CoS_STATE_BUSY;
	if (hObj->Instance->CC.enable==0)
		CoS_BIO_disable(hObj);
	else 
		CoS_BIO_enable(hObj);
	
//	if (hObj->CC->reinit)
//		CoS_BIO_Init(hObj);
	hObj->iState = CoS_STATE_OK;
	return CoS_STATE_OK;
}

 __weak void CoS_TIM3_IRQ_Task(CoS_BIO_HandleTypeDef* hObj)
{
	//CoS_LED_Update(hCoS->hHCLED);
	//UT_TermDeal();
}


static u8* const BIO_modeStr[4]={"Input","Output","Alternate","Analog"};
static u8* const BIO_spedStr[4]={"Low speed (2MHz)","Med speed(12.5-50MHz)","High speed(25-100MHz)","Very high speed(50-200MHz)"};
static u8* const BIO_pullStr[4]={"No pull","Pull up","Pull down","Researved"};
static u8* const BIO_ptypeStr[2]={"PP","OD"};
u8 CoS_Pin_StateReport(BIO_IO_TypdeDef* hPin,u8 trans_mode_blocking)
{
	u8 portname,pin,mode,pp=0,spd,pull,lk,af;
	//uint16_t sta,npin=0;
	uint32_t sta,npin,temp = 0;
	switch ((u32)hPin->Port)
	{
		case (u32)GPIOA:
		{
			portname='A';
			break;
		}
		case (u32)GPIOB:
		{
			portname='B';
			break;
		}
		case (u32)GPIOC:
		{
			portname='C';
			break;
		}
		case (u32)GPIOD:
		{
			portname='D';
			break;
		}
		case (u32)GPIOE:
		{
			portname='E';
			break;
		}
		case (u32)GPIOF:
		{
			portname='F';
			break;
		}
		case (u32)GPIOG:
		{
			portname='G';
			break;
		}
		case (u32)GPIOH:
		{
			portname='H';
			break;
		}
		case (u32)GPIOI:
		{
			portname='I';
			break;
		}
		case (u32)GPIOJ:
		{
			portname='J';
			break;
		}
		case (u32)GPIOK:
		{
			portname='K';
			break;
		}
		default:
		{
			portname='O';
			break;
		}
	}
	
	npin=hPin->Pin;
	pin=0;
	while (npin>>pin)
	{
		pin++;
	}
	pin--;
	//printf("Pin:%d\r\n",pin);
	temp=hPin->Port->MODER;
	mode=(temp>>(pin<<1))&0x03;
	temp=hPin->Port->OSPEEDR;
	spd=(temp>>(pin<<1))&0x03;
	temp=hPin->Port->PUPDR;
	pull=(temp>>(pin<<1))&0x03;
	temp=hPin->Port->OTYPER;
	if (temp&npin)
		pp=1;
	if (pin>7)
	{
		temp=hPin->Port->AFR[1];
		af=(temp>>((pin-8)<<2))&0x0f;
	}
	else
	{
		temp=hPin->Port->AFR[0];
		af=(temp>>((pin)<<2))&0x0f;
	}
	printf("{Pin:P%c%d, MODE:%s, Speed:%s, Pull:%s(mode:%s), AF:%d}\r\n",
	portname,pin,BIO_modeStr[mode],BIO_spedStr[spd],BIO_pullStr[pull],BIO_ptypeStr[pp],af);
	return 0;
}

static uint32_t const BIO_modes[12]={GPIO_MODE_INPUT,GPIO_MODE_OUTPUT_PP,GPIO_MODE_OUTPUT_OD,\
GPIO_MODE_AF_PP,GPIO_MODE_AF_OD,GPIO_MODE_ANALOG,GPIO_MODE_IT_RISING,GPIO_MODE_IT_FALLING,GPIO_MODE_IT_RISING_FALLING,\
GPIO_MODE_EVT_RISING,GPIO_MODE_EVT_FALLING,GPIO_MODE_EVT_RISING_FALLING};
static uint32_t const BIO_speeds[4]={GPIO_SPEED_FREQ_LOW,GPIO_SPEED_FREQ_MEDIUM,\
GPIO_SPEED_FREQ_HIGH,GPIO_SPEED_FREQ_VERY_HIGH};
static uint32_t const BIO_pulls[3]={GPIO_NOPULL,GPIO_PULLUP,GPIO_PULLDOWN};


u8 CoS_Pin_StateSave(BIO_IO_TypdeDef* hPin,GPIO_InitTypeDef *hStruct)
{
//	u8 portname,pin,mode,pp=0,spd,pull,lk,af;
//	//uint16_t sta,npin=0;
//	uint32_t sta,npin,temp = 0;
//	
//	npin=hPin->Pin;
//	pin=0;
//	while (npin>>pin)
//	{
//		pin++;
//	}
//	pin--;
//	temp=hPin->Port->MODER;
//	mode=(temp>>(pin<<1))&0x03;
//	temp=hPin->Port->OSPEEDR;
//	spd=(temp>>(pin<<1))&0x03;
//	temp=hPin->Port->PUPDR;
//	pull=(temp>>(pin<<1))&0x03;
//	temp=hPin->Port->OTYPER;
//	if (temp&npin)
//		pp=1;
//	if (pin>7)
//	{
//		temp=hPin->Port->AFR[1];
//		af=(temp>>((pin-8)<<2))&0x0f;
//	}
//	else
//	{
//		temp=hPin->Port->AFR[0];
//		af=(temp>>((pin)<<2))&0x0f;
//	}
//	printf("{Pin:P%c%d, MODE:%s, Speed:%s, Pull:%s(mode:%s), AF:%d}\r\n",
//	portname,pin,BIO_modeStr[mode],BIO_spedStr[spd],BIO_pullStr[pull],BIO_ptypeStr[pp],af);
	return 0;
}
















uint8_t CoS_Restart(u8 keynum)
{
	static uint8_t status=0;
	if (keynum==82)
	{
		if (status==0)
		{
			status=82;
			return 1;
		}
		if (status==82)
		{
			NVIC_SystemReset();//Reset all
			return 2;
		}
	}
	else
	{
		status=0;
		return 0;
	}
	return 0;
}



uint32_t CoS_MCUID(char i)
{
	switch (i)
	{
		case 0:
			return *(__IO u32*)(0x1FFF7A10);
		case 1:
			return *(__IO u32*)(0x1FFF7A14);
		case 2:
			return *(__IO u32*)(0x1FFF7A18);
		default :
			return *(__IO u32*)(0xE0042000);
	}
}