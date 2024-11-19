#include "coliasSense_Camera.h"
#include "coliasSense_board.h"
#include "delay.h"


/* original ones
static const u8 OV7670InitReg[]={
	0x11,0x80,//CLKRC
	0x12,0x08,//color bar YUV
	//0x13,0x8b,//Disable Auto Gain Control
	0x13,0x8f,//Enable Auto Gain Control
	0x0C,0x0C,//COM3 OUTPUT/SCALE
	0x3E,0x12,//COM14  PCLK,scale,DCW
	0x70,0x3A,//SCALING_XSC testpattern 
	0x71,0x35,//35 SCALING_YSC testpattern 
	0x72,0x22,//SCALING_DCWCTR
	0x73,0xF2,//SCALING_PCLK_DIV
	0xa2,0x2a,//SCALING_PCLK_DELAY
	0x1e,0x11,//MVFP mirror/flip
	0x15,0x00,//COM10 HREF/PCLK/VSYNC pattern
	0x3A,0x15,//TSLB Line Buffer Test Option
	0x67,0xaa,//MANU
	0x68,0x55,//MANV
	
};*/


static const u8 OV7670InitReg[]={
	0x11,0x80,//CLKRC
	0x12,0x08,//color bar YUV
	0x13,0x8b,//Disable Auto Gain Control
	//0x13,0x8f,//Enable Auto Gain Control
	0x14,0x6a,//Set auto gain ceiling to 128x(was 32x for 0x4a)
	0x0C,0x0C,//COM3 OUTPUT/SCALE
	0x3E,0x12,//COM14  PCLK,scale,DCW
	0x70,0x3A,//SCALING_XSC testpattern 
	0x71,0x35,//35 SCALING_YSC testpattern 
	0x72,0x22,//SCALING_DCWCTR
	0x73,0xF2,//SCALING_PCLK_DIV
	0xa2,0x2a,//SCALING_PCLK_DELAY
	0x1e,0x11,//MVFP mirror/flip
	0x15,0x00,//COM10 HREF/PCLK/VSYNC pattern
	0x3A,0x05,//TSLB Line Buffer Test Option,disabled MANUV
	//0x3A,0x15,//TSLB Line Buffer Test Option,enabled MANUV
	0x67,0xaa,//MANU
	0x68,0x55,//MANV
	0x00,0x08, //gain
	
};


static void CoS_DCMI_DMAXferCplt(DMA_HandleTypeDef *hdma);
//Rewrite HAL DCMI_DMAXferCplt, because the HAL version aotumatically
//enables the frame IT, which is not wanted in our application.
static void CoS_DCMI_DMAXferCplt(DMA_HandleTypeDef *hdma)
{
  uint32_t tmp = 0U;
 
  DCMI_HandleTypeDef* hdcmi = ( DCMI_HandleTypeDef* )((DMA_HandleTypeDef* )hdma)->Parent;
  
  if(hdcmi->XferCount != 0)
  {
    /* Update memory 0 address location */
    tmp = ((hdcmi->DMA_Handle->Instance->CR) & DMA_SxCR_CT);
    if(((hdcmi->XferCount % 2U) == 0U) && (tmp != 0U))
    {
      tmp = hdcmi->DMA_Handle->Instance->M0AR;
      HAL_DMAEx_ChangeMemory(hdcmi->DMA_Handle, (tmp + (8U*hdcmi->XferSize)), MEMORY0);
      hdcmi->XferCount--;
    }
    /* Update memory 1 address location */
    else if((hdcmi->DMA_Handle->Instance->CR & DMA_SxCR_CT) == 0U)
    {
      tmp = hdcmi->DMA_Handle->Instance->M1AR;
      HAL_DMAEx_ChangeMemory(hdcmi->DMA_Handle, (tmp + (8U*hdcmi->XferSize)), MEMORY1);
      hdcmi->XferCount--;
    }
  }
  /* Update memory 0 address location */
  else if((hdcmi->DMA_Handle->Instance->CR & DMA_SxCR_CT) != 0U)
  {
    hdcmi->DMA_Handle->Instance->M0AR = hdcmi->pBuffPtr;
  }
  /* Update memory 1 address location */
  else if((hdcmi->DMA_Handle->Instance->CR & DMA_SxCR_CT) == 0U)
  {
    tmp = hdcmi->pBuffPtr;
    hdcmi->DMA_Handle->Instance->M1AR = (tmp + (4U*hdcmi->XferSize));
    hdcmi->XferCount = hdcmi->XferTransferNumber;
  }
  
  /* Check if the frame is transferred */
  if(hdcmi->XferCount == hdcmi->XferTransferNumber)
  {
    /* Enable the Frame interrupt */
    //__HAL_DCMI_ENABLE_IT(hdcmi, DCMI_IT_FRAME);
    
    /* When snapshot mode, set dcmi state to ready */
    if((hdcmi->Instance->CR & DCMI_CR_CM) == DCMI_MODE_SNAPSHOT)
    {  
      hdcmi->State= HAL_DCMI_STATE_READY;
    }
  }
}
static void CAM_IT_Init(void) //interrupt
{	
	
}



static void CAM_DCMI_IO_Init(void)
{	
	GPIO_InitTypeDef GPIO_InitStruct;
	
	/**DCMI GPIO Configuration    
	PE4     ------> DCMI_D4
	PE5     ------> DCMI_D6
	PE6     ------> DCMI_D7
	PA4     ------> DCMI_HSYNC
	PA6     ------> DCMI_PIXCK
	PC6     ------> DCMI_D0
	PC7     ------> DCMI_D1
	PB6     ------> DCMI_D5
	PB7     ------> DCMI_VSYNC
	PE0     ------> DCMI_D2
	PE1     ------> DCMI_D3 
	*/
	GPIO_InitStruct.Pin = 	GPIO_PIN_4|GPIO_PIN_5|
							GPIO_PIN_6|GPIO_PIN_0|GPIO_PIN_1;
	GPIO_InitStruct.Mode = GPIO_MODE_AF_PP;
	GPIO_InitStruct.Pull = GPIO_PULLUP;
	GPIO_InitStruct.Speed = GPIO_SPEED_FREQ_VERY_HIGH;
	GPIO_InitStruct.Alternate = GPIO_AF13_DCMI;
	HAL_GPIO_Init(GPIOE, &GPIO_InitStruct);

	GPIO_InitStruct.Pin = GPIO_PIN_4|GPIO_PIN_6;
	GPIO_InitStruct.Mode = GPIO_MODE_AF_PP;
	GPIO_InitStruct.Pull = GPIO_PULLUP;
	GPIO_InitStruct.Speed = GPIO_SPEED_FREQ_VERY_HIGH;
	GPIO_InitStruct.Alternate = GPIO_AF13_DCMI;
	HAL_GPIO_Init(GPIOA, &GPIO_InitStruct);

	GPIO_InitStruct.Pin = GPIO_PIN_6|GPIO_PIN_7;
	GPIO_InitStruct.Mode = GPIO_MODE_AF_PP;
	GPIO_InitStruct.Pull = GPIO_PULLUP;
	GPIO_InitStruct.Speed = GPIO_SPEED_FREQ_VERY_HIGH;
	GPIO_InitStruct.Alternate = GPIO_AF13_DCMI;
	HAL_GPIO_Init(GPIOC, &GPIO_InitStruct);

	GPIO_InitStruct.Pin = GPIO_PIN_6|GPIO_PIN_7;
	GPIO_InitStruct.Mode = GPIO_MODE_AF_PP;
	GPIO_InitStruct.Pull = GPIO_PULLUP;
	GPIO_InitStruct.Speed = GPIO_SPEED_FREQ_VERY_HIGH;
	GPIO_InitStruct.Alternate = GPIO_AF13_DCMI;
	HAL_GPIO_Init(GPIOB, &GPIO_InitStruct);
	
	
}

static void CAM_DCMI_DMA_Init(DCMI_HandleTypeDef* hdcmi,DMA_HandleTypeDef* hdma_dcmi)
{	
	//DMA_HandleTypeDef* hdma_dcmi=hObj->Instance->hdcmi->DMA_Handle;
	__HAL_RCC_DMA2_CLK_ENABLE();
	
	/* Peripheral DMA init*/

	hdma_dcmi->Instance = DMA2_Stream1;
	hdma_dcmi->Init.Channel = DMA_CHANNEL_1;
	hdma_dcmi->Init.Direction = DMA_PERIPH_TO_MEMORY;
	hdma_dcmi->Init.PeriphInc = DMA_PINC_DISABLE;
	hdma_dcmi->Init.MemInc = DMA_MINC_ENABLE;
	hdma_dcmi->Init.PeriphDataAlignment = DMA_PDATAALIGN_WORD;
	hdma_dcmi->Init.MemDataAlignment = DMA_MDATAALIGN_WORD;
	hdma_dcmi->Init.Mode = DMA_CIRCULAR;
	hdma_dcmi->Init.Priority = DMA_PRIORITY_HIGH;
	hdma_dcmi->Init.FIFOMode = DMA_FIFOMODE_ENABLE;
	hdma_dcmi->Init.FIFOThreshold = DMA_FIFO_THRESHOLD_FULL;
	hdma_dcmi->Init.MemBurst = DMA_MBURST_SINGLE;
	hdma_dcmi->Init.PeriphBurst = DMA_PBURST_SINGLE;
	if (HAL_DMA_Init(hdma_dcmi) != HAL_OK)
	{
		Error_Handler();
	}

	__HAL_LINKDMA(hdcmi,DMA_Handle,*hdma_dcmi);
	HAL_NVIC_SetPriority(DMA2_Stream1_IRQn, 2, 0);
	//HAL_NVIC_EnableIRQ(DMA2_Stream1_IRQn);
}

void CAM_DCMI_Init(CoS_CAM_HandleTypeDef* hObj,CoS_CAM_InitTypeDef* hCInit) 
{	
	DCMI_HandleTypeDef* hdcmi;
	hObj->hdcmi=hCInit->hdcmi;
	hObj->hdcmi->DMA_Handle=hCInit->hdma_dcmi;
	hdcmi=hObj->hdcmi;
	__HAL_RCC_DCMI_CLK_ENABLE();
	CAM_DCMI_IO_Init();
	CAM_DCMI_DMA_Init(hdcmi,hdcmi->DMA_Handle);
	
	hdcmi->Instance = DCMI;
	hdcmi->Init.SynchroMode = DCMI_SYNCHRO_HARDWARE;
	hdcmi->Init.PCKPolarity = DCMI_PCKPOLARITY_RISING;
	hdcmi->Init.VSPolarity = DCMI_VSPOLARITY_HIGH;
	hdcmi->Init.HSPolarity = DCMI_HSPOLARITY_LOW;
	hdcmi->Init.CaptureRate = DCMI_CR_ALL_FRAME;
	hdcmi->Init.ExtendedDataMode = DCMI_EXTEND_DATA_8B;
	hdcmi->Init.JPEGMode = DCMI_JPEG_DISABLE;
	if (HAL_DCMI_Init(hdcmi) != HAL_OK)
	{
		Error_Handler();
	}
	__HAL_DCMI_DISABLE_IT(hdcmi, DCMI_IT_LINE | DCMI_IT_FRAME | DCMI_IT_ERR | DCMI_IT_OVR);
	HAL_NVIC_SetPriority(DCMI_IRQn, 2, 1);
	//HAL_NVIC_EnableIRQ(DCMI_IRQn);
	
}

static CoS_StatusTypeDef CAM_SCCB_Init(CoS_CAM_HandleTypeDef* hObj,CoS_CAM_InitTypeDef* hCInit)
{
	CoS_sI2C_InitTypeDef Init;
	Init.Dev_Addr=A_OV7670W;
	Init.Instance=&CoS_sI2C_Camera;
	Init.mode=sI2CMode_SCCB;
	Init.nDelay=0x0e;
	Init.PinSCL.Port=SCCB_SCL_GPIO_Port;
	Init.PinSCL.Pin=SCCB_SCL_Pin;
	Init.PinSDA.Port=SCCB_SDA_GPIO_Port;
	Init.PinSDA.Pin=SCCB_SDA_Pin;
	hObj->hSSCB=hCInit->hSSCB;
	return CoS_sI2C_Init(hObj->hSSCB,&Init);
}


CoS_StatusTypeDef CoS_CAM_Init(CoS_CAM_HandleTypeDef *hObj,CoS_CAM_InitTypeDef* hInit)
{
	uint8_t id,n,ad,val;
	/* Check the handle allocation */
	if(hObj == NULL)
	{
		return CoS_STATE_ERROR;
	}
	if( 
		hInit->Instance==NULL)
	{
		return CoS_STATE_ERROR;
	}
	__CoS_LOCK(hObj);
	if(hObj->iState == CoS_STATE_UNSET)
	{
		hObj->Lock = CoS_UNLOCKED;
		hObj->Instance=hInit->Instance;
	}/* Change peripheral state */
	hObj->iState = CoS_STATE_BUSY;
	/* Init the low level hardware */
	//CAM_IO_Init(hObj,Init);
	
	printf("{}Camera init...\r\n");
	CAM_SCCB_Init(hObj,hInit);
	HAL_RCC_MCOConfig(RCC_MCO1,RCC_MCO1SOURCE_HSE,RCC_MCODIV_1);
	delay_ms(5);
	
	CAM_IT_Init();
	CoS_CAM_cameraReset(hObj);
	delay_ms(5);
	id =CoS_CAM_ReadID(hObj);
	if (id!=ID_OV7670)
	{
		printf("{}Got camera ID 0x%02x. Camera OV7670 not found.\r\n",id);
		hObj->iState = CoS_STATE_UNSET;
		__CoS_UNLOCK(hObj);
		return CoS_STATE_ERROR;
	}
	else
	{
		printf("{}Got camera ID 0x%02x. Camera OV7670 found.\r\n",id);
		printf("{}writting camera init registers\r\n");
		for(n=0;n<sizeof(OV7670InitReg)/2;n++)
		{
			ad = OV7670InitReg[n<<1];
			val= OV7670InitReg[1+(n<<1)];
			CoS_sI2C_WriteReg(hObj->hSSCB,ad,val);
			printf("{}A:0x%02X\tD:0x%02X\r\n",ad,val);
			delay_ms(20);
		}
		CAM_DCMI_Init(hObj,hInit);
		hObj->Instance->cam_DMACount=1;
		hObj->Instance->cam_hFrameCountount=1;
		/* Initialize the Instance state */
		hObj->iState = CoS_STATE_OK;
		printf("{}Camera init OK.\r\n");
		CoS_camera_StateReport(hObj,1);
		//hObj->Instance->CC.state=CoS_STATE_OK;
		__CoS_UNLOCK(hObj);
		/* Return function status */
		return CoS_STATE_OK;
	}
}




void HAL_DCMI_VsyncEventCallback(DCMI_HandleTypeDef *hdcmi)
{
	//CoS_BIO_Toggle(&(hBIO1->LEDs[LED3]));
}


CoS_StatusTypeDef CoS_CAM_cameraReset(CoS_CAM_HandleTypeDef *hObj)
{
	
	return CoS_sI2C_WriteReg(hObj->hSSCB,0x12,0x80);
}

CoS_StatusTypeDef CoS_CAM_regReset(CoS_CAM_HandleTypeDef *hObj)
{
	
	return CoS_sI2C_WriteReg(hObj->hSSCB,0x12,0x80);
}

uint8_t CoS_CAM_ReadID(CoS_CAM_HandleTypeDef *hObj)
{
	__IO u8 temp=0,sta=0;
	sta=CoS_sI2C_ReadReg(hObj->hSSCB,0x0b,&temp);
	if (sta==CoS_STATE_OK)
		return temp;
	else
		return 0;

}

CoS_StatusTypeDef CoS_CAM_Start(CoS_CAM_HandleTypeDef *hObj)
{
	if(hObj == NULL)
	{
		return CoS_STATE_ERROR;
	}
	__CoS_LOCK(hObj);
	
	HAL_NVIC_EnableIRQ(DCMI_IRQn);
	__HAL_DCMI_ENABLE_IT(hCoS->hdcmi, DCMI_IT_VSYNC);
	HAL_NVIC_EnableIRQ(DMA2_Stream1_IRQn);
	__CoS_UNLOCK(hObj);
	/* Return function status */
	return CoS_STATE_OK;
}

CoS_StatusTypeDef CoS_CAM_DMASet(CoS_CAM_HandleTypeDef *hObj,uint32_t dmabase,uint32_t length)
{
	if(hObj == NULL)
	{
		return CoS_STATE_ERROR;
	}
	__CoS_LOCK(hObj);
	//hdcmi->DMA_Handle->XferM1CpltCallback = DCMI_DMAXferCplt;
	HAL_DCMI_Start_DMA(hObj->hdcmi,DCMI_MODE_CONTINUOUS,dmabase,length);
	hObj->hdcmi->DMA_Handle->XferCpltCallback = CoS_DCMI_DMAXferCplt;
	hObj->hdcmi->DMA_Handle->XferM1CpltCallback=CoS_DCMI_DMAXferCplt;
	__CoS_UNLOCK(hObj);
	/* Return function status */
	return CoS_STATE_OK;
}



CoS_StatusTypeDef CoS_CAM_Stop(CoS_CAM_HandleTypeDef *hObj)
{
	if(hObj == NULL)
	{
		return CoS_STATE_ERROR;
	}
	__CoS_LOCK(hObj);
	HAL_DCMI_Stop(hObj->hdcmi);
	HAL_NVIC_DisableIRQ(DMA2_Stream1_IRQn);
	HAL_NVIC_DisableIRQ(DCMI_IRQn);
	__CoS_UNLOCK(hObj);
	/* Return function status */
	return CoS_STATE_OK;
}




void CoS_camera_StateReport(CoS_CAM_HandleTypeDef *hObj,u8 trans_mode_blocking)
{
	printf("{  camera Initialize state:%d}\r\n",hObj->iState);
	if (hObj->iState!=CoS_STATE_OK)
		return;
	printf("{  SCCB Initialize state:%d}\r\n",hObj->hSSCB->iState);
	printf("{  DCMI Initialize state:%d}\r\n",hObj->hdcmi->State);
	printf("{  DMA of DCMI Initialize state:%d}\r\n",hObj->hdcmi->DMA_Handle->State);
	printf("{  DMA transfer count:%d}\r\n",hObj->Instance->cam_DMACount);
	printf("{  image captured count:%d}\r\n",hObj->Instance->cam_hFrameCountount);
}



CoS_StatusTypeDef CoS_camera_MANUV_Enable(CoS_CAM_HandleTypeDef *hObj)
{
	CoS_StatusTypeDef a;
	if(hObj == NULL)
	{
		return CoS_STATE_ERROR;
	}
	__CoS_LOCK(hObj);
	a= CoS_sI2C_WriteReg(hObj->hSSCB,0x3A,0x15);
	__CoS_UNLOCK(hObj);
	/* Return function status */
	return a;
}

CoS_StatusTypeDef CoS_camera_MANUV_Disable(CoS_CAM_HandleTypeDef *hObj)
{
	CoS_StatusTypeDef a;
	if(hObj == NULL)
	{
		return CoS_STATE_ERROR;
	}
	__CoS_LOCK(hObj);
	a= CoS_sI2C_WriteReg(hObj->hSSCB,0x3A,0x05);
	__CoS_UNLOCK(hObj);
	/* Return function status */
	return a;
}

































