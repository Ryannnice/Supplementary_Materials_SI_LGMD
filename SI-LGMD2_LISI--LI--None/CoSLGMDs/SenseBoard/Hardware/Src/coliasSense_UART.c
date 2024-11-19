#include "coliasSense_UART.h"
#include "coliasSense_board.h"
#include "delay.h"
#include "stm32f4xx_hal.h"

#define PHASE_HEAD 1U
#define PHASE_END 2U
#define PHASE_DATA 4U
#define PHASE_DATAE 8U
#define PHASE_IDLE 0U

#define PHASE_RHEAD 0x10U
#define PHASE_REND 0x20U
#define PHASE_RDATA 0x40U
#define PHASE_RDATAE 0x80U
#define PHASE_RIDLE 0x00U


static u8* mode_nameStr[]={"UART_MODE_UART","UART_MODE_CPtP","UART_MODE_PPP"};
static void UART_CPtP_frameStart(CoS_UART_TypeDef* hObj);
static void UART_CPtP_frameEnd(  CoS_UART_TypeDef* hObj);
static void UART_transmit(       CoS_UART_TypeDef *hObj,uint8_t ch);

static void CoS_UART_ByteDeal(CoS_UART_HandleTypeDef *hObj,u8 rx);
static void CoS_UART_PeripheralInit(CoS_UART_HandleTypeDef *hObj,
						CoS_UART_InitTypeDef* Init);
static u8 CoS_UART_RegInit(CoS_UART_HandleTypeDef *hObj,
						CoS_UART_InitTypeDef* Init);

static void CoS_UART_TxStart(CoS_UART_HandleTypeDef *hObj);

#ifdef USE_printf
	volatile uint8_t PRINTF_INIT_OK=0;
	#ifdef USE_printf_USART1
		#pragma import(__use_no_semihosting)															
		struct __FILE 
		{ 
			int handle; 
		}; 
		FILE __stdout;		
		/*
		_sys_exit(int x) 
		{ 
			x = x; 
		} */
		
		int fputc(int ch, FILE *f)
		{
			if (PRINTF_INIT_OK)
			{
				while((USART1->SR&USART_SR_TC)==RESET);
				USART1->DR = (ch & (uint16_t)0x01FF);
					return ch;
			}
			else
				return ch;
		}
	#else 
		#ifdef USE_printf_UART4
		#pragma import(__use_no_semihosting)															
		struct __FILE 
		{ 
			int handle; 
		}; 
		FILE __stdout;		
		/*
		_sys_exit(int x) 
		{ 
			x = x; 
		} */
		int fputc(int ch, FILE *f)
		{
			while((UART4->SR&USART_SR_TC)==RESET);
			UART4->DR = (ch & (uint16_t)0x01FF);
				return ch;
		}
		#endif
	#endif
#endif
void Uputc1(char ch)
{
	while((USART1->SR&UART_FLAG_TC)==RESET);
		USART1->DR =ch;
}	

void UARTn_putc(USART_TypeDef* UARTn,char ch)
{
	while((UARTn->SR&UART_FLAG_TC)==RESET);
		UARTn->DR =ch;
}

static void UART_CPtP_frameStart(CoS_UART_TypeDef* hObj)
{
//	if (hObj->CC.userdata==PHASE_HEAD_W)
//		UARTn_putc(hObj->huart->Instance,hObj->escape2);
//	else
//		UARTn_putc(hObj->huart->Instance,(hObj->escape2)+1);
	UARTn_putc(hObj->USARTn,hObj->escape1);
}

static void UART_CPtP_frameEnd(CoS_UART_TypeDef* hObj)
{
	UARTn_putc(hObj->USARTn,hObj->escape2);
}

static void UART_transmit(CoS_UART_TypeDef *hObj,uint8_t ch)
{
	if ((ch==hObj->escape1 || ch==hObj->escape2)&& \
		hObj->CC.mode==UART_MODE_CPtP && hObj->CC.userdata==PHASE_DATA)
	{//CPtP mode, phase data, escape hit
		UARTn_putc(hObj->USARTn,hObj->escape1);
		UARTn_putc(hObj->USARTn,ch-0x20);
	}
	else
	{
		UARTn_putc(hObj->USARTn,ch);
	}
}



void UART_transmit_run(CoS_UART_HandleTypeDef *hObj)
{
	CoS_UART_TypeDef* hIns=hObj->Instance;
	u16 n=hIns->Tx_PacketSize;
	u8 tr_end=0;
	u8 sta=0;
	if (hObj->iState!=CoS_STATE_OK)
		return;
	if (!hIns->CC.enable )
		return;
	sta=(hIns->CC.state)&CoS_STATE_BUSY_MASKT;
	if (sta!=CoS_STATE_BUSY_MASKT)
		return;
	if (hObj->Instance->CC.mode==UART_MODE_CPtP)
	{
		for (n=hObj->Instance->Tx_PacketSize;n!=0;n--)
		{
			switch (hIns->CC.userdata&0x0f)//PHASE
			{
				case PHASE_HEAD:
				{
					UART_CPtP_frameStart(hIns);
					hIns->CC.userdata=(hIns->CC.userdata&0xf0)+PHASE_DATA;
					break;
				}
				case PHASE_END:
				{
					UART_CPtP_frameEnd(hIns);
					if (n)
						n=1;
					hIns->CC.userdata=(hIns->CC.userdata&0xf0)+PHASE_IDLE;
					tr_end=1;
					break;
				}
				case PHASE_IDLE:
				{
//					if (*hIns->Tx_Base==((hObj->Instance->escape1)+1))
//						hIns->CC.userdata=(hIns->CC.userdata&0xf0)+PHASE_HEAD;
//					else
						hIns->CC.userdata=(hIns->CC.userdata&0xf0)+PHASE_HEAD;
					n++;
					break;
				}
				case PHASE_DATA:
				default:
				{
					UART_transmit(hIns,hIns->Tx_Base[hIns->Tx_Count++]);
					if (hIns->Tx_Count>=hIns->Tx_Length)
					{
						hIns->CC.userdata=(hIns->CC.userdata&0xf0)+PHASE_END;
					}
					break;
				}
			}
		}
	}
	else//plain
	{
		for (n=hObj->Instance->Tx_PacketSize;n!=0;n--)
		{
			UART_transmit(hIns,hIns->Tx_Base[hIns->Tx_Count++]);
			if (hIns->Tx_Count>=hIns->Tx_Length)
			{
				n=1;
				tr_end=1;
				break;
			}
		}
	}
	if (tr_end)
	{
		hIns->CC.state&=(!CoS_STATE_BUSY_MASKT);
	}
}

CoS_StatusTypeDef CoS_UART_sendFrame(CoS_UART_HandleTypeDef *hObj,uint8_t* Tx_Base,uint16_t Tx_Length)
{
 	CoS_StatusTypeDef sta;
	if(hObj == NULL)
		return CoS_STATE_ERROR;
	if (hObj->iState!=CoS_STATE_OK)
		return CoS_STATE_ERROR;
	if ((hObj->Instance->CC.state)&CoS_STATE_BUSY_MASKT==CoS_STATE_BUSY_MASKT)
		return CoS_STATE_BUSY;
	__CoS_LOCK(hObj);
	hObj->Instance->Tx_Base=Tx_Base;
	hObj->Instance->Tx_Length=Tx_Length;
	hObj->Instance->Tx_Count=0;
	hObj->Instance->CC.state |=CoS_STATE_BUSY_MASKT;
	__CoS_UNLOCK(hObj);
	return CoS_STATE_OK;
}

CoS_StatusTypeDef CoS_UART_TxRST(CoS_UART_HandleTypeDef *hObj)
{
	if(hObj == NULL)
		return CoS_STATE_ERROR;
	__CoS_LOCK(hObj);
	hObj->Instance->Tx_Count=0;
	hObj->Instance->Tx_Length=0;
	hObj->Instance->CC.state &=!CoS_STATE_BUSY_MASKT;
	__CoS_UNLOCK(hObj);
	return CoS_STATE_OK;
}


CoS_StatusTypeDef CoS_UART_Enable(CoS_UART_HandleTypeDef *hObj,uint8_t sta)
{
	if(hObj == NULL)
	{
		return CoS_STATE_ERROR;
	}
	__CoS_LOCK(hObj);
	if (sta)
		hObj->Instance->CC.enable=1;
	else
		hObj->Instance->CC.enable=0;
	__CoS_UNLOCK(hObj);
	return CoS_STATE_OK;
}

void HAL_UART_MspInit(UART_HandleTypeDef* uartHandle)
{

	GPIO_InitTypeDef GPIO_InitStruct;
	if(uartHandle->Instance==UART4)
	{
		/* Peripheral clock enable */
		__HAL_RCC_UART4_CLK_ENABLE();
	
		/**UART4 GPIO Configuration		
		PC10		 ------> UART4_TX
		PC11		 ------> UART4_RX 
		*/
		GPIO_InitStruct.Pin = GPIO_PIN_10|GPIO_PIN_11;
		GPIO_InitStruct.Mode = GPIO_MODE_AF_PP;
		GPIO_InitStruct.Pull = GPIO_PULLUP;
		GPIO_InitStruct.Speed = GPIO_SPEED_FREQ_MEDIUM;
		GPIO_InitStruct.Alternate = GPIO_AF8_UART4;
		HAL_GPIO_Init(GPIOC, &GPIO_InitStruct);

	}
	else if(uartHandle->Instance==USART1)
	{
		/* Peripheral clock enable */
		__HAL_RCC_USART1_CLK_ENABLE();
	
		/**USART1 GPIO Configuration		
		PA9		 ------> USART1_TX
		PA10		 ------> USART1_RX 
		*/
		GPIO_InitStruct.Pin = GPIO_PIN_9|GPIO_PIN_10;
		GPIO_InitStruct.Mode = GPIO_MODE_AF_PP;
		GPIO_InitStruct.Pull = GPIO_PULLUP;
		GPIO_InitStruct.Speed = GPIO_SPEED_FREQ_HIGH;
		GPIO_InitStruct.Alternate = GPIO_AF7_USART1;
		HAL_GPIO_Init(GPIOA, &GPIO_InitStruct);
	}
}



static void CoS_UART_PeripheralInit(CoS_UART_HandleTypeDef *hObj,
						CoS_UART_InitTypeDef* Init)
{
	HAL_UART_DeInit(hObj->huart);
	HAL_UART_MspInit(hObj->huart);
	if (HAL_UART_Init(hObj->huart) != HAL_OK)
		Error_Handler();
	if (hObj->huart->Instance==USART1)
	{
		HAL_NVIC_SetPriority(USART1_IRQn, 2, 1);
		__HAL_UART_ENABLE_IT(hObj->huart,UART_IT_RXNE);
		HAL_NVIC_EnableIRQ(USART1_IRQn);
	}
	else if (hObj->huart->Instance==UART4)
	{
		HAL_NVIC_SetPriority(UART4_IRQn, 2, 1);
		__HAL_UART_ENABLE_IT(hObj->huart,UART_IT_RXNE);
		HAL_NVIC_EnableIRQ(UART4_IRQn);
	}
	
}

static u8 CoS_UART_RegInit(CoS_UART_HandleTypeDef *hObj,
						CoS_UART_InitTypeDef* Init)
{
	hObj->Instance->Rx_Base=Init->RxBufferBase;
	hObj->Instance->Rx_BufSize=Init->RxBufferSize;
	hObj->Instance->Tx_PacketSize=Init->Tx_PacketSize;
	hObj->Instance->USARTn=Init->huart->Instance;
	hObj->Instance->escape1=Init->escape1;
	hObj->Instance->escape2=Init->escape2;
	hObj->Instance->CC.mode=Init->mode;
	hObj->Instance->CC.state=CoS_STATE_OK;
	hObj->Instance->CC.enable=0;
	hObj->Instance->Rx_Base=Init->RxBufferBase;
	hObj->Instance->Rx_BufSize=Init->RxBufferSize;
	if (Init->Readback)
		hObj->Instance->CCByte|=_CoS_UART_RxReadBack;
	if (Init->RxICOv)
		hObj->Instance->CCByte|=_CoS_UART_RxICOv;
	if (Init->RxRAWF)
		hObj->Instance->CCByte|=_CoS_UART_RxRAWF;
	switch ((u32)hObj->Instance->USARTn)
	{
		case (u32)USART1:
			return 1;
		case (u32)USART2:
			return 2;
		case (u32)USART3:
			return 3;
		case (u32)UART4:
			return 4;
		case (u32)UART5:
			return 5;
		case (u32)USART6:
			return 6;
		case (u32)UART7:
			return 7;
		case (u32)UART8:
			return 8;
	}
	return CoS_STATE_OK;
}


CoS_StatusTypeDef CoS_UART_Init(CoS_UART_HandleTypeDef *hObj,
								CoS_UART_InitTypeDef* Init)
{
	u8 Uid;
	/* Check the handle allocation */
	if(hObj == NULL)
	{
		return CoS_STATE_ERROR;
	}
	if( Init->Instance == NULL ||
		Init->huart == NULL ||
		Init->RxBufferBase==NULL||
		Init->RxBufferSize==0||
		Init->Tx_PacketSize==0)
	{
		return CoS_STATE_ERROR;
	}
	__CoS_LOCK(hObj);
	if(hObj->iState == CoS_STATE_UNSET)
	{  
		/* Allocate lock resource and initialize it */
		hObj->Lock = CoS_UNLOCKED;
		hObj->Instance=Init->Instance;
		hObj->huart=Init->huart;
		//hObj->Instance->hUART=Init->hUART;
	}
	/* Change peripheral state */
	hObj->iState = CoS_STATE_BUSY;
	/* Init the low level hardware */
	CoS_UART_PeripheralInit(hObj,Init);
	Uid=CoS_UART_RegInit(hObj,Init);
	if (INIT_INFO_printf)
		printf("{}UART %d Initialize OK.\r\n",Uid);
	if (hObj->Instance->USARTn==USART1)
		PRINTF_INIT_OK=1;
	/* Initialize the Instance state */
	hObj->iState = CoS_STATE_OK;
	CoS_UART_StateReport(hObj,Uid,1);
	__CoS_UNLOCK(hObj);
	/* Return function status */
	return CoS_STATE_OK;
}

void CoS_UART_StateReport(CoS_UART_HandleTypeDef *hObj,u8 Uid,u8 trans_mode_blocking)
{
	
	if (trans_mode_blocking)
	{
		printf("{  UART %d Initialize state:%d}\r\n",Uid,hObj->iState);
		if (hObj->iState!=CoS_STATE_OK)
			return;
		printf("{  mode: %s}\r\n",mode_nameStr[hObj->Instance->CC.mode]);
		printf("{  RxBufferSize: %d}\r\n",hObj->Instance->Rx_BufSize);
		printf("{  TxPacket Size: %d}\r\n",hObj->Instance->Tx_PacketSize);
		printf("{  escape 1: 0x%X}\r\n",hObj->Instance->escape1);
		printf("{  escape 2: 0x%X}\r\n",hObj->Instance->escape2);
	}
	else
		return;
}

void Flush_UART_BUFF(CoS_UART_TypeDef* hInstance)
{
	uint8_t n;
	for (n=0;n<hInstance->Rx_BufSize;n++)
	{
		*(hInstance->Rx_Base+n)=0;
	}
}


CoS_StatusTypeDef CoS_UART_RxRST(CoS_UART_HandleTypeDef *hObj)
{
	CoS_UART_TypeDef* Instance=hObj->Instance;
	__CoS_LOCK(hObj);
	Flush_UART_BUFF(Instance);
	//Instance->RxBusy=0;
	Instance->Rx_Count=0;
	Instance->Rx_OverrunCount=0;
	Instance->RxOK=0;
	__CoS_UNLOCK(hObj);
	return CoS_STATE_OK;
}

CoS_StatusTypeDef CoS_UART_ReceiveByte(CoS_UART_HandleTypeDef *hObj)
{
	__IO u8 rx;
	UART_HandleTypeDef* huart=hObj->huart;
	USART_TypeDef * USARTn=huart->Instance;
	uint32_t isrflags   = READ_REG(USARTn->SR);
	uint32_t cr1its     = READ_REG(USARTn->CR1);
	uint32_t cr3its     = READ_REG(USARTn->CR3);
	uint32_t errorflags = 0x00U;
	uint32_t dmarequest = 0x00U;
	/* If no error occurs */
	//__CoS_LOCK(hObj);
	rx = (uint16_t)(USARTn->DR & (uint16_t)0x01FF);
	errorflags = (isrflags & (uint32_t)(USART_SR_PE | USART_SR_FE | USART_SR_ORE | USART_SR_NE));
	if(errorflags == RESET)
	{
		/* UART in mode Receiver -------------------------------------------------*/
		if(((isrflags & USART_SR_RXNE) != RESET) && ((cr1its & USART_CR1_RXNEIE) != RESET))
		{
			CoS_UART_ByteDeal(hObj,rx);
//			__CoS_UNLOCK(hObj);
//			return CoS_STATE_OK;
		}
	}  
	//if((errorflags != RESET) && ((cr3its & (USART_CR3_EIE | USART_CR1_PEIE)) != RESET))
	else if(((cr3its & (USART_CR3_EIE | USART_CR1_PEIE)) != RESET))
	{
		/* UART parity error interrupt occurred ----------------------------------*/
		if(((isrflags & USART_SR_PE) != RESET) && ((cr1its & USART_CR1_PEIE) != RESET))
			huart->ErrorCode |= HAL_UART_ERROR_PE;
		/* UART noise error interrupt occurred -----------------------------------*/
		if(((isrflags & USART_SR_NE) != RESET) && ((cr3its & USART_CR3_EIE) != RESET))
			huart->ErrorCode |= HAL_UART_ERROR_NE;
		/* UART frame error interrupt occurred -----------------------------------*/
		if(((isrflags & USART_SR_FE) != RESET) && ((cr3its & USART_CR3_EIE) != RESET))
			huart->ErrorCode |= HAL_UART_ERROR_FE;
		/* UART Over-Run interrupt occurred --------------------------------------*/
		if(((isrflags & USART_SR_ORE) != RESET) && ((cr3its & USART_CR3_EIE) != RESET))
			huart->ErrorCode |= HAL_UART_ERROR_ORE;
		/* Call UART Error Call back function if need be --------------------------*/    
		if(huart->ErrorCode != HAL_UART_ERROR_NONE)
		{
			/* UART in mode Receiver -----------------------------------------------*/
			if(((isrflags & USART_SR_RXNE) != RESET) && ((cr1its & USART_CR1_RXNEIE) != RESET))
			{
				//rx = (uint16_t)(huart->Instance->DR & (uint16_t)0x01FF);
				CoS_UART_ByteDeal(hObj,rx);
			}
			/* If Overrun error occurs, or if any error occurs in DMA mode reception,
			consider error as blocking */
			dmarequest = HAL_IS_BIT_SET(huart->Instance->CR3, USART_CR3_DMAR);
			if(((huart->ErrorCode & HAL_UART_ERROR_ORE) != RESET) || dmarequest)
			{
				/* Blocking error : transfer is aborted
				Set the UART state ready to be able to start again the process,
				Disable Rx Interrupts, and disable Rx DMA request, if ongoing */
				CLEAR_BIT(huart->Instance->CR1, (USART_CR1_RXNEIE | USART_CR1_PEIE));
				CLEAR_BIT(huart->Instance->CR3, USART_CR3_EIE);
				/* At end of Rx process, restore huart->RxState to Ready */
				 huart->RxState= HAL_UART_STATE_READY;
			}
			else
			{
				/* Non Blocking error : transfer could go on. 
				Error is notified to user through user error callback */
				HAL_UART_ErrorCallback(huart);
				huart->ErrorCode = HAL_UART_ERROR_NONE;
			}
		}
		
	} /* End if some error occurs */
	//__CoS_UNLOCK(hObj);
		return CoS_STATE_OK;
}

static void CoS_UART_ByteDeal(CoS_UART_HandleTypeDef *hObj,u8 rx)
{
	//USART_TypeDef* hUARTn=hObj->huart->Instance;
	CoS_UART_TypeDef* hInstance=hObj->Instance;
	u16 rx_counter=hInstance->Rx_Count;
	u8 n;
	if (hInstance->CCByte & _CoS_UART_RxReadBack)
	{
		UARTn_putc(hObj->Instance->USARTn,rx);
	}
	if(!hInstance->RxOK)
	{
		if(rx_counter>=hInstance->Rx_BufSize-3)
		{
			if (hInstance->CCByte & _CoS_UART_RxICOv)
			{
				hInstance->RxOK=2;
				CoS_UART_RxRST(hObj);
				return;
			}
			rx_counter=hInstance->Rx_BufSize-10;
			hInstance->Rx_OverrunCount++;
			if (hInstance->Rx_OverrunCount>5)
			{
				CoS_UART_RxRST(hObj);
				return;
			}
		}
		hInstance->Rx_Base[rx_counter]=rx;// hSB_USART1.Rx_Buff[rx_counter++]=rx;
		if (hInstance->CC.mode==UART_MODE_UART)
		{
			if(rx==hInstance->escape2 && rx_counter>2)
			{
				if(hInstance->Rx_Base[rx_counter-1]==hInstance->escape1)
				{
					if(!(hInstance->CCByte & _CoS_UART_RxRAWF))
					{
						hInstance->Rx_Base[rx_counter-1]=0;
						hInstance->Rx_Base[rx_counter]=0;
						hInstance->Rx_Count=rx_counter-1;
					}
					hInstance->RxOK=1;
				}
			}
			hInstance->Rx_Count++;
		}
		else if (hInstance->CC.mode==UART_MODE_CPtP)
		{
			switch (hInstance->CC.userdata>>4)
			{
				case PHASE_IDLE:
				{
					if (rx==hInstance->escape2)//head == end is not allowed
						return;
					hInstance->Rx_Base[0]=rx;
					hInstance->escapeR=rx;
					hInstance->Rx_Count++;
					hInstance->CC.userdata = (hInstance->CC.userdata&0x0f) + PHASE_RHEAD;
					break;
				}
				case PHASE_HEAD:
				{
					hInstance->Rx_Base[rx_counter]=rx;
					hInstance->Rx_Count++;
					if (rx==hInstance->escape2)//end
					{
						hInstance->CC.userdata = (hInstance->CC.userdata&0x0f) + PHASE_IDLE;
						hInstance->RxOK=1;
						if(!(hInstance->CCByte & _CoS_UART_RxRAWF))
						{
							hInstance->Rx_Base[rx_counter]=0;
							hInstance->Rx_Count=rx_counter;
						}
						return;
					}
					else if (rx==hInstance->escapeR)
						hInstance->CC.userdata = (hInstance->CC.userdata&0x0f) + PHASE_RDATAE;
					else 
						hInstance->CC.userdata = (hInstance->CC.userdata&0x0f) + PHASE_RDATA;
					break;
				}
				case PHASE_DATAE:
				{
					hInstance->Rx_Base[rx_counter-1]=rx+0x20;
					hInstance->CC.userdata = (hInstance->CC.userdata&0x0f) + PHASE_RDATA;
					break;
				}
				case PHASE_DATA:
				{
					hInstance->Rx_Base[rx_counter]=rx;
					hInstance->Rx_Count++;
					if (rx==hInstance->escape2)//end
					{
						hInstance->CC.userdata = (hInstance->CC.userdata&0x0f) + PHASE_RIDLE;
						hInstance->RxOK=1;
						if(!(hInstance->CCByte & _CoS_UART_RxRAWF))
						{
							hInstance->Rx_Base[rx_counter]=0;
							hInstance->Rx_Count=rx_counter;
						}
						return;
					}
					else if (rx==hInstance->escapeR)
						hInstance->CC.userdata = (hInstance->CC.userdata&0x0f) + PHASE_RDATAE;
					else 
						hInstance->CC.userdata = (hInstance->CC.userdata&0x0f) + PHASE_RDATA;
					break;
				}
				
				case PHASE_END:
				default:
				{
					hInstance->CC.userdata = (hInstance->CC.userdata&0x0f) + PHASE_RIDLE;
					hInstance->RxOK=1;
					return;
				}
			}
		}
	}
}
//#define PHASE_HEAD 1U
//#define PHASE_END 2U
//#define PHASE_DATA 4U
//#define PHASE_IDLE 0U

uint8_t CoS_UART_RxState(CoS_UART_HandleTypeDef *hObj)
{
	return !!((hObj->Instance->CC.state)&CoS_STATE_BUSY_MASKR);
}

uint8_t CoS_UART_TxState(CoS_UART_HandleTypeDef *hObj)
{
	return !!((hObj->Instance->CC.state)&CoS_STATE_BUSY_MASKT);
}

CoS_StatusTypeDef CoS_UART_SetBaud(CoS_UART_HandleTypeDef *hObj,uint32_t newBaud)
{
	if(hObj == NULL)
	{
		return CoS_STATE_ERROR;
	}
	__CoS_LOCK(hObj);
	hObj->huart->Init.BaudRate=newBaud;
	if (HAL_UART_Init(hObj->huart) != HAL_OK)
		return CoS_STATE_ERROR;
	__CoS_UNLOCK(hObj);
	return CoS_STATE_OK;
}
//static void CoS_UART_TxEnd(CoS_UART_HandleTypeDef *hObj)
//{
//	
//}

//static void CoS_UART_TxStart(CoS_UART_HandleTypeDef *hObj)
//{
//	hObj->Instance->TxBusy=1;
//}

//CoS_StatusTypeDef CoS_UART_TxSend(CoS_UART_HandleTypeDef *hObj,uint8_t* chbase, uint16_t length)
//{
//	CoS_UART_TypeDef* hInstance=hObj->Instance;
//	if (!hInstance->TxBusy && !(hInstance->Tx_Length))
//	{
//		hInstance->Tx_Base=chbase;
//		hInstance->Tx_Length=length;
//		hInstance->Tx_Count=0;
//		CoS_UART_TxStart(hObj);
//		return CoS_STATE_OK;
//	}
//	else
//		return CoS_STATE_BUSY_T;
//}

//void CoS_UART_TxRun(CoS_UART_HandleTypeDef *hObj)
//{
//	u16 n=0;
//	CoS_UART_TypeDef* hInstance=hObj->Instance;
//	if (hInstance->TxBusy)
//	{
//		if (!hInstance->Tx_Length)
//		{
//			CoS_UART_TxRST(hObj);
//			return;
//		}
//		while(n++<(hInstance->Tx_PacketSize) && hInstance->Tx_Count < hInstance->Tx_Length)
//			UARTn_putc(hObj->huart->Instance,hInstance->Tx_Base[hInstance->Tx_Count++]);
//		if (n==hInstance->Tx_PacketSize)
//			hInstance->Tx_Count--;
//		if ( hInstance->Tx_Count >= hInstance->Tx_Length)
//			CoS_UART_TxRST(hObj);
//	}
//}