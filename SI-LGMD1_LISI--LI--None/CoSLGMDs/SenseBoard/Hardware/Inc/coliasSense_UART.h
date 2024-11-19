#ifndef __CoSUART_H
#define __CoSUART_H

#include "stm32f4xx_hal.h"
#include "coliasSense_def.h"


//#ifndef CUSTOM_UART_ESCAPE
//	#define USART1_Escape1 0x0DU
//	#define USART1_Escape2 0x0AU
//	#define UART4_Escape1 0x0DU
//	#define UART4_Escape2 0x0AU
//	#define UART4_Escape1Send 0x55U
//	#define UART4_Escape2Send 0xAAU
//#endif

#define UART_MODE_UART 0
#define UART_MODE_CPtP 1
#define UART_MODE_PPP 2 //to be added later

#define _CoS_UART_RxReadBack 0x01
#define _CoS_UART_RxICOv 0x02 //clear overrun immediatly
#define _CoS_UART_RxRAWF 0x04 //store raw frame on receiving
typedef struct{
	CoS_CCSTypeDef CC;
	USART_TypeDef * USARTn; 
	uint8_t* Rx_Base;
	uint8_t* Tx_Base;
	
	uint16_t Rx_BufSize;
	uint16_t Rx_Count;
	uint16_t Rx_OverrunCount;

	uint16_t Tx_Length;
	uint16_t Tx_Count;
	uint16_t Tx_PacketSize;
	
	uint8_t RxOK;
	//uint8_t RxBusy;
	uint8_t CCByte;
	//uint8_t TxBusy;
	uint8_t escape1;
	uint8_t escape2;
	uint8_t escapeR;
}CoS_UART_TypeDef;


typedef struct 
{
	CoS_UART_TypeDef * Instance; 
	CoS_LockTypeDef Lock; 
	__IO CoS_StatusTypeDef iState;
	UART_HandleTypeDef* huart;
}CoS_UART_HandleTypeDef;

typedef struct{
	CoS_UART_TypeDef* Instance;
	UART_HandleTypeDef* huart;
	uint8_t* RxBufferBase;
	uint16_t RxBufferSize;
	uint16_t Tx_PacketSize;
	uint8_t escape1;
	uint8_t escape2;
	uint8_t mode;
	uint8_t Readback;
	uint8_t RxICOv;
	uint8_t RxRAWF;
}CoS_UART_InitTypeDef;


CoS_StatusTypeDef CoS_UART_Init(CoS_UART_HandleTypeDef *hObj,
								CoS_UART_InitTypeDef* Init);
void Flush_UART_BUFF(CoS_UART_TypeDef* hInstance);
CoS_StatusTypeDef CoS_UART_RxRST(CoS_UART_HandleTypeDef *hObj);
CoS_StatusTypeDef CoS_UART_ReceiveByte(CoS_UART_HandleTypeDef *hObj);
CoS_StatusTypeDef CoS_UART_TxRST(CoS_UART_HandleTypeDef *hObj);
//CoS_StatusTypeDef CoS_UART_TxSend(CoS_UART_HandleTypeDef *hObj,uint8_t* chbase, uint16_t length);
void CoS_UART_TxRun(CoS_UART_HandleTypeDef *hObj);
void UART_transmit_run(CoS_UART_HandleTypeDef *hObj);
CoS_StatusTypeDef CoS_UART_sendFrame(CoS_UART_HandleTypeDef *hObj,uint8_t* Tx_Base,uint16_t Tx_Length);
CoS_StatusTypeDef CoS_UART_Enable(CoS_UART_HandleTypeDef *hObj,uint8_t sta);

uint8_t CoS_UART_TxState(CoS_UART_HandleTypeDef *hObj);
void CoS_UART_StateReport(CoS_UART_HandleTypeDef *hObj,u8 Uid,u8 trans_mode_blocking);
CoS_StatusTypeDef CoS_UART_SetBaud(CoS_UART_HandleTypeDef *hObj,uint32_t newBaud);
#endif
