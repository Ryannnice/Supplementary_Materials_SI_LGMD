#ifndef __SENSEBDCOM_H
#define __SENSEBDCOM_H

#include "stm32f4xx_hal.h"
#include "coliasSense_def.h"
#include "CCMRAM_DEF.h"
#include "coliasSense_interfaces.h"
#include "coliasSense_BIO.h"
#include "coliasSense_CLED.h"
#include "coliasSense_UART.h"
#include "coliasSense_Camera.h"
#include "coliasSense_Exbd.h"
#include "coliasSense_Motion.h"
#include "coliasSense_LGMD.h"
#include "coliasSense_MPU9250.h"
#define UART1_Buff_size 50
#define UART4_Buff_size 50
#define INIT_INFO_printf true
#define USE_printf_USART1


//#define A_OV7670W 42//0x42

typedef struct{
	u8 Enabled_printf;
	TIM_HandleTypeDef*       htim1;
	TIM_HandleTypeDef*       htim2;
	TIM_HandleTypeDef*       htim3;
	TIM_HandleTypeDef*       htim5;
	RNG_HandleTypeDef*       hrng;
	UART_HandleTypeDef*      huart1;
	UART_HandleTypeDef*      huart4;
	DMA_HandleTypeDef*       hdma_dcmi;
	DCMI_HandleTypeDef*      hdcmi;
	CoS_BIO_HandleTypeDef*   hHBIO;
	CoS_CLED_HandleTypeDef*  hHCLED;
	CoS_UART_HandleTypeDef*  hHUART1;
	CoS_UART_HandleTypeDef*  hHUART4;
	CoS_sI2C_HandleTypeDef*  hHsI2C_Camera;
	CoS_sI2C_HandleTypeDef*  hHsI2C_MPU;
	CoS_sI2C_HandleTypeDef*  hHsI2C_MAG;
	CoS_CAM_HandleTypeDef*   hHCoS_Camera;
	CoS_Motion_HandleTypeDef* hHCos_Motion;
	CoS_EXBD_HandleTypeDef*  hHCoS_EXBD;
//CoS_LGMD1_HandleTypeDef* hHLGMD1;
}SenseBoard_HandleTypeDef;

extern SenseBoard_HandleTypeDef* const  hCoS;

extern CoS_BIO_TypeDef BIO1;
extern CoS_CLED_TypeDef CLED1;
extern CoS_UART_TypeDef CoS_UART1;
extern CoS_UART_TypeDef CoS_UART4;
extern CoS_sI2C_TypeDef CoS_sI2C_Camera;
extern CoS_sI2C_TypeDef CoS_sI2C_MPU;
//extern CoS_LGMD1_TypeDef CoS_LGMD1;

CoS_StatusTypeDef CoS_SenseBoard_Init(SenseBoard_HandleTypeDef *hObj);
void CoS_TIM3_IRQ_Task(CoS_BIO_HandleTypeDef* hObj);
#endif

