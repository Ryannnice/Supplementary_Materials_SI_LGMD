#ifndef __CoSIT_H
#define __CoSIT_H

#include "stm32f4xx_hal.h"
#include "coliasSense_def.h"
#include "coliasSense_BIO.h"

#define sI2CMode_I2C 0
#define sI2CMode_SCCB 1

typedef struct 
{
	CoS_CCSTypeDef CC;		//common controller
	BIO_IO_TypdeDef PinSCL;
	BIO_IO_TypdeDef PinSDA;
	uint8_t* Rx_Base; //read regs
	uint8_t* Tx_Base; //write regs
	uint16_t nDelay;
	uint8_t  Rx_Length;
	uint8_t  Rx_Count;
	uint8_t  Tx_Length;
	uint8_t  Tx_Count;
	uint8_t  Dev_Addr; // private, better not modified by user
	uint8_t  DR;
}CoS_sI2C_TypeDef;

typedef struct{
	CoS_sI2C_TypeDef* Instance;
	BIO_IO_TypdeDef PinSCL;
	BIO_IO_TypdeDef PinSDA;
	uint16_t nDelay;
	uint8_t  Dev_Addr;
	uint8_t mode;
}CoS_sI2C_InitTypeDef;

typedef struct 
{
	CoS_sI2C_TypeDef* Instance; //instance address
	CoS_LockTypeDef Lock; 	//Locking object
	__IO CoS_StatusTypeDef iState; //init state
	uint8_t Dev_Addr; //[public, for modification by user
}CoS_sI2C_HandleTypeDef;

CoS_StatusTypeDef CoS_sI2C_Init(CoS_sI2C_HandleTypeDef *hObj,CoS_sI2C_InitTypeDef* Init);
CoS_StatusTypeDef CoS_sI2C_WriteReg(CoS_sI2C_HandleTypeDef *hObj,uint8_t regID, uint8_t regData);
CoS_StatusTypeDef CoS_sI2C_ReadReg(CoS_sI2C_HandleTypeDef *hObj,uint8_t regID,uint8_t* pRx);

CoS_StatusTypeDef CoS_sI2C_TxRST(CoS_sI2C_HandleTypeDef *hObj);
CoS_StatusTypeDef CoS_sI2C_WriteRegs(CoS_sI2C_HandleTypeDef *hObj,uint8_t regID,uint8_t* Tx_Base,uint8_t  Tx_Length);
CoS_StatusTypeDef CoS_sI2C_ReadRegs(CoS_sI2C_HandleTypeDef *hObj,uint8_t regID,uint8_t* Rx_Base,uint8_t  Rx_Length);
void CoS_sI2C_StateReport(CoS_sI2C_HandleTypeDef *hObj,u8 trans_mode_blocking);
#endif
