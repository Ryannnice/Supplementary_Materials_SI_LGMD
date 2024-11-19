#ifndef __CoSREG_H
#define __CoSREG_H

#include "stm32f4xx_hal.h"
#include "coliasSense_def.h"
#include "coliasSense_board.h"


typedef struct
{
	CoS_CCSTypeDef CC;		//common controller
	u32 offset;
	
}CoS_REGS_TypeDef;

#define regCAM_START  0x0000
#define regCAM_END    0x00ff
#define regCAM_offset 0x60000000

#define regMPU_START  0x0100
#define regMPU_END    0x01ff
#define regMPU_offset 0xA1000000

#define regMAG_START  0x0200
#define regMAG_END    0x021F
#define regMAG_offset 0x0200

#define regCOB_START  0x0300
#define regCOB_END    0x03ff
#define regCOB_offset 0x0300

#define regCOB_START  0x0300
#define regCOB_END    0x03ff
#define regCOB_offset 0x0300

CoS_StatusTypeDef CoS_REG_init(u16 addr,u8* pData,u16 length);
CoS_StatusTypeDef CoS_REG_read(u16 addr,u8* pData,u16 length);
CoS_StatusTypeDef CoS_REG_write_nack(u16 addr,u8* pData,u16 length);
CoS_StatusTypeDef CoS_REG_write_ack(void);
#endif