#ifndef __DELAY_H
#define __DELAY_H

#include "stm32f4xx_hal.h"
#include "coliasSense_def.h"
//#include "main.h"
void Delay_pre(__IO uint32_t nTime);
void TimingDelay_Decrement(void);
void delay_ms(uint16_t nms);
void delay_us(uint32_t nus);
void delay_short(uint32_t count);
void delay_long(uint32_t count);
#endif

