#include "delay.h"

__IO u16 ntime;

static __IO uint32_t TimingDelay;

void Delay_pre(__IO uint32_t nTime)
{ 
	  TimingDelay = nTime;
	  while(TimingDelay != 0);
}

/**
  * @brief  Decrements the TimingDelay variable.
  * @param  None
  * @retval None
  */
void TimingDelay_Decrement(void)
{
	if (TimingDelay != 0x00)
		TimingDelay--;
}

void delay_ms(u16 nms)
{	 		  	   
//	ntime=nms;
	SysTick_Config(SystemCoreClock / 1000);//1ms
	HAL_NVIC_SetPriority(SysTick_IRQn, 0 ,1);
	TimingDelay = nms;
	while(TimingDelay != 0);
//	while(ntime);
//	SysTick->CTRL=0x00;			  	    
}   
		    								   
void delay_us(u32 nus)
{		
	ntime=nus;
	SysTick_Config(SystemCoreClock / 1000000);//1us
	HAL_NVIC_SetPriority(SysTick_IRQn, 0 ,1);
	TimingDelay = nus;
	while(TimingDelay != 0);
}

void delay_short(u32 count)
{
	while(count--);
}

void delay_long(u32 count)
{
	u16 n;
	u32 thecount;
	for(n=0;n<1000;n++)
	{
		thecount = count;
		while(thecount--);
	}
}
