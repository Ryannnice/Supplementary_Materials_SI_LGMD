/**
  ******************************************************************************
  * File Name          : main.c
  * Description        : Main program body
  ******************************************************************************
  *
  * COPYRIGHT(c) 2016 STMicroelectronics
  *
  * Redistribution and use in source and binary forms, with or without modification,
  * are permitted provided that the following conditions are met:
  *   1. Redistributions of source code must retain the above copyright notice,
  *      this list of conditions and the following disclaimer.
  *   2. Redistributions in binary form must reproduce the above copyright notice,
  *      this list of conditions and the following disclaimer in the documentation
  *      and/or other materials provided with the distribution.
  *   3. Neither the name of STMicroelectronics nor the names of its contributors
  *      may be used to endorse or promote products derived from this software
  *      without specific prior written permission.
  *
  * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
  * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
  * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
  * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
  * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
  * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
  * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
  * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
  * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
  *
  ******************************************************************************
  */
/* Includes ------------------------------------------------------------------*/

#include "main.h"
//#include "LCD18TFT.h"
#include "user_terminal.h"
#include "coliasSense_board.h"

#include "delay.h"
#include <stdio.h>
#include <string.h>
/* USER CODE BEGIN Includes */

/* USER CODE END Includes */

/* Private variables ---------------------------------------------------------*/

/* USER CODE BEGIN PV */
/* Private variables ---------------------------------------------------------*/

/* USER CODE END PV */


/* USER CODE BEGIN PFP */
/* Private function prototypes -----------------------------------------------*/

/* USER CODE END PFP */

/* USER CODE BEGIN 0 */
u8 id[]=      "123456798012345678901234657980123456798012345678\r\n";
volatile u8 testdt[4]={0x3a,0x35,0x22,0xf2};
uint8_t Key2_psta=0;
u8 testStr[]={"                                                                    Hello Wrold!\r\n"};
volatile uint8_t enable_motion=0;
//uint32_t gCounter=0;
uint8_t Key2_virtual;


#pragma pack(1)
/*
u8 testRead(void)
{
	u8 dataStr[4];
	u8 sta;
	sta= ReadBytesWithACKFromColiasBasic(hCoS->hHUART4, 0xA3, dataStr, 1);
}
*/
int main(void)
{
	volatile u8 a=0,b,c,d=1,e;
	u8 f[10]={0},n;
	u8 tstr[]={0x11,0xa1};
	ColiasBasicID_TypeDef idn;
	CoS_SenseBoard_Init(hCoS);
	CoS_UART_Enable(hCoS->hHUART1,ENABLE);
	UT_Init();
	CoS_CAM_DMASet(hCoS->hHCoS_Camera,(uint32_t)Image,10692);
	CoS_CAM_Start(hCoS->hHCoS_Camera);
	BackgroundTask_Enable();
	//TFT_LCD_Init();
	//delay_ms(1800);
	//LCD_Clear(BLACK);
	//TiltSense.CC.enable=0;
	Key2_virtual=0;
	//tilt detecting
	//tilt_detection_init(&TiltSense);
	//LED1_OFF;
	//LED2_OFF;
	//LED3_OFF;
	//hCoS->hHCLED->Instance->CC.mode=CLED_Freeze;
	//CoS_CLED_setRGB(0,0,0);
	//infinite loop
	while (1)
	{
		
		/*
		For Qinbing Attention
		chage the function in while(1) for your own purpose.
		*/
		
		//neuron model processing
		LGMD_demo(&hLGMD);
		//motion
  	Decision_making(&hLGMD,enable_motion);
		//LCD display
		//LCD_imgshow(0,0,Owidth,Oheight,&Image[0][0][0],disp_mode);
		//check buttom board IT
		//ColiasBasicITEvents_Check((hCoS->hHCos_Motion->Instance),&(hCoS->hHCos_Motion->Instance->IT_Events));
		//check TCRTs
		/*
		if (enable_motion)
		{
			//uint16_t ADC_Value[3] = {0};
			
			//#define TCRT_Th 700
			
			//left bumper IT
			//if (ADC_Value[2] == 0x0000 || ADC_Value[2] == 0x0200)
			//	ADC_Value[0] = 1000;
				
			//ADC_Value[0] = ColiasBasicTCRTADC_Get(hCoS->hHCos_Motion->Instance, 1);
			//ADC_Value[1] = ColiasBasicTCRTADC_Get(hCoS->hHCos_Motion->Instance, 2);
			//ADC_Value[2] = ColiasBasicTCRTADC_Get(hCoS->hHCos_Motion->Instance, 3);
			
			//if (ADC_Value[0] == 0x0000 || ADC_Value[0] == 0x0200)
			//	ADC_Value[0] = 1000;
			
			//if (ADC_Value[1] == 0x0000 || ADC_Value[1] == 0x0200)
			//	ADC_Value[1] = 1000;
			
			//if (ADC_Value[2] == 0x0000 || ADC_Value[2] == 0x0200)
			//	ADC_Value[2] = 1000;
			
			//if( ADC_Value[0] < TCRT_Th)
			if(hCoS->hHCos_Motion->Instance->IT_Events.SensorIT & TCRT_Left_IT)
			{
				fast_motion_control(hCoS->hHCos_Motion->Instance,'R',0);
			}
			//middle bumper IT
			if(hCoS->hHCos_Motion->Instance->IT_Events.SensorIT & TCRT_Middle_IT)
			{
				fast_motion_control(hCoS->hHCos_Motion->Instance,'L',0);
			}
			//right bumper IT
			//if(ADC_Value[2] < TCRT_Th)
			if(hCoS->hHCos_Motion->Instance->IT_Events.SensorIT & TCRT_Right_IT)
			{
				fast_motion_control(hCoS->hHCos_Motion->Instance,'L',0);
			}
		}
		printf("%.2f %.2f %.2f\r\n",LGMD.Results.LGMD_out[hLGMD.currentDiffImage],LGMD.Results.FFI_out[hLGMD.currentDiffImage],LGMD.Results.spiFrequency);
		*/
		//button functions
		if (Key1State||Key2State) //button 1 pushed to move
		{
			enable_motion=1;
			//change default motion to Front
			uint8_t changeDefaultMotion[2]={0x90,0xa0};
			WriteBytesWithACKToColiasBasic(hCoS->hHCos_Motion->Instance->hHUART,changeDefaultMotion,2);
		}
		/*
		if(Key2State)	//button 2 pushed to retrieve neuron model outputs
		{
			printf("%.2f %.2f %.2f\r\n",LGMD.Results.LGMD_out[hLGMD.currentDiffImage],LGMD.Results.FFI_out[hLGMD.currentDiffImage],LGMD.Results.spiFrequency);
		}
		if(Key2_virtual)
		{
			printf("%.2f %.2f %.2f\r\n",LGMD.Results.LGMD_out[hLGMD.currentDiffImage],LGMD.Results.FFI_out[hLGMD.currentDiffImage],LGMD.Results.spiFrequency);
		}
		*/
		//tilt on
		//tilt_detection_run();
	}
	
}


/**
  * @brief  This function is executed in case of error occurrence.
  * @param  None
  * @retval None
  */
void Error_Handler2(void)
{
  /* USER CODE BEGIN Error_Handler */
  /* User can add his own implementation to report the HAL error return state */
  while(1) 
  {
  }
  /* USER CODE END Error_Handler */ 
}

#ifdef USE_FULL_ASSERT

/**
   * @brief Reports the name of the source file and the source line number
   * where the assert_param error has occurred.
   * @param file: pointer to the source file name
   * @param line: assert_param error line source number
   * @retval None
   */
void assert_failed(uint8_t* file, uint32_t line)
{
  /* USER CODE BEGIN 6 */
  /* User can add his own implementation to report the file name and line number,
    ex: printf("Wrong parameters value: file %s on line %d\r\n", file, line) */
  /* USER CODE END 6 */

}

#endif

/**
  * @}
  */ 

/**
  * @}
*/ 

/************************ (C) COPYRIGHT STMicroelectronics *****END OF FILE****/