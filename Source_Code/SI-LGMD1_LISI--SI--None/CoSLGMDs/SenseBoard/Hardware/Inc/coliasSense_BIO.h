#ifndef __CoSBIO_H
#define __CoSBIO_H

#include "stm32f4xx_hal.h"
#include "coliasSense_def.h"

#define Singleton_BIO


#define TICin2 (TIM2->CNT=0)
#define TOCin2 (TIM2->CNT)


#define Heartbeat_us 50000U //used for background task (TIM3)

#define LEDn 	3
#define LED_ON GPIO_PIN_RESET
#define LED_OFF GPIO_PIN_SET


#define LED1_ON (CoS_BIO_write(&(hBIO1->LEDs[LED1]),GPIO_PIN_RESET))
#define LED2_ON (CoS_BIO_write(&(hBIO1->LEDs[LED2]),GPIO_PIN_RESET))
#define LED3_ON (CoS_BIO_write(&(hBIO1->LEDs[LED3]),GPIO_PIN_RESET))

#define LED1_OFF (CoS_BIO_write(&(hBIO1->LEDs[LED1]),GPIO_PIN_SET))
#define LED2_OFF (CoS_BIO_write(&(hBIO1->LEDs[LED2]),GPIO_PIN_SET))
#define LED3_OFF (CoS_BIO_write(&(hBIO1->LEDs[LED3]),GPIO_PIN_SET))

#define LED1_Toggle (CoS_BIO_Toggle(&(hBIO1->LEDs[LED1])))
#define LED2_Toggle (CoS_BIO_Toggle(&(hBIO1->LEDs[LED2])))
#define LED3_Toggle (CoS_BIO_Toggle(&(hBIO1->LEDs[LED3])))

#define Key1State (!CoS_BIO_read(&(hBIO1->keys[Key1])))
#define Key2State (!CoS_BIO_read(&(hBIO1->keys[Key2])))
typedef enum
{
	 LED1=0U,
	 LED2=1U,
	 LED3=2U,
}CoS_BIO_LEDType;

typedef enum
{
	 Key1=0U,
	 Key2=1U,
}CoS_BIO_KeyType;

/**********************************/
typedef enum
{
	BIO_Key_Ready=0U,
	BIO_Key_Pressing=0x01U,
	BIO_Key_Pressing_Long=0x02U,
	BIO_Key_Released=0x10,
	BIO_Key_Released_Long=0x20U,
}CoS_BIO_KeyStateType;
/**********************************/

typedef struct{
	GPIO_TypeDef* Port;
	uint16_t Pin;
}BIO_IO_TypdeDef;


typedef struct 
{
	CoS_CCSTypeDef CC;		//common controller
	BIO_IO_TypdeDef keys[2];
	BIO_IO_TypdeDef LEDs[LEDn];
	uint32_t timerlog[10];
	uint64_t TIM3TickCount;
	CoS_BIO_KeyStateType keyState[2];
	uint32_t presstime[2];
	uint32_t presstime_threshold;
	IRQn_Type BackgroundIRQn;
}CoS_BIO_TypeDef;

typedef struct{
	CoS_BIO_TypeDef* Instance;
	TIM_HandleTypeDef* htim2;
	TIM_HandleTypeDef* htim3;
	TIM_HandleTypeDef* htim5;
	IRQn_Type BackgroundIRQn;
}CoS_BIO_InitTypeDef;



typedef struct strCoS_BIO_HandleTypeDef
{
	CoS_BIO_TypeDef* Instance; //instance address
	CoS_LockTypeDef Lock; 	//Locking object
	__IO CoS_StatusTypeDef iState; //init state
	TIM_HandleTypeDef* htim2;
	TIM_HandleTypeDef* htim3;
	TIM_HandleTypeDef* htim5;
}CoS_BIO_HandleTypeDef;










//extern SenseBoard_HandleTypeDef hSenseBoard;
extern uint64_t* pTIM3Tick;
extern CoS_BIO_TypeDef* hBIO1;


 __weak void Error_Handler(void);
void BackgroundTask_Enable(void) ;
void BackgroundTask_Disable(void) ;
CoS_StatusTypeDef CoS_BIO_Init(CoS_BIO_HandleTypeDef *hObj,CoS_BIO_InitTypeDef* initer);
CoS_StatusTypeDef CoS_BIO_DeInit(CoS_BIO_HandleTypeDef *hObj);
void CoS_BIO_write(BIO_IO_TypdeDef* hIO,GPIO_PinState sta);
void CoS_BIO_Toggle(BIO_IO_TypdeDef* hIO);
GPIO_PinState CoS_BIO_read(BIO_IO_TypdeDef* hIO);
CoS_StatusTypeDef CoS_BIO_enable(CoS_BIO_HandleTypeDef* hObj);
CoS_StatusTypeDef CoS_BIO_disable(CoS_BIO_HandleTypeDef* hObj);
CoS_StatusTypeDef CoS_BIO_control(CoS_BIO_HandleTypeDef* hObj);
u8 CoS_Pin_StateReport(BIO_IO_TypdeDef* pin,u8 trans_mode_blocking);
//__weak void CoS_TIM3_IRQ_Task(CoS_BIO_HandleTypeDef* hObj);
uint8_t CoS_Restart(u8 keynum);
uint32_t CoS_MCUID(char i);

#endif