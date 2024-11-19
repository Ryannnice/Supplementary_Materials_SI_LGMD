#ifndef __SENSEBDDEF_H
#define __SENSEBDDEF_H

#include "stm32f4xx.h"

#define IDLE 0
#define BUSY 1
#ifndef true
	#define true 1
#endif
#ifndef false
	#define false 0
#endif
#define ARM_MATH_CM4
#define TMS_Pin GPIO_PIN_13
#define TMS_GPIO_Port GPIOA

#define TCK_Pin GPIO_PIN_14
#define TCK_GPIO_Port GPIOA

#define LED1_Pin GPIO_PIN_2
#define LED1_GPIO_Port GPIOA

#define LED2_Pin GPIO_PIN_12
#define LED2_GPIO_Port GPIOC

#define LED3_Pin GPIO_PIN_13
#define LED3_GPIO_Port GPIOD

#define LED4_Pin GPIO_PIN_9
#define LED4_GPIO_Port GPIOE

#define LED5_Pin GPIO_PIN_11
#define LED5_GPIO_Port GPIOE

#define LED6_Pin GPIO_PIN_13
#define LED6_GPIO_Port GPIOE

#define S1_Pin GPIO_PIN_0
#define S1_GPIO_Port GPIOA

#define S2_Pin GPIO_PIN_1
#define S2_GPIO_Port GPIOA

#define GP1_Pin GPIO_PIN_15
#define GP1_GPIO_Port GPIOA

#define GP2_Pin GPIO_PIN_4
#define GP2_GPIO_Port GPIOC

#define FLASH_CS_Pin GPIO_PIN_13
#define FLASH_CS_GPIO_Port GPIOC

#define FLASH_D2_Pin GPIO_PIN_14
#define FLASH_D2_GPIO_Port GPIOC

#define FLASH_D3_Pin GPIO_PIN_15
#define FLASH_D3_GPIO_Port GPIOC

#define eSPI_P2_Pin GPIO_PIN_0
#define eSPI_P2_GPIO_Port GPIOC

#define eSPI_P3_Pin GPIO_PIN_1
#define eSPI_P3_GPIO_Port GPIOC

#define eSPI_MOSI_Pin GPIO_PIN_3
#define eSPI_MOSI_GPIO_Port GPIOC

#define sI2C_SCL_Pin GPIO_PIN_8
#define sI2C_SCL_GPIO_Port GPIOB

#define sI2C_SDA_Pin GPIO_PIN_11
#define sI2C_SDA_GPIO_Port GPIOB

#define ACC_INT_Pin GPIO_PIN_2
#define ACC_INT_GPIO_Port GPIOC

#define ACC_AD0_Pin GPIO_PIN_15
#define ACC_AD0_GPIO_Port GPIOE

#define ACC_nCS_Pin GPIO_PIN_8
#define ACC_nCS_GPIO_Port GPIOC

#define eSPI_P1_Pin GPIO_PIN_12
#define eSPI_P1_GPIO_Port GPIOB

#define eSPI_SCK_Pin GPIO_PIN_13
#define eSPI_SCK_GPIO_Port GPIOB

#define eSPI_MISO_Pin GPIO_PIN_14
#define eSPI_MISO_GPIO_Port GPIOB

#define MCO_CAM_Pin GPIO_PIN_8
#define MCO_CAM_GPIO_Port GPIOA

#define SCCB_SDA_Pin GPIO_PIN_2
#define SCCB_SDA_GPIO_Port GPIOE

#define SCCB_SCL_Pin GPIO_PIN_3
#define SCCB_SCL_GPIO_Port GPIOE

#define EXBD_TEST_Pin GPIO_PIN_12
#define EXBD_TEST_GPIO_Port GPIOD

typedef enum
{
	CoS_STATE_UNSET     = 0x00U,
	CoS_STATE_ERROR     = 0x01U,
	CoS_STATE_BUSY      = 0x02U,
	CoS_STATE_TIMEOUT   = 0x03U,
	CoS_STATE_OK        = 0x04U,
	CoS_STATE_BUSY_T    = 0x12U,
	CoS_STATE_BUSY_R    = 0x22U,
	CoS_STATE_BUSY_TR   = 0x32U,
	CoS_STATE_BUSY_MASK = 0x30U,
	CoS_STATE_BUSY_MASKT= 0x10U,
	CoS_STATE_BUSY_MASKR= 0x20U,
	//CoS_STATE_DATA_OK   = 0x40U,
	CoS_STATE_SUSPENDED = 0x05U
}CoS_StatusTypeDef;

typedef enum
{
	CoS_UNLOCKED = 0x00U,
	CoS_LOCKED   = 0x01U 
}CoS_LockTypeDef;

typedef struct
{
	uint8_t enable;  //flag:enable
	uint8_t userdata; //customized data
	uint8_t mode; //a flag determine the mode
	CoS_StatusTypeDef state; 
}CoS_CCSTypeDef;	

typedef struct
{
	CoS_StatusTypeDef id;
	uint8_t* str;
}CoS_Str_IDType;


#define CoS_IS_BIT_SET(REG, BIT)           (((REG) & (BIT)) != RESET)
#define CoS_IS_BIT_CLR(REG, BIT)           (((REG) & (BIT)) == RESET)
#define CoS_RESET_HANDLE_STATE(__HANDLE__) ((__HANDLE__)->State = 0U)
#define __CoS_HANDLE_ENABLE(__HANDLE__)    ((__HANDLE__)->CC->enable =  1)
#define __CoS_HANDLE_DISABLE(__HANDLE__)   ((__HANDLE__)->CC->enable =  0)
#define __CoS_LOCK(__HANDLE__)    do{                                      \
                                    if((__HANDLE__)->Lock == CoS_LOCKED)   \
                                    {                                      \
                                       return CoS_STATE_BUSY;              \
                                    }                                      \
                                    else                                   \
                                    {                                      \
                                       (__HANDLE__)->Lock = CoS_LOCKED;    \
                                    }                                      \
                                  }while (0)


#define __CoS_UNLOCK(__HANDLE__)  do{                                       \
                                      (__HANDLE__)->Lock = CoS_UNLOCKED;    \
                                    }while (0)          
#define OFFSET_OF(type, member) (uint32_t)(&(((type *)0)->member))
#define STRUCT_ENTRY(ptr, type, member) (type *)((uint32_t)(ptr)-OFFSET_OF(type, member))

	 /*!< STM32F40x Standard Peripheral Library old types (maintained for legacy purpose) */
typedef int32_t  s32;
typedef int16_t s16;
typedef int8_t  s8;

typedef const int32_t sc32;  /*!< Read Only */
typedef const int16_t sc16;  /*!< Read Only */
typedef const int8_t sc8;   /*!< Read Only */

typedef __IO int32_t  vs32;
typedef __IO int16_t  vs16;
typedef __IO int8_t   vs8;

typedef __I int32_t vsc32;  /*!< Read Only */
typedef __I int16_t vsc16;  /*!< Read Only */
typedef __I int8_t vsc8;   /*!< Read Only */

typedef uint32_t  u32;
typedef uint16_t u16;
typedef uint8_t  u8;

typedef const uint32_t uc32;  /*!< Read Only */
typedef const uint16_t uc16;  /*!< Read Only */
typedef const uint8_t uc8;   /*!< Read Only */

typedef __IO uint32_t  vu32;
typedef __IO uint16_t vu16;
typedef __IO uint8_t  vu8;

typedef __I uint32_t vuc32;  /*!< Read Only */
typedef __I uint16_t vuc16;  /*!< Read Only */
typedef __I uint8_t vuc8;   /*!< Read Only */

#pragma diag_suppress=1
#pragma diag_suppress=68 //#68-D: integer conversion resulted in a change of sign
#pragma diag_suppress=144// #144-D: a value of type "volatile u8 *" cannot be used to initialize an entity of type "void *"
#pragma diag_suppress=167//#167-D: argument of type "volatile u8 *" is incompatible with parameter of type "u8 *"
#pragma diag_suppress=177
#pragma diag_suppress=188// #188-D: enumerated type mixed with another type
#pragma diag_suppress=513 // #513-D: a value of type "u8 *" cannot be assigned to an entity of type "char *"
#pragma diag_suppress=550 // #550-D: variable "t" was set but never used
#pragma diag_suppress=1035//#1035-D: single-precision operand implicitly converted to double-precision

#define USE_printf


#endif

