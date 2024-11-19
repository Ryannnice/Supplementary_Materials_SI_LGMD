#include "coliasSense_Exbd.h"
#include "coliasSense_board.h"
#include "delay.h"
static CoS_EXBD_NUMTypeDef Exboard_idStr[5]={0,1,2,4,8};
static u8* const Exboard_nameStr[5]={"none","USB","Bluetooth","FPV","SWARM"};
static u8* const Exboard_descriptionStr[5]={
	"no extend board has been detected",
	"extend board with USB-serial port and SD card slot",
	"extend board with Bluetooth-serial port and SD card slot",
	"extend board with an extra camera, LED and wireless transmitter",
	"extend board with basic SWARM supportting hardware"};


u8 CoS_EXBD_checkPin(CoS_EXBD_HandleTypeDef *hObj,u8 npin,u8 data)
{
	GPIO_InitTypeDef GPIO_InitStruct;
	u8 n,sta=0;
	//GPIO_InitTypeDef GPIO_InitStruct2;
	if (npin>7)
		return 0;
	GPIO_InitStruct.Pin = EXBD_TEST_Pin;
	GPIO_InitStruct.Mode = GPIO_MODE_INPUT;
	GPIO_InitStruct.Pull = GPIO_NOPULL;
	HAL_GPIO_Init(EXBD_TEST_GPIO_Port, &GPIO_InitStruct);
	delay_ms(2);
	GPIO_InitStruct.Pin = 1<<npin;
	GPIO_InitStruct.Mode = GPIO_MODE_OUTPUT_PP;
	GPIO_InitStruct.Pull = GPIO_PULLUP;
	GPIO_InitStruct.Speed = GPIO_SPEED_FREQ_MEDIUM;
	HAL_GPIO_Init(GPIOD, &GPIO_InitStruct);
	delay_ms(2);
	
	for (n=0;n<8;n++)
	{
		HAL_GPIO_WritePin(GPIOD, 1<<npin, (data>>n)&0x01);
		delay_ms(2);
		sta=sta+(HAL_GPIO_ReadPin(EXBD_TEST_GPIO_Port, EXBD_TEST_Pin)<<n);
	}
	
	GPIO_InitStruct.Mode = GPIO_MODE_INPUT;
	GPIO_InitStruct.Pull = GPIO_NOPULL;
	HAL_GPIO_Init(GPIOD, &GPIO_InitStruct);
	delay_ms(2);
	if (sta==data)
		return 1;
	else
		return 0;
}

u8 CoS_EXBD_check(CoS_EXBD_HandleTypeDef *hObj)
{
	u8 n,tmp=0;
	for (n=0;n<8;n++)
		tmp=tmp+(CoS_EXBD_checkPin(hObj,n,0x6d)<<n);
	for (n=0;n<5;n++)
	{
		if (tmp==Exboard_idStr[n])
			break;
	}
	if (n==5)
		return 0;
	else
		return n;
	return 0;
}
	
static CoS_StatusTypeDef CoS_ExBD_IO_Init(CoS_EXBD_HandleTypeDef *hObj)
{
	switch (hObj->Instance->current_extend)
	{
		case Exbd_USB:
		case Exbd_BT:
		{
			GPIO_InitTypeDef GPIO_InitStruct;
			GPIO_InitStruct.Pin = 0xff;
			GPIO_InitStruct.Mode = GPIO_MODE_OUTPUT_PP;
			GPIO_InitStruct.Pull = GPIO_PULLUP;
			HAL_GPIO_Init(GPIOD, &GPIO_InitStruct);
			HAL_GPIO_WritePin(GPIOD,0xff,GPIO_PIN_SET);
			break;
		}
		default:
			break;
	}
	return CoS_STATE_OK;
}


void CoS_EXBD_GPIOD_Init(CoS_EXBD_HandleTypeDef *hObj,u8 mode)
{
	GPIO_InitTypeDef GPIO_InitStruct;
	GPIO_InitStruct.Pin = EXBD_TEST_Pin &&0xff;
	GPIO_InitStruct.Mode = GPIO_MODE_INPUT;
	GPIO_InitStruct.Pull = GPIO_NOPULL;
	HAL_GPIO_Init(EXBD_TEST_GPIO_Port, &GPIO_InitStruct);
}
CoS_StatusTypeDef CoS_EXBD_Init(CoS_EXBD_HandleTypeDef *hObj,
								CoS_EXBD_InitTypeDef* Init)
{
	/* Check the handle allocation */
	if(hObj == NULL)
	{
		return CoS_STATE_ERROR;
	}
	if( Init->Instance == NULL )
	{
		return CoS_STATE_ERROR;
	}
	__CoS_LOCK(hObj);
	if(hObj->iState == CoS_STATE_UNSET)
	{  
		/* Allocate lock resource and initialize it */
		hObj->Lock = CoS_UNLOCKED;
		hObj->Instance=Init->Instance;
		//hObj->Instance->hUART=Init->hUART;
	}
	/* Change peripheral state */
	hObj->iState = CoS_STATE_BUSY;
	
	hObj->Instance->current_extend= (CoS_EXBD_NUMTypeDef)CoS_EXBD_check(hObj);
	CoS_ExBD_IO_Init(hObj);
	/* Init the low level hardware */
	if (INIT_INFO_printf)
		printf("{}Found extend board %d(%s).\r\n",\
	hObj->Instance->current_extend,Exboard_nameStr[hObj->Instance->current_extend]);
	/* Initialize the Instance state */
	hObj->iState = CoS_STATE_OK;
	__CoS_UNLOCK(hObj);
	/* Return function status */
	return CoS_STATE_OK;
}
CoS_StatusTypeDef CoS_EXBD_PDControl(CoS_EXBD_HandleTypeDef *hObj,u8 pins,CoS_EXBD_PDCONTROL ctrl)
{
	/* Check the handle allocation */
	if(hObj == NULL || hObj->iState!=CoS_STATE_OK)
	{
		return CoS_STATE_ERROR;
	}
	u8 tmp;
	if (hObj->Instance->current_extend==Exbd_USB || hObj->Instance->current_extend==Exbd_BT)
	{
		switch (ctrl)
		{
			case ExPD_VALUE:
			{
				GPIOD->ODR=(GPIOD->ODR&0xffffff00)+(~(pins));
				break;
			}
			case ExPD_ON:
			{
				tmp=~(GPIOD->ODR&0x000000ff);
				tmp|=pins;
				GPIOD->ODR=(GPIOD->ODR&0xffffff00)+(~(tmp));
				break;
			}
			case ExPD_OFF:
			{
				tmp=~(GPIOD->ODR&0x000000ff);
				tmp&=(~pins);
				GPIOD->ODR=(GPIOD->ODR&0xffffff00)+(~(tmp));
				break;
			}
			case ExPD_TOGGLE:
			{
				tmp=~(GPIOD->ODR&0x000000ff);
				tmp^=pins;
				GPIOD->ODR=(GPIOD->ODR&0xffffff00)+(~(tmp));
				break;
			}
			default:
				break;
		}
//	
	
	}
	/* Return function status */
	return CoS_STATE_OK;
}