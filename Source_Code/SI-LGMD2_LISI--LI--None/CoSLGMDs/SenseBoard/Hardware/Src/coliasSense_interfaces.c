#include "coliasSense_interfaces.h"
#include "coliasSense_board.h"
#include "delay.h"

#define SCL_H     (HAL_GPIO_WritePin(hObj->PinSCL.Port,hObj->PinSCL.Pin,GPIO_PIN_SET))
#define SCL_L     (HAL_GPIO_WritePin(hObj->PinSCL.Port,hObj->PinSCL.Pin,GPIO_PIN_RESET))

#define SDA_H     (HAL_GPIO_WritePin(hObj->PinSDA.Port,hObj->PinSDA.Pin,GPIO_PIN_SET))
#define SDA_L     (HAL_GPIO_WritePin(hObj->PinSDA.Port,hObj->PinSDA.Pin,GPIO_PIN_RESET))

#define SDA_STATE (HAL_GPIO_ReadPin(hObj->PinSDA.Port,hObj->PinSDA.Pin))

static uint8_t           sI2C_Start(CoS_sI2C_TypeDef* hObj);
static void              sI2C_Stop(CoS_sI2C_TypeDef* hObj);
static void              sI2C_Ack(CoS_sI2C_TypeDef* hObj);
static void              sI2C_NoAck(CoS_sI2C_TypeDef* hObj);
static uint8_t           sI2C_WaitAck(CoS_sI2C_TypeDef* hObj);
static uint8_t           sI2C_Read(CoS_sI2C_TypeDef* hObj);
static void              sI2C_Write(CoS_sI2C_TypeDef* hObj,uint8_t data);
static CoS_StatusTypeDef sI2C_WriteReg(CoS_sI2C_TypeDef* hObj);
static CoS_StatusTypeDef sI2C_ReadReg(CoS_sI2C_TypeDef* hObj);
static CoS_StatusTypeDef sI2C_WriteRegs(CoS_sI2C_TypeDef* hObj);
static CoS_StatusTypeDef sI2C_ReadRegs(CoS_sI2C_TypeDef* hObj);
static u8* const mode_nameStr[2]={"I2C","SCCB"};

static uint8_t sI2C_Start(CoS_sI2C_TypeDef* hObj)
{

	__IO u8 i=0;
	SDA_H;
	SCL_H;
	delay_short(hObj->nDelay);
	while((!SDA_STATE))
	{
		if(i++>=249)
			return (1);
		delay_short(hObj->nDelay);
	}
	SDA_L;
	delay_short(hObj->nDelay);
	i=0;
	while((SDA_STATE == 1) && (i<250))
	{
		if(i++>=249)
			return (1);
		delay_short(hObj->nDelay);
	}        
	SDA_L;
	SCL_L;
	delay_short(hObj->nDelay);
	return 0;                //SUCCESS;
}

static void sI2C_Stop(CoS_sI2C_TypeDef* hObj)
{
	SDA_L;
	delay_short(hObj->nDelay);
	SCL_H;
	delay_short(hObj->nDelay);
	SDA_H;
	delay_short(hObj->nDelay);
}

static void sI2C_Ack(CoS_sI2C_TypeDef* hObj)
{
//	SCL_L;
//	delay_short(hObj->nDelay);
	SDA_L;
	delay_short(hObj->nDelay);
	SCL_H;
	delay_short(hObj->nDelay);
	SCL_L;
	delay_short(hObj->nDelay);
	SDA_L;
	delay_short(hObj->nDelay);
}

static void sI2C_NoAck(CoS_sI2C_TypeDef* hObj)
{
	SDA_H;
	delay_short(hObj->nDelay);
	SCL_H;
	delay_short(hObj->nDelay);
	SCL_L;
	delay_short(hObj->nDelay);
	SDA_L;
	delay_short(hObj->nDelay);
}

static uint8_t sI2C_WaitAck(CoS_sI2C_TypeDef* hObj)
{
	//SCL_L;
	//delay_short(hObj->nDelay);
	SDA_H;
	delay_short(hObj->nDelay);
	SCL_H;
	delay_short(hObj->nDelay);
	if ((SDA_STATE))
	{
		SCL_L;
		delay_short(hObj->nDelay);
			return 1; //got no ack
	}
	SCL_L;
	//delay_short(hObj->nDelay);
	return 0;  //got ack
}

static uint8_t sI2C_Read(CoS_sI2C_TypeDef* hObj)
{
	u8 read=0x00,j;
	//SDA_H; 
	delay_short(hObj->nDelay);
	for(j=8;j>0;j--) 
	{		     
		SCL_H;
		delay_short(hObj->nDelay);
		read=SDA_STATE+(read<<1); 
		SCL_L;
		delay_short(hObj->nDelay);
	}
	return read;
}

static void sI2C_Write(CoS_sI2C_TypeDef* hObj,uint8_t data) 
{
	u8 i=8;
	while(i--)
	{
		if(data&0x80)
			SDA_H;  
		else 
			SDA_L;   
		data<<=1;
		delay_short(hObj->nDelay);
		SCL_H;
		delay_short(hObj->nDelay);
		SCL_L;
	}
//	SCL_L;
//	delay_short(hObj->nDelay);
}  

static CoS_StatusTypeDef sI2C_WriteReg(CoS_sI2C_TypeDef* hObj)
{
	hObj->CC.state=CoS_STATE_BUSY_T;
	if(sI2C_Start(hObj))
	{
		sI2C_Stop(hObj);
		hObj->CC.state=CoS_STATE_OK;
		return CoS_STATE_ERROR;
	}
	sI2C_Write(hObj,hObj->Dev_Addr);   //I2C_SendByte(((REG_Address & 0x0700) >>7) | SlaveAddress & 0xFFFE);
	if(sI2C_WaitAck(hObj))
	{
		sI2C_Stop(hObj); //got no ack, busy or error
		hObj->CC.state=CoS_STATE_OK;
		return CoS_STATE_BUSY;
	}
	sI2C_Write(hObj,hObj->DR);
	if(sI2C_WaitAck(hObj))
	{
		sI2C_Stop(hObj); //got no ack, busy or error
		hObj->CC.state=CoS_STATE_OK;
		return CoS_STATE_BUSY;
	}
	sI2C_Write(hObj,*(hObj->Tx_Base));
	if(sI2C_WaitAck(hObj))
	{
		sI2C_Stop(hObj); //got no ack, busy or error
		hObj->CC.state=CoS_STATE_OK;
		return CoS_STATE_BUSY;
	}
	sI2C_Stop(hObj);
	hObj->Tx_Count=1;
	hObj->CC.state=CoS_STATE_OK;
	return CoS_STATE_OK;
}

static CoS_StatusTypeDef sI2C_ReadReg(CoS_sI2C_TypeDef* hObj)
{
	hObj->CC.state=CoS_STATE_BUSY_R;
	if(sI2C_Start(hObj))
	{
		sI2C_Stop(hObj);
		hObj->CC.state=CoS_STATE_OK;
		return CoS_STATE_ERROR;
	}
	sI2C_Write(hObj,hObj->Dev_Addr);   //I2C_SendByte(((REG_Address & 0x0700) >>7) | SlaveAddress & 0xFFFE);
	if(sI2C_WaitAck(hObj))
	{
		sI2C_Stop(hObj); //got no ack, busy or error
		hObj->CC.state=CoS_STATE_OK;
		return CoS_STATE_BUSY;
	}
	sI2C_Write(hObj,hObj->DR);
	if(sI2C_WaitAck(hObj))
	{
		sI2C_Stop(hObj); //got no ack, busy or error
		hObj->CC.state=CoS_STATE_OK;
		return CoS_STATE_BUSY;
	}
	
	if (hObj->CC.mode==sI2CMode_SCCB)
		sI2C_Stop(hObj);
	if(sI2C_Start(hObj))
	{
		sI2C_Stop(hObj);
		hObj->CC.state=CoS_STATE_OK;
		return CoS_STATE_ERROR;
	}
	sI2C_Write(hObj,(hObj->Dev_Addr)+1);
	if(sI2C_WaitAck(hObj))
	{
		sI2C_Stop(hObj); //got no ack, busy or error
		hObj->CC.state=CoS_STATE_OK;
		return CoS_STATE_BUSY;
	}
	*(hObj->Rx_Base)=sI2C_Read(hObj);
	sI2C_NoAck(hObj);
	sI2C_Stop(hObj);
	hObj->Rx_Count=1;
	hObj->CC.state=CoS_STATE_OK;
	return CoS_STATE_OK;
}

static CoS_StatusTypeDef sI2C_WriteRegs(CoS_sI2C_TypeDef* hObj)
{
	u8 count=hObj->Tx_Count;
	u8 nums;
	hObj->CC.state=CoS_STATE_BUSY_T;
	if (hObj->CC.mode==sI2CMode_I2C)//Continuous write is supported only in I2C mode 
	{
		if(sI2C_Start(hObj))
		{
			sI2C_Stop(hObj);
			hObj->CC.state=CoS_STATE_OK;
			return CoS_STATE_ERROR;
		}
		sI2C_Write(hObj,hObj->Dev_Addr);   //I2C_SendByte(((REG_Address & 0x0700) >>7) | SlaveAddress & 0xFFFE);
		if(sI2C_WaitAck(hObj))
		{
			sI2C_Stop(hObj); //got no ack, busy or error
			hObj->CC.state=CoS_STATE_OK;
			return CoS_STATE_BUSY;
		}
		sI2C_Write(hObj,hObj->DR);
		if(sI2C_WaitAck(hObj))
		{
			sI2C_Stop(hObj); //got no ack, busy or error
			hObj->CC.state=CoS_STATE_OK;
			return CoS_STATE_BUSY;
		}
		nums=hObj->Tx_Length;
		//for(count=hObj->Tx_Count;count<hObj->Tx_Length;count++)
		while (count<nums)
		{
			sI2C_Write(hObj,*((hObj->Tx_Base)+count));
			hObj->Tx_Count=count++;
			if(sI2C_WaitAck(hObj))
			{
				sI2C_Stop(hObj); //got no ack, busy or error
				hObj->CC.state=CoS_STATE_OK;
				return CoS_STATE_BUSY;
			}
		}
		sI2C_Stop(hObj);
		//hObj->Tx_Count=1;
		hObj->CC.state=CoS_STATE_OK;
		return CoS_STATE_OK;
	}
	else if (hObj->CC.mode==sI2CMode_SCCB)
	{
		for(count=hObj->Tx_Count;count<hObj->Tx_Length;count++)
		{
			if(sI2C_Start(hObj))
			{
				sI2C_Stop(hObj);
				hObj->CC.state=CoS_STATE_OK;
				return CoS_STATE_ERROR;
			}
			sI2C_Write(hObj,hObj->Dev_Addr);   //I2C_SendByte(((REG_Address & 0x0700) >>7) | SlaveAddress & 0xFFFE);
			if(sI2C_WaitAck(hObj))
			{
				sI2C_Stop(hObj); //got no ack, busy or error
				hObj->CC.state=CoS_STATE_OK;
				return CoS_STATE_BUSY;
			}
			sI2C_Write(hObj,(hObj->DR)+count);
			if(sI2C_WaitAck(hObj))
			{
				sI2C_Stop(hObj); //got no ack, busy or error.
				hObj->CC.state=CoS_STATE_OK;
				return CoS_STATE_BUSY;
			}
			sI2C_Write(hObj,*((hObj->Tx_Base)+count));
			if(sI2C_WaitAck(hObj))
			{
				sI2C_Stop(hObj); //got no ack, busy or error
				hObj->CC.state=CoS_STATE_OK;
				return CoS_STATE_BUSY;
			}
			sI2C_Stop(hObj);
			hObj->Tx_Count=count;
		}
		hObj->CC.state=CoS_STATE_OK;
		return CoS_STATE_OK;
	}
	else
	{
		hObj->CC.state=CoS_STATE_ERROR;
		return CoS_STATE_ERROR;
	}
}

static CoS_StatusTypeDef sI2C_ReadRegs(CoS_sI2C_TypeDef* hObj)
{
	u8 count = 0;
	u8 temp,nums;
	hObj->CC.state=CoS_STATE_BUSY_T;
	if (hObj->CC.mode==sI2CMode_I2C)//Continuous write is supported only in I2C mode 
	{
		if(sI2C_Start(hObj))
		{
			sI2C_Stop(hObj);
			hObj->CC.state=CoS_STATE_OK;
			return CoS_STATE_ERROR;
		}
		sI2C_Write(hObj,hObj->Dev_Addr);   //I2C_SendByte(((REG_Address & 0x0700) >>7) | SlaveAddress & 0xFFFE);
		if(sI2C_WaitAck(hObj))
		{
			sI2C_Stop(hObj); //got no ack, busy or error
			hObj->CC.state=CoS_STATE_OK;
			return CoS_STATE_BUSY;
		}
		sI2C_Write(hObj,hObj->DR);
		if(sI2C_WaitAck(hObj))
		{
			sI2C_Stop(hObj); //got no ack, busy or error
			hObj->CC.state=CoS_STATE_OK;
			return CoS_STATE_BUSY;
		}
		if(sI2C_Start(hObj))
		{
			sI2C_Stop(hObj);
			hObj->CC.state=CoS_STATE_OK;
			return CoS_STATE_ERROR;
		}
		sI2C_Write(hObj,(hObj->Dev_Addr)+1);   //I2C_SendByte(((REG_Address & 0x0700) >>7) | SlaveAddress & 0xFFFE);
		if(sI2C_WaitAck(hObj))
		{
			sI2C_Stop(hObj); //got no ack, busy or error
			hObj->CC.state=CoS_STATE_OK;
			return CoS_STATE_BUSY;
		}
		nums=hObj->Rx_Length;
		count=hObj->Rx_Count;
		if (nums>1)
		{
			while(count<nums)
			{
				temp=sI2C_Read(hObj);
				if (count==nums-1)
					sI2C_NoAck(hObj);
				else
					sI2C_Ack(hObj);
				*((hObj->Rx_Base)+count)=temp;
				hObj->Rx_Count=++count;
			}
		}
		sI2C_Stop(hObj);
		//hObj->Tx_Count=1;
		hObj->CC.state=CoS_STATE_OK;
		return CoS_STATE_OK;
	}
	else
	{
		for (count=hObj->Rx_Count;count<hObj->Rx_Length;count++)
		{
			if(sI2C_Start(hObj))
			{
				sI2C_Stop(hObj);
				hObj->CC.state=CoS_STATE_OK;
				return CoS_STATE_ERROR;
			}
			sI2C_Write(hObj,hObj->Dev_Addr);   //I2C_SendByte(((REG_Address & 0x0700) >>7) | SlaveAddress & 0xFFFE);
			if(sI2C_WaitAck(hObj))
			{
				sI2C_Stop(hObj); //got no ack, busy or error
				hObj->CC.state=CoS_STATE_OK;
				return CoS_STATE_BUSY;
			}
			sI2C_Write(hObj,(hObj->DR)+count);
			if(sI2C_WaitAck(hObj))
			{
				sI2C_Stop(hObj); //got no ack, busy or error
				hObj->CC.state=CoS_STATE_OK;
				return CoS_STATE_BUSY;
			}
			sI2C_Stop(hObj);
			if(sI2C_Start(hObj))
			{
				sI2C_Stop(hObj);
				hObj->CC.state=CoS_STATE_OK;
				return CoS_STATE_ERROR;
			}
			sI2C_Write(hObj,(hObj->Dev_Addr)+1);
			if(sI2C_WaitAck(hObj))
			{
				sI2C_Stop(hObj); //got no ack, busy or error
				hObj->CC.state=CoS_STATE_OK;
				return CoS_STATE_BUSY;
			}
			*((hObj->Rx_Base)+count)=sI2C_Read(hObj);
			sI2C_NoAck(hObj);
			sI2C_Stop(hObj);
			hObj->Rx_Count=count;
		}
		hObj->CC.state=CoS_STATE_OK;
		return CoS_STATE_OK;
	}
	
}


static void sI2C_IO_Init(CoS_sI2C_HandleTypeDef* hObj,CoS_sI2C_InitTypeDef* Init)
{
	GPIO_InitTypeDef PInit;
	hObj->Instance->PinSCL=Init->PinSCL;
	hObj->Instance->PinSDA=Init->PinSDA;
	PInit.Pin=Init->PinSCL.Pin;
	PInit.Mode=GPIO_MODE_OUTPUT_OD;
	PInit.Pull=GPIO_PULLUP;
	PInit.Speed=GPIO_SPEED_FREQ_MEDIUM;
	HAL_GPIO_Init(Init->PinSCL.Port, &PInit);
	PInit.Pin=Init->PinSDA.Pin;
	HAL_GPIO_Init(Init->PinSDA.Port, &PInit);
	
	hObj->Instance->CC.mode=Init->mode;
	hObj->Dev_Addr=Init->Dev_Addr;
	hObj->Instance->nDelay=Init->nDelay;
}



CoS_StatusTypeDef CoS_sI2C_Init(CoS_sI2C_HandleTypeDef *hObj,CoS_sI2C_InitTypeDef* Init)
{
	/* Check the handle allocation */
	if(hObj == NULL)
	{
		return CoS_STATE_ERROR;
	}
	if( Init->PinSCL.Port == NULL ||
		Init->PinSDA.Port == NULL ||
		Init->nDelay == 0|| 
		Init->Dev_Addr==0||
		Init->Instance==NULL)
	{
		return CoS_STATE_ERROR;
	}
	__CoS_LOCK(hObj);
	if(hObj->iState == CoS_STATE_UNSET)
	{
		hObj->Lock = CoS_UNLOCKED;
		hObj->Instance=Init->Instance;
	}/* Change peripheral state */
	hObj->iState = CoS_STATE_BUSY;
	/* Init the low level hardware */
	sI2C_IO_Init(hObj,Init);
	/* Initialize the Instance state */
	
	
	sI2C_Stop(hObj->Instance);
	CoS_sI2C_TxRST(hObj);
	hObj->Instance->CC.state=CoS_STATE_OK;
	hObj->iState = CoS_STATE_OK;
	if (INIT_INFO_printf)
		printf("{sI2C Initialize OK.}\r\n");
	CoS_sI2C_StateReport(hObj,1);
	__CoS_UNLOCK(hObj);
	/* Return function status */
	return CoS_STATE_OK;
}

void CoS_sI2C_StateReport(CoS_sI2C_HandleTypeDef *hObj,u8 trans_mode_blocking)
{
	if (trans_mode_blocking)
	{
		printf("{  sI2C Initialize state:%d}\r\n",hObj->iState);
		if (hObj->iState!=CoS_STATE_OK)
			return;
		printf("{  interface mode: %s}\r\n",mode_nameStr[hObj->Instance->CC.mode]);
		printf("{  device address: %d(0x%2x)}\r\n",hObj->Dev_Addr,hObj->Dev_Addr);
		printf("{  pin SCL:");
		CoS_Pin_StateReport(&hObj->Instance->PinSCL,trans_mode_blocking);
		printf("{  pin SDA:");
		CoS_Pin_StateReport(&hObj->Instance->PinSDA,trans_mode_blocking);
		printf("{  pulse delay: %d}\r\n",hObj->Instance->nDelay);
	}
	else
		return;
}
CoS_StatusTypeDef CoS_sI2C_setAddr(CoS_sI2C_HandleTypeDef *hObj,uint8_t Dev_Addr)
{
	if(hObj == NULL)
	{
		return CoS_STATE_ERROR;
	}
	__CoS_LOCK(hObj);
	hObj->Instance->Dev_Addr=Dev_Addr;
	__CoS_UNLOCK(hObj);
	return CoS_STATE_OK;
}
CoS_StatusTypeDef CoS_sI2C_TxRST(CoS_sI2C_HandleTypeDef *hObj)
{
	if(hObj == NULL)
	{
		return CoS_STATE_ERROR;
	}
	__CoS_LOCK(hObj);
	//TODO: check transmit status
	hObj->Instance->Rx_Count=0;
	hObj->Instance->Rx_Length=0;
	hObj->Instance->Tx_Count=0;
	hObj->Instance->Tx_Length=0;
	hObj->Instance->Rx_Base=&(hObj->Instance->CC.userdata);
	hObj->Instance->Tx_Base=&(hObj->Instance->CC.userdata);
	hObj->Instance->CC.state=CoS_STATE_OK;
	__CoS_UNLOCK(hObj);
	return CoS_STATE_OK;
}

CoS_StatusTypeDef CoS_sI2C_WriteReg(CoS_sI2C_HandleTypeDef *hObj,uint8_t regID,uint8_t regData)
{
	CoS_StatusTypeDef sta;
	if(hObj == NULL)
		return CoS_STATE_ERROR;
	if (hObj->iState!=CoS_STATE_OK)
		return CoS_STATE_ERROR;
	if (hObj->Instance->CC.state!=CoS_STATE_OK)
		return CoS_STATE_BUSY;
	__CoS_LOCK(hObj);
	hObj->Instance->Dev_Addr=hObj->Dev_Addr;
	hObj->Instance->DR=regID;
	hObj->Instance->Tx_Base=&regData;
	hObj->Instance->Tx_Length=1;
	sta=sI2C_WriteReg(hObj->Instance);
	hObj->Instance->CC.state=CoS_STATE_OK;
	__CoS_UNLOCK(hObj);
	return sta;
}

CoS_StatusTypeDef CoS_sI2C_ReadReg(CoS_sI2C_HandleTypeDef *hObj,uint8_t regID,uint8_t* pData)
{
	CoS_StatusTypeDef sta;
	if(hObj == NULL)
		return CoS_STATE_ERROR;
	if (hObj->iState!=CoS_STATE_OK)
		return CoS_STATE_ERROR;
	if (hObj->Instance->CC.state!=CoS_STATE_OK)
		return CoS_STATE_BUSY;
	__CoS_LOCK(hObj);
	hObj->Instance->Dev_Addr=hObj->Dev_Addr;
	hObj->Instance->DR=regID;
	hObj->Instance->Rx_Base=pData;
	hObj->Instance->Rx_Length=1;
	sta=sI2C_ReadReg(hObj->Instance);
	//hObj->Instance->CC.state=CoS_STATE_OK;
	__CoS_UNLOCK(hObj);
	return sta;
}

CoS_StatusTypeDef CoS_sI2C_WriteRegs(CoS_sI2C_HandleTypeDef *hObj,uint8_t regID,uint8_t* Tx_Base,uint8_t  Tx_Length)
{
	CoS_StatusTypeDef sta;
	if(hObj == NULL)
		return CoS_STATE_ERROR;
	if (hObj->iState!=CoS_STATE_OK)
		return CoS_STATE_ERROR;
	if (hObj->Instance->CC.state!=CoS_STATE_OK)
		return CoS_STATE_BUSY;
	__CoS_LOCK(hObj);
	hObj->Instance->Dev_Addr=hObj->Dev_Addr;
	hObj->Instance->DR=regID;
	hObj->Instance->Tx_Base=Tx_Base;
	hObj->Instance->Tx_Length=Tx_Length;
	hObj->Instance->Tx_Count=0;
	sta=sI2C_WriteRegs(hObj->Instance);
	//hObj->Instance->CC.state=CoS_STATE_OK;
	__CoS_UNLOCK(hObj);
	return sta;
}
CoS_StatusTypeDef CoS_sI2C_ReadRegs(CoS_sI2C_HandleTypeDef *hObj,uint8_t regID,uint8_t* Rx_Base,uint8_t  Rx_Length)
{
	CoS_StatusTypeDef sta;
	if(hObj == NULL)
		return CoS_STATE_ERROR;
	if (hObj->iState!=CoS_STATE_OK)
		return CoS_STATE_ERROR;
	if (hObj->Instance->CC.state!=CoS_STATE_OK)
		return CoS_STATE_BUSY;
	__CoS_LOCK(hObj);
	hObj->Instance->Dev_Addr=hObj->Dev_Addr;
	hObj->Instance->DR=regID;
	hObj->Instance->Rx_Base=Rx_Base;
	hObj->Instance->Rx_Length=Rx_Length;
	hObj->Instance->Rx_Count=0;
	sta=sI2C_ReadRegs(hObj->Instance);
	//hObj->Instance->CC.state=CoS_STATE_OK;
	__CoS_UNLOCK(hObj);
	return sta;
}









