#include "coliasSense_board.h"

 CoS_BIO_TypeDef BIO1;
static CoS_BIO_HandleTypeDef HBIO1;

 CoS_CLED_TypeDef CLED1;
static CoS_CLED_HandleTypeDef HCLED1;
//CoS_BIO_TypeDef BIO2;

static uint8_t UART1_Buff[UART1_Buff_size];
 CoS_UART_TypeDef CoS_UART1;
static CoS_UART_HandleTypeDef HCoS_UART1;

static uint8_t UART4_Buff[UART4_Buff_size];
 CoS_UART_TypeDef CoS_UART4;
static CoS_UART_HandleTypeDef HCoS_UART4;

 CoS_sI2C_TypeDef CoS_sI2C_Camera;
static CoS_sI2C_HandleTypeDef HCoS_sI2C_Camera;

 CoS_sI2C_TypeDef CoS_sI2C_MPU;
static CoS_sI2C_HandleTypeDef HCoS_sI2C_MPU;
static CoS_sI2C_HandleTypeDef HCoS_sI2C_MAG;


CoS_CAM_TypeDef CoS_Camera;
static CoS_CAM_HandleTypeDef HCoS_Camera;

CoS_Motion_TypeDef Cos_Motion;
static CoS_Motion_HandleTypeDef HCos_Motion;

CoS_EXBD_TypeDef CoS_EXBD;
static CoS_EXBD_HandleTypeDef HCoS_EXBD;

//CoS_LGMD1_TypeDef CoS_LGMD1;
//static CoS_LGMD1_HandleTypeDef HCoS_LGMD1;
////handle of Colias Sense Board definition////
//handle of Colias Sense Board control struct, local
//static CoS_CCSTypeDef CC_SenseBoard;
//handle of Colias Sense Board, GLOBAL

static TIM_HandleTypeDef htim1;
static TIM_HandleTypeDef htim2;
static TIM_HandleTypeDef htim3;
static TIM_HandleTypeDef htim5;
static RNG_HandleTypeDef hrng;
static UART_HandleTypeDef huart1;
static UART_HandleTypeDef huart4;

static DMA_HandleTypeDef hdma_dcmi;
static DCMI_HandleTypeDef hdcmi;

SenseBoard_HandleTypeDef HSenseBoard;
SenseBoard_HandleTypeDef* const  hCoS=&HSenseBoard;



CoS_StatusTypeDef CoS_BIO_structInit(SenseBoard_HandleTypeDef *hObj)
{
	CoS_BIO_InitTypeDef BIO_Init;
	BIO_Init.htim2=&htim2;
	BIO_Init.htim3=&htim3;
	BIO_Init.htim5=&htim5;
	BIO_Init.Instance=&BIO1;
	BIO_Init.BackgroundIRQn=TIM3_IRQn;
	CoS_BIO_Init(hObj->hHBIO,&BIO_Init);
	hObj->hHBIO->Instance->CC.mode=CLED_Freeze;
	return CoS_STATE_OK;
}

CoS_StatusTypeDef CoS_CLED_structInit(SenseBoard_HandleTypeDef *hObj)
{
	CoS_CLED_InitTypeDef Init;
	Init.htim1=&htim1;
	Init.hrng=&hrng;
	Init.Instance= &CLED1;
	CoS_CLED_Init(hObj->hHCLED,&Init);
	return CoS_STATE_OK;
}



CoS_StatusTypeDef CoS_UART1_structInit(SenseBoard_HandleTypeDef *hObj)
{
	CoS_UART_InitTypeDef CInit;
	uint8_t sta;
	UART_InitTypeDef HInit;
	HInit.BaudRate=115200;
	HInit.HwFlowCtl=UART_HWCONTROL_NONE;
	HInit.Mode=UART_MODE_TX_RX;
	HInit.OverSampling=UART_OVERSAMPLING_8;
	HInit.Parity=UART_PARITY_NONE;
	HInit.StopBits=UART_STOPBITS_1;
	HInit.WordLength=UART_WORDLENGTH_8B;
	huart1.Init=HInit;
	huart1.Instance=USART1;
	CInit.huart=&huart1;
	CInit.Instance=&CoS_UART1;
	CInit.RxBufferBase=UART1_Buff;
	CInit.RxBufferSize=UART1_Buff_size;
	CInit.Tx_PacketSize=10;
	CInit.escape1=0x0d;
	CInit.escape2=0x0a;
	CInit.mode=UART_MODE_UART;
	CInit.Readback=0;
	CInit.RxICOv=1;
	CInit.RxRAWF=1;
//	CInit.escape1=0x7d;
//	CInit.escape2=0x7e;
//	CInit.mode=UART_MODE_CPtP;
	sta=CoS_UART_Init(hObj->hHUART1,&CInit);
	return sta;
//	Init.hUART=&huart1;
//	Init.Instance=&CoS_UART1;
//	Init.RxBufferBase=UART1_Buff;
//	Init.RxBufferSize=UART1_Buff_size;
//	Init.Baud=115200U;
//	return CoS_STATE_OK;
}

//CoS_StatusTypeDef CoS_UART4_structInit(SenseBoard_HandleTypeDef *hObj)
//{
//	CoS_UART_InitTypeDef CInit;
//	uint8_t sta;
//	UART_InitTypeDef HInit;
//	HInit.BaudRate = 38400;
//	HInit.HwFlowCtl = UART_HWCONTROL_NONE;
//	HInit.Mode = UART_MODE_TX_RX;
//	HInit.OverSampling = UART_OVERSAMPLING_16;
//	HInit.Parity = UART_PARITY_NONE;
//	HInit.StopBits = UART_STOPBITS_1;
//	HInit.WordLength = UART_WORDLENGTH_8B;
//	huart4.Init = HInit;
//	huart4.Instance = UART4;
//	CInit.huart = &huart4;
//	CInit.Instance = &CoS_UART4;
//	CInit.RxBufferBase = UART4_Buff;
//	CInit.RxBufferSize = UART4_Buff_size;
//	CInit.Tx_PacketSize = 10;
//	//CInit.escape1 = 0x0d;
//	//CInit.escape2 = 0x0a;
//	//CInit.mode = UART_MODE_UART;
//	CInit.Readback = 0;
//	CInit.RxICOv = 1;	//???
//	CInit.RxRAWF = 0;
//	CInit.escape1=0x7d;
//	CInit.escape2=0x7F;
//	CInit.mode=UART_MODE_CPtP;
//	sta = CoS_UART_Init(hObj->hHUART4, &CInit);
//	if (sta == CoS_STATE_OK)
//		printf("USART4 Init OK.\r\n");
//	return sta;

//}

CoS_StatusTypeDef CoS_sI2C_Camera_structInit(SenseBoard_HandleTypeDef *hObj)
{
	CoS_sI2C_InitTypeDef Init;
	Init.Dev_Addr=0x42;
	Init.Instance=&CoS_sI2C_Camera;
	Init.mode=sI2CMode_SCCB;
	Init.nDelay=0x0e;
	Init.PinSCL.Port=SCCB_SCL_GPIO_Port;
	Init.PinSCL.Pin=SCCB_SCL_Pin;
	Init.PinSDA.Port=SCCB_SDA_GPIO_Port;
	Init.PinSDA.Pin=SCCB_SDA_Pin;

	return CoS_sI2C_Init(hObj->hHsI2C_Camera,&Init);
}

CoS_StatusTypeDef CoS_sI2C_MPU_structInit(SenseBoard_HandleTypeDef *hObj)
{
	CoS_sI2C_InitTypeDef Init;
	u8 sta;
	Init.Dev_Addr=0xD0;
	Init.Instance=&CoS_sI2C_MPU;
	Init.mode=sI2CMode_I2C;
	Init.nDelay=0x0e;
	Init.PinSCL.Port=sI2C_SCL_GPIO_Port;
	Init.PinSCL.Pin=sI2C_SCL_Pin;
	Init.PinSDA.Port=sI2C_SDA_GPIO_Port;
	Init.PinSDA.Pin=sI2C_SDA_Pin;
	
	sta=CoS_sI2C_Init(hObj->hHsI2C_MPU,&Init);
	if (sta!=CoS_STATE_OK)
		return sta;
	Init.Dev_Addr=0x18;
	return CoS_sI2C_Init(hObj->hHsI2C_MAG,&Init);
}


CoS_StatusTypeDef CoS_cam_structInit(SenseBoard_HandleTypeDef *hObj)
{
	CoS_CAM_InitTypeDef Init;
	Init.hdcmi=hObj->hdcmi;
	Init.hdma_dcmi=hObj->hdma_dcmi;
	Init.hSSCB=hObj->hHsI2C_Camera;
	Init.Instance=&CoS_Camera;
	return CoS_CAM_Init(hObj->hHCoS_Camera,&Init);
}

CoS_StatusTypeDef Cos_Motion_structInit(SenseBoard_HandleTypeDef *hObj)
{
	CoS_Motion_InitTypeDef MInit;
	CoS_UART_InitTypeDef UInit;
	UART_InitTypeDef HInit;
	uint8_t sta;
	/*HAL huart Init structure */
	HInit.BaudRate = 38400;
	HInit.HwFlowCtl = UART_HWCONTROL_NONE;
	HInit.Mode = UART_MODE_TX_RX;
	HInit.OverSampling = UART_OVERSAMPLING_16;
	HInit.Parity = UART_PARITY_NONE;
	HInit.StopBits = UART_STOPBITS_1;
	HInit.WordLength = UART_WORDLENGTH_8B;
	hObj->huart4->Init= HInit;
	hObj->huart4->Instance = UART4;
	
	/*CoS HUART Init structure*/
	UInit.huart = (hObj->huart4);
	UInit.Instance = &CoS_UART4;
	UInit.RxBufferBase = UART4_Buff;
	UInit.RxBufferSize = UART4_Buff_size;
	UInit.Tx_PacketSize = 10;
	UInit.Readback = 0;
	UInit.RxICOv = 1;	//???
	UInit.RxRAWF = 0;
	UInit.escape1=0x78;
	UInit.escape2=0x7F;
	UInit.mode=UART_MODE_CPtP;
	sta=CoS_UART_Init(hObj->hHUART4,&UInit);
	if (sta!=CoS_STATE_OK)
		return sta;
	
	CoS_UART_Enable(hCoS->hHUART4,ENABLE);
	CoS_UART_RxRST(hObj->hHUART4);
	MInit.hHUART = hObj->hHUART4;
	MInit.Instance = &Cos_Motion;
	delay_ms(300);
	return CoS_Motion_Init(hObj->hHCos_Motion, &MInit);

}

CoS_StatusTypeDef CoS_EXBD_structInit(SenseBoard_HandleTypeDef *hObj)
{
	CoS_EXBD_InitTypeDef Init;

	Init.Instance=&CoS_EXBD;
	CoS_EXBD_Init(hObj->hHCoS_EXBD,&Init);
	
	if (Key2State)
	{
		printf("Key1 pushed, entering high speed UART (921600)\r\n");
		if (CoS_UART_SetBaud(hCoS->hHUART1,921600)==CoS_STATE_OK)
			hCoS->hHUART1->Instance->Tx_PacketSize=800;
	}
	return CoS_STATE_OK;

}
//CoS_StatusTypeDef CoS_LGMD1_structInit(SenseBoard_HandleTypeDef *hObj)
//{
//	CoS_LGMD1_InitTypeDef Init;
//	Init.Instance=&CoS_LGMD1;
//	Init.hCam=&HCoS_Camera;
//	Init.Image_BufBase=&Image[0][0][0];
//	Init.Diff_Image_BufBase=&Diff_Image[0][0][0];
//	Init.processRate=30;
//	Init.xstart=0;
//	Init.ystart=0;
//	Init.layer_Width=99;
//	Init.layer_Height=72;
//	Init.img_Width=99;
//	Init.img_Height=72;
//	Init.TdeCde=30;
//	Init.TdeCdeomega=500;
//	Init.Kf_threshold=128;
//	Init.Ts=128;
//	Init.TFFI=256;
//	Init.nsp=4;
//	Init.delta_C=0.01;
//	Init.Cw=4;
//	Init.Wi=0.4;
//	Init.alpha_ffi=0.02;
//	Init.C1=80000;
//	Init.C2=70000;
//	return CoS_LGMD1_Init(hObj->hHLGMD1,&Init);
//}

CoS_StatusTypeDef CoS_SenseBoard_Init(SenseBoard_HandleTypeDef *hObj)
{
//	hObj->CC=&CC_SenseBoard;
//	
	hObj->htim1=&htim1;
	hObj->htim2=&htim2;
	hObj->htim3=&htim3;
	hObj->htim5=&htim5;
	hObj->hrng=&hrng;
	hObj->huart1=&huart1;
	hObj->huart4=&huart4;
	hObj->hdcmi=&hdcmi;
	hObj->hdma_dcmi=&hdma_dcmi;
	
	
	hObj->hHBIO=&HBIO1;
	hObj->hHUART1=&HCoS_UART1;
	hObj->hHCLED=&HCLED1;
	hObj->hHsI2C_Camera=&HCoS_sI2C_Camera;
	hObj->hHsI2C_MPU=&HCoS_sI2C_MPU;
	hObj->hHsI2C_MAG=&HCoS_sI2C_MAG;
	hObj->hHCoS_Camera=&HCoS_Camera;
	hObj->hHUART4=&HCoS_UART4;
	hObj->hHCos_Motion = &HCos_Motion;
	hObj->hHCoS_EXBD=&HCoS_EXBD;
	
	
	CoS_BIO_structInit(hObj);
	CoS_UART1_structInit(hObj);
	CoS_CLED_structInit(hObj);
	CoS_sI2C_MPU_structInit(hObj);
	CoS_cam_structInit(hObj);

	Cos_Motion_structInit(hObj);
	CoS_EXBD_structInit(hObj);
	CoS_BIO_enable(hObj->hHBIO);
	CoS_CLED_Enable(hObj->hHCLED);
	
	
	
//	hObj->hHLGMD1=&HCoS_LGMD1;
//	CoS_LGMD1_structInit(hObj);
	
//	LGMDCC.hModel=&LGMDs[0];
	hLGMD.Model=&LGMD;
	hLGMD.hFrameCount=&(hObj->hHCoS_Camera->Instance->cam_hFrameCountount);
//	LGMDCC.hFrameCount=&(hObj->hHCoS_Camera->Instance->cam_hFrameCountount);
//	LGMDCC.imgTC;
//	LGMD_Param_Init(&LGMDCC);
	LGMD_Param_Init(&hLGMD);
	
//	LGMD2CC.hModel=&LGMDs[1];
//	LGMD2CC.hFrameCount=&(hObj->hHCoS_Camera->Instance->cam_hFrameCountount);
//	LGMDCC.imgTC;
//	LGMD_Param_Init(&LGMD2CC);
	
	printf("{}CPUID: %08X-%08X-%08X",CoS_MCUID(0),CoS_MCUID(1),CoS_MCUID(2));
	printf("-%08X\r\n",CoS_MCUID(3));
	
	return CoS_STATE_OK;
}



void CoS_TIM3_IRQ_Task(CoS_BIO_HandleTypeDef* hObj)
{
	CoS_StatusTypeDef status;
	CoS_LED_Update(hCoS->hHCLED);
	//CoS_UART_TxRun(hCoS->hHUART1);
	UART_transmit_run(hCoS->hHUART1);
	UART_transmit_run(hCoS->hHUART4);

//	status = ColiasBasicITEvents_Check(hCoS->hHCos_Motion->Instance, &(hCoS->hHCos_Motion->Instance->IT_Events));			
//			while (status != CoS_STATE_OK)
//	{
//		LED1_Toggle;
//		status = ColiasBasicITEvents_Check(hCoS->hHCos_Motion->Instance, &(hCoS->hHCos_Motion->Instance->IT_Events));			
//	}
	//MotionTask(hCoS->hHCos_Motion->Instance);
}