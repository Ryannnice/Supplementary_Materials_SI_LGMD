#ifndef __USERTER_H
#define __USERTER_H

#include "stm32f4xx_hal.h"
#include "coliasSense_def.h"



void UT_TermDeal(void);

#define CMBufferSize 50 //keep it greater than RXBufferSize
#define UT_MaxArgs 16

typedef union
{
	float Dfloat;
	uint32_t DU32;
	uint16_t DU16;
	uint8_t DU8;
	int32_t DS32;
	int16_t DS16;
	int8_t DS8;
}B32DataUType;

typedef enum
{
	UTypeID_float=7U,
	UTypeID_U32=1U,
	UTypeID_S32=2U,
	UTypeID_U16=3U,
	UTypeID_S16=4U,
	UTypeID_U8=5U,
	UTypeID_S8=6U,
}B32TypeIDdef;

typedef struct 
{
	B32DataUType data;
	B32TypeIDdef type;
}B32DataType;

typedef struct
{
	
	int32_t minnum;
	int32_t maxnum;
	u8 ctrl;
}UT_SvalRangeTypedef;
typedef struct
{
	
	uint32_t minnum;
	uint32_t maxnum;
	u8 ctrl;
}UT_UvalRangeTypedef;
typedef struct
{
	
	float minnum;
	float maxnum;
	u8 ctrl;
}UT_FvalRangeTypedef;

typedef struct
{
	B32TypeIDdef type;
	const char* name;
	const char* helpwd;
}CMD_NodeNameStr_Typedef;


typedef union
{
	UT_SvalRangeTypedef SRange;
	UT_UvalRangeTypedef URange;
	UT_FvalRangeTypedef FRange;
}UT_ValRangeUniondef;

typedef struct 
{
	void* addr;
	B32TypeIDdef type;
	char* name;
	char* helpwd;
	UT_ValRangeUniondef Cr;
}CMD_NodeTypedef_old;

typedef struct 
{
	void* hData;
	UT_ValRangeUniondef * hLim;
	 CMD_NodeNameStr_Typedef const*  hDetails;
}CMD_NodeTypedef;

typedef enum
{
	UT_ArgErr_OK =0x00,
	UT_ArgErr_overflow=0x01,
	UT_ArgErr_insufficient=0x02,
	UT_ArgErr_longopt=0x04,
	UT_ArgErr_unhandled=0x08,
	UT_ArgErr_unknown=0x80,
	UT_ArgErr_help=0x90,
	UT_ArgErr_OKmask =0x0f,
}CMD_ArgErrEnum;

typedef enum
{
	UT_ValErr_OK =0x00,
	UT_ValErr_noarg=0x01,
	UT_ValErr_noaddr=0x02,
	UT_ValErr_nodatatype=0x03,
	UT_ValErr_nan=0x04,
	UT_ValErr_empty=0x05,
	UT_ValErr_outrange=0x06,
	UT_ValErr_help=0x07,
	UT_ValErr_unknown=0x08,
}UT_ValErrEnum;


typedef struct
{
	char Cargc; //args count
	char Nargc; //last valid arg number
	char ErrID;
	char bzero[2];
	char* Cargv[UT_MaxArgs];
}CMDargTypedef;


typedef struct
{
	__IO u8    CM_Buff[CMBufferSize];
	__IO u8    CM_Buffcount;
	__IO u8    F_NotEmpty;
	__IO u8    F_DealBusy;
	__IO u8    C_Readback;
	__IO u8    C_TimeReadback;
	__IO u8    C_mode_CPtP;
	__IO u8    C_do_dataTranfer;
	__IO u8    S_cat;
	__IO u8    S_subcat;
	__IO u8    ErrID;
	__IO u8    userdata;
	CMDargTypedef ArgStruct;
}CMD_StructTypedef;

typedef struct
{
	 u8   start;
	 u8   period;
	 u16   tasks_to_exceute;
	 u16   tasks_remain;
}CMD_DataTransferTypedef;


extern CMD_StructTypedef UT_CMD;
/******************/
extern uint8_t Key2_virtual;
/******************/

//u8 UT_get_args(char * str, B32ValueStructTypedef* dest,u8 maxnum);
void UT_FlushBuff(CMD_StructTypedef* hUT);
//void UT_commandcpyUSART1(void);
u8 UT_commandcpy(void);
void UT_TermDeal(void);
void UT_Init(void);
void timed_data_init(void);
char timed_data_transmit(void);
char timed_data_set(u8 addr,u16 data);
void SendFrame(void);
















#endif