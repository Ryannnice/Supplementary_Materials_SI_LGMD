


#include "user_terminal.h"
#include "user_terminal_private.h"
#include "coliasSense_board.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


CMD_StructTypedef UT_CMD;
static CMD_DataTransferTypedef UT_Dtimer;

char UT_valSet(CMDargTypedef* Sarg,CMD_NodeTypedef* hTarget);
static int UT_Help(const CMD_NodeTypedef* pNode);
/***********************
Actions 
************************/



/*
Flush User terminal buffer and init it.
*/
void UT_FlushBuff(CMD_StructTypedef* hUT)
{
	u16 n;
	for (n=0;n<CMBufferSize;n++)
	{
		hUT->CM_Buff[n]=0;
		hUT->CM_Buffcount=0;
		hUT->F_NotEmpty=0;
	}
}
/*
Copy a command from UART receive buffer.
The number of characters should be less than UART buffer size
*/
//void UT_commandcpyUSART1(CMD_StructTypedef* hUT,SB_UartStructTypedef* hSB_UARTn)
//{
//	u16 n;
//	for (n=0;n<hSB_UARTn->Rx_count;n++)
//	{
//		
//		hUT->CM_Buff[n]=hSB_UARTn->Rx_Buff[n];
//	}
//	hUT->CM_Buffcount=n;
//	
//}
/*
Copy a command from selected buffer.
*/
u8 UT_commandcpy(void)
{
	if (hCoS->hHUART1->Instance->RxOK) //got a command from USART, which has high priority
	{
		UT_CMD.F_NotEmpty=1;
		memcpy(UT_CMD.CM_Buff,hCoS->hHUART1->Instance->Rx_Base,hCoS->hHUART1->Instance->Rx_Count);
		CoS_UART_RxRST(hCoS->hHUART1);
//		UT_commandcpyUSART1(&UT_CMD,&hSB_USART1);
//		Flush_UART_BUFF(&hSB_USART1);
//		hSB_USART1.F_RXOK=0;
		return 1;
	}
	else//from SD
	{
		
	}
	return 0;
}

/*
Parse buffer string into( int argc,char* argv ) format.
argument1: the source string buffer
argument2: the destination structure of Parsed arguments property
*/
int UT_CMD_StringParse(char* commandLine,CMDargTypedef* Sarg)
{
	char* p2;
	int n;
	Sarg->Cargc=0;
	Sarg->Nargc=0;
	Sarg->ErrID=0;
	for (n=0;n<UT_MaxArgs;n++)
		Sarg->Cargv[n]=&(Sarg->bzero[0]);
	p2=strtok(commandLine, " ");
	while (p2 && (Sarg->Cargc) < UT_MaxArgs - 1)
	{
		Sarg->Cargv[(Sarg->Cargc)++] = p2;
		p2 = strtok(0, " ");
	}
	Sarg->Cargv[Sarg->Cargc] = 0;
	return Sarg->Cargc;
}


/*
find a match command from optstr in given arguments
An approach of the GNU based getopt.
argument1: the source parsed arguments property structure
argument2: the index of input arguments to be parsed
argument3: the availible option string. Format should be :
	"[character1][min left args][character2][min left args]"
output: status.
	-1: error happened. the detail of error are stored in the input structure
	[character1]: 
*/
char UT_getopt( CMDargTypedef* Sarg, char idx, char* optstr)
{
	char* opt;
	signed char leastnum,n,sel;
	Sarg->ErrID=UT_ArgErr_unknown;
	if (idx>=Sarg->Cargc)//idx exceed avalable Arg counts
	{
		Sarg->ErrID=UT_ArgErr_overflow;
		return -1; 
	}
	sel=*(Sarg->Cargv[idx]);
	if (strlen(Sarg->Cargv[idx])>1)//selected string is not a single letter
	{
		if (sel=='-' && *(Sarg->Cargv[idx+1])=='h') //and it is a help command
		{	
			Sarg->ErrID=UT_ArgErr_help;
			return 'H';
		}
		else
		Sarg->ErrID=UT_ArgErr_longopt;
		//return -1;  //do not return, continue to use the first instead
	}
	for (n=0;n<strlen(optstr);n+=2)
	{
		if (sel==*(optstr+n))
		{
			leastnum=*(optstr+n+1)-48; //least number count
			if ((Sarg->Cargc)<(leastnum+idx+1))
			{
				Sarg->ErrID=UT_ArgErr_insufficient;
				return -1; //not enough args
			}
			else
			{				
				Sarg->ErrID&=UT_ArgErr_OKmask;
				Sarg->Nargc=idx;
				return sel;
			}
		}
	}
	Sarg->ErrID=UT_ArgErr_unhandled;
	return -1;
}
	
void Valprint_new(CMD_NodeTypedef* hTarget)
{
	switch (hTarget->hDetails->type)
	{
		case UTypeID_float:
			printf("%s: %f\r\n",hTarget->hDetails->name,*(float*)(hTarget->hData));
		break;
		case UTypeID_U32:
			printf("%s: %d\r\n",hTarget->hDetails->name,*(uint32_t*)(hTarget->hData));
		break;
		case UTypeID_S32:
			printf("%s: %d\r\n",hTarget->hDetails->name,*(int32_t*)(hTarget->hData));
		break;
		case UTypeID_U16:
			printf("%s: %d\r\n",hTarget->hDetails->name,*(uint16_t*)(hTarget->hData));
		break;
		case UTypeID_S16:
			printf("%s: %d\r\n",hTarget->hDetails->name,*(int16_t*)(hTarget->hData));
		break;
		case UTypeID_U8:
			printf("%s: %d\r\n",hTarget->hDetails->name,*(uint8_t*)(hTarget->hData));
		break;
		case UTypeID_S8:
			printf("%s: %d\r\n",hTarget->hDetails->name,*(int8_t*)(hTarget->hData));
		break;
		default:
			break;
		
	}
}

void Valprint(const CMD_NodeTypedef_old* destVar)
{
	switch (destVar->type)
	{
		case UTypeID_float:
			printf("%s: %f\r\n",destVar->name,*(float*)(destVar->addr));
		break;
		case UTypeID_U32:
			printf("%s: %d\r\n",destVar->name,*(uint32_t*)(destVar->addr));
		break;
		case UTypeID_S32:
			printf("%s: %d\r\n",destVar->name,*(int32_t*)(destVar->addr));
		break;
		case UTypeID_U16:
			printf("%s: %d\r\n",destVar->name,*(uint16_t*)(destVar->addr));
		break;
		case UTypeID_S16:
			printf("%s: %d\r\n",destVar->name,*(int16_t*)(destVar->addr));
		break;
		case UTypeID_U8:
			printf("%s: %d\r\n",destVar->name,*(uint8_t*)(destVar->addr));
		break;
		case UTypeID_S8:
			printf("%s: %d\r\n",destVar->name,*(int8_t*)(destVar->addr));
		break;
		default:
			break;
		
	}
}


/*
get a unsigned value from argument string.
argument1: input string
argument2: destination data address
argument3: data range check struct
*/
char UT_val_getu(char* str,uint32_t* pdata,UT_UvalRangeTypedef* Cr)
{
	uint32_t val;
	char* tail;
	val=strtoul(str,&tail,0); //convert string to ulong
	if (str==tail)
	{
		if (*str=='-' && *(str+1)=='h')
			return UT_ValErr_help;//need help
		else if (*str==0)
			return UT_ValErr_empty;
		else
			return UT_ValErr_nan;//cannot get a number
	}
	else
	{
		if ((Cr->ctrl)&&0xf0) //check with max
		{
			if (val>(Cr->maxnum))
				return UT_ValErr_outrange;
		}
		if ((Cr->ctrl)&&0x0f) //check with min
		{
			if (val<(Cr->minnum))
				return UT_ValErr_outrange;
		}
		*pdata=val;
		return UT_ValErr_OK;
			
	}
}
/*
get a signed value from argument string.
argument1: input string
argument2: destination data address
argument3: data range check struct
*/
char UT_val_gets(char* str,int32_t* pdata, UT_SvalRangeTypedef* Cr)
{
	int32_t val;
	char* tail;
	val=strtol(str,&tail,0); //convert string to long
	if (str==tail)
	{
		if (*str=='-' && *(str+1)=='h')
			return UT_ValErr_help;//need help
		else if (*str==0)
			return UT_ValErr_empty;
		else
			return UT_ValErr_nan;//cannot get a number
	}
	else
	{
		if ((Cr->ctrl)&&0xf0) //check with max
		{
			if (val>(Cr->maxnum))
				return UT_ValErr_outrange;
		}
		if ((Cr->ctrl)&&0x0f) //check with min
		{
			if (val<(Cr->minnum))
				return UT_ValErr_outrange;
		}
		*pdata=val;
		return UT_ValErr_OK;
			
	}
}

/*
get a float value from argument string.
argument1: input string
argument2: destination data address
argument3: data range check struct
*/
char UT_val_getf(char* str,float* pdata,  const  UT_FvalRangeTypedef* Cr)
{
	float val;
	char* tail;
	val=strtof(str,&tail); //convert string to float
	if (str==tail)
	{
		if (*str=='-' && *(str+1)=='h')
			return UT_ValErr_help;//need help
		else if (*str==0)
			return UT_ValErr_empty;
		else
			return UT_ValErr_nan;//cannot get a number
	}
	else
	{
		if ((Cr->ctrl)&&0xf0) //check with max
		{
			if (val>(Cr->maxnum))
				return UT_ValErr_outrange;
		}
		if ((Cr->ctrl)&&0x0f) //check with min
		{
			if (val<(Cr->minnum))
				return UT_ValErr_outrange;
		}
		*pdata=val;
		return UT_ValErr_OK;
			
	}
}



char UT_outrangeReport(CMD_NodeTypedef* hTarget)
{
	UT_ValRangeUniondef* hRange=hTarget->hLim;
	switch (hTarget->hDetails->type)
	{
		case UTypeID_float:
		{
			float fmax,fmin;
			fmax=hRange->FRange.maxnum;
			fmin=hRange->FRange.minnum;
			if (hRange->FRange.ctrl&&0xf0)
			{
				if (hRange->FRange.ctrl&&0x0f)
					printf("error:arg of %s out of range (max %f, min %f).\r\n",hTarget->hDetails->name,fmax,fmin);
				else
					printf("error:arg of %s out of range (max %f).\r\n",hTarget->hDetails->name,fmax);
			}
			else
			{
				if (hRange->FRange.ctrl&&0x0f)
					printf("error:arg of %s out of range (min %f).\r\n",hTarget->hDetails->name,fmin);
				else
					printf("error:arg of %s out of range.\r\n",hTarget->hDetails->name);
			}
			return 0;
		}
		case UTypeID_U32:
		case UTypeID_U16:
		case UTypeID_U8:
		{
			uint32_t umax,umin;
			umax=hRange->URange.maxnum;
			umin=hRange->URange.minnum;
			if (hRange->URange.ctrl&&0xf0)
			{
				if (hRange->URange.ctrl&&0x0f)
					printf("error:arg of %s out of range (max %d, min %d).\r\n",hTarget->hDetails->name,umax,umin);
				else
					printf("error:arg of %s out of range (max %d).\r\n",hTarget->hDetails->name,umax);
			}
			else
			{
				if (hRange->URange.ctrl&&0x0f)
					printf("error:arg of %s out of range (min %d).\r\n",hTarget->hDetails->name,umin);
				else
					printf("error:arg of %s out of range.\r\n",hTarget->hDetails->name);
			}
			return 0;
		}
		case UTypeID_S32:
		case UTypeID_S16:
		case UTypeID_S8:
		{
			int32_t smax,smin;
			smax=hRange->SRange.maxnum;
			smin=hRange->SRange.minnum;
			if (hRange->SRange.ctrl&&0xf0)
			{
				if (hRange->SRange.ctrl&&0x0f)
					printf("error:arg of %s out of range (max %d, min %d).\r\n",hTarget->hDetails->name,smax,smin);
				else
					printf("error:arg of %s out of range (max %d).\r\n",hTarget->hDetails->name,smax);
			}
			else
			{
				if (hRange->SRange.ctrl&&0x0f)
					printf("error:arg of %s out of range (min %d).\r\n",hTarget->hDetails->name,smin);
				else
					printf("error:arg of %s out of range.\r\n",hTarget->hDetails->name);
			}
			return 0;
		}
		default:
		{
			printf("unknown error happend with arg of %s(outrange).\r\n",hTarget->hDetails->name);
			return 0;
		}
	}
}

/*
check input command arg, which could be no arg, help string.
Try to convert into number if has valid instruction and limit
Check following flags: 
1. has following argument. if yes, translate it 
//(not here)2. has valid target data address. if no, this is not a setting command
during translate stage, check the following

3. requires help string. if yesm tell the srting
4. has valid data type. if no, is a commnad string
5. within range. if no, tell the data range
******************************************************
                   has arg:       yes          no
-h?    data type:                            
 yes    yes                       -h        do nothing 
 yes    no                        -h        do nothing 
 no     yes                    translate    do nothing 
 no     no                      command     do nothing 
******************************************************
translate status
OK
outrange
not a number
unknown
*/
char UT_argCheck(CMDargTypedef* Sarg,CMD_NodeTypedef* hTarget,B32DataType* hData,u8 silent)
{
	int errn=UT_ArgErr_unknown;
	char* str;
	if ((Sarg->Nargc)+1>=(Sarg->Cargc)) //no arg
	{
		errn=UT_ValErr_noarg;
		return errn;
	}
	else //has arg, try to transfer arg into number
	{
		str=(Sarg->Cargv[(Sarg->Nargc)+1]);
		if (*str=='-' &&*(str+1)=='h')
		{
			UT_Help(hTarget);
			return UT_ValErr_help;
		}
		hData->type=hTarget->hDetails->type;
		switch (hTarget->hDetails->type)
		{
			case UTypeID_float:
			{
				float fval;
				//UT_FvalRangeTypedef Cr={0x11,0.8,-0.8};
				errn=UT_val_getf(str,&fval,&(hTarget->hLim->FRange));
				if (!errn)
				{
					hData->data.Dfloat=fval;
					errn=UT_ValErr_OK;
				}
				break;
			}
			case UTypeID_S32:
			{
				int32_t val;
				errn=UT_val_gets(str,&val,&(hTarget->hLim->SRange));
				if (!errn)
				{
					hData->data.DS32=val;
					errn=UT_ValErr_OK;
				}
				break;
			}
			case UTypeID_U32:
			{
				uint32_t val;
				errn=UT_val_getu(str,&val,&(hTarget->hLim->URange));
				if (!errn)
				{
					hData->data.DU32=val;
					errn=UT_ValErr_OK;
				}
				break;
			}
			case UTypeID_S16:
			{
				int val;
				errn=UT_val_gets(str,&val,&(hTarget->hLim->SRange));
				if (!errn)
				{
					hData->data.DS16=val;
					errn=UT_ValErr_OK;
				}
				break;
			}
			case UTypeID_U16:
			{
				uint32_t val;
				errn=UT_val_getu(str,&val,&(hTarget->hLim->URange));
				if (!errn)
				{
					hData->data.DU16=val;
					errn=UT_ValErr_OK;
				}
				break;
			}
			case UTypeID_S8:
			{
				int val;
				errn=UT_val_gets(str,&val,&(hTarget->hLim->SRange));
				if (!errn)
				{
					hData->data.DS8=val;
					errn=UT_ValErr_OK;
				}
				break;
			}
			case UTypeID_U8:
			{
				int val;
				errn=UT_val_getu(str,&val,&(hTarget->hLim->URange));
				if (!errn)
				{
					hData->data.DU8=val;
					errn=UT_ValErr_OK;
				}
				break;
			}
			default:
			{
				errn=UT_ValErr_nodatatype;//unknown data type
				break;
			}
		}
	}
	switch (errn)
	{
		case UT_ValErr_OK:
			return UT_ValErr_OK;
		case UT_ValErr_outrange:
		{
			UT_outrangeReport(hTarget);
			return UT_ValErr_outrange;
		}
		case UT_ValErr_nan:
		{
			if (!silent)
				printf("argument of %s :(%s) is not a valid number.\r\n",hTarget->hDetails->name,str);
			return UT_ValErr_nan;
		}
		default:
		{
			if (!silent)
				printf("error happened with argument of %s:(%s).\r\n",hTarget->hDetails->name,str);
			return UT_ValErr_unknown;
		}
	}
	
}
/*
set a data by specified address.
input: target data, translated 
*/
static char* const str_verb_is="is";
static char* const str_verb_has="has been changed to";
char UT_valSet(CMDargTypedef* Sarg,CMD_NodeTypedef* hTarget)
{
	B32DataType data;
	char errn;
	char* str_verb=str_verb_is;
	errn=UT_argCheck(Sarg,hTarget,&data,0);
	Sarg->Nargc++;
	switch (errn)
	{
		case UT_ValErr_OK: //set
		{
			if (hTarget->hData==0)
				return UT_ValErr_noaddr;
			str_verb=str_verb_has;
			switch (hTarget->hDetails->type)
			{
				case UTypeID_float:
				{
					*(float*)(hTarget->hData)=data.data.Dfloat;
					break;
				}
				case UTypeID_U32:
				{
					*(uint32_t*)(hTarget->hData)=data.data.DU32;
					break;
				}
				case UTypeID_U16:
				{
					*(uint16_t*)(hTarget->hData)=data.data.DU16;
					break;
				}
				case UTypeID_U8:
				{
					*(uint8_t*)(hTarget->hData)=data.data.DU8;
					break;
				}
				case UTypeID_S32:
				{
					*(int32_t*)(hTarget->hData)=data.data.DS32;
					break;
				}
				case UTypeID_S16:
				{
					*(int16_t*)(hTarget->hData)=data.data.DS16;
					break;
				}
				case UTypeID_S8:
				{
					*(int8_t*)(hTarget->hData)=data.data.DS8;
					break;
				}
			}
		}
		case UT_ValErr_noarg: //read
		{
			if (hTarget->hData==0)
				return UT_ValErr_noaddr;
			switch (hTarget->hDetails->type)
			{
				case UTypeID_float:
				{	printf("Value of %s %s %f\r\n",hTarget->hDetails->name,str_verb,*(float*)(hTarget->hData));
					break;}
				case UTypeID_U32:
				{	printf("Value of %s %s %d\r\n",hTarget->hDetails->name,str_verb,*(uint32_t*)(hTarget->hData));
					break;}
				case UTypeID_U16:
				{	printf("Value of %s %s %d\r\n",hTarget->hDetails->name,str_verb,*(uint16_t*)(hTarget->hData));
					break;}
				case UTypeID_U8:
				{	printf("Value of %s %s %d\r\n",hTarget->hDetails->name,str_verb,*(uint8_t*)(hTarget->hData));
					break;}
				case UTypeID_S32:
				{	printf("Value of %s %s %d\r\n",hTarget->hDetails->name,str_verb,*(int32_t*)(hTarget->hData));
					break;}
				case UTypeID_S16:
				{	printf("Value of %s %s %d\r\n",hTarget->hDetails->name,str_verb,*(int16_t*)(hTarget->hData));
					break;}
				case UTypeID_S8:
				{	printf("Value of %s %s %d\r\n",hTarget->hDetails->name,str_verb,*(int8_t*)(hTarget->hData));
					break;
				}
			}
		}
		default :
			return errn;
	}
}


/*
set a data by specified address.
argument1: input string arguments. The Nargc should be the last valid string
argument2: destination data address, defined by node
*/
char UT_SB_valSet_old(CMDargTypedef* Sarg,CMD_NodeTypedef_old* destVar)
{
	int errn=UT_ArgErr_unknown;
	char do_read=0, do_help=0;
	char* str;
	if (Sarg->Cargc-Sarg->Nargc<2) //no more availible args,read
		errn=UT_ValErr_noarg;
	else if (destVar->addr) //has valid destination value address
		{
			str=(Sarg->Cargv[(Sarg->Nargc)+1]);
			switch (destVar->type)
			{
				case UTypeID_float:
				{
					float fval;
					//UT_FvalRangeTypedef Cr={0x11,0.8,-0.8};
					errn=UT_val_getf(str,&fval,&(destVar->Cr.FRange));
					if (!errn)
					{
						*(float*)(destVar->addr)=fval;
						errn=UT_ValErr_OK;
					}
					break;
				}
				case UTypeID_S32:
				{
					int val;
					errn=UT_val_gets(str,&val,&(destVar->Cr.SRange));
					if (!errn)
					{
						*(int32_t*)(destVar->addr)=val;
						errn=UT_ValErr_OK;
					}
					break;
				}
				case UTypeID_U32:
				{
					int val;
					errn=UT_val_getu(str,&val,&(destVar->Cr.URange));
					if (!errn)
					{
						*(uint32_t*)(destVar->addr)=val;
						errn=UT_ValErr_OK;
					}
					break;
				}
				case UTypeID_S16:
				{
					int val;
					errn=UT_val_gets(str,&val,&(destVar->Cr.SRange));
					if (!errn)
					{
						*(int16_t*)(destVar->addr)=val;
						errn=UT_ValErr_OK;
					}
					break;
				}
				case UTypeID_U16:
				{
					int val;
					errn=UT_val_getu(str,&val,&(destVar->Cr.URange));
					if (!errn)
					{
						*(uint16_t*)(destVar->addr)=val;
						errn=UT_ValErr_OK;
					}
					break;
				}
				case UTypeID_S8:
				{
					int val;
					errn=UT_val_gets(str,&val,&(destVar->Cr.SRange));
					if (!errn)
					{
						*(int8_t*)(destVar->addr)=val;
						errn=UT_ValErr_OK;
					}
					break;
				}
				case UTypeID_U8:
				{
					int val;
					errn=UT_val_getu(str,&val,&(destVar->Cr.URange));
					if (!errn)
					{
						*(uint8_t*)(destVar->addr)=val;
						errn=UT_ValErr_OK;
					}
					break;
				}
				default:
				{
					errn=UT_ValErr_nodatatype;//unknown data type
					break;
				}
			}
		}
	else 
		errn=UT_ValErr_noaddr;
	//deal with errs
	switch (errn)
	{
		case UT_ValErr_OK: //no error
		{
			printf("Set %s OK.\r\n",destVar->name);
			
		}
		case UT_ValErr_noarg: //just need to read the value
		{
			Valprint(destVar);
			return UT_ValErr_OK;
		}
		case UT_ValErr_noaddr: //not valid address
		{
			printf("invalid address \r\n");
			break;
		}
		case UT_ValErr_nodatatype: //not valid data type definition
		{
			printf("invalid data type \r\n");
			break;
		}
		case UT_ValErr_nan: //not a valid number string
		{
			if (*str!='-' ||*(str+1)!='h')
			{
				printf("argument of %s is not a number\r\n",destVar->name);
				break;
			}
			else
			{
				//CoS_UART_sendFrame(hCoS->hHUART1,destVar->helpwd,strlen(destVar->helpwd));
				return UT_ValErr_OK;
			}
			
		}
		case UT_ValErr_outrange: //not a valid number string
		{
			printf("data out of range when setting %s \r\n",destVar->name);
			break;
		}
		case UT_ValErr_unknown: //not a valid number string
		default:
			{
			printf("unknown error when setting a data \r\n");
			break;
		}
		
	}
	return errn;
	
}




/***********************
Command actions 
************************/



static int UT_Help(const CMD_NodeTypedef* pNode)
{
	CoS_UART_sendFrame(hCoS->hHUART1,(uint8_t*)pNode->hDetails->helpwd,strlen(pNode->hDetails->helpwd));
	return 0;
}
int UT_Help_old(const CMD_NodeTypedef_old* pNode)
{
	CoS_UART_sendFrame(hCoS->hHUART1,(uint8_t*)pNode->helpwd,strlen(pNode->helpwd));
	return 0;
}
u8 st=0;
CMD_NodeTypedef_old CMD_I_a={&st,UTypeID_U8,"MPU9250 write address","I2C write address(mapped with MPU9250)\r\n"};

char UT_CM_IMU(CMDargTypedef* Sarg,char arg)
{
	char d;
	CMD_NodeTypedef currentNode;
	B32DataType data;
	currentNode.hData=0;
	currentNode.hDetails=&Str_I_node;
	currentNode.hLim=0;
	UT_argCheck(Sarg,&currentNode,&data,1);
	d=UT_getopt(Sarg,arg,"r1w2W2a0b0e0u0x1");
	switch (d)
	{
		case 'a':
		{	
			
			return (UT_SB_valSet_old(Sarg,&CMD_I_a));
		}
		case 'e':
		{	
			currentNode.hData=&TiltSense.CC.enable;
			currentNode.hDetails=&Str_I_MPU_enable;
			currentNode.hLim=(UT_ValRangeUniondef*)&Lim_I_MPU_enable;
			return UT_valSet(Sarg,&currentNode);
		}
		case 'u':
		{	
			currentNode.hData=&TiltSense.CC.userdata;
			currentNode.hDetails=&Str_I_MPU_userdata;
			currentNode.hLim=(UT_ValRangeUniondef*)&Lim_I_MPU_userdata;
			return UT_valSet(Sarg,&currentNode);
		}
		case 'x':
		{	
			u8 sta;
			currentNode.hDetails=&Str_I_MPU_reset;
			currentNode.hLim=(UT_ValRangeUniondef*)&Lim_I_MPU_reset;
			sta=UT_argCheck(Sarg,&currentNode,&data,0);
			if (sta==UT_ValErr_help)
				return UT_ValErr_help;
			if (sta==UT_ValErr_OK)
			{
				if (data.data.DU8==0x82)
					tilt_detection_init(&TiltSense);
				if (data.data.DU8==0x54)
					Tilt_detect_reset(&TiltSense);
				else
					return UT_ValErr_outrange;
			}
			return sta;
		}
//		case 'b':
//		{
//			CMD_NodeTypedef Lr_Node={&hCoS->hHCLED->Instance->Rs.freq,(UT_ValRangeUniondef*)&Lr_Cr,&Lr_Str};
//			return UT_valSet(Sarg,&Lr_Node,1);
//		}
		case 'r':
		{	
			uint32_t addr,sta1,sta2;
			uint8_t data;
			UT_UvalRangeTypedef VCaddr={0,0,0};
			sta1=UT_val_getu(Sarg->Cargv[arg+1],&addr,&VCaddr);
			if (sta1 !=UT_ValErr_OK ) //got no valid address
			{
				printf("Unrecognized argument in MPU9250 register address.\r\n");
				return UT_ValErr_unknown;
			}
				//sta1 = SB_sI2C_ReadReg(ADR_MPUI2C,addr&0xff,&data);
				sta1 = CoS_sI2C_ReadReg(hCoS->hHsI2C_MPU,addr&0xff,&data);
			if (sta1 ==CoS_STATE_OK)
				printf("Read MPU9250(0x%02X) register %d (0x%02X): %d (0x%02X)\r\n",hCoS->hHsI2C_MPU->Dev_Addr,addr,addr,data,data);
			else 
				printf("Read MPU9250 register ERROR%d\r\n",STRUCT_ENTRY(hCoS->hHBIO, struct strCoS_BIO_HandleTypeDef, htim3));
			return UT_ValErr_OK;
		}
//		case 'g':
//		{
//			uint8_t data;
//			IICreadBytes(A_MPU9250W,MPU9250_RA_WHO_AM_I,1,&data);
//			printf("Read MPU9250 register 0X%02X: 0X%02X\r\n",MPU9250_RA_WHO_AM_I,data);
//			return UT_ValErr_OK;
//		}
		case 'w':
		case 'W':
		{	
			uint32_t addr,alt,sta1,sta2;
			uint8_t data;
			UT_UvalRangeTypedef VCaddr={0,0,0};
			sta1=UT_val_getu(Sarg->Cargv[arg+1],&addr,&VCaddr);
			sta2=UT_val_getu(Sarg->Cargv[arg+2],&alt,&VCaddr);
			if (sta1 !=UT_ValErr_OK ||sta2 !=UT_ValErr_OK) //got no valid address
			{
				printf("Unrecognized argument in MPU9250 register address.\r\n");
				return UT_ValErr_unknown;
			}
			//sta1 = SB_sI2C_WriteReg(ADR_MPUI2C,addr&0xff,alt);
			sta1 = CoS_sI2C_WriteReg(hCoS->hHsI2C_MPU,addr&0xff,alt);
			if (sta1 ==CoS_STATE_OK)
				printf("Write MPU9250(0x%02X) register  %d (0x%02X): %d (0x%02X)\r\n",hCoS->hHsI2C_MPU->Dev_Addr,addr,addr,alt,alt);
			else 
				printf("Write MPU9250 register ERROR\r\n");
			return UT_ValErr_OK;
		}
		default: //error happen
		{
			printf("argument error\r\n");
			return d;
		}
		
	}
}

char UT_CM_Camera(CMDargTypedef* Sarg,char arg)
{
	char d;
	d=UT_getopt(Sarg,arg,"r1w2p0");
	switch (d)
	{
		case 'r':
		{	
//			uint32_t addr,sta1;
//			uint8_t data;
//			const UT_UvalRangeTypedef VCaddr={0,0xc9,0x11};
//			sta1=UT_val_getu(Sarg->Cargv[arg+1],&addr,&VCaddr);
//			if (sta1 !=UT_ValErr_OK) //got no valid address
//			{
//				printf("Unrecognized argument in Camera register address.\r\n");
//				return UT_ValErr_unknown;
//			}
//			sta1 = SB_SCCB_ReadReg(A_OV7670W,addr&0xff,&data);
//			printf("Read Camera register 0X%02X: 0X%02X\r\n",addr,data);
//			return UT_ValErr_OK;
			
			uint32_t addr,sta1,sta2;
			uint8_t data;
			UT_UvalRangeTypedef VCaddr={0,0xc9,0xff};
			sta1=UT_val_getu(Sarg->Cargv[arg+1],&addr,&VCaddr);
			if (sta1 !=UT_ValErr_OK ) //got no valid address
			{
				printf("Unrecognized argument in camera register address.\r\n");
				return UT_ValErr_unknown;
			}
			sta1 = CoS_sI2C_ReadReg(hCoS->hHsI2C_Camera,addr&0xff,&data);
			if (sta1 ==CoS_STATE_OK)
				printf("Read OV7670(0x%02X) register %d (0x%02X): %d (0x%02X)\r\n",hCoS->hHsI2C_Camera->Dev_Addr,addr,addr,data,data);
			else 
				printf("Read OV7670 register ERROR%d\r\n",STRUCT_ENTRY(hCoS->hHBIO, struct strCoS_BIO_HandleTypeDef, htim3));
			return UT_ValErr_OK;
		}
		case 'w':
		{	
//			uint32_t addr,data,sta1,sta2;
//			//uint8_t data;
//			const UT_UvalRangeTypedef VCaddr={0,0xc9,0x11};
//			const UT_UvalRangeTypedef VCdata={0,0xca,0};
//			sta1=UT_val_getu(Sarg->Cargv[arg+1],&addr,&VCaddr);
//			if (sta1 !=UT_ValErr_OK) //got no valid address
//			{
//				printf("Unrecognized argument in Camera register address.\r\n");
//				return UT_ValErr_unknown;
//			}
//			sta2=UT_val_getu(Sarg->Cargv[arg+2],&data,&VCdata);
//			if (sta2!=UT_ValErr_OK) //got no valid address
//			{
//				printf("Unrecognized argument in Camera register data.\r\n");
//				return UT_ValErr_unknown;
//			}
//			sta2 = SB_SCCB_WriteReg(A_OV7670W,addr&0xff,data&0xff);
//			printf("Write Camera register 0X%02X with 0X%02X\r\n",addr,data);
//			return UT_ValErr_OK;
			
			uint32_t addr,alt,sta1,sta2;
			uint8_t data;
			UT_UvalRangeTypedef VCaddr={0,0xc9,0xff};
			sta1=UT_val_getu(Sarg->Cargv[arg+1],&addr,&VCaddr);
			VCaddr.ctrl=0;
			sta2=UT_val_getu(Sarg->Cargv[arg+2],&alt,&VCaddr);
			if (sta1 !=UT_ValErr_OK ||sta2 !=UT_ValErr_OK) //got no valid address
			{
				printf("Unrecognized argument in OV7670 register address.\r\n");
				return UT_ValErr_unknown;
			}
			sta1 = CoS_sI2C_WriteReg(hCoS->hHsI2C_Camera,addr&0xff,alt);
			if (sta1 ==CoS_STATE_OK)
				printf("Write OV7670(0x%02X) register  %d (0x%02X): %d (0x%02X)\r\n",hCoS->hHsI2C_Camera->Dev_Addr,addr,addr,alt,alt);
			else 
				printf("Write OV7670 register ERROR\r\n");
			return UT_ValErr_OK;
		}
//		case 's':
//		{
//			SB_CAM_Start();
//			return UT_ValErr_OK;
//		}
//		case 'S':
//		{
//			SB_CAM_Stop();
//			return UT_ValErr_OK;
//		}
		//attention
		case 'p':
		{
			SendFrame();
			//u16* p;
			//p=&Image[0][0][0];
			//CoS_UART_sendFrame(hCoS->hHUART1,(uint8_t*)p,Image_Height*Image_Width*2);
			return UT_ValErr_OK;
		}
//		case 'n':
//		{
//			SendFrames(0);
//			return UT_ValErr_OK;
//		}
		default: //error happen
		{
			printf("argument error\r\n");
			return d;
		}
		
	}
}

//Send frame to buffer
void SendFrame(void)
{
	uint8_t i,i1,j;
	uint8_t tmp;
	for(i=0;i<Image_Height;i++)
		{
			i1=i+1;
			putchar('$');
			putchar(i+11);
			putchar(Image_Height+11);
			putchar('[');
			for(j=0;j<Image_Width;j++)
			{
				tmp=Image[0][i][j]&0x00ff;
				if (tmp==10)
					tmp=11;
				putchar(tmp);
			}
			if (i1!=Image_Height)
				printf("z\t\t");
			else
			{
				delay_ms(1);
				printf(" \r\n");
			}
		}
}

//char UT_CM_Camera(CMDargTypedef* Sarg,char arg)
//{
//	char d;
//	d=UT_getopt(Sarg,arg,"R0r1w2W2s0S0p0n0n0");
//	switch (d)
//	{
//		case 'r':
//		{	
//			uint32_t addr,sta1;
//			uint8_t data;
//			const UT_UvalRangeTypedef VCaddr={0,0xc9,0x11};
//			sta1=UT_val_getu(Sarg->Cargv[arg+1],&addr,&VCaddr);
//			if (sta1 !=UT_ValErr_OK) //got no valid address
//			{
//				printf("Unrecognized argument in Camera register address.\r\n");
//				return UT_ValErr_unknown;
//			}
//			sta1 = SB_SCCB_ReadReg(A_OV7670W,addr&0xff,&data);
//			printf("Read Camera register 0X%02X: 0X%02X\r\n",addr,data);
//			return UT_ValErr_OK;
//		}
//		case 'w':
//		case 'W':
//		{	
//			uint32_t addr,data,sta1,sta2;
//			//uint8_t data;
//			const UT_UvalRangeTypedef VCaddr={0,0xc9,0x11};
//			const UT_UvalRangeTypedef VCdata={0,0xca,0};
//			sta1=UT_val_getu(Sarg->Cargv[arg+1],&addr,&VCaddr);
//			if (sta1 !=UT_ValErr_OK) //got no valid address
//			{
//				printf("Unrecognized argument in Camera register address.\r\n");
//				return UT_ValErr_unknown;
//			}
//			sta2=UT_val_getu(Sarg->Cargv[arg+2],&data,&VCdata);
//			if (sta2!=UT_ValErr_OK) //got no valid address
//			{
//				printf("Unrecognized argument in Camera register data.\r\n");
//				return UT_ValErr_unknown;
//			}
//			sta2 = SB_SCCB_WriteReg(A_OV7670W,addr&0xff,data&0xff);
//			printf("Write Camera register 0X%02X with 0X%02X\r\n",addr,data);
//			return UT_ValErr_OK;
//		}
//		case 's':
//		{
//			SB_CAM_Start();
//			return UT_ValErr_OK;
//		}
//		case 'S':
//		{
//			SB_CAM_Stop();
//			return UT_ValErr_OK;
//		}
//		case 'p':
//		{
//			u16* p;
//			p=&Image[0][0][0];
//			CoS_UART_sendFrame(hCoS->hHUART1,(const char*)p,14256);
//			return UT_ValErr_OK;
//		}
//		case 'n':
//		{
//			SendFrames(0);
//			return UT_ValErr_OK;
//		}
//		default: //error happen
//		{
//			printf("argument error\r\n");
//			return d;
//		}
//		
//	}
//}
char UT_CM_SYS(CMDargTypedef* Sarg,char arg)
{
	char d;
	UT_ValErrEnum sta;
	B32DataType data;
	CMD_NodeTypedef currentNode;
	currentNode.hData=0;
	currentNode.hDetails=&Str_Y_node;
	currentNode.hLim=0;
	UT_argCheck(Sarg,&currentNode,&data,1);
	d=UT_getopt(Sarg,arg,"?0M0X1r0T0D0p0n0");
	switch (d)
	{
		case '?':
		{
			currentNode.hDetails=&Str_Y_report;
			if (UT_argCheck(Sarg,&currentNode,&data,1)==UT_ValErr_help)
				return UT_ValErr_help;
			printf("System status:\r\n");
			return UT_ValErr_OK;
		}
		case 'M':
		{	
			currentNode.hDetails=&Str_Y_MCUID;
			if (UT_argCheck(Sarg,&currentNode,&data,1)==UT_ValErr_help)
				return UT_ValErr_help;
			printf("CPUID: %08X-%08X-%08X-%08X\r\n",CoS_MCUID(0),CoS_MCUID(1),CoS_MCUID(2),CoS_MCUID(3));
			return UT_ValErr_OK;
		}
		case 'X':
		{	
			currentNode.hDetails=&Str_Y_RST;
			currentNode.hLim=(UT_ValRangeUniondef*)&Lim_Y_RST;
			sta=UT_argCheck(Sarg,&currentNode,&data,0);
			if (sta==UT_ValErr_help)
				return UT_ValErr_help;
			if (sta==UT_ValErr_OK)
			{
				switch (CoS_Restart(data.data.DU8))
				{
					case 1:
					{
						printf("write 82(0x52) one more time to restart\r\n");
						return UT_ValErr_OK;
					}
					case 0:
					{
						printf("Cancelled. Write 82(0x52) twice to restart\r\n");
						return UT_ValErr_OK;
					}
					default:
					{
						printf("Write 82(0x52) twice to restart\r\n");
						return UT_ValErr_OK;
					}
				}
			}
			return sta;
		}
		case 'r':

		{	
			currentNode.hData=&UT_CMD.C_Readback;
			currentNode.hDetails=&Str_Y_readback;
			currentNode.hLim=(UT_ValRangeUniondef*)&Lim_Y_readback;
			return UT_valSet(Sarg,&currentNode);
		}
		case 'n':
		{	
			currentNode.hData=&hCoS->hHUART1->Instance->Tx_PacketSize;
			currentNode.hDetails=&Str_Y_packet;
			currentNode.hLim=(UT_ValRangeUniondef*)&Lim_Y_packet;
			return UT_valSet(Sarg,&currentNode);
		}
		case 'p':
		{	
			currentNode.hData=&UT_CMD.C_TimeReadback;
			currentNode.hDetails=&Str_Y_time;
			currentNode.hLim=(UT_ValRangeUniondef*)&Lim_Y_time;
			return UT_valSet(Sarg,&currentNode);
		}
		case 'D':
		{	
			currentNode.hData=&UT_CMD.C_do_dataTranfer;
			currentNode.hDetails=&Str_Y_data;
			currentNode.hLim=(UT_ValRangeUniondef*)&Lim_Y_data;
			return UT_valSet(Sarg,&currentNode);
		}
		//attention
		case 'T':	//timer
		{	
			
			uint32_t addr,alt,sta1,sta2;
			uint16_t data[4]={0},n=0;
			UT_UvalRangeTypedef VCaddr={1,3,0xff};
			sta1=UT_val_getu(Sarg->Cargv[arg+1],&addr,&VCaddr);
			VCaddr.ctrl=0;
			sta2=UT_val_getu(Sarg->Cargv[arg+2],&alt,&VCaddr);

			sta1=timed_data_set(addr,alt);
			if (sta1 ==CoS_STATE_OK)
			{
				
				printf("Write OK\r\n");
				return UT_ValErr_OK;
			}
			
			printf("Write register ERROR\r\n");
			return sta1;
		}
		default: //error happen
		{
			printf("argument error\r\n");
			return d;
		}
		
	}
}
//char UT_CM_Sense_LED(CMDargTypedef* Sarg,char arg)
//{
//	uint32_t addr,data,sta1,sta2;
//	const UT_UvalRangeTypedef VCaddr={1,3,0x11};
//	const UT_UvalRangeTypedef VCdata={0,3,0x11};
//	const char* Str[2]={"OFF/0","ON/1"};
//	sta1=UT_val_getu(Sarg->Cargv[arg],&addr,&VCaddr);
//	if (sta1 !=UT_ValErr_OK) //got no valid address
//	{		printf("Unrecognized argument in LED setting.\r\n");
//		return UT_ValErr_unknown;	}
//	if (arg+1==Sarg->Cargc)  //read the LED value if no data given
//	{		printf("LED %d is in state %s.\r\n",addr,Str[SB_LED_Get((SB_LEDType)addr-1)]);
//		return UT_ValErr_OK;	}
//	sta2=UT_val_getu(Sarg->Cargv[arg+1],&data,&VCdata);
//	if (sta2!=UT_ValErr_OK) //got no valid setting number
//	{		printf("Unrecognized argument in LED setting.\r\n");
//		return UT_ValErr_outrange;	}
//	switch (data)
//	{	case 0:
//		case 1:
//		{	SB_LED_Set((SB_LEDType)addr-1,(GPIO_PinState)data);
//			printf("LED %d set OK\r\n",addr);
//			break;		}
//		case 2:
//		{	SB_LED_Toggle((SB_LEDType)addr-1);
//			printf("LED %d toggle OK\r\n", addr);
//			break;	}
//		case 3:
//		{	printf("Flashing LED function to be added(LED setting).\r\n");
//			break;	}
//		default:
//		{	printf("Unrecognized commmand in LED setting.\r\n");
//			return UT_ValErr_outrange;	}
//	}
//	return UT_ValErr_OK;
//}

//char UT_CM_Sense_Key(CMDargTypedef* Sarg,char arg)
//{
//	uint32_t addr,data,sta1,sta2;
//	const UT_UvalRangeTypedef VCaddr={1,2,0x11};
//	const char* Str[2]={"0 (pressed)","1 (unpressed)"};
//	sta1=UT_val_getu(Sarg->Cargv[arg],&addr,&VCaddr);
//	if (sta1 !=UT_ValErr_OK) //got no valid address
//	{
//		printf("Unrecognized argument in key.\r\n");
//		return UT_ValErr_unknown;
//	}
//	printf("Key %d is in state %s.\r\n",addr,Str[SB_Keys_GetState((SB_KeyType)addr-1)]);
//	return UT_ValErr_OK;
//}


//char UT_CM_Sense(CMDargTypedef* Sarg,char arg)
//{
//	char d;
//	d=UT_getopt(Sarg,arg,"L0C0K0H0");
//	switch (d)
//	{
//		case 'L':
//		{	return UT_CM_Sense_LED(Sarg,arg+1);}
//		case 'K':
//		{	return UT_CM_Sense_Key(Sarg,arg+1);}
//		case 'C': 
//		{	return (UT_SB_valSet_old(Sarg,&CMD_r3));}
//		case 'H':
//		{	UT_Help(&CMD_Sense);
//			return 0;
//		}
//		default: //error happen
//		{
//			printf("argument error\r\n");
//			return d;
//		}
//		
//	}
//}

char UT_CM_Motion(CMDargTypedef* Sarg,char arg)
{
	char d;
	uint8_t *pStr;
	uint8_t motionStr_Front[3]={0,0x2f,0};
	uint8_t motionStr_Left_turn[4]={0,0x03,0x49,0};
	uint8_t motionStr_Right_turn[4]={0,0x03,0x69,0};
	uint8_t motionStr_back[3]={0,0x88,0};
	uint8_t motionStr_stop[3]={0,0x08,0};
	d=UT_getopt(Sarg,arg,"w2c1r1H0t3F0R0L0S0B0A0");
	switch (d)
	{
		case 'w':
		{
			uint32_t addr,alt,sta1,sta2;
			uint8_t data[20]={0},n=0;
			UT_UvalRangeTypedef VCaddr={0,0xff,0xf0};
			do{
				sta1=UT_val_getu(Sarg->Cargv[arg+n],(uint32_t*)&data[n++],&VCaddr);
			}while (sta1==UT_ValErr_OK);
			if (n<2) //got no valid address
			{
				printf("not enough arguments.\r\n");
				return UT_ValErr_unknown;
			}
			sta1=WriteBytesWithACKToColiasBasic(hCoS->hHUART4,data,n-1);
			if (sta1 ==CoS_STATE_OK)
				printf("Write OK\r\n");
			else if (sta1 ==CoS_STATE_TIMEOUT)
				printf("Write time out\r\n");
			else
				printf("Write register ERROR\r\n");
			return UT_ValErr_OK;
		}
//		case 'W':
//		{	return UT_CM_Sense_Key(Sarg,arg+1);}
		case 'c':
		case 'r':
		{
			uint32_t addr,alt,sta1,sta2;
			uint8_t paddr[5]={0};
			uint8_t data[20]={0},n=0,k;
			UT_UvalRangeTypedef VCaddr={0,0xff,0xf0};
			do{
				sta1=UT_val_getu(Sarg->Cargv[arg+n],(uint32_t*)&paddr[n++],&VCaddr);
			}while (sta1==UT_ValErr_OK);
			n--;
			if (n<1) //got no valid address
			{
				printf("not enough arguments.\r\n");
				return UT_ValErr_unknown;
			}
			if (n==1)
				paddr[1]=1;
			if (paddr[1]<16)
			{
				if (n>3)
					printf("argument more than two. Take first two as valid.\r\n");
				sta1=ReadBytesWithACKFromColiasBasic(hCoS->hHUART4,paddr[0],data,paddr[1]);
			}
			else
			{
				printf("too many regs to read. keep it less than 15.\r\n");
				return UT_ValErr_outrange;
			}
			if (sta1 ==CoS_STATE_OK)
			{
				for(k=0;k<paddr[1];k++)
				{
					printf("r%d: 0x%02X, 0x%02X\r\n",k,paddr[0]+k,data[k+1]);
				}
			}
			else if (sta1 ==CoS_STATE_TIMEOUT)
				printf("Read time out\r\n");
			else
				printf("Read register ERROR\r\n");
			return UT_ValErr_OK;
		}
		case 't':
		{
			uint32_t addr,alt,sta1,sta2;
			uint8_t data[20]={0},n=0;
			UT_UvalRangeTypedef VCaddr={0,0xff,0xf0};
			do{
				sta1=UT_val_getu(Sarg->Cargv[arg+n],(uint32_t*)&data[n++],&VCaddr);
			}while (sta1==UT_ValErr_OK);
			if (n<3) //got no valid address
			{
				printf("not enough arguments.\r\n");
				return UT_ValErr_unknown;
			}
			sta1=BitSetColiasBasicReg(hCoS->hHUART4,data[0],data[1],data[2]);
			if (sta1 ==CoS_STATE_OK)
				printf("Set OK\r\n");
			else if (sta1 ==CoS_STATE_TIMEOUT)
				printf("Set time out\r\n");
			else
				printf("Set register ERROR\r\n");
			return UT_ValErr_OK;
			
		}
		case 'F':
		{
			//WriteBytesWithACKToColiasBasic(hCoS->hHUART4,data,n-1);
			pStr=motionStr_Front;
			*pStr=0xc0;
			return WriteBytesWithACKToColiasBasic(hCoS->hHUART4,pStr,strlen(pStr));
			//return 0;
		}
		case 'S':
		{
			//WriteBytesWithACKToColiasBasic(hCoS->hHUART4,data,n-1);
			pStr=motionStr_stop;
			*pStr=0xc0;
			return WriteBytesWithACKToColiasBasic(hCoS->hHUART4,pStr,strlen(pStr));
			//return 0;
		}
		case 'B':
		{
			//WriteBytesWithACKToColiasBasic(hCoS->hHUART4,data,n-1);
			pStr=motionStr_back;
			*pStr=0xc0;
			return WriteBytesWithACKToColiasBasic(hCoS->hHUART4,pStr,strlen(pStr));
			//return 0;
		}
		case 'R':
		{
			//WriteBytesWithACKToColiasBasic(hCoS->hHUART4,data,n-1);
			pStr=motionStr_Right_turn;
			*pStr=0xc0;
			return WriteBytesWithACKToColiasBasic(hCoS->hHUART4,pStr,strlen(pStr));
			//return 0;
		}
		case 'L':
		{
			//WriteBytesWithACKToColiasBasic(hCoS->hHUART4,data,n-1);
			pStr=motionStr_Left_turn;
			*pStr=0xc0;
			return WriteBytesWithACKToColiasBasic(hCoS->hHUART4,pStr,strlen(pStr));
			//return 0;
		}
		case 'A':
		{
			//WriteBytesWithACKToColiasBasic(hCoS->hHUART4,data,n-1);
			uint8_t changeDefaultMotion[2]={0x90,0xa0};
			WriteBytesWithACKToColiasBasic(hCoS->hHCos_Motion->Instance->hHUART,changeDefaultMotion,2);
			return 0;
		}
		case 'H':
		{	UT_Help_old(&CMD_ColiasBasic);
			return 0;
		}
		default: //error happen
		{
			printf("argument error\r\n");
			return d;
		}
		
	}
}

/*
For Qinbing Attention
The features related to serial port are realized here
*/
void UT_TermDeal(void)
{
	u8 errcode;
	if (!UT_CMD.F_NotEmpty)
		UT_commandcpy();
	if (UT_CMD.C_mode_CPtP==0)
	{
		if (UT_CMD.F_NotEmpty)
		{
			u32 t;
			char result;
			
			UT_CMD_StringParse(UT_CMD.CM_Buff,&UT_CMD.ArgStruct);
			if (UT_CMD.C_Readback)
			{
				u16 n;
				printf(">>%s\r\n",UT_CMD.CM_Buff);
				for (n=0;n<UT_CMD.ArgStruct.Cargc;n++)
					printf("arg %d: %s\r\n",n,UT_CMD.ArgStruct.Cargv[n]);
			}
			result=UT_getopt(&(UT_CMD.ArgStruct),0,"H0I1M1C1Y1");
			switch (result)
			{
				case 'Y':
				{
					errcode=UT_CM_SYS(&UT_CMD.ArgStruct,1);
					break;
				}
	//			case 'S':
	//			{
	//				errcode=UT_CM_Sense(&UT_CMD.ArgStruct,1);
	//				break;
	//			}
				case 'C':
				{
					errcode=UT_CM_Camera(&UT_CMD.ArgStruct,1);
					break;
				}
				case 'I':
				{
					errcode=UT_CM_IMU(&UT_CMD.ArgStruct,1);
					break;
				}
				case 'M':
				{
					errcode=UT_CM_Motion(&UT_CMD.ArgStruct,1);
					break;
				}
				default:
					printf("Undefined command\r\n");
				case 'H':
					UT_Help_old(&CMD_root);
					break;
				
			}
			//printf(styggh.description->name);
	//		if (errcode!=UT_ErrID_OK);
	//			printf(Str_Errs[errcode]);
			UT_FlushBuff(&UT_CMD);
			if (UT_CMD.C_TimeReadback)
				printf("command time occupation:%8.1f us.\r\n",(float)t/10);
			UT_CMD.F_NotEmpty=0;
		}
	}
	/*********************/
	timed_data_transmit();
	/*********************/
}

void UT_Init(void)
{
	UT_FlushBuff(&UT_CMD);
	UT_CMD.C_Readback=0;
	UT_CMD.C_TimeReadback=0;
	UT_CMD.C_mode_CPtP=0;
	UT_CMD.C_do_dataTranfer=0;
	hCoS->hHUART1->Instance->Tx_PacketSize=100;
	timed_data_init();
}
/****************************/
void timed_data_init(void)
{

	UT_Dtimer.period=0;
	//attention
	/*********************/
	UT_Dtimer.tasks_to_exceute=100;
	/*********************/
}
/****************************/
char timed_data_transmit(void)
{
	static u8 counter=0;
	u8 did=0;
	if (UT_Dtimer.start==0)
	{
		return 0;
	}
	if (UT_Dtimer.tasks_remain==0)
	{
		//printf("timer start\r\n");
		UT_Dtimer.tasks_remain=UT_Dtimer.tasks_to_exceute;
	}

	counter++;
	if (counter>=UT_Dtimer.period)//excecute
	{
		counter=0;
		
		if (UT_Dtimer.tasks_remain)
		{
			//attention
			/****************************/
			printf("%d %.2f\r\n",100-UT_Dtimer.tasks_remain+1, LGMD.Results.LGMD_out[hLGMD.currentDiffImage]);
			//printf("%d \r\n",hCoS->hHBIO->Instance->timerlog[2]-hCoS->hHBIO->Instance->timerlog[0]);
			Key2_virtual=1;
			/****************************/
			//printf("task%d\r\n",UT_Dtimer.tasks_remain);
			did=1;
		}
		UT_Dtimer.tasks_remain--;
		if(UT_Dtimer.tasks_remain>65520)
			UT_Dtimer.tasks_remain=65533;
	}
	if (UT_Dtimer.tasks_remain==0)
	{
		//printf("timer stop\r\n");
		Key2_virtual=0;
		UT_Dtimer.start=0;
	}
	return 0;
}

static char timed_data_set_start(u8 data)
{
	UT_Dtimer.start=!!data;
	if (data)
		UT_Dtimer.tasks_remain=0;
	return CoS_STATE_OK;
}

static char timed_data_set_period(u8 data)
{
	UT_Dtimer.period=data;
	return CoS_STATE_OK;
}

static char timed_data_set_task(u16 data)
{
	UT_Dtimer.tasks_to_exceute=data;
	return CoS_STATE_OK;
}

char timed_data_set(u8 addr,u16 data)
{
	switch (addr)
	{
		case 1:
			return timed_data_set_start(data&0xff);
		case 2:
			return timed_data_set_period(data&0xff);
		case 3:
			return timed_data_set_task(data);
		default:
			return UT_ValErr_outrange;
	}
}


