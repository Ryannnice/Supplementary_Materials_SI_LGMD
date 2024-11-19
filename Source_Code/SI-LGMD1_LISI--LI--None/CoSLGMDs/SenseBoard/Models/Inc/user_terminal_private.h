#ifndef __USERTERPRT_H
#define __USERTERPRT_H


#define UT_ErrID_OK 0x00
#define UT_ErrID_unknown 0x01
#define UT_ErrID_insufArg 0x02
#define UT_ErrID_unrecog 0x03
#define UT_ErrID_outRangeArg 0x04
#define UT_ErrID_UndefinedArg 0x05

CMD_StructTypedef UT_CMD;
float tflo1=0.2;
int tint=-3;
u8 tu8=5;

CMD_NodeTypedef_old CMD_root={0,UTypeID_U8,"root","Available commands: s[Y]stem, [S]ense board, [C]amera, [A]udio, [H]elp"};
CMD_NodeTypedef_old CMD_sYstem={0,UTypeID_U8,"sYstem","Available commands: \r\n [R]:Readback  control.\r\n [T]: Timing function\r\n [H]: Help\r\n"};
CMD_NodeTypedef_old CMD_Y_r={&UT_CMD.C_Readback,UTypeID_U8,"Readback","Readback command string or not\r\n"};
CMD_NodeTypedef_old CMD_Y_t={&UT_CMD.C_TimeReadback,UTypeID_U8,"Time count","Readback the time count of each command\r\n"};

CMD_NodeTypedef_old CMD_Sense={0,UTypeID_U8,"Sense board","Available commands: \r\n [L]:LED control.\r\n [C]: Multi color LED\r\n [K]: Keys\r\n [H]: Help\r\n"};
CMD_NodeTypedef_old CMD_S_L={0,UTypeID_U8,"Sense board","Available commands: \r\n [L]:LED control.\r\n [C]: Multi color LED\r\n [K]: Keys\r\n [H]: Help\r\n"};
CMD_NodeTypedef_old CMD_ColiasBasic={0,UTypeID_U8,"Colias Basic board control",\
"Available commands: \r\n [w]:write one register.\r\n [W]: write multiple registers\r\n [r]: read one register\r\n [R]: Read multiple registers\r\n"};



CMD_NodeTypedef_old CMD_r={.addr=&tflo1,UTypeID_float,"name1","test float number",.Cr.FRange={-4,5.2,0x11}};
CMD_NodeTypedef_old CMD_r2={&tint,UTypeID_S16,"name2","test int number",.Cr.SRange={-80,70,0x11}};
CMD_NodeTypedef_old CMD_r3={&tu8,UTypeID_U8,"name3","test u8 number.",.Cr.URange={0,9,0x11}};
CMD_NodeTypedef_old CMD_r4={.addr=&tu8,UTypeID_S8,"name3","test u8 number.",.Cr.SRange={0,5,0}};

char  Str_default123[]="NULL";

 char Str_default[]="NULL";
 char* Str_pool_default[]={"null","null2"};

 char*  Str_sYsR_name="Uart Command Readback";
 char*  Str_sYsT_name="Uart Command Time Count";

 char*  Str_SET_RST[]={"RESET","SET"};


 char* const Str_Errs[]={
	"Command OK.",
	"Error: not enough argument number. Expected at least %d (%d given) ",
	"Error: unrecognized argument ",
	"Error: argument <%d> out of range. Expected %d to %d ",
	"Error: argument out of range. Expected %f to %f ",
	"Error: unknown error"};

 char*const Str_echo="ECHO>>";

 char*const Str_err="Error occored dealing with commands.\r\n";

 char*const Str_unrecog="Unrecognized input command.\r\n";

 char*const Str_insuf_arg="Insufficient input Arguments. Required at least %d but got %d.\r\n";

 char*const Str_toomany_arg="Too many input Arguments. Required at most %d but got %d.\r\n";

 char*const Str_outrangeU="Argument out of range. Expected %d to %d .\r\n";

 char*const Str_outrangeF="Argument out of range. Expected %f to %f .\r\n";

 char*const Str_tobeadded="Function not finished, to be added later.\r\n";

 char* const Str_help1=
"\r\n****Help****\r\n\
The Colias Sense board accept command start with letter:\r\n\
Y  for sYstem commands.\r\n\
S  for Senseboard basic commands.\r\n\
C  for Camera commands.\r\n\
A  for Audio commands.\r\n\
M  for Motion commands.\r\n\
H  for help commands.\r\n";


 char* const Str_help_sb=
"\r\n****Help****Sense Board****\r\n\
Colias Sense Basic device commands\r\n\
LED, Color LED, Keys and USART.\r\n\
SL [LED_ID](1-3) [cmd(optional)]\r\n\
\t [cmd] = null: read [LED_ID]=1-3 status;\r\n \
\t [cmd] = 0/1: [LED_ID] OFF/ON;\r\n\
\t [cmd] = 2: [LED_ID] toggle;\r\n\
\t [cmd] = 3: [LED_ID] flash;\r\n\
SC [Channel_ID](0-3) [cmd0(optional)] [cmd1(optional)] [cmd2-4(optional)]\r\n\
\t [cmd0] = null: read Color LED settings.\r\n\
\t [cmd0] = 0: mode settings.\r\n\
\t\t [cmd1] = null: read mode settings\r\n\
\t\t [cmd1] = 0: Color LED OFF\r\n\
\t\t [cmd1] = 1: Set a random color\r\n\
\t\t [cmd1] = 2: Set to breath effect (by default)\r\n\
\t\t [cmd1] = 3: Set to RED\r\n\
\t\t [cmd1] = 4: Set to GREEN\r\n\
\t\t [cmd1] = 5: Set to BLUE\r\n\
\t\t [cmd1] = 6: Set to manual and freeze\r\n\
\t\t manual color according to [cmd2-4] (0-65535) in R,G and B respectively.\r\n\
\t [cmd0] = 1/2/3: Breath generator R/G/B settings.\r\n\
\t\t [cmd1] = null: read selected settings\r\n\
\t\t [cmd1] = 0: reset to default settings.\r\n\
\t\t3:channel Blue setting;\r\n\
**Help end.\r\n";


CMD_NodeNameStr_Typedef const Yr_Str_Readback={UTypeID_U32,"C_Readback","Readback command or not"};
UT_ValRangeUniondef     const Yr_Cr={.URange={0,9,0x11}};
CMD_NodeTypedef Yr_Node={&UT_CMD.C_Readback,(UT_ValRangeUniondef*)&Yr_Cr,&Yr_Str_Readback};

CMD_NodeNameStr_Typedef const Lr_Str={UTypeID_float,"Red Freq","Red channel frequency"};
UT_ValRangeUniondef     const Lr_Cr={.FRange={0,100,0x11}};
//CMD_NodeTypedef Lr_Node={&hCoS-,(UT_ValRangeUniondef*)&Yr_Cr,&Yr_Str_Readback};

CMD_NodeNameStr_Typedef const Str_Y_node={0,"system config","The system configs. Avlailable commands:\r\n\
\t[?]: report\r\n\t[M]: MCU_ID\r\n\t[X]: reset\r\n\t[r]: read back control\r\n\t[H]:help\r\n"};
CMD_NodeNameStr_Typedef const Str_Y_report={0,"system report","The brief status of system\r\n"};
CMD_NodeNameStr_Typedef const Str_Y_MCUID={0,"MCUID","The ID of MCU chip\r\n"};
CMD_NodeNameStr_Typedef const Str_Y_RST={UTypeID_U8,"RESET","Reset command. enabled by Setting it twice with arg 82. Other numbers to cancel.\r\n"};
UT_ValRangeUniondef     const Lim_Y_RST={.URange={0,0xff,0x00}};
CMD_NodeNameStr_Typedef const Str_Y_readback={UTypeID_U8,"UT readback","To enable read back input command.\r\n"};
UT_ValRangeUniondef     const Lim_Y_readback={.URange={0,1,0x11}};
CMD_NodeNameStr_Typedef const Str_Y_time={UTypeID_U8,"UT time calc","To calculate the time used by the command.\r\n"};
UT_ValRangeUniondef     const Lim_Y_time={.URange={0,1,0x11}};
CMD_NodeNameStr_Typedef const Str_Y_CPtP={UTypeID_U8,"UT time calc","To calculate the time used by the command.\r\n"};
UT_ValRangeUniondef     const Lim_Y_CPtP={.URange={0,1,0x11}};

CMD_NodeNameStr_Typedef const Str_Y_data={UTypeID_U8,"UT data transfer","To transfer data periodically back by serial port.\r\n"};
UT_ValRangeUniondef     const Lim_Y_data={.URange={0,1,0x11}};

CMD_NodeNameStr_Typedef const Str_Y_T_statr={UTypeID_U8,"UT data transfer","To transfer data periodically back by serial port.\r\n"};
UT_ValRangeUniondef     const Lim_Y_T_start={.URange={0,1,0x11}};

CMD_NodeNameStr_Typedef const Str_Y_T_p={UTypeID_U8,"UT data transfer","To transfer data periodically back by serial port.\r\n"};
UT_ValRangeUniondef     const Lim_Y_T_p={.URange={0,1,0x11}};

CMD_NodeNameStr_Typedef const Str_Y_T_t={UTypeID_U8,"UT data transfer","To transfer data periodically back by serial port.\r\n"};
UT_ValRangeUniondef     const Lim_Y_T_t={.URange={0,1,0x11}};
CMD_NodeNameStr_Typedef const Str_Y_packet={UTypeID_U16,"UART1 transfer packet size","Bytes transferred in each uart packet. 100 is good for normal transmission, 1k for large data\r\n"};
UT_ValRangeUniondef     const Lim_Y_packet={.URange={1,2000,0x11}};

CMD_NodeNameStr_Typedef const Str_Y_LCFEnable={UTypeID_U8,"LCF Model Enable","Enable/Disable auto run\r\n"};
UT_ValRangeUniondef     const Lim_Y_LCFEnable={.URange={0,1,0x11}};

CMD_NodeNameStr_Typedef const Str_I_node={0,"MPU9250 access","The Access of MPU9250 sensor. \r\n"};

CMD_NodeNameStr_Typedef const Str_I_MPU_enable={UTypeID_U8,"MPU9250 Enable","Enable/Disable MPU9250 auto updata\r\n"};
UT_ValRangeUniondef     const Lim_I_MPU_enable={.URange={0,1,0x11}};

CMD_NodeNameStr_Typedef const Str_I_MPU_userdata={UTypeID_U8,"MPU9250 readback","Enable/Disable MPU9250 read data back\r\n"};
UT_ValRangeUniondef     const Lim_I_MPU_userdata={.URange={0,1,0x11}};

CMD_NodeNameStr_Typedef const Str_I_MPU_reset={UTypeID_U8,"MPU9250 RESET","enabled by Setting it twice with arg 82. Other numbers to cancel.\r\n"};
UT_ValRangeUniondef     const Lim_I_MPU_reset={.URange={0,1,0}};
/*
Y commands
?0M0X1R0T0
MCUID

Readback //read back command or not
reset

S commands
S? //brief report
SL //LED
Sk //key
SC//multi color led

M commands
M? //brief report
Mr
Mw
MF
MR
ML
MB
MS
MX
MY
MZ
MD //default motion setup
	MD 
MC //custom motion
MQ //motion queue

B commands //basic board
B? //brief report of basic board
Br //read a basic board register
Bw //write a basic board register
BH //help

C commands //camera
C?
Cr //read a camera register
Cw //write a camera register
Cc //camera control (stop/ start)
CC //capture image


I commands
I?
Ir//read
Iw//write
Id//reset
IMr//Mag write
IMw//Mag read
IMd//Mag reset
IH
E commands
E?
ED //GPIOD control
ES //SD control
EH

argcheck 
if help

if 


*/

/**********************
function declerations
**********************/

//void UT_commandcpyUSART1(CMD_StructTypedef* hUT,SB_UartStructTypedef* hSB_UARTn);


#endif