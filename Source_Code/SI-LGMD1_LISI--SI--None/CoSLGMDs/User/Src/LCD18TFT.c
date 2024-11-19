//#include "coliasSense_BIO.h"
//#include "coliasSense_board.h"
#include "LCD18TFT.h"
#include "delay.h"
#include "fontLCD.h"
static void TFTLCD_IOInit(void);
static void LCD_Writ_Bus(char data);	
static void TFT_LCD_REG_Init(void);

#define SSCKD 2
u16 BACK_COLOR, POINT_COLOR;   //背景色，画笔色	

u16 CBuff[LCD_H*LCD_W]={0};

const char welcome_psr[] = "Colias Vision   ";
const char welcome_psr2[]="            V4.5";
const char welcome_psr3[]="  Visomorphic   ";
const char welcome_psr4[]="      25/11/2018";

u16 COLORTABLE[]={
 WHITE,         	 
 BLACK,        	 
 BLUE,        	 
 BRED,           
 GRED,			 
 GBLUE	,		 
 RED,        	 
 MAGENTA,     	 
 GREEN,     	 
 CYAN,     	 
 YELLOW,    	 
 BROWN,	 
 BRRED,	 
 GRAY,		 
 DARKBLUE,   	 
 LIGHTBLUE,    
 GRAYBLUE,   
 LIGHTGREEN ,   
 LGRAY,	 
 LGRAYBLUE,   
 LBBLUE          
};
static void TFTLCD_IOInit(void)
{
	GPIO_InitTypeDef PInit;

	PInit.Pin=TFT_CS_Pin;
	PInit.Mode=GPIO_MODE_OUTPUT_PP;
	PInit.Pull=GPIO_PULLUP;
	PInit.Speed=GPIO_SPEED_FREQ_LOW;
	HAL_GPIO_Init(TFT_CS_GPIO_Port, &PInit);
	
	PInit.Pin=TFT_RST_Pin;
	HAL_GPIO_Init(TFT_RST_GPIO_Port, &PInit);
	
	PInit.Pin=TFT_RST_Pin;
	HAL_GPIO_Init(TFT_RST_GPIO_Port, &PInit);
	
	PInit.Pin=TFT_A0_Pin;
	HAL_GPIO_Init(TFT_A0_GPIO_Port, &PInit);
	
	PInit.Pin=TFT_SDA_Pin;
	HAL_GPIO_Init(TFT_SDA_GPIO_Port, &PInit);
	
	PInit.Pin=TFT_SCK_Pin;
	HAL_GPIO_Init(TFT_SCK_GPIO_Port, &PInit);
	
	PInit.Pin=TFT_LED_Pin;
	HAL_GPIO_Init(TFT_LED_GPIO_Port, &PInit);
	
	TFT_CS_ON;
	TFT_LED_ON;
	TFT_RST_OFF;
	
	
}

void TFT_LCD_Init(void)
{
	TFTLCD_IOInit();
	
	TFT_CS_ON;  //unselect
	TFT_LED_ON; //initiallt start backlight
	TFT_RST_OFF; // rest chip
	
	
	TFT_LED_ON;
	delay_ms(100);
	TFT_LED_OFF;
	delay_ms(100);
	TFT_LED_ON;
	
	TFT_RST_OFF;
	delay_ms(100);
	TFT_RST_ON;
	delay_ms(100);
	TFT_CS_OFF;  //select on
	
	TFT_LCD_REG_Init();
	
	LCD_Clear(DARKBLUE);
	BACK_COLOR=DARKBLUE;
	POINT_COLOR=WHITE; 
	//LCD_DISP_PIC((u8*)picc1,160*128*2-1);
	LCD_ShowString(10,10,(const u8*)welcome_psr);
	LCD_ShowString(10,35,(const u8*)welcome_psr2);
	POINT_COLOR=LIGHTGREEN; 
	LCD_ShowString(10,60,(const u8*)welcome_psr3);
	LCD_ShowString(10,85,(const u8*)welcome_psr4);
}

static void TFT_LCD_REG_Init(void)
{
//	while (1)
//	{
//		TFT_CS_ON;  //unselect
//		LCD_WR_REG(0x05); 
//		delay_short(500);
//		LCD_WR_REG(0x0A); 
//		delay_short(1500);
//		LCD_WR_DATA8(0x05);
//		delay_short(500);
//		LCD_WR_DATA8(0x0A);
//		delay_short(500);
//		TFT_CS_OFF;  
//	
//	}
	LCD_WR_REG(0x11); //Sleep out
	delay_ms(120); //Delay 120ms
	//------------------------------------ST7735S Frame Rate-----------------------------------------//
	LCD_WR_REG(0xB1);
	LCD_WR_DATA8(0x05);
	LCD_WR_DATA8(0x3C);
	LCD_WR_DATA8(0x3C);
	LCD_WR_REG(0xB2);
	LCD_WR_DATA8(0x05);
	LCD_WR_DATA8(0x3C);
	LCD_WR_DATA8(0x3C);
	LCD_WR_REG(0xB3);
	LCD_WR_DATA8(0x05);
	LCD_WR_DATA8(0x3C);
	LCD_WR_DATA8(0x3C);
	LCD_WR_DATA8(0x05);
	LCD_WR_DATA8(0x3C);
	LCD_WR_DATA8(0x3C);
	//------------------------------------End ST7735S Frame Rate-----------------------------------------//
	LCD_WR_REG(0xB4); //Dot inversion
	LCD_WR_DATA8(0x03);
	LCD_WR_REG(0xC0);
	LCD_WR_DATA8(0x28);
	LCD_WR_DATA8(0x08);
	LCD_WR_DATA8(0x04);
	LCD_WR_REG(0xC1);
	LCD_WR_DATA8(0XC0);
	LCD_WR_REG(0xC2);
	LCD_WR_DATA8(0x0D);
	LCD_WR_DATA8(0x00);
	LCD_WR_REG(0xC3);
	LCD_WR_DATA8(0x8D);
	LCD_WR_DATA8(0x2A);
	LCD_WR_REG(0xC4);
	LCD_WR_DATA8(0x8D);
	LCD_WR_DATA8(0xEE);
	//---------------------------------End ST7735S Power Sequence-------------------------------------//
	LCD_WR_REG(0xC5); //VCOM
	LCD_WR_DATA8(0x1A);
	LCD_WR_REG(0x36); //MX, MY, RGB mode
	LCD_WR_DATA8(0xA0);
	//------------------------------------ST7735S Gamma Sequence-----------------------------------------//
	LCD_WR_REG(0xE0);
	LCD_WR_DATA8(0x04);
	LCD_WR_DATA8(0x22);
	LCD_WR_DATA8(0x07);
	LCD_WR_DATA8(0x0A);
	LCD_WR_DATA8(0x2E);
	LCD_WR_DATA8(0x30);
	LCD_WR_DATA8(0x25);
	LCD_WR_DATA8(0x2A);
	LCD_WR_DATA8(0x28);
	LCD_WR_DATA8(0x26);
	LCD_WR_DATA8(0x2E);
	LCD_WR_DATA8(0x3A);
	LCD_WR_DATA8(0x00);
	LCD_WR_DATA8(0x01);
	LCD_WR_DATA8(0x03);
	LCD_WR_DATA8(0x13);
	LCD_WR_REG(0xE1);
	LCD_WR_DATA8(0x04);
	LCD_WR_DATA8(0x16);
	LCD_WR_DATA8(0x06);
	LCD_WR_DATA8(0x0D);
	LCD_WR_DATA8(0x2D);
	LCD_WR_DATA8(0x26);
	LCD_WR_DATA8(0x23);
	LCD_WR_DATA8(0x27);
	LCD_WR_DATA8(0x27);
	LCD_WR_DATA8(0x25);
	LCD_WR_DATA8(0x2D);
	LCD_WR_DATA8(0x3B);
	LCD_WR_DATA8(0x00);
	LCD_WR_DATA8(0x01);
	LCD_WR_DATA8(0x04);
	LCD_WR_DATA8(0x13);
	//------------------------------------End ST7735S Gamma Sequence-----------------------------------------//
	LCD_WR_REG(0x3A); //65k mode
	LCD_WR_DATA8(0x05);
	LCD_WR_REG(0x29); //Display on
}



static void LCD_Writ_Bus(char data)   
{	
	volatile char n,ddata=data;
	for (n=0;n<8;n++)
	{
		if ((ddata & 0x80)==0)
			TFT_SDA_OFF;
		else
			TFT_SDA_ON;
		ddata=ddata<<1;
		//LCDELAY;
		TFT_SCK_OFF;
		//LCDELAY;
		TFT_SCK_ON;
		//LCDELAY;
	}
	
} 

void LCD_WR_DATA8(char da) 
{
    TFT_A0_ON;
	//LCDELAY;
	LCD_Writ_Bus(da);
	//LCDELAY;
}  
 void LCD_WR_DATA(int da)
{
    TFT_A0_ON;
	//LCDELAY;
	LCD_Writ_Bus(da>>8);
	//LCDELAY;
	LCD_Writ_Bus(da);
	//LCDELAY;
}	  
void LCD_WR_REG(char da)	 
{
    TFT_A0_OFF;
	//LCDELAY;
	LCD_Writ_Bus(da);
	//LCDELAY;
}
 void LCD_WR_REG_DATA(int reg,int da)
{
    LCD_WR_REG(reg);
	LCD_WR_DATA(da);
}
void Address_set(u16 x1,u16 y1,u16 x2,u16 y2)
{  
   LCD_WR_REG(0x2a);
   LCD_WR_DATA8(x1>>8);
   LCD_WR_DATA8(x1);
   LCD_WR_DATA8(x2>>8);
   LCD_WR_DATA8(x2);
  
   LCD_WR_REG(0x2b);
   LCD_WR_DATA8(y1>>8);
   LCD_WR_DATA8(y1);
   LCD_WR_DATA8(y2>>8);
   LCD_WR_DATA8(y2);

   LCD_WR_REG(0x2C);

}


//清屏函数
//Color:要清屏的填充色
void LCD_Clear(u16 Color)
{
	u8 VH,VL;
	u16 i,j;
	VH=Color>>8;
	VL=Color;	
	Address_set(0,0,LCD_W-1,LCD_H-1);
    for(i=0;i<LCD_W;i++)
	 {
	  for (j=0;j<LCD_H;j++)
	   	{
        	 LCD_WR_DATA8(VH);
			 LCD_WR_DATA8(VL);	

	    }

	  }
}
//画点
//POINT_COLOR:此点的颜色
void LCD_DrawPoint(u16 x,u16 y)
{
	Address_set(x,y,x,y);//设置光标位置 
	LCD_WR_DATA(POINT_COLOR); 	    
} 	 
//画一个大点
//POINT_COLOR:此点的颜色
void LCD_DrawPoint_big(u16 x,u16 y)
{
	LCD_Fill(x-1,y-1,x+1,y+1,POINT_COLOR);
} 
//在指定区域内填充指定颜色
//区域大小:
//  (xend-xsta)*(yend-ysta)
void LCD_Fill(u16 xsta,u16 ysta,u16 xend,u16 yend,u16 color)
{          
	u16 i,j; 
	Address_set(xsta,ysta,xend,yend);      //设置光标位置 
	for(i=ysta;i<=yend;i++)
	{													   	 	
		for(j=xsta;j<=xend;j++)LCD_WR_DATA(color);//设置光标位置 	    
	} 					  	    
}  
//画线
//x1,y1:起点坐标
//x2,y2:终点坐标  
void LCD_DrawLine(u16 x1, u16 y1, u16 x2, u16 y2)
{
	u16 t; 
	int xerr=0,yerr=0,delta_x,delta_y,distance; 
	int incx,incy,uRow,uCol; 

	delta_x=x2-x1; //计算坐标增量 
	delta_y=y2-y1; 
	uRow=x1; 
	uCol=y1; 
	if(delta_x>0)incx=1; //设置单步方向 
	else if(delta_x==0)incx=0;//垂直线 
	else {incx=-1;delta_x=-delta_x;} 
	if(delta_y>0)incy=1; 
	else if(delta_y==0)incy=0;//水平线 
	else{incy=-1;delta_y=-delta_y;} 
	if( delta_x>delta_y)distance=delta_x; //选取基本增量坐标轴 
	else distance=delta_y; 
	for(t=0;t<=distance+1;t++ )//画线输出 
	{  
		LCD_DrawPoint(uRow,uCol);//画点 
		xerr+=delta_x ; 
		yerr+=delta_y ; 
		if(xerr>distance) 
		{ 
			xerr-=distance; 
			uRow+=incx; 
		} 
		if(yerr>distance) 
		{ 
			yerr-=distance; 
			uCol+=incy; 
		} 
	}  
}    
//画矩形
void LCD_DrawRectangle(u16 x1, u16 y1, u16 x2, u16 y2)
{
	LCD_DrawLine(x1,y1,x2,y1);
	LCD_DrawLine(x1,y1,x1,y2);
	LCD_DrawLine(x1,y2,x2,y2);
	LCD_DrawLine(x2,y1,x2,y2);
}
//在指定位置画一个指定大小的圆
//(x,y):中心点
//r    :半径
void Draw_Circle(u16 x0,u16 y0,u8 r)
{
	int a,b;
	int di;
	a=0;b=r;	  
	di=3-(r<<1);             //判断下个点位置的标志
	while(a<=b)
	{
		LCD_DrawPoint(x0-b,y0-a);             //3           
		LCD_DrawPoint(x0+b,y0-a);             //0           
		LCD_DrawPoint(x0-a,y0+b);             //1       
		LCD_DrawPoint(x0-b,y0-a);             //7           
		LCD_DrawPoint(x0-a,y0-b);             //2             
		LCD_DrawPoint(x0+b,y0+a);             //4               
		LCD_DrawPoint(x0+a,y0-b);             //5
		LCD_DrawPoint(x0+a,y0+b);             //6 
		LCD_DrawPoint(x0-b,y0+a);             
		a++;
		//使用Bresenham算法画圆     
		if(di<0)di +=4*a+6;	  
		else
		{
			di+=10+4*(a-b);   
			b--;
		} 
		LCD_DrawPoint(x0+a,y0+b);
	}
} 
//在指定位置显示一个字符

//num:要显示的字符:" "--->"~"
//mode:叠加方式(1)还是非叠加方式(0)
//在指定位置显示一个字符

//num:要显示的字符:" "--->"~"

//mode:叠加方式(1)还是非叠加方式(0)
void LCD_ShowChar(u16 x,u16 y,u8 num,u8 mode)
{
    u8 temp;
    u8 pos,t;
	u16 x0=x;
	u16 colortemp=POINT_COLOR;      
    if(x>LCD_W-16||y>LCD_H-16)return;	    
	//设置窗口		   
	num=num-' ';//得到偏移后的值
	Address_set(x,y,x+8-1,y+16-1);      //设置光标位置 
	if(!mode) //非叠加方式
	{
		for(pos=0;pos<16;pos++)
		{ 
			temp=asc2_1608[(u16)num*16+pos];		 //调用1608字体
			for(t=0;t<8;t++)
		    {                 
		        if(temp&0x01)POINT_COLOR=colortemp;
				else POINT_COLOR=BACK_COLOR;
				LCD_WR_DATA(POINT_COLOR);	
				temp>>=1; 
				x++;
		    }
			x=x0;
			y++;
		}	
	}else//叠加方式
	{
		for(pos=0;pos<16;pos++)
		{
		    temp=asc2_1608[(u16)num*16+pos];		 //调用1608字体
			for(t=0;t<8;t++)
		    {                 
		        if(temp&0x01)LCD_DrawPoint(x+t,y+pos);//画一个点     
		        temp>>=1; 
		    }
		}
	}
	POINT_COLOR=colortemp;	    	   	 	  
}   
//m^n函数
u32 mypow(u8 m,u8 n)
{
	u32 result=1;	 
	while(n--)result*=m;    
	return result;
}			 
//显示2个数字
//x,y :起点坐标	 
//len :数字的位数
//color:颜色
//num:数值(0~4294967295);	
void LCD_ShowNum(u16 x,u16 y,u32 num,u8 len)
{         	
	u8 t,temp;
	u8 enshow=0;
	num=(u16)num;
	for(t=0;t<len;t++)
	{
		temp=(num/mypow(10,len-t-1))%10;
		if(enshow==0&&t<(len-1))
		{
			if(temp==0)
			{
				LCD_ShowChar(x+8*t,y,' ',0);
				continue;
			}else enshow=1; 
		 	 
		}
	 	LCD_ShowChar(x+8*t,y,temp+48,0); 
	}
} 

//显示2个数字
//x,y:起点坐标
//num:数值(0~99);	 
void LCD_Show2Num(u16 x,u16 y,u16 num,u8 len)
{         	
	u8 t,temp;						   
	for(t=0;t<len;t++)
	{
		temp=(num/mypow(10,len-t-1))%10;
	 	LCD_ShowChar(x+8*t,y,temp+'0',0); 
	}
} 
//显示字符串
//x,y:起点坐标  
//*p:字符串起始地址
//用16字体
void LCD_ShowString(u16 x,u16 y,const u8 *p)
{         
    while(*p!='\0')
    {       
        if(x>LCD_W-16){x=0;y+=16;}
        if(y>LCD_H-16){y=x=0;}
        LCD_ShowChar(x,y,*p,0);
        x+=8;
        p++;
    }  
}

void LCD_DISP_PIC(u8* pstart,u16 size_char)
{
	u16 n;
	Address_set(0,0,LCD_W-1,LCD_H-1);
	for (n=0;n<size_char;n++)
	LCD_WR_DATA8(*(pstart+n));
}


void LCD_imgshow_cam_raw_data(u16 x,u16 y, u16 width, u16 height, u16* pstart, u8 mode)
{
//In this function, the input image is YUV mode used in most cases.
//This function transfers YUV data into RGB565 first, then draw it on LCD.
	
}


u16 pixel_YUV2RGB565(u8 Y,u8 U,u8 V)
{
	s16 R,G,B,u,v;
	u16 P;
	u=U-128;
	v=V-128;
	R = Y + v + ((103*v)>>8);
	G = Y - ((u*88)>>8) - ((v*183)>>8);
	B = Y + u + ((u*198)>>8);
	if (R<0)
		R=0;
	else if (R>255)
		R=255;
	
	if (G<0)
		G=0;
	else if (G>255)
		G=255;
	if (B<0)
		B=0;
	if (B>255)
		B=255;
	P= ((R & 0xF8)<<8) | ((G & 0xFC) <<3 ) | ((B & 0xF8 ) >> 3);
	return P;
}

u32 img_YUV2RGB565(u32 d1, u8 mode)
{
	u8 Y1,U,V,Y2;
	s16 R1,G1,B1,R2,G2,B2;
	u16 P1,P2;
	u8* pd;
	u32* pd1;
	mode=(mode&0x0c)>>2;
	pd1=&d1; //treat d1 as four u8 data array
	pd=(u8*)pd1;
	if (mode==0)
	{
		Y1 =*pd;
		U  =*(pd+1);
		Y2 =*(pd+2);
		V  =*(pd+3);
	}
	else if (mode==1)
	{
		Y1 =*pd;
		V  =*(pd+1);
		Y2 =*(pd+2);
		U  =*(pd+3);
	}
	else if (mode ==2)
	{
		U  =*pd;
		Y1 =*(pd+1);
		V  =*(pd+2);
		Y2 =*(pd+3);
	}
	else if (mode ==3)
	{
		V  =*pd;
		Y1=*(pd+1);
		U =*(pd+2);
		Y2=*(pd+3);
	}
	
	P1= pixel_YUV2RGB565(Y1,U,V);
	P2= pixel_YUV2RGB565(Y2,U,V);
	
	return (P1<<16)|P2;
}




void LCD_imgshow(u16 x,u16 y, u16 width, u16 height, u16* pstart, u8 mode)
{
	/*mode: bit[0] 	=0 direct map,treat source data as 1d array
					=1 grid map,treat source data as 2d array
	bit [1] =0: source is RGB565 
			=1: source is YUV422
 	bit[3:2] YUV swap pattern
		00 YUYV
		01 YVYU
		10 UYVY
		11 VYUY
		
	*/
	u16 n;
	u16 data;
	if ((mode&0x01)==0) //direct map mode, treat source data as 1d array
	{
		if (width>LCD_W ||height>LCD_H)
			return;
		Address_set(x,y,x+width-1,y+height-1);
		for (n=0;n<width*height;n++)
		{
			if (mode &0x02==0)// RGB565 data
			{
				data=*(pstart+n);
				LCD_WR_DATA(data);
			}
			else
			{
				u32 Tdata,Rdata;
				Tdata=*(u32*)(pstart+n);
				Rdata=img_YUV2RGB565(Tdata,mode);
				LCD_WR_DATA((Rdata>>16)&0xffff);
				LCD_WR_DATA(Rdata&0xffff);
				n=n+1;
			}
		}
	}
	else if (mode&0x01) 
	{
		u16 i,j;
		if (width>LCD_W )
			width = LCD_W;
		if (height>LCD_H)
			height=LCD_H;
		Address_set(x,y,x+width-1,y+height-1);
		for (i=0;i<height;i++)
			for (j=0;j<width;j++)
			{
				data=*(pstart+i*width+j);
				LCD_WR_DATA(data);

			}
	}
}


/*
void LCD_imgshow(u16 x,u16 y, u16 width, u16 height, u8* pstart)
{
	u16 n,i;
	u8* p;
	if (width>LCD_W ||height>LCD_H)
		return;
	Address_set(x,y,x+width-1,y+height-1);
	for (n=0;n<height;n++)
		for (i=0;i<width*2;i++)
		{
			p=pstart+height*n+i;
			LCD_WR_DATA8(*(p));
			//LCD_WR_DATA8(*(p+1));
		}
}
*/





void YUV2RGB(ColordataTypedef* CYUV,ColordataTypedef* CRGB)
{
	float r,g,b,y,u,v;
	y=CYUV->C1;
	u=CYUV->C2;
	v=CYUV->C3;
	r=y + 1.4075* (v-0.5);
	g=y - 0.3455* (u-0.5)-(0.7169*(v-0.5));
	b=y + 1.7790* (u-0.5);
	if (r<0)
		r=0;
	if (g<0)
		g=0;
	if (b<0)
		b=0;
	CRGB->C1=r;
	CRGB->C2=g;
	CRGB->C3=b;
}

void RGB2HSV(ColordataTypedef* CRGB,ColordataTypedef* CHSV)
{
	float cmax,cmin,delta,r,g,b,h,s,v;
	int hsel;
	r=CRGB->C1;
	g=CRGB->C2;
	b=CRGB->C3;
	
	cmin = nmin(nmin(r, g), b);
	cmax = nmax(nmax(r, g), b);
	
	delta=cmax-cmin;
	if (delta<=0.01)
	{
		h=0;
		s=0;
	}
	else
	{
		s=delta/cmax;
		if (r>=cmax)
				h=(g-b)/delta;
		else if (g>=cmax)
				h=2+((b-r)/delta);
		else
			h=4+((r-g)/delta);
		h*=60;
		if (h<0)
			h+=360;
	}
	CHSV->C1=h;
	CHSV->C2=s;
	CHSV->C3=cmax;
}






