#ifndef __LCD18TFT_H
#define __LCD18TFT_H

#include "stm32f4xx_hal.h"
#include "coliasSense_def.h"

#define LCD_H 128
#define LCD_W 160


#define TFT_CS_Pin GPIO_PIN_15 //e7
#define TFT_CS_GPIO_Port GPIOD

#define TFT_RST_Pin GPIO_PIN_14 //e6
#define TFT_RST_GPIO_Port GPIOD

#define TFT_A0_Pin GPIO_PIN_13 //e5
#define TFT_A0_GPIO_Port GPIOD

#define TFT_SDA_Pin GPIO_PIN_12 //e4
#define TFT_SDA_GPIO_Port GPIOD

#define TFT_SCK_Pin GPIO_PIN_11 //e3
#define TFT_SCK_GPIO_Port GPIOD

#define TFT_LED_Pin GPIO_PIN_10 //e2
#define TFT_LED_GPIO_Port GPIOD


#define TFT_CS_ON   (HAL_GPIO_WritePin(TFT_CS_GPIO_Port,  TFT_CS_Pin, 1))
#define TFT_CS_OFF  (HAL_GPIO_WritePin(TFT_CS_GPIO_Port,  TFT_CS_Pin, 0))

#define TFT_RST_ON  (HAL_GPIO_WritePin(TFT_RST_GPIO_Port, TFT_RST_Pin, 1))
#define TFT_RST_OFF (HAL_GPIO_WritePin(TFT_RST_GPIO_Port, TFT_RST_Pin, 0))

#define TFT_A0_ON   (HAL_GPIO_WritePin(TFT_A0_GPIO_Port,  TFT_A0_Pin, 1))
#define TFT_A0_OFF  (HAL_GPIO_WritePin(TFT_A0_GPIO_Port,  TFT_A0_Pin, 0))

#define TFT_SDA_ON  (HAL_GPIO_WritePin(TFT_SDA_GPIO_Port, TFT_SDA_Pin, 1))
#define TFT_SDA_OFF (HAL_GPIO_WritePin(TFT_SDA_GPIO_Port, TFT_SDA_Pin, 0))

#define TFT_SCK_ON  (HAL_GPIO_WritePin(TFT_SCK_GPIO_Port, TFT_SCK_Pin, 1))
#define TFT_SCK_OFF (HAL_GPIO_WritePin(TFT_SCK_GPIO_Port, TFT_SCK_Pin, 0))

#define TFT_LED_ON  (HAL_GPIO_WritePin(TFT_LED_GPIO_Port, TFT_LED_Pin, 1))
#define TFT_LED_OFF (HAL_GPIO_WritePin(TFT_LED_GPIO_Port, TFT_LED_Pin, 0))

#define LCDELAY delay_short(SSCKD)

extern u16 BACK_COLOR, POINT_COLOR;   //背景色，画笔色	
extern u16 COLORTABLE[];
extern const u8 picc1[];
extern u16 CBuff[];


typedef struct
{
	float C1;//YUV:Y  RGB:R  HSV:H
	float C2;//YUV:U  RGB:G  HSV:S
	float C3;//YUV:V  RGB:B  HSV:V
}ColordataTypedef;//




/*
const char welcome_psr[] ;
const char welcome_psr2[];
const char welcome_psr3[];
const char welcome_psr4[];
	*/
void LCD_WR_DATA8(char da);
 void LCD_WR_DATA(int da);
 void LCD_WR_REG(char da);
 void LCD_WR_REG_DATA(int reg,int da);
 void Address_set(u16 x1,u16 y1,u16 x2,u16 y2);
 
 void TFT_LCD_Init(void);
 void LCD_Clear(u16 Color);
 void LCD_DrawPoint(u16 x,u16 y);//画点
void LCD_DrawPoint_big(u16 x,u16 y);//画一个大点
u16  LCD_ReadPoint(u16 x,u16 y); //读点
void Draw_Circle(u16 x0,u16 y0,u8 r);
void LCD_DrawLine(u16 x1, u16 y1, u16 x2, u16 y2);
void LCD_DrawRectangle(u16 x1, u16 y1, u16 x2, u16 y2);		   
void LCD_Fill(u16 xsta,u16 ysta,u16 xend,u16 yend,u16 color);
void LCD_ShowChar(u16 x,u16 y,u8 num,u8 mode);//显示一个字符
void LCD_ShowNum(u16 x,u16 y,u32 num,u8 len);//显示数字
void LCD_Show2Num(u16 x,u16 y,u16 num,u8 len);//显示2个数字
void LCD_ShowString(u16 x,u16 y,const u8 *p);		 //显示一个字符串,16字体
 void LCD_DISP_PIC(u8* pstart,u16 size_char);
 void LCD_imgshow(u16 x,u16 y, u16 width, u16 height, u16* pstart, u8 mode);
 //画笔颜色
#define WHITE         	 0xFFFF
#define BLACK         	 0x0000	  
#define BLUE         	 0x001F  
#define BRED             0XF81F
#define GRED 			 0XFFE0
#define GBLUE			 0X07FF
#define RED           	 0xF800
#define MAGENTA       	 0xF81F
#define GREEN         	 0x07E0
#define CYAN          	 0x7FFF
#define YELLOW        	 0xFFE0
#define BROWN 			 0XBC40 //棕色
#define BRRED 			 0XFC07 //棕红色
#define GRAY  			 0X8430 //灰色
//GUI颜色

#define DARKBLUE      	 0X01CF	//深蓝色
#define LIGHTBLUE      	 0X7D7C	//浅蓝色  
#define GRAYBLUE       	 0X5458 //灰蓝色
//以上三色为PANEL的颜色 
 
#define LIGHTGREEN     	 0X841F //浅绿色
#define LGRAY 			 0XC618 //浅灰色(PANNEL),窗体背景色

#define LGRAYBLUE        0XA651 //浅灰蓝色(中间层颜色)
#define LBBLUE           0X2B12 //浅棕蓝色(选择条目的反色)
 
#endif