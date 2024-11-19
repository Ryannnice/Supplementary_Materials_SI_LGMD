/*
 * Description: LGMD-Plus
 * File: coliasSence_LGMD.c
 * Author: Qinbing Fu
 * Date: 2019 Aug
*/


#include "coliasSense_LGMD.h"
#include "coliasSense_board.h"
#include "delay.h"
#include <math.h>
#include <arm_math.h>
#include <stdlib.h>


uint16_t Image[3][Image_Height][Image_Width];   //store three frames of the initial image
int8_t Diff_Image[2][Image_Height][Image_Width];
LGMD_pControlTypedef hLGMD;

uint8_t LGMD_Param_Init(LGMD_pControlTypedef* hLGMD)
{
	LGMDstruct_Params *Params = &hLGMD->Model->Params;
	
	Params->Tffi = 10;   // 100
	Params->Threshold = 3.0; // threshold in S layer
	Params->nsp = 4;
	Params->nts = 4;
	Params->n = 0;
	Params->flag = 0;
	Params->tau_P = 4;     //4
	Params->tau_E = 6;     //6
	Params->beta1 = 0.3;   //0.3
	Params->beta2 = 0.04;  //0.04
	Params->beta3 = 5.0;     //5
	Params->Taz = 3.0;       //3
	Params->Ts = 0.72;     //0.78           // in arena : 0.65
	Params->Kspi = 1.0;
	hLGMD->Enable = 1;
	hLGMD->status = 0;
	hLGMD->processRate = 30;
	hLGMD->currentImage = 1;
	hLGMD->currentDiffImage = 0;
	hLGMD->AGC_enable_period = 0;
	hLGMD->processCount = 1;
	hLGMD->processCountLast = 1;
	Params->interval = 1000/hLGMD->processRate;
	Params->alpha1 = Params->interval/(Params->interval+Params->tau_P);
	Params->alpha2 = Params->interval/(Params->interval+Params->tau_E);
	printf(" LGMD's parameter initialization is complete");
	return 0;
}

uint8_t Calc_LGMDs_OutOfChannel(LGMD_pControlTypedef* hLGMD)
{
	uint8_t If_m[3] = {2,0,1};   // current image
	uint8_t If_m1[3] = {1,2,0};  // prior image
	int8_t *diff_p;    // point the current diffimage
	uint16_t *cur_image;  // point the current image
	uint16_t *pre_image; // point the prior image
	uint16_t ncell, n;
	int16_t temp;
	uint8_t *cur_E_hat,*pre_E_hat,*cur_P_hat,*pre_P_hat;
	float sum_ffi = 0;  // store ffi	
	float *cur_ffi;
	float alpha1,alpha2;
	
	alpha1 = hLGMD->Model->Params.alpha1;
	alpha2 = hLGMD->Model->Params.alpha2;
	ncell = Image_Height * Image_Width;
	hLGMD->currentDiffImage = !(hLGMD->currentDiffImage);
	cur_E_hat = &hLGMD->Model->Layers.E_hat[hLGMD->currentDiffImage][0][0];
	pre_E_hat = &hLGMD->Model->Layers.E_hat[!hLGMD->currentDiffImage][0][0];
	cur_P_hat = &hLGMD->Model->Layers.P_hat[hLGMD->currentDiffImage][0][0];
	pre_P_hat = &hLGMD->Model->Layers.P_hat[!hLGMD->currentDiffImage][0][0];
	cur_ffi = &hLGMD->Model->Results.FFI_out[hLGMD->currentDiffImage];
	cur_image = (uint16_t*)(((uint8_t*)&Image[If_m[(*(hLGMD->hFrameCount))%3]][0][0])+1);
	pre_image = (uint16_t*)(((uint8_t*)&Image[If_m1[(*(hLGMD->hFrameCount))%3]][0][0])+1);
	diff_p = &Diff_Image[hLGMD->currentDiffImage][0][0];
	for(n=0;n<ncell;n++)
	{
		temp = abs((*cur_image - *pre_image)>>9);
		if(temp > 127)
			temp = 127;
		if(temp < -127)
			temp = -127; 
		*diff_p = temp;
		
		//calculate FFI
		sum_ffi += abs(*diff_p);
		
		//calculate E_hat and P_hat
		*cur_E_hat = alpha2*(*diff_p)+(1-alpha2)*(*pre_E_hat);
		*cur_P_hat = alpha1*(*diff_p)+(1-alpha1)*(*pre_P_hat);
	  diff_p++;
		cur_image++;
		pre_image++;
		cur_P_hat++;
		pre_P_hat++;
		cur_E_hat++;
		pre_E_hat++;
	}
	*cur_ffi = sum_ffi / ncell;
	return 0;
}

uint8_t Fast_LGMD(LGMD_pControlTypedef* hLGMD)
{
	int8_t *cur_diff,*pre_diff;  //point the Diff_Image
	uint8_t Tffi,Taz,nsp,nts,cur_image;
	uint8_t *cur_E_hat,*pre_E_hat,*cur_P_hat,*pre_P_hat;
	uint8_t *layS,*SI_lay;
	uint8_t *flag,*num,*spike;
	uint16_t ncell,n,t;
	uint8_t i,j,height,width;
	uint16_t temg;
	int16_t S_inhibition,L_inhibition;
	float alpha1,alpha2,*beta1,beta2,*beta3;
	float s,ffi,Ts,Kspi,threshold;
	float *mp,*cur_smp;
	LGMDstruct_Layer *Layers = &hLGMD->Model->Layers;
	LGMDstruct_Params *Params = &hLGMD->Model->Params;
	LGMDstruct_Result *Results = &hLGMD->Model->Results;
	
	threshold = Params->Threshold;
	alpha1 = Params->alpha1;
	alpha2 = Params->alpha2;
	beta1 = &Params->beta1;
	beta2 = Params->beta2;
	beta3 = &Params->beta3;
	height = Image_Height;
	Kspi = Params->Kspi;
	width = Image_Width;
	ncell = height*width;
	cur_image = hLGMD->currentDiffImage;
	cur_diff = &Diff_Image[cur_image][0][0];  //P(x,y,t) = E(x,y,t)
  pre_diff = &Diff_Image[!cur_image][0][0]; 
	cur_E_hat = &Layers->E_hat[cur_image][0][0];
	cur_P_hat = &Layers->P_hat[cur_image][0][0];
	Tffi = Params->Tffi;
	Taz = Params->Taz;
	num = &Params->n;
	flag = &Params->flag;
	Ts = Params->Ts;
	ffi = Results->FFI_out[!cur_image];   // FFI has delay
	mp = &Results->MP;
	layS = &Layers->S[0][0];
	SI_lay = &Layers->SI[0][0];
	cur_smp = &Results->SMP;
	spike=&Results->SPIKE;
	*mp = 0;
  *flag = 0;
	
	//1800frames : 1min.    9000frames:  5min

	// default : LI + SI 


	// only LI
	if(*hLGMD->hFrameCount == 5 )
	{
		LED1_Toggle;
	}
	if(*(hLGMD->hFrameCount) >= 9000 && *(hLGMD->hFrameCount) < 18000 ){
		//there is no Self-Inhibition
		if(*hLGMD->hFrameCount == 9000 )
		{
			LED1_Toggle;
		}

		*beta3 = 0;
	}
	
	// no LI, no SI
	if(*(hLGMD->hFrameCount) >= 18000 && *(hLGMD->hFrameCount) < 27000 ){
		// NO SI, NO LI
		LED1_Toggle;
		*beta1 = 0;
		*beta3 = 0;
	}

	if(*(hLGMD->hFrameCount) == 27000 ){
		// NO SI, NO LI
		LED1_Toggle;
		LED1_Toggle;
		*beta1 = 0;
		*beta3 = 0;
	}



	
  //calculate L_inhibition
	for(i=1;i<height-1;i++)
	{
		layS = &Layers->S[i][0];
		cur_E_hat = &Layers->E_hat[cur_image][i][0];
		cur_diff = &Diff_Image[cur_image][i][0];        //Excition
		for(j=1;j<width-1;j++)
		{
			L_inhibition = (*(cur_E_hat+j-width)>>2)   //adjacent
			           +(*(cur_E_hat+j-1)>>2)
			           +(*(cur_E_hat+j+1)>>2)
			           +(*(cur_E_hat+j+width)>>2)
			           +(*(cur_E_hat+j-width-1)>>3)  //diagonal
			           +(*(cur_E_hat+j-width+1)>>3)
			           +(*(cur_E_hat+j+width-1)>>3)
			           +(*(cur_E_hat+j+width+1)>>3);
		  if((*(cur_diff+j))*L_inhibition>=0)
			{
				s = *(cur_diff+j)-L_inhibition*(*beta1);
				if(s<0)
					s=0;
				*(layS+j) = (uint8_t)s;
			}
			else
				*(layS+j) = 0;
		}
	}
	
	//calculate S_inhibition
	for(i=1;i<height-1;i=i+3)
	{
		layS = &Layers->S[i][0];
		cur_P_hat = &Layers->P_hat[cur_image][i][0];
		SI_lay = &Layers->SI[i][0];
		cur_diff = &Diff_Image[cur_image][i][0];
		for(j=1;j<width-1;j=j+3)
		{
			t = *(cur_diff+j)
			+*(cur_diff+j-width)+*(cur_diff+j-1)+*(cur_diff+j+1)+*(cur_diff+j+width)
			+*(cur_diff+j-width-1)+*(cur_diff+j-width+1)+*(cur_diff+j+width-1)+*(cur_diff+j+width+1);
			
			*(SI_lay+j) = *(cur_P_hat+j)
			             +(*(cur_P_hat+j-width)>>2)+(*(cur_P_hat+j-1)>>2)+(*(cur_P_hat+j+1)>>2)+(*(cur_P_hat+j+width)>>2)
			             +(*(cur_P_hat+j-width-1)>>3)+(*(cur_P_hat+j-width+1)>>3)+(*(cur_P_hat+j+width-1)>>2)+(*(cur_P_hat+j+width+1)>>2);
			*(SI_lay+j-width-1) = *(cur_P_hat+j-width-1)+(*(cur_P_hat+j-width)>>2)+(*(cur_P_hat+j-1)>>2)+(*(cur_P_hat+j)>>3);
			*(SI_lay+j-width) = (*(cur_P_hat+j-width-1)>>2)+*(cur_P_hat+j-width)+(*(cur_P_hat+j-width+1)>>2)
			                   +(*(cur_P_hat+j-1)>>3)+(*(cur_P_hat+j)>>2)+(*(cur_P_hat+j+1)>>3);
			*(SI_lay+j-width+1) = (*(cur_P_hat+j-width)>>2)+*(cur_P_hat+j-width+1)+(*(cur_P_hat+j)>>3)+(*(cur_P_hat+j+1)>>2);
			*(SI_lay+j-1) = (*(cur_P_hat+j-width-1)>>2)+(*(cur_P_hat+j-width)>>3)+*(cur_P_hat+j-1)
			               +(*(cur_P_hat+j)>>2)+(*(cur_P_hat+j+width-1)>>2)+(*(cur_P_hat+j+width)>>3);
			*(SI_lay+j+1) = (*(cur_P_hat+j-width)>>3)+(*(cur_P_hat+j-width+1)>>2)+*(cur_P_hat+j+1)
			               +(*(cur_P_hat+j)>>2)+(*(cur_P_hat+j+width)>>3)+(*(cur_P_hat+j+width+1)>>2);
			*(SI_lay+j+width-1) = (*(cur_P_hat+j-1)>>2)+(*(cur_P_hat+j)>>3)+*(cur_P_hat+j+width-1)+(*(cur_P_hat+j+width)>>2);
			*(SI_lay+j+width) = (*(cur_P_hat+j-1)>>3)+(*(cur_P_hat+j)>>2)+(*(cur_P_hat+j+1)>>3)
			                   +(*(cur_P_hat+j+width-1)>>2)+*(cur_P_hat+j+width)+(*(cur_P_hat+j+width+1)>>2);
			*(SI_lay+j+width+1) = (*(cur_P_hat+j)>>3)+(*(cur_P_hat+j+1)>>2)+(*(cur_P_hat+j+width)>>2)+*(cur_P_hat+j+width+1);
			if(t>Taz)
			{
				*(SI_lay+j) = (uint8_t)beta2*(*(SI_lay+j));
				*(SI_lay+j-width-1) = (uint8_t)beta2*(*(SI_lay+j-width-1));
				*(SI_lay+j-width) = (uint8_t)beta2*(*(SI_lay+j-width));
				*(SI_lay+j-width+1) = (uint8_t)beta2*(*(SI_lay+j-width+1));
				*(SI_lay+j-1) = (uint8_t)beta2*(*(SI_lay+j-1));
				*(SI_lay+j+1) = (uint8_t)beta2*(*(SI_lay+j+1));
				*(SI_lay+j+width-1) = (uint8_t)beta2*(*(SI_lay+j+width-1));
				*(SI_lay+j+width) = (uint8_t)beta2*(*(SI_lay+j+width));
				*(SI_lay+j+width+1) = (uint8_t)beta2*(*(SI_lay+j+width+1));
			}
		}
	}

	// gain the layS
	for(i=1;i<height-1;i++)
	{
		layS = &Layers->S[i][0];
		SI_lay = &Layers->SI[i][0];
		for(j=1;j<width-1;j++)
		{
			if((*(layS+j))*(*(SI_lay+j))>=0)
			{
				s = *(layS+j)-(*beta3)*(*(SI_lay+j));
				if(s<threshold)    // the activation threshold of S layer
					s=0;
				*(layS+j) = (uint8_t)s;
			}
			else
				*(layS+j) = 0;
			*mp += *(layS+j);
		}
	}
	
	if(ffi>Tffi)
	{
		*mp=0;
		//LED1_Toggle;
	}
	*cur_smp = 1/(1+exp(-(*mp)/(ncell*Kspi)));
	
	//spike manchism
	if((*cur_smp)>=Ts)
	{
		*spike=1;
		*num += *spike;
		if((*num)>=4)
		{
			*flag = 1;
			*num = 0;
		}
	}
	else
	{
		*spike=0;
		*num = 0;
	}
	return 0;
}

uint8_t Decision_making(LGMD_pControlTypedef* hLGMD,uint8_t allow_motion)
{
	uint8_t command = 0;
	uint8_t collision;
	collision = hLGMD->Model->Params.flag;
	uint8_t fCount=*hLGMD->hFrameCount%3;

	if(*hLGMD->hFrameCount >= 27000 ) // 27000 --> 15 min
	{
		command='S';
	}else{
		if (collision)	//Collision recognition
		{
			LED2_Toggle;

			if(fCount==0)
			{
				command='L';		//long right
			}
			else
			{
				command='R';		//long left
			}
		}
	}

	/*
	if (collision)	//Collision recognition
	{
		LED2_Toggle;
	
		if(fCount==0)
		{
			command='L';		//long right
		}
		else
		{
			command='R';		//long left
		}
	}
	*/

		//attention
		
	


	if(allow_motion)
	{
		/************* Pay Attention! ******************/
		//if motion sequence is not empty, do not push motions
		uint8_t readMotionQueueStatus[10]={0};
		ReadBytesWithACKFromColiasBasic(hCoS->hHCos_Motion->Instance->hHUART,0x97,readMotionQueueStatus,1);
		if((readMotionQueueStatus[1]&0x0f)==0x09)
		{
			fast_motion_control(hCoS->hHCos_Motion->Instance,command,0);	
		}
	}
	
	return 0;
}

uint8_t LGMD_demo(LGMD_pControlTypedef* hLGMD)
{
	uint8_t tmpOT;
	float SMP = hLGMD->Model->Results.SMP;
	if(hLGMD->Enable==0)
		return 0;
	while (*(hLGMD ->hFrameCount) == hLGMD->processCountLast)
		hLGMD->status=0;	//wait until next frame if processing duration is less than one frame
	/*
	if (hLGMD->processRate)
		while ((*(hLGMD->hFrameCount))*(hLGMD->processRate) < (hLGMD->processCount)*OV7670_FPS)	//wait until expected processing rate	
			hLGMD->status=0;	//wait until next frame if processing duration is less than one frame
	*/
	
	hLGMD->status=1;
	TICin2; //start clocking,reset timer to 0
	hCoS->hHBIO->Instance->timerlog[0]=TOCin2; //save timer tick
	Calc_LGMDs_OutOfChannel(hLGMD);
	hCoS->hHBIO->Instance->timerlog[1]=TOCin2; //save timer tick
	Fast_LGMD(hLGMD);
	hCoS->hHBIO->Instance->timerlog[2]=TOCin2; //save timer tick
	//delay_ms(10);
	
	//print elapsed time
	//printf("%d \r\n",hCoS->hHBIO->Instance->timerlog[2]-hCoS->hHBIO->Instance->timerlog[0]);
	
	hLGMD->processCountLast = *(hLGMD ->hFrameCount);
	(hLGMD->processCount)++;
	//hCoS->hHBIO->Instance->timerlog[2]=TOCin2; //save timer tick
	tmpOT=256*(SMP-0.5);
	
	GPIOD->ODR=(GPIOD->ODR&0xffffff00)+(~tmpOT);
	
	return 0;
}
