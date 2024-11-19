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
uint8_t S[Image_Height][Image_Width];
uint8_t G[Image_Height][Image_Width];
LGMD_pControlTypedef hLGMD;

uint8_t LGMD_Param_Init(LGMD_pControlTypedef* hLGMD)
{
	float on_base_tau = 15;
	float off_base_tau = 60;
	LGMDstruct_Params *Params = &hLGMD->Model->Params;
	Params->Tpm = 20;
	Params->n = 0;
	Params->flag = 0;
	Params->Cw = 4;
	Params->CdeTde = 35;
	Params->a4 = 4;
	Params->clip = 0;
	Params->a1 = 0.1;
	Params->ffi_tau = 90;
	Params->Won_base = 1;
	Params->Woff_base = 0.5;
	Params->o1 = 0.5;
	Params->o2 = 1;
	Params->o3 = 1;
	Params->Delatc = 0.01;
	Params->Kspi = 0.5;  	//0.5~1
	Params->Tspi = 0.64;      //0.65~0.78
	Params->Tsfa = 0.003;
	Params->sfa_tau = 500;   //500~1000
	Params->on_tau = 5.0;
	Params->off_tau = 20.0;
	Params->beta1 = 0.3;
	Params->beta2 = 0.1;          
	Params->beta3 = 3.0;
	Params->beta4 = 4.0;
	Params->r1 = 3.0;
	Params->r2 = 0.005;
	
	hLGMD->currentDiffImage=0;
	hLGMD->currentImage=1;
	hLGMD->Enable=1;
	hLGMD->Model=&LGMD;
	hLGMD->processCount=1;
	hLGMD->processCountLast=1;
	hLGMD->status=0;
	hLGMD->processRate=30;
	hLGMD->AGC_enable_period=0;
	
	Params->interval = 1000/hLGMD->processRate;
	Params->a2 = Params->interval/(Params->interval+Params->ffi_tau);  //delay coefficient in PM pathway
	Params->a3 = Params->sfa_tau/(Params->sfa_tau+Params->interval);   //delay coefficient in SFA mechanism
	Params->on_center_tau = on_base_tau;
	Params->on_near_tau = 2*on_base_tau;
	Params->on_diag_tau = 3*on_base_tau;
	Params->on_delay_center = Params->interval/(Params->interval+Params->on_center_tau);
	Params->on_delay_near = Params->interval/(Params->interval+Params->on_near_tau);
	Params->on_delay_diag = Params->interval/(Params->interval+Params->on_diag_tau);
	Params->off_center_tau = off_base_tau;
	Params->off_near_tau = 2*off_base_tau;
	Params->off_diag_tau = 3*off_base_tau;
	Params->off_delay_center = Params->interval/(Params->interval+Params->off_center_tau);
	Params->off_delay_near = Params->interval/(Params->interval+Params->off_near_tau);
	Params->off_delay_diag = Params->interval/(Params->interval+Params->off_diag_tau);
	Params->on_delay = Params->interval/(Params->interval+Params->on_tau);
	Params->off_delay = Params->interval/(Params->interval+Params->off_tau);
	
	printf(" LGMD's parameter initialization is complete");
	return 0;
}

uint8_t Calc_LGMDs_OutOfChannel(LGMD_pControlTypedef* hLGMD)
{
	uint8_t Im[3] = {2,0,1};
	uint8_t Im1[3] = {1,2,0};
	int8_t *diff_image;
	uint8_t Tpm,clip,i,j;
	uint16_t ncell,n;
	uint16_t *cur_image,*pre_image;
	int16_t temg;
	uint8_t *cur_Pon,*pre_Pon,*cur_Poff,*pre_Poff,*on_hat,*off_hat;
	float a1,a2,a5,a6,*Won_base,*Woff_base,sum_ffi,ffi;
	float *w1,*w2,*cur_ffi,*pre_ffi;

	clip = hLGMD->Model->Params.clip;
  hLGMD->currentDiffImage = !hLGMD->currentDiffImage;
	Tpm = hLGMD->Model->Params.Tpm;
	ncell = Image_Width*Image_Height;
	a1 = hLGMD->Model->Params.a1;
	a2 = hLGMD->Model->Params.a2;
	a5 = hLGMD->Model->Params.on_delay;
	a6 = hLGMD->Model->Params.off_delay;
	Won_base = &hLGMD->Model->Params.Won_base;
	Woff_base = &hLGMD->Model->Params.Woff_base;
	w1 = &hLGMD->Model->Params.w1;
	w2 = &hLGMD->Model->Params.w2;
	cur_ffi = &hLGMD->Model->Results.FFI_out[hLGMD->currentDiffImage];
	pre_ffi = &hLGMD->Model->Results.FFI_out[!hLGMD->currentDiffImage];
	cur_image = (uint16_t*)(((uint8_t*)&Image[Im[(*(hLGMD->hFrameCount))%3]][0][0])+1);
  pre_image = (uint16_t*)(((uint8_t*)&Image[Im1[(*(hLGMD->hFrameCount))%3]][0][0])+1);
	diff_image = &Diff_Image[hLGMD->currentDiffImage][0][0];
	cur_Pon = &hLGMD->Model->Layers.ON[hLGMD->currentDiffImage][0][0];
	pre_Pon = &hLGMD->Model->Layers.ON[!hLGMD->currentDiffImage][0][0];
	cur_Poff = &hLGMD->Model->Layers.OFF[hLGMD->currentDiffImage][0][0];
	pre_Poff = &hLGMD->Model->Layers.OFF[!hLGMD->currentDiffImage][0][0];
	on_hat = &hLGMD->Model->Layers.ON_hat[0][0];
	off_hat = &hLGMD->Model->Layers.OFF_hat[0][0];
	
	sum_ffi = 0;
	for(n=0;n<ncell;n++)
	{
		temg = (*cur_image - *pre_image)>>9;
		if(temg>127)
			temg = 127;
		if(temg<-127)
			temg = -127;
		*diff_image = temg;
		sum_ffi += abs(*diff_image);
		
		//half-waving rectify mechanism
		if(*diff_image>=clip)
		{
			*cur_Pon = *(diff_image);
			*cur_Poff = 0;
		}
		else
		{
			*cur_Pon = 0;
			*cur_Poff = abs(*diff_image);
	  }
		*on_hat = a5*(*cur_Pon)+(1-a5)*(*pre_Pon);
		*off_hat = a6*(*cur_Poff)+(1-a6)*(*pre_Poff);
		
		diff_image++;
		cur_image++;
		pre_image++;
		cur_Pon++;
		cur_Poff++;
		pre_Pon++;
		pre_Poff++;
		on_hat++;
		off_hat++;
	}
	
	*cur_ffi = sum_ffi/ncell;
  ffi = a2*(*cur_ffi)+(1-a2)*(*pre_ffi);
	*w1 = (ffi)/Tpm;
	if(*Won_base>*w1)
		*w1 = *Won_base;
	*w2 = (ffi)/Tpm;
	if(*Woff_base>*w2)
		*w2 = *Woff_base; 
	
	return 0;
}

uint8_t Fast_LGMD(LGMD_pControlTypedef* hLGMD)
{
	LGMDstruct_Layer *Layers = &hLGMD->Model->Layers;
	LGMDstruct_Params *Params = &hLGMD->Model->Params;
	LGMDstruct_Result *Results = &hLGMD->Model->Results;
	int8_t *diff_image;
	uint8_t cur_image,clip;
	uint8_t width,height,i,j;
	uint8_t *cur_Pon,*pre_Pon,*cur_Poff,*pre_Poff,*on_hat,*off_hat;
	uint8_t *layS,*layG,*SI_on,*SI_off;
	uint8_t t_max,Cw,CdeTde,a4;
	uint8_t *flag,*num,*spike;
	uint16_t ncell,n,center,near,diag,inhibition1,inhibition2,t1,t2;
	uint16_t m1,m2,m3,m4,m5,m6,m7,m8,m9;
	uint16_t n1,n2,n3,n4,n5,n6,n7,n8,n9;
	uint16_t pre_center,pre_near,pre_diag,sum;
	float a1,*w1,*w2,s,Son,Soff,o1,o2,o3,Delatc,scale,a3,beta1,beta2,*beta3,*beta4;
	float on_center,on_near,on_diag,off_center,off_near,off_diag;
	float *mp,*cur_smp,*pre_smp,*cur_lgmd_out,*pre_lgmd_out;
	float Kspi,Tsfa,Tspi,r1,r2,a5,a6;
	
	t_max = 0;
	clip = hLGMD->Model->Params.clip;
	height = Image_Height;
	width = Image_Width;
	ncell = (height-H_Start*2)*(width-W_Start*2);
	a1 = Params->a1;
	w1 = &Params->w1;
	w2 = &Params->w2;
	o1 = Params->o1;
	o2 = Params->o2;
	o3 = Params->o3;
	a3 = Params->a3;
	a4 = Params->a4;
	a5 = hLGMD->Model->Params.on_delay;
	a6 = hLGMD->Model->Params.off_delay;
	Cw = Params->Cw;
	CdeTde = Params->CdeTde;
	Delatc = Params->Delatc;
	beta1 = Params->beta1;
	beta2 = Params->beta2;
	beta3 = &Params->beta3;
	beta4 = &Params->beta4;
	Kspi = Params->Kspi;
	Tsfa = Params->Tsfa;
	Tspi = Params->Tspi;
	r1 = Params->r1;
	r2 = Params->r2;
	mp = &Results->MP;
	flag = &Params->flag;
	num = &Params->n;
	spike = &Results->SPIKE;
	cur_image = hLGMD->currentDiffImage;
	on_center = Params->on_delay_center;
	on_near = Params->on_delay_near;
	on_diag = Params->on_delay_diag;
	off_center = Params->off_delay_center;
	off_near = Params->off_delay_near;
	off_diag = Params->off_delay_diag;
	diff_image = &Diff_Image[cur_image][0][0];
	cur_Pon = &Layers->ON[cur_image][0][0];
	pre_Pon = &Layers->ON[!cur_image][0][0];
	cur_Poff = &Layers->OFF[cur_image][0][0];
	pre_Poff = &Layers->OFF[!cur_image][0][0];
	on_hat = &Layers->ON_hat[0][0];
	off_hat = &Layers->OFF_hat[0][0];
	cur_smp = &Results->SMP[cur_image];
	pre_smp = &Results->SMP[!cur_image];
	cur_lgmd_out = &Results->LGMD_out[cur_image];
	pre_lgmd_out = &Results->LGMD_out[!cur_image];
	layS = &S[0][0];
	layG = &G[0][0];
	*mp = 0;
	*flag = 0;
	

	//1800frames : 1min.    9000frames:  5min

	// default : LI + SI 
	if(*hLGMD->hFrameCount == 5 )
	{
		LED1_Toggle;
	}

	// only LI
	if(*(hLGMD->hFrameCount) >= 9000 && *(hLGMD->hFrameCount) < 18000 ){
		//there is no Self-Inhibition
		if(*hLGMD->hFrameCount == 9000 )
		{
			LED1_Toggle;
		}

		*beta3 = 0;
		*beta4 = 0;
	}
	
	// no LI, no SI
	if(*(hLGMD->hFrameCount) >= 18000 && *(hLGMD->hFrameCount) < 27000 ){
		// NO SI, NO LI
		LED1_Toggle;
		*beta4 = 0;
		*beta3 = 0;
		*w1 = 0;
		*w2 = 0;
	}
	
	/*
	for(i=H_Start;i<Image_Height-H_Start;i++)
  {
	  diff_image = &Diff_Image[cur_image][i][0];
		cur_Pon = &hLGMD->Model->Layers.ON[cur_image][i-H_Start][0];
		cur_Poff = &hLGMD->Model->Layers.OFF[cur_image][i-H_Start][0];
		pre_Pon = &hLGMD->Model->Layers.ON[!cur_image][i-H_Start][0];
		pre_Poff = &hLGMD->Model->Layers.OFF[!cur_image][i-H_Start][0];
		on_hat = &hLGMD->Model->Layers.ON_hat[i-H_Start][0];
	  off_hat = &hLGMD->Model->Layers.OFF_hat[i-H_Start][0];
		for(j=W_Start;j<Image_Width-W_Start;j++)
		{
			//half-waving rectify mechanism
			if(*(diff_image+j)>=clip)
			{
				*(cur_Pon-W_Start+j) = *(diff_image+j);
			  *(cur_Poff-W_Start+j) = 0;
		  }
		  else
		  {
			  *(cur_Pon-W_Start+j) = 0;
			  *(cur_Poff-W_Start+j) = abs(*(diff_image+j));
		  }
			
			//calculate Pon_hat and Poff_hat to gain S_inhibition
		  *(on_hat-W_Start+j) = a5*(*(cur_Pon-W_Start+j))+(1-a5)*(*(pre_Pon-W_Start+j));
		  *(off_hat-W_Start+j) = a6*(*(cur_Poff-W_Start+j))+(1-a6)*(*(pre_Poff-W_Start+j));
		}
	}	
	*/
	
	//calculate L_inhibition
	for(i=H_Start;i<Image_Height-H_Start;i++)
	{
		cur_Pon = &Layers->ON[cur_image][i][0];
		pre_Pon = &Layers->ON[!cur_image][i][0];
		cur_Poff = &Layers->OFF[cur_image][i][0];
		pre_Poff = &Layers->OFF[!cur_image][i][0];
		layS = &S[i][0];
		
		for(j=W_Start;j<Image_Width-W_Start;j++)
		{
			//calculate the inhibition of ON channel
			center = (*(cur_Pon+j))<<1;
			near = (*(cur_Pon+j-width)>>1)+(*(cur_Pon+j-1)>>1)+(*(cur_Pon+j+1)>>1)+(*(cur_Pon+j+width)>>1);
			diag = (*(cur_Pon+j-width-1)>>2)+(*(cur_Pon+j-width+1)>>2)+(*(cur_Pon+j+width-1)>>2)+(*(cur_Pon+j+width+1)>>2);
			pre_center = (*(pre_Pon+j))<<1;
			pre_near = (*(pre_Pon+j-width)>>1)+(*(pre_Pon+j-1)>>1)+(*(pre_Pon+j+1)>>1)+(*(pre_Pon+j+width)>>1);
			pre_diag = (*(pre_Pon+j-width-1)>>2)+(*(pre_Pon+j-width+1)>>2)+(*(pre_Pon+j+width-1)>>2)+(*(pre_Pon+j+width+1)>>2);
			inhibition1 = (uint16_t)(center*on_center+near*on_near+diag*on_diag+pre_center*(1-on_center)+pre_near*(1-on_near)+pre_diag*(1-on_diag));
			
			//calculate the inhibition of OFF chanel
			center = *(cur_Poff+j);
			near = (*(cur_Poff+j-width)>>2)+(*(cur_Poff+j-1)>>2)+(*(cur_Poff+j+1)>>2)+(*(cur_Poff+j+width)>>2);
			diag = (*(cur_Poff+j-width-1)>>3)+(*(cur_Poff+j-width+1)>>3)+(*(cur_Poff+j+width-1)>>3)+(*(cur_Poff+j+width+1)>>3);
      pre_center = *(pre_Poff+j);
			pre_near = (*(pre_Poff+j-width)>>2)+(*(pre_Poff+j-1)>>2)+(*(pre_Poff+j+1)>>2)+(*(pre_Poff+j+width)>>2);
			pre_diag = (*(pre_Poff+j-width-1)>>3)+(*(pre_Poff+j-width+1)>>3)+(*(pre_Poff+j+width-1)>>3)+(*(pre_Poff+j+width+1)>>3);
			inhibition2 = (uint16_t)(center*off_center+near*off_near+diag*off_diag+pre_center*(1-off_center)+pre_near*(1-off_near)+pre_diag*(1-off_diag));
			
			Son = *(cur_Pon+j)-(*w1)*inhibition1;
			if(Son<0)
				Son = 0;
			Soff = *(cur_Poff+j)-(*w2)*inhibition2;
			if(Soff<0)
				Soff = 0;
			*(layS+j) = (uint8_t)(o1*Son+o2*Soff);
		}
	}
	
	//calculate S_inhibition
	for(i=H_Start;i<Image_Height-H_Start;i=i+3)
	{
		on_hat = &Layers->ON_hat[i][0];
		off_hat = &Layers->OFF_hat[i][0];
		//SI_on = &Layers->SI_ON[i][0];
		//SI_off = &Layers->SI_OFF[i][0];	
		cur_Pon = &Layers->ON[cur_image][i][0];
		cur_Poff = &Layers->OFF[cur_image][i][0];
		layS = &S[i][0];
		for(j=W_Start;j<Image_Width-W_Start;j=j+3)
		{
			t1 = *(cur_Pon+j)
			+*(cur_Pon+j-width)+*(cur_Pon+j-1)+*(cur_Pon+j+1)+*(cur_Pon+j+width)
			+*(cur_Pon+j-width-1)+*(cur_Pon+j-width+1)+*(cur_Pon+j+width-1)+*(cur_Pon+j+width+1);
			
			t2 = *(cur_Poff+j)
			+*(cur_Poff+j-width)+*(cur_Poff+j-1)+*(cur_Poff+j+1)+*(cur_Poff+j+width)
			+*(cur_Poff+j-width-1)+*(cur_Poff+j-width+1)+*(cur_Poff+j+width-1)+*(cur_Poff+j+width+1);
			
			//ON channel
			/*
			*(SI_on+j) = (*(on_hat+j)<<2)
			             +*(on_hat+j-width)+*(on_hat+j-1)+*(on_hat+j+1)+*(on_hat+j+width)
			             +(*(on_hat+j-width-1)>>1)+(*(on_hat+j-width+1)>>1)+*(on_hat+j+width-1)+*(on_hat+j+width+1);
			*(SI_on+j-width-1) = (*(on_hat+j-width-1)<<2)+*(on_hat+j-width)+*(on_hat+j-1)+(*(on_hat+j)>>1);
			*(SI_on+j-width) = *(on_hat+j-width-1)+(*(on_hat+j-width)<<2)+*(on_hat+j-width+1)
			                   +(*(on_hat+j-1)>>1)+*(on_hat+j)+(*(on_hat+j+1)>>1);
			*(SI_on+j-width+1) = *(on_hat+j-width)+(*(on_hat+j-width+1)<<2)+(*(on_hat+j)>>1)+*(on_hat+j+1);
			*(SI_on+j-1) = *(on_hat+j-width-1)+(*(on_hat+j-width)>>1)+(*(on_hat+j-1)<<2)
			               +*(on_hat+j)+*(on_hat+j+width-1)+(*(on_hat+j+width)>>1);
			*(SI_on+j+1) = (*(on_hat+j-width)>>1)+*(on_hat+j-width+1)+(*(on_hat+j+1)<<2)
			               +*(on_hat+j)+(*(on_hat+j+width)>>1)+*(on_hat+j+width+1);
			*(SI_on+j+width-1) = *(on_hat+j-1)+(*(on_hat+j)>>1)+(*(on_hat+j+width-1)<<2)+*(on_hat+j+width);
			*(SI_on+j+width) = (*(on_hat+j-1)>>1)+*(on_hat+j)+(*(on_hat+j+1)>>1)
			                   +*(on_hat+j+width-1)+(*(on_hat+j+width)<<2)+*(on_hat+j+width+1);
			*(SI_on+j+width+1) = (*(on_hat+j)>>1)+*(on_hat+j+1)+*(on_hat+j+width)+(*(on_hat+j+width+1)<<2);
			*/
			m1 = (*(on_hat+j-width-1)<<2)+*(on_hat+j-width)+*(on_hat+j-1)+(*(on_hat+j)>>1);
			m2 = *(on_hat+j-width-1)+(*(on_hat+j-width)<<2)+*(on_hat+j-width+1)
			     +(*(on_hat+j-1)>>1)+*(on_hat+j)+(*(on_hat+j+1)>>1);
			m3 = *(on_hat+j-width)+(*(on_hat+j-width+1)<<2)+(*(on_hat+j)>>1)+*(on_hat+j+1);
			m4 = *(on_hat+j-width-1)+(*(on_hat+j-width)>>1)+(*(on_hat+j-1)<<2)
			     +*(on_hat+j)+*(on_hat+j+width-1)+(*(on_hat+j+width)>>1);
			m5 = (*(on_hat+j)<<2)
			     +*(on_hat+j-width)+*(on_hat+j-1)+*(on_hat+j+1)+*(on_hat+j+width)
			     +(*(on_hat+j-width-1)>>1)+(*(on_hat+j-width+1)>>1)+*(on_hat+j+width-1)+*(on_hat+j+width+1);	
      m6 = (*(on_hat+j-width)>>1)+*(on_hat+j-width+1)+(*(on_hat+j+1)<<2)
			     +*(on_hat+j)+(*(on_hat+j+width)>>1)+*(on_hat+j+width+1);
			m7 = *(on_hat+j-1)+(*(on_hat+j)>>1)+(*(on_hat+j+width-1)<<2)+*(on_hat+j+width);
      m8 = (*(on_hat+j-1)>>1)+*(on_hat+j)+(*(on_hat+j+1)>>1)
			     +*(on_hat+j+width-1)+(*(on_hat+j+width)<<2)+*(on_hat+j+width+1);		
      m9 = (*(on_hat+j)>>1)+*(on_hat+j+1)+*(on_hat+j+width)+(*(on_hat+j+width+1)<<2);					 
			
			//OFF channel
			/*
			*(SI_off+j) = (*(off_hat+j)<<1)
			             +(*(off_hat+j-width)>>1)+(*(off_hat+j-1)>>1)+(*(off_hat+j+1)>>1)+(*(off_hat+j+width)>>1)
			             +(*(off_hat+j-width-1)>>2)+(*(off_hat+j-width+1)>>2)+(*(off_hat+j+width-1)>>1)+(*(off_hat+j+width+1)>>1);
			*(SI_off+j-width-1) = (*(off_hat+j-width-1)<<1)+(*(off_hat+j-width)>>1)+(*(off_hat+j-1)>>1)+(*(off_hat+j)>>2);
			*(SI_off+j-width) = (*(off_hat+j-width-1)>>1)+(*(off_hat+j-width)<<1)+(*(off_hat+j-width+1)>>1)
			                   +(*(off_hat+j-1)>>2)+(*(off_hat+j)>>1)+(*(off_hat+j+1)>>2);
			*(SI_off+j-width+1) = (*(off_hat+j-width)>>1)+(*(off_hat+j-width+1)<<1)+(*(off_hat+j)>>2)+(*(off_hat+j+1)>>1);
			*(SI_off+j-1) = (*(off_hat+j-width-1)>>1)+(*(off_hat+j-width)>>2)+(*(off_hat+j-1)<<1)
			               +(*(off_hat+j)>>1)+(*(off_hat+j+width-1)>>1)+(*(off_hat+j+width)>>2);
			*(SI_off+j+1) = (*(off_hat+j-width)>>2)+(*(off_hat+j-width+1)>>1)+(*(off_hat+j+1)<<1)
			               +(*(off_hat+j)>>1)+(*(off_hat+j+width)>>2)+(*(off_hat+j+width+1)>>1);
			*(SI_off+j+width-1) = (*(off_hat+j-1)>>1)+(*(off_hat+j)>>2)+(*(off_hat+j+width-1)<<1)+(*(off_hat+j+width)>>1);
			*(SI_off+j+width) = (*(off_hat+j-1)>>2)+(*(off_hat+j)>>1)+(*(off_hat+j+1)>>2)
			                   +(*(off_hat+j+width-1)>>1)+(*(off_hat+j+width)<<1)+(*(off_hat+j+width+1)>>1);
			*(SI_off+j+width+1) = (*(off_hat+j)>>2)+(*(off_hat+j+1)>>1)+(*(off_hat+j+width)>>1)+(*(off_hat+j+width+1)<<1);
			*/
			n1 = (*(off_hat+j-width-1)<<1)+(*(off_hat+j-width)>>1)+(*(off_hat+j-1)>>1)+(*(off_hat+j)>>2);
			n2 = (*(off_hat+j-width-1)>>1)+(*(off_hat+j-width)<<1)+(*(off_hat+j-width+1)>>1)
			     +(*(off_hat+j-1)>>2)+(*(off_hat+j)>>1)+(*(off_hat+j+1)>>2);
			n3 = (*(off_hat+j-width)>>1)+(*(off_hat+j-width+1)<<1)+(*(off_hat+j)>>2)+(*(off_hat+j+1)>>1);
			n4 = (*(off_hat+j-width-1)>>1)+(*(off_hat+j-width)>>2)+(*(off_hat+j-1)<<1)
			     +(*(off_hat+j)>>1)+(*(off_hat+j+width-1)>>1)+(*(off_hat+j+width)>>2);
		  n5 = (*(off_hat+j)<<1)
			     +(*(off_hat+j-width)>>1)+(*(off_hat+j-1)>>1)+(*(off_hat+j+1)>>1)+(*(off_hat+j+width)>>1)
			     +(*(off_hat+j-width-1)>>2)+(*(off_hat+j-width+1)>>2)+(*(off_hat+j+width-1)>>1)+(*(off_hat+j+width+1)>>1);
			n6 = (*(off_hat+j-width)>>2)+(*(off_hat+j-width+1)>>1)+(*(off_hat+j+1)<<1)
			     +(*(off_hat+j)>>1)+(*(off_hat+j+width)>>2)+(*(off_hat+j+width+1)>>1);
			n7 = (*(off_hat+j-1)>>1)+(*(off_hat+j)>>2)+(*(off_hat+j+width-1)<<1)+(*(off_hat+j+width)>>1);
			n8 = (*(off_hat+j-1)>>2)+(*(off_hat+j)>>1)+(*(off_hat+j+1)>>2)
			     +(*(off_hat+j+width-1)>>1)+(*(off_hat+j+width)<<1)+(*(off_hat+j+width+1)>>1);
			n9 = (*(off_hat+j)>>2)+(*(off_hat+j+1)>>1)+(*(off_hat+j+width)>>1)+(*(off_hat+j+width+1)<<1);
			
			if(t1>r1)
			{
				/*
				*(SI_on+j) = (uint8_t)beta1*(*(SI_on+j));
				*(SI_on+j-width-1) = (uint8_t)beta1*(*(SI_on+j-width-1));
				*(SI_on+j-width) = (uint8_t)beta1*(*(SI_on+j-width));
				*(SI_on+j-width+1) = (uint8_t)beta1*(*(SI_on+j-width+1));
				*(SI_on+j-1) = (uint8_t)beta1*(*(SI_on+j-1));
				*(SI_on+j+1) = (uint8_t)beta1*(*(SI_on+j+1));
				*(SI_on+j+width-1) = (uint8_t)beta1*(*(SI_on+j+width-1));
				*(SI_on+j+width) = (uint8_t)beta1*(*(SI_on+j+width));
				*(SI_on+j+width+1) = (uint8_t)beta1*(*(SI_on+j+width+1));
				*/
				m1 = (uint16_t)beta1*m1;
				m2 = (uint16_t)beta1*m2;
				m3 = (uint16_t)beta1*m3;
				m4 = (uint16_t)beta1*m4;
				m5 = (uint16_t)beta1*m5;
				m6 = (uint16_t)beta1*m6;
				m7 = (uint16_t)beta1*m7;
				m8 = (uint16_t)beta1*m8;
				m9 = (uint16_t)beta1*m9;
			}
			
			if(t2>r2)
			{
				/*
				*(SI_off+j) = (uint8_t)beta2*(*(SI_off+j));
				*(SI_off+j-width-1) = (uint8_t)beta2*(*(SI_off+j-width-1));
				*(SI_off+j-width) = (uint8_t)beta2*(*(SI_off+j-width));
				*(SI_off+j-width+1) = (uint8_t)beta2*(*(SI_off+j-width+1));
				*(SI_off+j-1) = (uint8_t)beta2*(*(SI_off+j-1));
				*(SI_off+j+1) = (uint8_t)beta2*(*(SI_off+j+1));
				*(SI_off+j+width-1) = (uint8_t)beta2*(*(SI_off+j+width-1));
				*(SI_off+j+width) = (uint8_t)beta2*(*(SI_off+j+width));
				*(SI_off+j+width+1) = (uint8_t)beta2*(*(SI_off+j+width+1));
				*/
				n1 = (uint16_t)beta1*n1;
				n2 = (uint16_t)beta1*n2;
				n3 = (uint16_t)beta1*n3;
				n4 = (uint16_t)beta1*n4;
				n5 = (uint16_t)beta1*n5;
				n6 = (uint16_t)beta1*n6;
				n7 = (uint16_t)beta1*n7;
				n8 = (uint16_t)beta1*n8;
				n9 = (uint16_t)beta1*n9;
			}
			*(layS+j-width-1) = (uint8_t)(*(layS+j-width-1)-o1*(*beta3)*m1-o2*(*beta4)*n1);
			*(layS+j-width) = (uint8_t)(*(layS+j-width)-o1*(*beta3)*m2-o2*(*beta4)*n2);
			*(layS+j-width+1) = (uint8_t)(*(layS+j-width+1)-o1*(*beta3)*m3-o2*(*beta4)*n3);
			*(layS+j-1) = (uint8_t)(*(layS+j-1)-o1*(*beta3)*m4-o2*(*beta4)*n4);
			*(layS+j) = (uint8_t)(*(layS+j)-o1*(*beta3)*m5-o2*(*beta4)*n5);
			*(layS+j+1) = (uint8_t)(*(layS+j+1)-o1*(*beta3)*m6-o2*(*beta4)*n6);
			*(layS+j+width-1) = (uint8_t)(*(layS+j+width-1)-o1*(*beta3)*m7-o2*(*beta4)*n7);
			*(layS+j+width)= (uint8_t)(*(layS+j+width)-o1*(*beta3)*m8-o2*(*beta4)*n8);
			*(layS+j+width+1)= (uint8_t)(*(layS+j+width+1)-o1*(*beta3)*m9-o2*(*beta4)*n9);
		}
	}
	
	//calculate S layer
	/*
	for(i=H_Start;i<Image_Height-H_Start;i++)
	{
		layS = &S[i][0];
		SI_on = &Layers->SI_ON[i][0];
		SI_off = &Layers->SI_OFF[i][0];	
		for(j=W_Start;j<Image_Width-W_Start;j++)
		{
			Son = beta3*(*(SI_on+j));
			Soff = beta4*(*(SI_off+j));
			s = o1*Son +o2*Soff;
			*(layS+j) = (uint8_t)(*(layS+j)-s);
		}
	}
	*/
	
	//calculate Ce
	for(i=H_Start;i<Image_Height-H_Start;i++)
	{
		layG = &G[i][0];
		layS = &S[i][0];
		for(j=W_Start;j<Image_Width-W_Start;j++)
		{
			sum = *(layS+j)+*(layS+j-width)+*(layS+j-1)+*(layS+j+1)+*(layS+j+width)+*(layS+j-width-1)+*(layS+j-width+1)+*(layS+j+width-1)+*(layS+j+width+1);
			*(layG+j) = (uint8_t)(sum/9);
			if((*(layG+j))>t_max)
			t_max = *(layG+j);
		}
	}

	//calculate G
	scale = Delatc+(t_max*1.0)/Cw;
	for(i=H_Start;i<Image_Height-H_Start;i++)
	{
		layS = &S[i][0];
		layG = &G[i][0];
		for(j=W_Start;j<Image_Width-W_Start;j++)
		{
			*(layG+j) = (uint8_t)((*(layS+j))*(*(layG+j))*(1/scale));
			if(*(layG+j)<CdeTde)
				*(layG+j) = 0;
			*mp += *(layG+j);
		}
	}
	
	// sigmoid membrane potential
	*cur_smp = 1/(1+exp(-(*mp)/(ncell*Kspi)));
	
	//SFA mechanism
	if((*cur_smp-*pre_smp)<=Tsfa)
		*cur_lgmd_out = a3*(*(pre_lgmd_out)+*(cur_smp)-*(pre_smp));
	else
		*cur_lgmd_out = a3*(*cur_smp);
	if(*cur_lgmd_out<0.5)
		*cur_lgmd_out = 0.5;
	
	//detect collision
	*spike = floor(exp(a4*(*cur_lgmd_out-Tspi)));
	if((*spike)!=0)
	{
		*num += *spike;
		if(*num>=4)
		{
			*flag = 1;
			*num = 0;
		}
	}
	else
	{
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
	tmpOT=256*((hLGMD->Model->Results.LGMD_out[hLGMD->currentDiffImage])-0.5);
	
	GPIOD->ODR=(GPIOD->ODR&0xffffff00)+(~tmpOT);
	
	return 0;
}
