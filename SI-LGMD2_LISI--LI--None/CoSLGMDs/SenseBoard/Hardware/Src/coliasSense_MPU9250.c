#include "coliasSense_MPU9250.h"
#include "coliasSense_interfaces.h"
#include "coliasSense_board.h"
#include "delay.h"
#include "math.h"

Tilt_detectTypedef TiltSense;

u8 IICwriteBits(CoS_sI2C_HandleTypeDef *hObj,u8 reg,u8 bitStart,u8 length,u8 data)
{
	u8 b;
	if (CoS_sI2C_ReadReg(hObj, reg, &b) != 0) 
	{
		u8 mask = (0xFF << (bitStart + 1)) | 0xFF >> ((8 - bitStart) + length - 1);
		data <<= (8 - length);
		data >>= (7 - bitStart);
		b &= mask;
		b |= data;
		return CoS_sI2C_WriteReg(hObj, reg, b);
	} else {
		return 0;
	}
}

u8 IICwriteBit(CoS_sI2C_HandleTypeDef *hObj, u8 reg, u8 bitNum, u8 data){
	u8 b;
	CoS_sI2C_ReadReg(hObj, reg, &b);
	b = (data != 0) ? (b | (1 << bitNum)) : (b & ~(1 << bitNum));

	return CoS_sI2C_WriteReg(hObj, reg, b);
}

int tilt_detection_init(Tilt_detectTypedef* hObj)
{
	printf("{initing MPU9250. Keep the board steady.}\r\n");
	hObj->hHMPU=hCoS->hHsI2C_MPU;
	//setClockSource(MPU9250_CLOCK_PLL_XGYRO)
	 IICwriteBits(hObj->hHMPU, MPU9250_RA_PWR_MGMT_1, MPU9250_PWR1_CLKSEL_BIT, MPU9250_PWR1_CLKSEL_LENGTH, MPU9250_CLOCK_PLL_XGYRO);
	//setFullScaleGyroRange(MPU9250_GYRO_FS_1000);
	 IICwriteBits(hObj->hHMPU, MPU9250_RA_GYRO_CONFIG, MPU9250_GCONFIG_FS_SEL_BIT, MPU9250_GCONFIG_FS_SEL_LENGTH, MPU9250_GYRO_FS_1000);
	//setFullScaleAccelRange(MPU9250_ACCEL_FS_2);
	IICwriteBits(hObj->hHMPU, MPU9250_RA_ACCEL_CONFIG, MPU9250_ACONFIG_AFS_SEL_BIT, MPU9250_ACONFIG_AFS_SEL_LENGTH, MPU9250_ACCEL_FS_2);
//setSleepEnabled(0);
	IICwriteBit(hObj->hHMPU, MPU9250_RA_PWR_MGMT_1, MPU9250_PWR1_SLEEP_BIT, 0);
	hObj->meanGy=0;
	hObj->meanGy_decay=0.9;
	hObj->meanGy_nGain=0.1;
	hObj->meanGy_threshold=-0.05;
	hObj->OT_threshold=10;
	hObj->do_LED_flash=1;
	MPU_selfCalibration_run(hObj);
	Tilt_detect_reset(hObj);
	hObj->CC.enable=1;
	
return 0;
}

float MPU9250_getAx(Tilt_detectTypedef* hObj)
{
	u8 accxouth,accxoutl;
	float x;
	short xx;
	CoS_sI2C_ReadReg(hObj->hHMPU,MPU9250_RA_ACCEL_XOUT_H,&accxouth);
	CoS_sI2C_ReadReg(hObj->hHMPU,MPU9250_RA_ACCEL_XOUT_L,&accxoutl);
	xx=((accxouth<<8)|accxoutl);
	x=xx/16384.0f;
	return x;
}
float MPU9250_getAy(Tilt_detectTypedef *hObj)
{
	u8 accxouth,accxoutl;
	float x;
	short xx;
	CoS_sI2C_ReadReg(hObj->hHMPU,MPU9250_RA_ACCEL_YOUT_H,&accxouth);
	CoS_sI2C_ReadReg(hObj->hHMPU,MPU9250_RA_ACCEL_YOUT_L,&accxoutl);
	xx=((accxouth<<8)|accxoutl);
	x=xx/16384.0f;
	return x;
}
float MPU9250_getAz(Tilt_detectTypedef* hObj)
{
	u8 accxouth,accxoutl;
	float x;
	short xx;
	CoS_sI2C_ReadReg(hObj->hHMPU,MPU9250_RA_ACCEL_ZOUT_H,&accxouth);
	CoS_sI2C_ReadReg(hObj->hHMPU,MPU9250_RA_ACCEL_ZOUT_L,&accxoutl);
	xx=((accxouth<<8)|accxoutl);
	x=xx/16384.0f;
	return x;
}

float MPU9250_getGx(Tilt_detectTypedef* hObj)
{
	u8 accxouth,accxoutl;
	float x;
	short xx;
	CoS_sI2C_ReadReg(hObj->hHMPU,MPU9250_RA_GYRO_XOUT_H,&accxouth);
	CoS_sI2C_ReadReg(hObj->hHMPU,MPU9250_RA_GYRO_XOUT_L,&accxoutl);
	xx=((accxouth<<8)|accxoutl);
	x=xx/656.0f;
	return x;
}
float MPU9250_getGy(Tilt_detectTypedef* hObj)
{
	u8 accxouth,accxoutl;
	float x;
	short xx;
	CoS_sI2C_ReadReg(hObj->hHMPU,MPU9250_RA_GYRO_YOUT_H,&accxouth);
	CoS_sI2C_ReadReg(hObj->hHMPU,MPU9250_RA_GYRO_YOUT_L,&accxoutl);
	xx=((accxouth<<8)|accxoutl);
	x=xx/656.0f;
	return x;
}
float MPU9250_getGz(Tilt_detectTypedef* hObj)
{
	u8 accxouth,accxoutl;
	float x;
	short xx;
	CoS_sI2C_ReadReg(hObj->hHMPU,MPU9250_RA_GYRO_ZOUT_H,&accxouth);
	CoS_sI2C_ReadReg(hObj->hHMPU,MPU9250_RA_GYRO_ZOUT_L,&accxoutl);
	xx=((accxouth<<8)|accxoutl);
	x=xx/656.0f;
	return x;
}
float MPU9250_getTemp(Tilt_detectTypedef* hObj)
{
	u8 accxouth,accxoutl;
	float x;
	short xx;
	CoS_sI2C_ReadReg(hObj->hHMPU,MPU9250_RA_TEMP_OUT_H,&accxouth);
	CoS_sI2C_ReadReg(hObj->hHMPU,MPU9250_RA_TEMP_OUT_L,&accxoutl);
	xx=((accxouth<<8)|accxoutl);
	x=xx/340.0f+36.53f;
	return x;
}

void IMU_getValues(Tilt_detectTypedef* hObj) 
{
	hObj->ax = MPU9250_getAx(hObj);
	hObj->ay = MPU9250_getAy(hObj);
	hObj->az = MPU9250_getAz(hObj);
	hObj->gx = MPU9250_getGx(hObj) - hObj->gxoffset;
	hObj->gy = MPU9250_getGy(hObj) - hObj->gyoffset;
	hObj->gz = MPU9250_getGz(hObj) - hObj->gzoffset;
	hObj->tempC = MPU9250_getTemp(hObj);
	hObj->gxsum+=hObj->gx;
	hObj->gysum+=hObj->gy;
	hObj->gzsum+=hObj->gz;
	hObj->sum_count++;
//    acc_gyro_mag[6] = 0;
//    acc_gyro_mag[7] = 0;
//    acc_gyro_mag[8] = 0;
}
int tilt_detection_run(void)
{
	u32 timer[2];
	float gnorm;
	if (!TiltSense.CC.enable)
		return 0;
	timer[0]=TOCin2;
	IMU_getValues(&TiltSense);
	timer[1]=TOCin2;
	gnorm=sqrt(TiltSense.ax*TiltSense.ax+TiltSense.ay*TiltSense.ay+TiltSense.az*TiltSense.az);
	
	TiltSense.meanGy=TiltSense.meanGy*TiltSense.meanGy_decay+TiltSense.meanGy_nGain*(TiltSense.ay-TiltSense.ayoffset);
	if (TiltSense.meanGy<TiltSense.meanGy_threshold)
		TiltSense.OT_count++;
	else
		TiltSense.OT_count=0;
	if (TiltSense.OT_count>TiltSense.OT_threshold && TiltSense.do_LED_flash)
		LED3_Toggle;
	if (TiltSense.CC.userdata)
	{

		printf("Gn%7.3f ",gnorm);
		printf("Ax%7.3f Ay%7.3f Az%7.3f ",TiltSense.ax,TiltSense.ay,TiltSense.az);
		printf("Gx%7.3f Gy%7.3f Gz%7.3f T %7.3f ",TiltSense.gx,TiltSense.gy,TiltSense.gz,TiltSense.tempC);
		printf("Sx%7.3f Sy%7.3f Sz%7.3f",TiltSense.gxsum,TiltSense.gysum,TiltSense.gzsum);

		printf("\r\n");
	}
	return 0;
}

void Tilt_detect_reset(Tilt_detectTypedef* hObj)
{
	hObj->gxsum=0;
	hObj->gysum=0;
	hObj->gzsum=0;
	hObj->sum_count=0;
}

int MPU_selfCalibration_run(Tilt_detectTypedef* hObj)
{
	u8 n,tasks=50,tmpEnable;
	float tmpResult[12]={0};
	
	tmpEnable=hObj->CC.enable;
	hObj->CC.enable=0;
	Tilt_detect_reset(hObj);
	
	hObj->gxoffset=0;
	hObj->gyoffset=0;
	hObj->gzoffset=0;
	
	for (n=0;n<tasks;n++)
	{
		delay_ms(10);
		IMU_getValues(hObj);
		tmpResult[0]+=hObj->ax;
		tmpResult[1]+=hObj->ay;
		tmpResult[2]+=hObj->az;
		tmpResult[3]+=hObj->gx;
		tmpResult[4]+=hObj->gy;
		tmpResult[5]+=hObj->gz;

	}
	hObj->axoffset=tmpResult[0]/tasks;
	hObj->ayoffset=tmpResult[1]/tasks;
	hObj->azoffset=tmpResult[2]/tasks;
	hObj->gxoffset=tmpResult[3]/tasks;
	hObj->gyoffset=tmpResult[4]/tasks;
	hObj->gzoffset=tmpResult[5]/tasks;
	printf("MPU9250 Calibrated. ");
	printf("Ax:%7.4f, Ay:%7.4f, Az:%7.4f, ",hObj->axoffset,hObj->ayoffset,hObj->azoffset);
	printf("Gx:%7.4f, Gy:%7.4f, Gz:%7.4f\r\n",hObj->gxoffset,hObj->gyoffset,hObj->gzoffset);
	hObj->CC.enable=tmpEnable;
	return tmpEnable;
}







