/* 
 * File:   dsPIC33FJ64GS606TIMER.h
 * Author: Yiou He
 * Comments: header file of the definitions in the main file
 * Revision history: v1
 */

// This is a guard condition so that contents of this file are not included
// more than once.  
#ifndef XC_HEADER_mainDefinition_H
#define	XC_HEADER_mainDefinition_H

#include <xc.h> // include processor files - each processor file is guarded.  
#include "p33FJ64GS606.h"
#include "stdio.h"
#include "string.h"

//define constants
#define initPWMTimerFreq 700.0 //unit in kHz
#define MinPWMTimerFreq 480.0 //unit in kHz
#define MaxPWMTimerFreq 700.0 //unit in kHz
#define initPWMdeadtime 45.0   //unit in ns
#define Timer12Freq 29.49/64.0*1000.0 //Timer 1 and 2's frequency in terms of 1 ms, times 2.0 adjusted for the 58 MHz
#define ADCtimerPeriod 1  //unit in ms, the period of ADC sampling is 1 ms
#define IC1timer2Period 2.5*Timer12Freq  //unit in ms, the period of the timer for RC signal 
#define HVOFFminperiod 0.7*Timer12Freq    //0.7,unit in ms, minimum pulse width of RC on signal 
#define HVOFFmaxperiod 1.1*Timer12Freq    //1.1,unit in ms, max pulse width of RC on signal
#define HVONminperiod 1.7*Timer12Freq   //1.7, unit in ms, minimum pulse width of RC off signal
#define HVONmaxperiod 2.2*Timer12Freq   //2.2,unit in ms, max pulse width of RC off signal
#define ADCScale 3.3/1023.0  //ADC scale
#define CalibratedScale 0.7  //calibrated input and output voltage scale factor
#define InputVScale 0.4327//0.4338//0.4596//0.4598 //0.3197   //,Calibrated on 06102017 //0.324    //Calibrated on March 2017
#define InputVOffset 2.0681//1.8397//1.0663//2.5139 //-0.0686  //,Calibrated on 06102017   //0.5296  //Calibrated on March 2017
#define InputIScale 0.0197//0.0213//0.021//0.0143    //,Calibrated on 06102017//0.0323 Calibrated on March 2017
#define InputIOffset -9.3941//-9.6692//-9.4371//-7.8522  //,Calibrated on 06102017//-17.031 Calibrated on March 2017
#define OutputVScale 76.574//81.363//82.66//82.936 //58.082  //,Calibrated on 06102017//56.577     //Calibrated on March 2017
#define OutputVOffset 534.23//367.65//225.48//-266.63 //342.97  //,Calibrated on 06102017//132.27    //Calibrated on March 2017
#define TempScale 1.0/(0.005*(100.0/24.9+1.0))*ADCScale
#define TempOffset 15.0  //change to 0 on June 3   //15 degree offset
#define AccScale 1.0/(100.0/24.9+1.0)*ADCScale/CalibratedScale
#define AccOffset 0.0
#define LogicBVScale 0.0059   //0.0065, Calibrated on Dec 04, 2016
#define LogicBVOffset -0.3132  //1053, Calibrated on Dec 05, 2016
//Sanity check value
#define MaxTemp 130.0
#define MaxInputV 250.0 //*CalibratedScale
#define MinInputV 150.0 //*CalibratedScale //145
#define MaxInputI 4.8 //*CalibratedScale  //if set as 4 A, can't even start
#define MaxOutputV 45000.0 //*CalibratedScale
#define MinOutputV 0 //30000.0 //
#define MinLogicBV 2.9 //*CalibratedScale
#define MaxLogicBV 4.5 //*CalibratedScale
#define ScaleFrom10To8 255.0/1023.0
#define MaxLoc 524287   //2^19 - 1, for FRAM address
//define desired control point
#define DesiredOutputV 40500.0 //*CalibratedScale //unit in Volts, the full desired output voltage,0.7 discount
//define timer
//#define EndTime 120 //unit in second

#endif	/* XC_HEADER_maindefinition_H */