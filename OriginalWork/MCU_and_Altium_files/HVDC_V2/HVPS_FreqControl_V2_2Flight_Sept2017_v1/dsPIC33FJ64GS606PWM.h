/* 
 * File:   dsPIC33FJ64GS606PWM.h
 * Author: Yiou He
 * Comments: header file of dsPIC33FJ64GS606 PWM unit
 * Revision history: v1
 */

// This is a guard condition so that contents of this file are not included
// more than once.  
#ifndef XC_HEADER_PWM_H
#define	XC_HEADER_PWM_H

#include <xc.h> // include processor files - each processor file is guarded.  
#include "p33FJ64GS606.h"
#include "stdio.h"
#include "string.h"
#include "math.h"
#include "main_PWMoutputDefinitions.h"

void initPWM(int period, int duty, int deadtime);
void EnablePWM();    //Start PWM in complement mode
void DisablePWM();    //Disable PWM
//void PrechargePWM();  //Precharge the high side 5V 
void PreparePWM();  //Assign pin ownership to IO pins, precharge the highside 5V 
void ChangePD(int period, int duty, int deadtime);   //Change period, duty, deadtime
void initChangePD(int period, int duty, int deadtime);

#endif	/* XC_HEADER_PWM_H */