/* 
 * File:   dsPIC33FJ64GS606TIMER.h
 * Author: Yiou He
 * Comments: header file of the dsPIC33FJ64GS606 TIMER unit
 * Revision history: v1
 */

// This is a guard condition so that contents of this file are not included
// more than once.  
#ifndef XC_HEADER_TIMER_H
#define	XC_HEADER_TIMER_H

#include <xc.h> // include processor files - each processor file is guarded.  
#include "p33FJ64GS606.h"
#include "stdio.h"
#include "string.h"
#include "main_PWMoutputDefinitions.h"

void initTimer1(int ADCperiod);
void EnableTimer1();
void initTimer2();
void EnableTimer2();
void initTimer45(int EndTimeMSB,int EndTimeLSB);
void EnableTimer45();
void DisableTimer45();
void initIC1();

#endif	/* XC_HEADER_TIMER_H */