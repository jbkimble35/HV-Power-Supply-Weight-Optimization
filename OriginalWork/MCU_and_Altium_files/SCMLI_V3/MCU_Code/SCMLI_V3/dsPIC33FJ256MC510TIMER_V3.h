/* 
 * File:   dsPIC33FJ64GS606TIMER_test.h
 * Author: Suzanne O'Meara
 * Comments: header file of the dsPIC33FJ64GS606 TIMER unit
 * Revision history: v1
 */

// This is a guard condition so that contents of this file are not included
// more than once.  
#ifndef XC_HEADER_TIMER_H
#define	XC_HEADER_TIMER_H

#include <xc.h> // include processor files - each processor file is guarded.  
#include "p33FJ256MC510.h"
#include "stdio.h"
#include "string.h"
#include "main_SCMLI_V3.h"

void initTimer1(int SwitchPeriod);
void EnableTimer1();
void DisableTimer1();
void initTimer2(int StartPeriod);
void EnableTimer2();
void DisableTimer2();

#endif	/* XC_HEADER_TIMER_H */