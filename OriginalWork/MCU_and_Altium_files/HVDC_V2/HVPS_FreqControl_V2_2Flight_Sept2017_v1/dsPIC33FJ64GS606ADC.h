/* 
 * File:   dsPIC33FJ64GS606ADC.h
 * Author: Yiou He
 * Comments: header file of the dsPIC33FJ64GS606 ADC unit
 * Revision history: v1
 */

// This is a guard condition so that contents of this file are not included
// more than once.  
#ifndef XC_HEADER_ADC_H
#define	XC_HEADER_ADC_H

#include <xc.h> // include processor files - each processor file is guarded.  
#include "p33FJ64GS606.h"
#include "stdio.h"
#include "string.h"

void initADC(int period);
void EnableADC();
void DisableADC();

#endif	/* XC_HEADER_ADC_H*/