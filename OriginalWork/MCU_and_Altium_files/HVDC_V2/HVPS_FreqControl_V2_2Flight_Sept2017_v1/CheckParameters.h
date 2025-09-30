/* 
 * File:   dsPIC33FJ64GS606TIMER.h
 * Author: Yiou He
 * Comments: header file of the check functions, including preflight check and during flight ADC value check
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
#include "main_PWMoutputDefinitions.h"

int CheckMechamical(double AccX,double AccY,double AccZ,double Temperature);  //Check to see if all values are OK
int CheckElectrical(double InputVoltage,double InputCurrent,double OutputVoltage,double LogicBV);  //Check to see if all values are OK

#endif	/* XC_HEADER_maindefinition_H */