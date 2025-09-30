/* 
 * File:   dsPIC33FJ64GS606IO.h
 * Author: Yiou He
 * Comments: header file of dsPIC33FJ64GS606 IO file
 * Revision history: v1
 */

// This is a guard condition so that contents of this file are not included
// more than once.  
#ifndef XC_HEADER_IO_H
#define	XC_HEADER_IO_H

#include <xc.h> // include processor files - each processor file is guarded.  
#include "p33FJ64GS606.h"
#include "stdio.h"
#include "string.h"

void initIO();
void DischargeCaps();    //Start to discharge the input caps
void Precharge();         //Precharge the input caps
void NormalOperation();  //Connect 200V to the switches
void HVon(); //HV on, LED light on.
#endif	/* XC_HEADER_IO_H */